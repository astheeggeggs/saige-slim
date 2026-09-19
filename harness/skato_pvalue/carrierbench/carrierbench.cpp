// carrierbench.cpp -- dense (as shipped) vs carrier-space assembly of the
// SKAT-O region variance matrix (VarMat; the q x q Phi is a weighted subset of it).
//
// Standalone. Touches nothing in src/ or R/. See NOTES.md for the census facts
// each modelling choice is anchored to.
//
// Arms:
//   A_dense  : exactly what src/Main.cpp does today -- per-marker dense gtilde via
//              getadjGFast's split form (carrier-space forward gather, dense N x p
//              back-projection, dense N subtraction), written as STRIDED rows of a
//              column-major (100 x N) P1 and CONTIGUOUS columns of an (N x 100) P2,
//              chunked at 100 markers, copied by value into a chunk store, then the
//              nchunks x nchunks dgemm triangle of Main.cpp:4160-4214.
//   A_syrk   : the same assembly but one (N x q) matrix A = sqrt(w).*gtilde and a
//              single cblas_dsyrk. Not what ships; it is the cheap dense lever that
//              any carrier-space prize must beat, not just dgemm.
//   B_merge  : the proposed change as specified -- markers held as (idx,val) over
//              carriers; Gram entry by sorted-list INTERSECTION merge, O(nnz_j+nnz_l);
//              U = X'W G by gather, O(p*nnz); Phi = A - U'(X'WX)^-1 U.
//   B_samp   : same maths, sample-major (SpGEMM) accumulation instead of pairwise
//              merges -- cost sum_i c_i^2 rather than q^2 * nnz. Included so that a
//              B_merge loss cannot be mistaken for a carrier-space loss.
//
// Identity both arms compute (W = diag(mu2*tau) = Sigma^-1 on the diagonal-Sigma path):
//   gtilde_j' W gtilde_l = g_j'W g_l - (X'W g_j)' (X'W X)^-1 (X'W g_l)

#include <cblas.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <random>
#include <chrono>
#include <algorithm>

extern "C" void openblas_set_num_threads(int);

static const int P = 7;            // census: p = 7 (minus1 age sex1 PC1..PC4)
static const int CHUNKCAP = 100;   // census: markers_per_chunk_in_groupTest default

typedef std::chrono::steady_clock clk;
static double secs(clk::time_point a, clk::time_point b) {
  return std::chrono::duration<double>(b - a).count();
}
static double median3(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  return v[v.size() / 2];
}

struct Data {
  long N; int q;
  std::vector<std::vector<int> >    idx;   // per marker, sorted carrier indices
  std::vector<std::vector<double> > val;   // per marker, genotype at those carriers
  std::vector<double> w;        // N     diag(Sigma^-1)
  std::vector<double> XV;       // P x N col-major  (X'W)
  std::vector<double> XXVXinv;  // N x P col-major  (X (X'WX)^-1)
  std::vector<double> XWXinv;   // P x P
  long totnnz;
};

static void invert_pp(const double* in, double* out) {
  double a[P][2 * P];
  for (int i = 0; i < P; i++) {
    for (int j = 0; j < P; j++) { a[i][j] = in[i + P * j]; a[i][P + j] = (i == j) ? 1.0 : 0.0; }
  }
  for (int c = 0; c < P; c++) {
    int piv = c;
    for (int r = c + 1; r < P; r++) if (fabs(a[r][c]) > fabs(a[piv][c])) piv = r;
    for (int j = 0; j < 2 * P; j++) std::swap(a[c][j], a[piv][j]);
    double d = a[c][c];
    for (int j = 0; j < 2 * P; j++) a[c][j] /= d;
    for (int r = 0; r < P; r++) {
      if (r == c) continue;
      double f = a[r][c];
      if (f == 0.0) continue;
      for (int j = 0; j < 2 * P; j++) a[r][j] -= f * a[c][j];
    }
  }
  for (int i = 0; i < P; i++) for (int j = 0; j < P; j++) out[i + P * j] = a[i][P + j];
}

static void makeData(Data& d, long N, int q, double frac, unsigned seed) {
  d.N = N; d.q = q;
  std::mt19937_64 rng(seed);
  std::normal_distribution<double> nrm(0.0, 1.0);
  std::uniform_real_distribution<double> uni01(0.0, 1.0);
  std::uniform_int_distribution<long> unii(0, N - 1);

  // W: binary trait, prevalence ~1%, w = mu(1-mu)*tau
  d.w.resize(N);
  for (long i = 0; i < N; i++) { double mu = 0.002 + 0.048 * uni01(rng); d.w[i] = mu * (1 - mu); }

  // X: N x P col-major, col 0 = -1 (the 'minus1' intercept), rest standard normal
  std::vector<double> X((size_t)N * P);
  for (long i = 0; i < N; i++) X[i] = -1.0;
  for (int k = 1; k < P; k++) for (long i = 0; i < N; i++) X[i + (size_t)N * k] = nrm(rng);

  double XWX[P * P];
  for (int a = 0; a < P; a++) for (int b = 0; b < P; b++) {
    double s = 0.0;
    for (long i = 0; i < N; i++) s += X[i + (size_t)N * a] * d.w[i] * X[i + (size_t)N * b];
    XWX[a + P * b] = s;
  }
  d.XWXinv.assign(P * P, 0.0);
  invert_pp(XWX, &d.XWXinv[0]);

  d.XV.assign((size_t)P * N, 0.0);
  for (long i = 0; i < N; i++) for (int k = 0; k < P; k++) d.XV[k + (size_t)P * i] = X[i + (size_t)N * k] * d.w[i];
  d.XXVXinv.assign((size_t)N * P, 0.0);
  for (int k = 0; k < P; k++) for (long i = 0; i < N; i++) {
    double s = 0.0;
    for (int m = 0; m < P; m++) s += X[i + (size_t)N * m] * d.XWXinv[m + P * k];
    d.XXVXinv[i + (size_t)N * k] = s;
  }

  long nnz = (long)llround(frac * (double)N);
  if (nnz < 1) nnz = 1;
  if (nnz > N) nnz = N;
  d.idx.assign(q, std::vector<int>());
  d.val.assign(q, std::vector<double>());
  std::vector<char> seen(N, 0);
  d.totnnz = 0;
  for (int j = 0; j < q; j++) {
    long got = 0;
    std::vector<int> picked; picked.reserve(nnz);
    while (got < nnz) { long r = unii(rng); if (!seen[r]) { seen[r] = 1; picked.push_back((int)r); got++; } }
    std::sort(picked.begin(), picked.end());
    d.idx[j] = picked;
    d.val[j].resize(nnz);
    for (long e = 0; e < nnz; e++) { d.val[j][e] = (uni01(rng) < 0.05) ? 2.0 : 1.0; seen[picked[e]] = 0; }
    d.totnnz += nnz;
  }
}

// ---------------------------------------------------------------- A_dense
struct DenseTimes { double asm_s, gram_s; double gram_flops; };

static void armDense(const Data& d, std::vector<double>& Phi, DenseTimes& t) {
  const long N = d.N; const int q = d.q;
  int nch = (q + CHUNKCAP - 1) / CHUNKCAP;
  std::vector<int> mrows(nch), off(nch);
  for (int c = 0; c < nch; c++) { off[c] = c * CHUNKCAP; mrows[c] = std::min(CHUNKCAP, q - off[c]); }

  std::vector<double> P1((size_t)CHUNKCAP * N), P2((size_t)N * CHUNKCAP), gt(N);
  std::vector<std::vector<double> > P1st(nch), P2st(nch);
  for (int c = 0; c < nch; c++) { P1st[c].resize((size_t)mrows[c] * N); P2st[c].resize((size_t)N * mrows[c]); }
  Phi.assign((size_t)q * q, 0.0);

  clk::time_point t0 = clk::now();
  for (int c = 0; c < nch; c++) {
    for (int r = 0; r < mrows[c]; r++) {
      int j = off[c] + r;
      const std::vector<int>& ix = d.idx[j]; const std::vector<double>& vv = d.val[j];
      double u[P] = {0, 0, 0, 0, 0, 0, 0};
      for (size_t e = 0; e < ix.size(); e++) {          // getadjGFast forward half: carrier space already
        const double* col = &d.XV[(size_t)P * ix[e]];
        double v = vv[e];
        for (int k = 0; k < P; k++) u[k] += col[k] * v;
      }
      memset(&gt[0], 0, sizeof(double) * (size_t)N);      // GVec.zeros() + scatter
      for (size_t e = 0; e < ix.size(); e++) gt[ix[e]] = vv[e];
      cblas_dgemv(CblasColMajor, CblasNoTrans, (int)N, P, -1.0, &d.XXVXinv[0], (int)N, u, 1, 1.0, &gt[0], 1);
      double* p1 = &P1[0]; double* p2 = &P2[(size_t)N * r];
      for (long i = 0; i < N; i++) p1[r + (size_t)CHUNKCAP * i] = gt[i];   // STRIDED row write
      for (long i = 0; i < N; i++) p2[i] = gt[i] * d.w[i];                 // contiguous col write
    }
    // RegionChunkStore::add takes the matrix BY VALUE -- copy, trimmed to mrows
    for (long i = 0; i < N; i++)
      for (int r = 0; r < mrows[c]; r++) P1st[c][r + (size_t)mrows[c] * i] = P1[r + (size_t)CHUNKCAP * i];
    memcpy(&P2st[c][0], &P2[0], sizeof(double) * (size_t)N * mrows[c]);
  }
  clk::time_point t1 = clk::now();

  double fl = 0.0;
  std::vector<double> blk;
  for (int ci = 0; ci < nch; ci++) {
    for (int cj = ci; cj < nch; cj++) {
      int mi = mrows[ci], mj = mrows[cj];
      blk.resize((size_t)mi * mj);
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, mi, mj, (int)N,
                  1.0, &P1st[ci][0], mi, &P2st[cj][0], (int)N, 0.0, &blk[0], mi);
      fl += 2.0 * mi * mj * (double)N;
      for (int b = 0; b < mj; b++) for (int a = 0; a < mi; a++) {
        double v = blk[a + (size_t)mi * b];
        Phi[(off[ci] + a) + (size_t)q * (off[cj] + b)] = v;
        if (ci != cj) Phi[(off[cj] + b) + (size_t)q * (off[ci] + a)] = v;   // the .t() write-back
      }
    }
  }
  clk::time_point t2 = clk::now();
  t.asm_s = secs(t0, t1); t.gram_s = secs(t1, t2); t.gram_flops = fl;
}

// ---------------------------------------------------------------- A_syrk
static void armSyrk(const Data& d, std::vector<double>& Phi, DenseTimes& t) {
  const long N = d.N; const int q = d.q;
  std::vector<double> A((size_t)N * q), gt(N), sw(N);
  for (long i = 0; i < N; i++) sw[i] = sqrt(d.w[i]);
  Phi.assign((size_t)q * q, 0.0);

  clk::time_point t0 = clk::now();
  for (int j = 0; j < q; j++) {
    const std::vector<int>& ix = d.idx[j]; const std::vector<double>& vv = d.val[j];
    double u[P] = {0, 0, 0, 0, 0, 0, 0};
    for (size_t e = 0; e < ix.size(); e++) {
      const double* col = &d.XV[(size_t)P * ix[e]]; double v = vv[e];
      for (int k = 0; k < P; k++) u[k] += col[k] * v;
    }
    memset(&gt[0], 0, sizeof(double) * (size_t)N);
    for (size_t e = 0; e < ix.size(); e++) gt[ix[e]] = vv[e];
    cblas_dgemv(CblasColMajor, CblasNoTrans, (int)N, P, -1.0, &d.XXVXinv[0], (int)N, u, 1, 1.0, &gt[0], 1);
    double* a = &A[(size_t)N * j];
    for (long i = 0; i < N; i++) a[i] = gt[i] * sw[i];
  }
  clk::time_point t1 = clk::now();
  cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, q, (int)N, 1.0, &A[0], (int)N, 0.0, &Phi[0], q);
  for (int a2 = 0; a2 < q; a2++) for (int b = 0; b < a2; b++) Phi[a2 + (size_t)q * b] = Phi[b + (size_t)q * a2];
  clk::time_point t2 = clk::now();
  t.asm_s = secs(t0, t1); t.gram_s = secs(t1, t2); t.gram_flops = (double)q * q * (double)N;
}

// ------------------------------------------------- carrier-space shared tail
static void carrierCorrect(const Data& d, const std::vector<double>& U, std::vector<double>& Phi) {
  const int q = d.q;
  std::vector<double> T((size_t)P * q, 0.0);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, P, q, P,
              1.0, &d.XWXinv[0], P, &U[0], P, 0.0, &T[0], P);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, q, q, P,
              -1.0, &U[0], P, &T[0], P, 1.0, &Phi[0], q);
}

// ---------------------------------------------------------------- B_merge
struct CarTimes { double u_s, gram_s, corr_s; double steps; };

static void armMerge(const Data& d, std::vector<double>& Phi, CarTimes& t) {
  const int q = d.q;
  std::vector<double> U((size_t)P * q, 0.0);
  Phi.assign((size_t)q * q, 0.0);

  clk::time_point t0 = clk::now();
  for (int j = 0; j < q; j++) {
    const std::vector<int>& ix = d.idx[j]; const std::vector<double>& vv = d.val[j];
    double* u = &U[(size_t)P * j];
    for (size_t e = 0; e < ix.size(); e++) {
      const double* col = &d.XV[(size_t)P * ix[e]]; double v = vv[e];
      for (int k = 0; k < P; k++) u[k] += col[k] * v;
    }
  }
  clk::time_point t1 = clk::now();
  double steps = 0.0;
  const double* wp = &d.w[0];
  for (int j = 0; j < q; j++) {
    const int* ia = &d.idx[j][0]; const double* va = &d.val[j][0]; size_t na = d.idx[j].size();
    for (int l = j; l < q; l++) {
      const int* ib = &d.idx[l][0]; const double* vb = &d.val[l][0]; size_t nb = d.idx[l].size();
      size_t a = 0, b = 0; double s = 0.0;
      while (a < na && b < nb) {
        int xa = ia[a], xb = ib[b];
        if (xa < xb) a++;
        else if (xb < xa) b++;
        else { s += wp[xa] * va[a] * vb[b]; a++; b++; }
      }
      steps += (double)(a + b);
      Phi[j + (size_t)q * l] = s;
      if (j != l) Phi[l + (size_t)q * j] = s;
    }
  }
  clk::time_point t2 = clk::now();
  carrierCorrect(d, U, Phi);
  clk::time_point t3 = clk::now();
  t.u_s = secs(t0, t1); t.gram_s = secs(t1, t2); t.corr_s = secs(t2, t3); t.steps = steps;
}

// ---------------------------------------------------------------- B_samp
static void armSamp(const Data& d, std::vector<double>& Phi, CarTimes& t) {
  const long N = d.N; const int q = d.q;
  std::vector<double> U((size_t)P * q, 0.0);
  Phi.assign((size_t)q * q, 0.0);

  clk::time_point t0 = clk::now();
  for (int j = 0; j < q; j++) {
    const std::vector<int>& ix = d.idx[j]; const std::vector<double>& vv = d.val[j];
    double* u = &U[(size_t)P * j];
    for (size_t e = 0; e < ix.size(); e++) {
      const double* col = &d.XV[(size_t)P * ix[e]]; double v = vv[e];
      for (int k = 0; k < P; k++) u[k] += col[k] * v;
    }
  }
  clk::time_point t1 = clk::now();
  // transpose carrier lists to sample-major CSR (work the change would also have to do)
  std::vector<long> cnt(N + 1, 0);
  for (int j = 0; j < q; j++) { const std::vector<int>& ix = d.idx[j]; for (size_t e = 0; e < ix.size(); e++) cnt[ix[e] + 1]++; }
  for (long i = 0; i < N; i++) cnt[i + 1] += cnt[i];
  std::vector<int> cm(d.totnnz); std::vector<double> cv(d.totnnz); std::vector<long> pos(cnt.begin(), cnt.end() - 1);
  for (int j = 0; j < q; j++) {
    const std::vector<int>& ix = d.idx[j]; const std::vector<double>& vv = d.val[j];
    for (size_t e = 0; e < ix.size(); e++) { long p = pos[ix[e]]++; cm[p] = j; cv[p] = vv[e]; }
  }
  double steps = 0.0;
  for (long i = 0; i < N; i++) {
    long s0 = cnt[i], s1 = cnt[i + 1];
    if (s1 == s0) continue;
    double wi = d.w[i];
    for (long a = s0; a < s1; a++) {
      double wa = wi * cv[a]; int ja = cm[a];
      for (long b = a; b < s1; b++) Phi[ja + (size_t)q * cm[b]] += wa * cv[b];  // upper triangle
    }
    double c = (double)(s1 - s0); steps += c * (c + 1) * 0.5;
  }
  for (int a = 0; a < q; a++) for (int b = 0; b < a; b++) Phi[a + (size_t)q * b] = Phi[b + (size_t)q * a];
  clk::time_point t2 = clk::now();
  carrierCorrect(d, U, Phi);
  clk::time_point t3 = clk::now();
  t.u_s = secs(t0, t1); t.gram_s = secs(t1, t2); t.corr_s = secs(t2, t3); t.steps = steps;
}

// ---------------------------------------------------------------- compare
static void cmpPhi(const std::vector<double>& A, const std::vector<double>& B, int q,
                   double& relmax, double& relfro) {
  double mx = 0.0, dmx = 0.0, fa = 0.0, fd = 0.0;
  for (size_t i = 0; i < (size_t)q * q; i++) {
    double a = A[i], b = B[i], dd = fabs(a - b);
    if (fabs(a) > mx) mx = fabs(a);
    if (dd > dmx) dmx = dd;
    fa += a * a; fd += dd * dd;
  }
  relmax = (mx > 0) ? dmx / mx : dmx;
  relfro = (fa > 0) ? sqrt(fd / fa) : sqrt(fd);
}

int main(int argc, char** argv) {
  openblas_set_num_threads(1);
  std::string mode = (argc > 1) ? argv[1] : "check";
  int reps = (argc > 2) ? atoi(argv[2]) : 3;

  if (mode == "check") {
    // CORRECTNESS CONTROL, before any timing.
    long Ns[3] = {3000, 3000, 50000};
    int qs[3] = {11, 30, 11};
    double fr[3] = {5e-3, 2e-2, 5e-4};
    printf("# correctness control: arms must agree on Phi to 1e-10 relative\n");
    printf("N\tq\tnnz/N\tpair\trel_max\trel_fro\n");
    int bad = 0;
    for (int c = 0; c < 3; c++) {
      Data d; makeData(d, Ns[c], qs[c], fr[c], 1234u + c);
      std::vector<double> Pa, Ps, Pm, Pv; DenseTimes ta, ts; CarTimes tm, tv;
      armDense(d, Pa, ta); armSyrk(d, Ps, ts); armMerge(d, Pm, tm); armSamp(d, Pv, tv);
      const char* nm[3] = {"A_syrk-vs-A_dense", "B_merge-vs-A_dense", "B_samp-vs-A_dense"};
      const std::vector<double>* ot[3] = {&Ps, &Pm, &Pv};
      for (int k = 0; k < 3; k++) {
        double rm, rf; cmpPhi(Pa, *ot[k], d.q, rm, rf);
        printf("%ld\t%d\t%.1e\t%s\t%.3e\t%.3e\n", Ns[c], qs[c], fr[c], nm[k], rm, rf);
        if (!(rm < 1e-10)) bad++;
      }
    }
    printf("# VERDICT: %s\n", bad ? "DISAGREE -- timing would be meaningless" : "AGREE (all < 1e-10 relative)");
    return bad ? 2 : 0;
  }

  if (mode == "detail") {
    // component breakdown + cost per elementary operation, for the realistic points
    long dN[6]   = {50000, 50000, 50000, 50000, 400000, 50000};
    int  dq[6]   = {11, 30, 134, 134, 30, 134};
    double df[6] = {5e-4, 5e-4, 5e-4, 5e-1, 5e-4, 2e-2};
    printf("N\tq\tfrac\tnnz\tdense_asm_ms\tdense_gram_ms\tgram_gflops\tns_per_dgemm_flop\t"
           "merge_u_ms\tmerge_gram_ms\tmerge_corr_ms\tmerge_steps\tns_per_merge_step\t"
           "samp_u_ms\tsamp_gram_ms\tsamp_steps\tsamp_floor_ns_per_sample\n");
    for (int c = 0; c < 6; c++) {
      Data d; makeData(d, dN[c], dq[c], df[c], 555u + c);
      std::vector<double> Pa, Pm, Pv; DenseTimes ta; CarTimes tm, tv;
      std::vector<double> va, vg, vmu, vmg, vmc, vvu, vvg;
      double steps_m = 0, steps_v = 0, fl = 0;
      for (int r = 0; r <= 3; r++) {
        armDense(d, Pa, ta); armMerge(d, Pm, tm); armSamp(d, Pv, tv);
        if (r == 0) continue;
        va.push_back(ta.asm_s); vg.push_back(ta.gram_s);
        vmu.push_back(tm.u_s); vmg.push_back(tm.gram_s); vmc.push_back(tm.corr_s);
        vvu.push_back(tv.u_s); vvg.push_back(tv.gram_s);
        steps_m = tm.steps; steps_v = tv.steps; fl = ta.gram_flops;
      }
      double mg = median3(vg), mmg = median3(vmg), mvg = median3(vvg);
      printf("%ld\t%d\t%.1e\t%ld\t%.3f\t%.3f\t%.2f\t%.4f\t%.3f\t%.3f\t%.3f\t%.4g\t%.2f\t%.3f\t%.3f\t%.4g\t%.2f\n",
             dN[c], dq[c], df[c], (long)(d.totnnz / d.q),
             median3(va) * 1e3, mg * 1e3, fl / mg / 1e9, mg * 1e9 / fl,
             median3(vmu) * 1e3, mmg * 1e3, median3(vmc) * 1e3, steps_m, mmg * 1e9 / steps_m,
             median3(vvu) * 1e3, mvg * 1e3, steps_v, mvg * 1e9 / (double)dN[c]);
      fflush(stdout);
    }
    return 0;
  }

  // sweep
  long Ns[3] = {3000, 50000, 400000};
  int qs[4] = {4, 11, 30, 134};
  double frs[6] = {5e-5, 5e-4, 1.2e-3, 5e-3, 2e-2, 0.5};
  printf("N\tq\tfrac\tnnz\tdense_asm_ms\tdense_gram_ms\tdense_ms\tsyrk_ms\tmerge_ms\tsamp_ms\tdense_gflops\tratio_merge\tratio_samp\tratio_merge_vs_syrk\tpred_flopratio\n");
  for (int a = 0; a < 3; a++) {
    for (int b = 0; b < 4; b++) {
      for (int c = 0; c < 6; c++) {
        long N = Ns[a]; int q = qs[b]; double frac = frs[c];
        long nnz = (long)llround(frac * (double)N); if (nnz < 1) nnz = 1;
        Data d; makeData(d, N, q, frac, 99u + 7 * a + 13 * b + 29 * c);
        std::vector<double> Pa, Ps, Pm, Pv; DenseTimes ta, ts; CarTimes tm, tv;
        std::vector<double> vda, vdg, vs, vm, vv;
        for (int r = 0; r <= reps; r++) {   // r == 0 is a warmup, discarded
          armDense(d, Pa, ta); armSyrk(d, Ps, ts); armMerge(d, Pm, tm); armSamp(d, Pv, tv);
          if (r == 0) continue;
          vda.push_back(ta.asm_s); vdg.push_back(ta.gram_s);
          vs.push_back(ts.asm_s + ts.gram_s);
          vm.push_back(tm.u_s + tm.gram_s + tm.corr_s);
          vv.push_back(tv.u_s + tv.gram_s + tv.corr_s);
        }
        double mda = median3(vda), mdg = median3(vdg), ms = median3(vs), mm = median3(vm), mv = median3(vv);
        double dense = mda + mdg;
        double gf = ta.gram_flops / mdg / 1e9;
        printf("%ld\t%d\t%.1e\t%ld\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t%.3f\t%.3f\t%.3f\t%.1f\n",
               N, q, frac, nnz, mda * 1e3, mdg * 1e3, dense * 1e3, ms * 1e3, mm * 1e3, mv * 1e3,
               gf, dense / mm, dense / mv, ms / mm, (double)N / (double)nnz);
        fflush(stdout);
      }
    }
  }
  return 0;
}
