// Phase G2 (2026-09-14): the types mainRegionInCPP used to declare at
// function scope, moved here VERBATIM (namespace Region). Nothing else moved:
// no Work struct over the frame's 89 locals, no phase lift, no new
// translation unit -- src/Main.cpp includes this file AFTER its file-scope
// helpers (MarkerStatePool, MarkerWorkerState, RefuseMainThreadGuard,
// RegionChunkStore, the SAIGE_TSAN_* macros), which these types use as
// complete types, so this header is a section of Main.cpp lifted out for
// readability and is not includable on its own. See
// docs/state/handoff-2026-09-14-codebase-audit.md section 6.2 for why the
// rest was refused (armadillo's partial-accumulator unrolling makes a TU
// split a 1-ULP hazard for the printed goldens).
#ifndef SAIGE_REGION_LOOP_HPP
#define SAIGE_REGION_LOOP_HPP

namespace Region {

// E4 (2026-09-14): the per-set table, ONE record per weight set. It used to
// be four parallel arrays over three index domains (setKindVec/setColVec by
// set, betaDistVec by the set's beta ROW, erAllowedVec by set with its
// vector<bool> constructor default standing in for "ER allowed" on beta and
// No_Weights sets). Here every field is explicit for every kind. The
// per-kind semantics are design doc D2 / the table in the comment above.
//   erAllowed: the legacy (!isWeightCustomized || isEqualWeights) predicate
//   per set -- true for beta and No_Weights, the all-ones test for custom,
//   FlexRV::isUnitRow for a transform set (only unw x unw is identically 1;
//   a zero weight is not a unit weight, so lof_unw is NOT exempt).
struct WeightSet {
  int kind = -1;
  unsigned int col = 0;                       // custom column / transform column / beta row
  boost::math::beta_distribution<> dist;      // beta sets only
  bool erAllowed = true;
};

// The per-cell test records (passes B and D) split the same way: the
// genotype side (URCell / BurdenCell below, persistent across regions) is
// shared, the model side lives here.
// E5 (2026-09-14): the 14 result fields a cell test produces, ONCE.
// URCellModel, BurdenCellModel and RegionPrepModel repeated them
// field-for-field; inheritance (not composition) so every existing
// `prm.Beta` / `cm.pval` spelling stays as it is. The G1 rows are NOT in
// the base on purpose: URCellModel::G1 holds the raw G1tilde_P_G2tilde_Vec
// while RegionPrepModel::G1row holds the already-weighted row -- one name
// for two quantities is the LEDGER #118/#133/#135 shape.
struct CellResult {
  double Beta = 0.0, seBeta = 0.0, Tstat = 0.0, varT = 0.0, gy = 0.0;
  std::string pval, pval_noSPA;
  bool isSPAConverge = false;
  double Beta_c = 0.0, seBeta_c = 0.0, Tstat_c = 0.0, varT_c = 0.0;
  std::string pval_c, pval_noSPA_c;
};

struct URCellModel : CellResult {
  SAIGE::TestContext ctx;
  bool isER = false;
  arma::rowvec G1;               // G1tilde_P_G2tilde_Vec, conditional runs only
};

struct BurdenCellModel : CellResult {
  SAIGE::TestContext ctx;
  bool isER = false;
};

struct RegionModelWork {
  const SAIGE::SAIGEClass * model = nullptr;
  int modelIndex = 0;
  std::string traitType;
  bool isSingleVarianceRatio = true;
  // the per-index outputs (q = q0 markers + the ultra-rare pseudo-markers)
  std::vector<double> BetaVec, seBetaVec, TstatVec, TstatVec_flip, gyVec, varTVec;
  std::vector<std::string> pvalVec, pvalNAVec;
  std::vector<unsigned char> isSPAConvergeVec;   // not vector<bool>: bit-packed writes race across markers
  std::vector<double> Beta_cVec, seBeta_cVec, Tstat_cVec, varT_cVec;
  std::vector<std::string> pval_cVec, pvalNA_cVec;
  arma::mat G1tilde_P_G2tilde_Weighted_Mat;
  arma::vec MACCase_GroupVec, MACControl_GroupVec;
  // binary models: sample -> its position in m_case_indices / m_ctrl_indices
  // (both sorted, from arma::find), -1 otherwise; CarrierStore::paritySumRanked
  std::vector<int32_t> caseRank, ctrlRank;
  std::vector<double> AF_caseVec, AF_ctrlVec, N_case_homVec, N_ctrl_hetVec, N_case_hetVec, N_ctrl_homVec;
  std::vector<uint32_t> N_caseVec, N_ctrlVec, N_Vec;
  arma::mat genoSumMat;        // SKAT-O only: this model's residualised copy of the shared matrix
  arma::mat blockP1, blockP2;  // n x blockSize: column k = sqrt(vr) * gtilde / sqrt(vr) * P2Vec of marker k
  arma::mat VarMat;
  std::unique_ptr<RegionChunkStore> chunkStore;
  std::vector<URCellModel> urOut;
  std::vector<BurdenCellModel> burdenOut;
  double cctpval = 1.0, cctpval_cond = 1.0;
  bool iswriteOutput = false;
  Rcpp::List OutList;
};

// Each marker's work is split the way the marker path's was in S4:
//   prepareAndTest(i, k, scratch, reader)  read + impute/flip + QC + weights,
//       then (rare markers) the test context and getMarkerPval. Everything
//       it produces lands in preps[k] (one record per marker in the block)
//       and in column k of blockP1/blockP2. It reads the model and never
//       writes it, touches no cross-marker state and no per-index vector, so
//       one marker per worker is safe by construction.
//   accumulate(i, k)  the serial part, in marker order: the slot writes into
//       the per-index vectors, the running chunk index, the P1/P2 row copy,
//       the float `+=` into genoSumMat (order-dependent), the max-collapse
//       into genoURMat, the group counters, the anno x MAF indicator row and
//       the chunk flush.
// Markers are processed in blocks of m1 (= markers_per_chunk): phase 1 over
// the block -- serially at t_nThreads == 1, on a parallelFor otherwise --
// then phase 2 over the same block in order. Every expression in phase 2 is
// the pre-restructure one over the values phase 1 stored, so the output is
// byte-identical whatever the thread count.
struct RegionPrepModel : CellResult {
  arma::rowvec G1row;       // G1tilde_P_G2tilde_Weighted_Mat.row(i), already weighted
  double MACcasegroup = 0.0;
  double MACcontrolgroup = 0.0;
  double AF_case = 0.0, AF_ctrl = 0.0;
  uint32_t N_case = 0, N_ctrl = 0;
  double N_case_hom = 0.0, N_case_het = 0.0, N_ctrl_hom = 0.0, N_ctrl_het = 0.0;
};

struct RegionPrep {
  bool readFailed = false;
  bool passQC = false;
  bool isRare = false;
  bool tested = false;      // test outputs below and blockP1/P2 column k are valid
  bool deferredToMainThread = false;   // the ER arm, main-thread only: re-run there
  bool hasError = false;
  std::string errMsg;
  bool flip = false;
  // marker fields (written for every marker that was read, QC or not)
  std::string chr, pos, ref, alt, marker, info;
  double altFreq = 0.0, altCounts = 0.0, missingRate = 0.0, imputeInfo = 0.0;
  double MAF = 0.0;
  double MAC = 0.0;
  arma::uvec nz;            // nonzero sample indices after impute/flip
  arma::vec gnz;            // their genotype values (all phase 2 needs of GVec)
  arma::vec w0;             // per weight set
  std::vector<RegionPrepModel> m;   // one per model: test outputs and case/control counts
};

// Phase 1 on a parallelFor. Each range checks out a worker state (one
// scratch + one reader clone) from the marker path's pool, refuses the R
// API for its duration, and records a marker's failure IN ITS OWN prep
// record: the main thread then acts on the lowest index, which is what a
// serial run would have stopped at (the marker path's post-passes, S6).
struct RegionBlockWorker : public RcppParallel::Worker {
  RegionBlockWorker(MarkerStatePool & t_pool, unsigned int t_iStart,
                    std::vector<RegionPrep> & t_preps,
                    const std::function<void(unsigned int, unsigned int, Marker::MarkerScratch&, Geno::GenoReader&)> & t_fn)
    : m_pool(t_pool), m_iStart(t_iStart), m_preps(t_preps), m_fn(t_fn) {}
  void operator()(std::size_t t_begin, std::size_t t_end) override {
    std::unique_ptr<MarkerWorkerState> state = m_pool.checkout();
    RefuseMainThreadGuard noRAPI;
    for(std::size_t k = t_begin; k < t_end; k++){
      try {
        m_fn(m_iStart + k, k, *state->scratch[0], *state->readers[0]);
      } catch (const SKATMainThread::NeedsMainThread&) {
        m_preps[k] = RegionPrep();
        m_preps[k].deferredToMainThread = true;
      } catch (const std::exception& e) {
        m_preps[k] = RegionPrep();
        m_preps[k].hasError = true;
        m_preps[k].errMsg = e.what();
      }
    }
    m_pool.checkin(std::move(state));
    SAIGE_TSAN_HAPPENS_BEFORE(&m_pool);
  }
  MarkerStatePool & m_pool;
  unsigned int m_iStart;
  std::vector<RegionPrep> & m_preps;
  const std::function<void(unsigned int, unsigned int, Marker::MarkerScratch&, Geno::GenoReader&)> & m_fn;
};

// The ultra-rare pass and the BURDEN emitter each end in one getMarkerPval
// per (anno x max_MAF x weight) cell. The cells are independent until their
// outputs are written -- each reads its own column of genoURMat/genoSumMat
// and the const model -- so the tests run on the pass-A pool, one cell per
// range, and everything order-dependent (the residualisation and genoSumMat
// updates that feed them, the running chunk index, the output slots and the
// BURDEN carry-forward) stays serial on either side. Same protocol as the
// marker loop: a failure is recorded in the cell's own record and the main
// thread acts on the lowest index; the ER arm is deferred and re-run on the
// main thread (ComputeExact re-seeds at entry, so the order is immaterial).
struct CellOutcome {
  bool deferredToMainThread = false;
  bool hasError = false;
  std::string errMsg;
};

struct RegionCellWorker : public RcppParallel::Worker {
  RegionCellWorker(MarkerStatePool & t_pool, std::vector<CellOutcome> & t_out,
                   const std::function<void(unsigned int, Marker::MarkerScratch&)> & t_fn)
    : m_pool(t_pool), m_out(t_out), m_fn(t_fn) {}
  void operator()(std::size_t t_begin, std::size_t t_end) override {
    // The FORK edge, twin of the join edge declared on the pool below.
    // Everything the main thread wrote before parallelFor -- this worker's
    // own members among them -- happens-before every worker read of it, but
    // the edge lives inside libtbb, which is NOT instrumented, so TSAN
    // reported the constructor's writes racing this loop's reads: 1 warning
    // at t=4 and 2 at t=8 on the region workload, in every build including
    // d36e7277's, which left the clean arm's baseline non-zero and gave a
    // future real race somewhere to hide. Declaring it is not a suppression:
    // it is main-vs-worker only. The worker-vs-worker class stays fully
    // visible -- the hasVarRatio race (selectVarianceRatio's out-flag
    // shared across workers; d36e7277 racing, e72e4186 clean, and since E1
    // the flag is a local so the class cannot recur) was DETECTED with the
    // join edge already declared, and went to zero only when that bug was
    // fixed (harness/work/tsan-prefix vs harness/work/tsan).
    SAIGE_TSAN_HAPPENS_AFTER(this);
    std::unique_ptr<MarkerWorkerState> state = m_pool.checkout();
    RefuseMainThreadGuard noRAPI;
    for(std::size_t k = t_begin; k < t_end; k++){
      try {
        m_fn((unsigned int)k, *state->scratch[0]);
      } catch (const SKATMainThread::NeedsMainThread&) {
        m_out[k] = CellOutcome();
        m_out[k].deferredToMainThread = true;
      } catch (const std::exception& e) {
        m_out[k] = CellOutcome();
        m_out[k].hasError = true;
        m_out[k].errMsg = e.what();
      }
    }
    m_pool.checkin(std::move(state));
    SAIGE_TSAN_HAPPENS_BEFORE(&m_pool);
  }
  MarkerStatePool & m_pool;
  std::vector<CellOutcome> & m_out;
  const std::function<void(unsigned int, Marker::MarkerScratch&)> & m_fn;
};

//the last chunk for UR
// Step 3b (2026-09-06): the per-cell getMarkerPval calls are split out of
// this loop. The loop below (B1) is the pre-existing per-cell body minus the
// test and the writes that consume its outputs; it fills one URCell per
// tested cell. The tests then run over those cells (B2: serial at one
// thread, on the pass-A pool otherwise) and a serial replay (B3) writes the
// outputs in cell order. Every expression is the pre-split one over the
// same values, so the output is byte-identical whatever the thread count.
struct URCell {
  bool tested = false;
  unsigned int i = 0;            // output slot q0 + jmr
  unsigned int rowInChunk = 0;   // P1Mat row / P2Mat column
  // the collapsed pseudo-marker is cs.ur's column jmr (+ cs.urFill, which is
  // always 0 since LEDGER #130), scattered at test time
  std::vector<arma::uword> nz;         // its nonzero sample indices (the zero list is derived at test time)
  double altFreq = 0.0;
  bool flip = false;
  double w0 = 0.0;
};

// Step 3b (2026-09-06): the per-cell getMarkerPval is split out of the
// emitter. D1 walks the cells in order and prepares one BurdenCell for
// every cell the emitter would test (the same predicate over the same
// operands; step 4: the test context per model beside it); D2 runs
// those tests, serially at one thread and one (cell, model) per range on
// the pass-A pool otherwise; D3 is the pre-split emitter loop per model,
// reading each cell's result instead of computing it, so the
// carry-forward cache sees exactly the sequence it always did.
struct BurdenCell {
  bool tested = false;
  // the burden vector is cs.sum's column jmr, scattered at test time
  std::vector<arma::uword> nz;          // zero list derived at test time (complementIndices)
  double altFreq = 0.0;
};

} // namespace Region

#endif // SAIGE_REGION_LOOP_HPP
