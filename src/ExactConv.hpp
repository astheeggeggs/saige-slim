// LEDGER #176 (2026-10-02). Exact two-sided tail of
//   S = sum_i g_i x_i + N(rm, rs^2),   x_i ~ Bernoulli(mu_i) independent,
// by convolution of the carriers' terms on a lattice: O(n x grid), not 2^n.
// It is the fast saddlepoint kernel's own model (the carriers exact, everyone
// else a Gaussian remainder) with the tail computed instead of approximated.
// Weights are rounded to the NEAREST lattice point and what the lattice does
// not carry goes to the Gaussian part (a term that rounds to zero entirely,
// otherwise the mean of its rounding error); DOWN and UP give a bracket that
// contains the exact tail. With no Gaussian part an atom at q counts half
// (mid-p). Reference implementation and its controls:
// harness/tools/spa_exact_ref.py. Used ONLY by the stretch gate's shadow mode
// (SAIGE_test.cpp); it reports nothing into any result file.
#ifndef EXACTCONV_HPP
#define EXACTCONV_HPP
#include <vector>
#include <cmath>
#include <cstddef>
namespace ExactConv {
struct Result { double p = NAN, upper = NAN, lower = NAN, lo = NAN, hi = NAN; std::size_t gridSize = 0; double step = 0.0; };
enum Rounding { NEAREST, DOWN, UP };
inline double normUpper(double z){ return 0.5 * std::erfc(z / M_SQRT2); }
// one pass at one rounding; returns upper + lower tail (mid-p at an atom when rs == 0)
inline double pass(const double* g, const double* mu, std::size_t n, double rm, double rs,
                   double q, double qinv, double step, Rounding mode, double & upper, double & lower, std::size_t & size){
  std::vector<long> k(n);
  long neg = 0, pos = 0; double var = rs * rs, mean = rm;
  for(std::size_t i = 0; i < n; i++){
    const double r = g[i] / step;
    k[i] = (long)(mode == NEAREST ? std::llround(r) : (mode == DOWN ? std::floor(r) : std::ceil(r)));
    if(mode == NEAREST){
      // what the lattice does not carry goes to the Gaussian part: the whole term
      // when it rounds to zero, otherwise the mean of its rounding error
      const double e = g[i] - (double)k[i] * step;
      mean += e * mu[i];
      if(k[i] == 0) var += g[i] * g[i] * mu[i] * (1.0 - mu[i]);
    }
    if(k[i] < 0) neg += k[i]; else pos += k[i];
  }
  size = (std::size_t)(pos - neg + 1);
  std::vector<double> pmf(size, 0.0);
  long a = -neg, b = -neg;                 // occupied range [a, b]
  pmf[(std::size_t)a] = 1.0;
  for(std::size_t i = 0; i < n; i++){
    const long ki = k[i]; if(ki == 0) continue;
    const double m = mu[i], c = 1.0 - m;
    if(ki > 0){
      for(long j = b + ki; j >= a; j--){
        const double from = (j - ki >= a && j - ki <= b) ? pmf[(std::size_t)(j - ki)] : 0.0;
        const double stay = (j <= b) ? pmf[(std::size_t)j] : 0.0;
        pmf[(std::size_t)j] = stay * c + from * m;
      }
      b += ki;
    }else{
      for(long j = a + ki; j <= b; j++){
        const double from = (j - ki >= a && j - ki <= b) ? pmf[(std::size_t)(j - ki)] : 0.0;
        const double stay = (j >= a) ? pmf[(std::size_t)j] : 0.0;
        pmf[(std::size_t)j] = stay * c + from * m;
      }
      a += ki;
    }
  }
  const double sd = std::sqrt(var);
  upper = 0.0; lower = 0.0;
  const double tol = 0.5 * step;
  for(long j = a; j <= b; j++){
    const double pj = pmf[(std::size_t)j]; if(pj <= 0.0) continue;
    const double v = (double)(j + neg) * step + mean;
    if(sd > 0.0){
      upper += pj * normUpper((q - v) / sd);
      lower += pj * normUpper((v - qinv) / sd);
    }else{
      if(v > q + tol) upper += pj; else if(v >= q - tol) upper += 0.5 * pj;
      if(v < qinv - tol) lower += pj; else if(v <= qinv + tol) lower += 0.5 * pj;
    }
  }
  return upper + lower;
}
// q above the mean and qinv its mirror, or the reverse; the two tails are summed.
inline Result twoSided(const double* g, const double* mu, std::size_t n, double rm, double rs,
                       double q, double qinv, std::size_t grid, bool bracket){
  Result r; double span = 0.0;
  for(std::size_t i = 0; i < n; i++) span += std::fabs(g[i]);
  if(!(span > 0.0)) return r;
  const double hiq = q > qinv ? q : qinv, loq = q > qinv ? qinv : q;
  r.step = span / (double)grid;
  r.p = pass(g, mu, n, rm, rs, hiq, loq, r.step, NEAREST, r.upper, r.lower, r.gridSize);
  if(bracket){
    double u, l; std::size_t s;
    const double d = pass(g, mu, n, rm, rs, hiq, loq, r.step, DOWN, u, l, s);
    const double up = pass(g, mu, n, rm, rs, hiq, loq, r.step, UP, u, l, s);
    r.lo = d < up ? d : up; r.hi = d < up ? up : d;
  }
  if(r.p > 1.0) r.p = 1.0;
  return r;
}
}
#endif
