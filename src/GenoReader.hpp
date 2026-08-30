#ifndef GENOREADER_HPP
#define GENOREADER_HPP

// Genotype-reader interface (Phase 3 refactor).
//
// One implementation exists today (PlinkReader, wrapping PLINK::PlinkClass).
// The registry below is the open registration point for future formats —
// PGEN/plink2 input returns as a new implementation registered here (see
// CLAUDE.md and the Phase-3 design notes). The R layer keeps addressing
// formats by the genoType string ("plink", later "pgen", ...).
//
// NOTE (transitional, until the cursor-internalization commit): getOneMarker
// still carries t_gIndex_prev, the PLINK sequential-read optimization leaked
// into the old call signature. It is slated to become reader-private cursor
// state; do not build new callers against it.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <memory>
#include <functional>
#include <string>
#include <vector>

namespace Geno {

class GenoReader {
public:
  virtual ~GenoReader() = default;

  virtual std::string formatName() const = 0;   // must equal the registered genoType string
  virtual uint32_t getN0() const = 0;           // samples in the genotype file
  virtual uint32_t getN() const = 0;            // samples in the analysis

  // gIndex addresses the 0-based row of the variant metadata file (bim).
  // Returns false when no marker could be read (end of stream).
  virtual bool getOneMarker(uint64_t& t_gIndex_prev,
                            uint64_t& t_gIndex,
                            std::string& t_ref,
                            std::string& t_alt,
                            std::string& t_marker,
                            uint32_t& t_pd,
                            std::string& t_chr,
                            double& t_altFreq,
                            double& t_altCounts,
                            double& t_missingRate,
                            double& t_imputeInfo,
                            bool t_isOutputIndexForMissing,
                            std::vector<uint>& t_indexForMissing,
                            bool t_isOnlyOutputNonZero,
                            std::vector<uint>& t_indexForNonZero,
                            arma::vec& t_GVec) = 0;

  virtual void close() = 0;

  // Phase-4 hook: a clone shares immutable metadata but owns private
  // stream/buffer state, enabling per-thread readers.
  virtual std::unique_ptr<GenoReader> clone() const = 0;
};

// Factory contract: `spec` is a named list of format-specific inputs
// (for "plink": bimFile, famFile, bedFile, AlleleOrder).
using ReaderFactory = std::function<std::unique_ptr<GenoReader>(
    const Rcpp::List& spec,
    const std::vector<std::string>& sampleInModel)>;

// The open registration point. Registration is explicit (no static-init
// order games): makeReader() lazily registers the built-ins on first use.
void registerReader(const std::string& genoType, ReaderFactory factory);
std::unique_ptr<GenoReader> makeReader(const std::string& genoType,
                                       const Rcpp::List& spec,
                                       const std::vector<std::string>& sampleInModel);

} // namespace Geno

#endif
