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
// Cursor/seek behavior is reader-private. gIndex addresses the 0-based row of
// the variant metadata file (bim/pvar); random access must be supported —
// every implementation must return the requested record regardless of prior
// read position (LEDGER.md #16 was exactly this contract being violated).

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <memory>
#include <functional>
#include <string>
#include <vector>

namespace Geno {

// Per-marker metadata filled by getOneMarker.
struct MarkerMeta {
  std::string chr, id, ref, alt;
  uint32_t pos = 0;
  double altFreq = 0, altCounts = 0, missingRate = 0, imputeInfo = 1;
};

class GenoReader {
public:
  virtual ~GenoReader() = default;

  virtual std::string formatName() const = 0;   // must equal the registered genoType string
  virtual uint32_t getN0() const = 0;           // samples in the genotype file
  virtual uint32_t getN() const = 0;            // samples in the analysis

  // gIndex addresses the 0-based row of the variant metadata file (bim).
  // Returns false when no marker could be read (end of stream).
  virtual bool getOneMarker(uint64_t t_gIndex,
                            MarkerMeta& t_meta,
                            bool t_isOutputIndexForMissing,
                            std::vector<uint>& t_indexForMissing,
                            bool t_isOnlyOutputNonZero,
                            std::vector<uint>& t_indexForNonZero,
                            arma::vec& t_GVec) = 0;

  virtual void close() = 0;

  // True when the reader carries no mutable per-read stream state that would
  // be shared across a fork(). A reader that seeks a FILE* is NOT fork-safe:
  // parent and children share one file description and therefore one offset,
  // so concurrent readers silently decode each other's records.
  virtual bool isForkSafe() const = 0;

  // A clone is an independent handle onto the SAME open genotype data, so a
  // worker can read without coordinating with any other. Implementations own
  // the underlying data jointly: whatever the last surviving handle releases is
  // what gets released, so close() on one handle must never invalidate another.
  // Implemented for PLINK on the mmap path (S5); PGEN still throws.
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
