#ifndef GENOREADER_HPP
#define GENOREADER_HPP

// Genotype-reader interface (Phase 3 refactor).
//
// Two implementations exist: PlinkReader (wrapping PLINK::PlinkClass) and
// PgenReader (wrapping PGEN::PgenClass, added in Phase 7). The registry below
// is the open registration point for further formats. The R layer addresses
// formats by the genoType string ("plink", "pgen", ...).
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

  // Item 3: the marker's metadata and QC-gate quantities (altFreq, altCounts,
  // missingRate, imputeInfo, and the missing count) from the genotype codes
  // alone, without writing t_GVec. Optional: a reader that cannot do this
  // cheaper than a full read returns false and the caller reads in full. A
  // reader that returns true MUST produce the same values getOneMarker would,
  // bit for bit -- the caller's QC decision is made from them.
  virtual bool countOneMarker(uint64_t /*t_gIndex*/, MarkerMeta& /*t_meta*/,
                              uint32_t& /*t_numMissing*/) { return false; }

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
  // How that joint ownership is arranged is the implementation's business, and
  // the two built-ins do it differently: PLINK clones SHARE one decoder because
  // its reads are const pointer arithmetic into an mmap (S5), while PGEN clones
  // each build their own pgenlib reader over shared metadata, because a pgenlib
  // read mutates the decoder (S8). Refusing is still allowed -- PLINK's pread
  // fallback does -- and the marker scan degrades to one thread when it happens.
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
