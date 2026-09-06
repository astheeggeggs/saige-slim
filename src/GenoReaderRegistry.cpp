// Reader registry + the PLINK/PGEN implementations (wrapping PLINK::PlinkClass
// and PGEN::PgenClass by composition — the decode logic itself is untouched;
// see GenoReader.hpp).

#include "GenoReader.hpp"
#include "PLINK.hpp"
#include "PGEN.hpp"
#include <map>
#include <stdexcept>

namespace Geno {

namespace {

std::map<std::string, ReaderFactory>& registry() {
  static std::map<std::string, ReaderFactory> r;
  return r;
}

// S5: clones SHARE one PlinkClass rather than each opening the .bed again.
// That is sound because PlinkClass::getOneMarker is const -- it writes no
// member on either path -- so concurrent readers cannot interfere, and it is
// cheaper than re-reading the bim/fam metadata per worker. Ownership is shared,
// so the map and the FILE* die with the LAST owner: closing the primary while a
// clone is still reading is safe by construction, not by convention.
class PlinkReader : public GenoReader {
public:
  PlinkReader(const Rcpp::List& spec, const std::vector<std::string>& sampleInModel)
      : m_impl(std::make_shared<PLINK::PlinkClass>(
                   Rcpp::as<std::string>(spec["bimFile"]),
                   Rcpp::as<std::string>(spec["famFile"]),
                   Rcpp::as<std::string>(spec["bedFile"]),
                   Rcpp::as<std::string>(spec["AlleleOrder"]))) {
    std::vector<std::string> ids = sampleInModel;
    m_impl->setPosSampleInPlink(ids);
  }

  std::string formatName() const override { return "plink"; }
  uint32_t getN0() const override { return impl().getN0(); }
  uint32_t getN() const override { return impl().getN(); }

  bool getOneMarker(uint64_t t_gIndex, MarkerMeta& t_meta,
                    bool t_isOutputIndexForMissing, std::vector<uint>& t_indexForMissing,
                    bool t_isOnlyOutputNonZero, std::vector<uint>& t_indexForNonZero,
                    arma::vec& t_GVec) override {
    bool isTrueGenotype = true;  // PLINK-private; always true on the step-2 path
    impl().getOneMarker(t_gIndex, t_meta.ref, t_meta.alt, t_meta.id,
                        t_meta.pos, t_meta.chr,
                        t_meta.altFreq, t_meta.altCounts, t_meta.missingRate,
                        t_meta.imputeInfo,
                        t_isOutputIndexForMissing, t_indexForMissing,
                        t_isOnlyOutputNonZero, t_indexForNonZero,
                        isTrueGenotype, t_GVec);
    return true;
  }

  bool countOneMarker(uint64_t t_gIndex, MarkerMeta& t_meta, uint32_t& t_numMissing) override {
    impl().countOneMarker(t_gIndex, t_meta.ref, t_meta.alt, t_meta.id, t_meta.pos, t_meta.chr,
                          t_meta.altFreq, t_meta.altCounts, t_meta.missingRate,
                          t_meta.imputeInfo, t_numMissing);
    return true;
  }

  // Releases THIS owner's reference. The .bed is unmapped and closed when the
  // last owner goes away, so a clone that is still reading keeps working.
  void close() override { m_impl.reset(); }

  // Fork-safe only on the mmap path: reads are pointer arithmetic into a
  // MAP_PRIVATE mapping, with no shared offset. The pread fallback carries no
  // shared offset either, but it has never been exercised at biobank scale, so
  // the conservative answer is kept.
  bool isForkSafe() const override { return impl().isMemMapped(); }

  std::unique_ptr<GenoReader> clone() const override {
    if(!impl().isMemMapped()){
      // The fallback decodes correctly (pread is offset-free) but has no
      // performance story and no scale testing; refuse rather than hand out a
      // reader whose behaviour under threads nothing has measured.
      throw std::logic_error(
          "PlinkReader::clone() requires the mmap path; this .bed fell back to pread");
    }
    return std::unique_ptr<GenoReader>(new PlinkReader(m_impl));
  }

private:
  explicit PlinkReader(std::shared_ptr<PLINK::PlinkClass> t_impl)
      : m_impl(std::move(t_impl)) {}

  PLINK::PlinkClass& impl() const {
    if(!m_impl)
      throw std::logic_error("PlinkReader used after close()");
    return *m_impl;
  }

  std::shared_ptr<PLINK::PlinkClass> m_impl;
};

// PGEN (Phase 7; cloneable since S8). Clones do NOT share one PgenClass the
// way PlinkReader's do -- a pgenlib read mutates the reader's own decode state
// (fread buffer, LD cache, workspaces), so sharing one handle would be a data
// race. Each clone gets its own plink2::PgenReader onto the SAME shared
// PgenShared (metadata + PgenFileInfo), which is pgenlib's documented mode-2
// arrangement; see the header comment in src/PGEN.hpp.
class PgenReader : public GenoReader {
public:
  PgenReader(const Rcpp::List& spec, const std::vector<std::string>& sampleInModel)
      : m_sampleInModel(sampleInModel),
        m_impl(std::make_shared<PGEN::PgenClass>(
                   Rcpp::as<std::string>(spec["pgenFile"]),
                   Rcpp::as<std::string>(spec["psamFile"]),
                   Rcpp::as<std::string>(spec["pvarFile"]),
                   m_sampleInModel)) {}

  std::string formatName() const override { return "pgen"; }
  uint32_t getN0() const override { return impl().getN0(); }
  uint32_t getN() const override { return impl().getN(); }

  bool getOneMarker(uint64_t t_gIndex, MarkerMeta& t_meta,
                    bool t_isOutputIndexForMissing, std::vector<uint>& t_indexForMissing,
                    bool t_isOnlyOutputNonZero, std::vector<uint>& t_indexForNonZero,
                    arma::vec& t_GVec) override {
    impl().getOneMarker(t_gIndex, t_meta.ref, t_meta.alt, t_meta.id,
                        t_meta.pos, t_meta.chr,
                        t_meta.altFreq, t_meta.altCounts, t_meta.missingRate,
                        t_meta.imputeInfo,
                        t_isOutputIndexForMissing, t_indexForMissing,
                        t_isOnlyOutputNonZero, t_indexForNonZero,
                        t_GVec);
    return true;
  }

  // Releases THIS handle's decoder and its reference to the shared file info.
  // A clone keeps reading: it holds its own decoder and its own reference, and
  // the shared state dies with the last of them.
  void close() override { m_impl.reset(); }

  // Not fork-safe, and this one is a live hazard rather than conservatism: a
  // pgenlib reader owns an open FILE* that it fseeks, so parent and child would
  // share one file description and one offset. (The forked scan is gone -- S7 --
  // so nothing in the pipeline asks.)
  bool isForkSafe() const override { return false; }

  std::unique_ptr<GenoReader> clone() const override {
    return std::unique_ptr<GenoReader>(new PgenReader(
        std::make_shared<PGEN::PgenClass>(impl().shared())));
  }

private:
  explicit PgenReader(std::shared_ptr<PGEN::PgenClass> t_impl)
      : m_impl(std::move(t_impl)) {}

  PGEN::PgenClass& impl() const {
    if(!m_impl)
      throw std::logic_error("PgenReader used after close()");
    return *m_impl;
  }

  std::vector<std::string> m_sampleInModel;  // PgenClass stores a reference during construction only
  std::shared_ptr<PGEN::PgenClass> m_impl;
};

void registerBuiltinReaders() {
  static bool done = false;
  if (done) return;
  done = true;
  registerReader("plink", [](const Rcpp::List& spec,
                             const std::vector<std::string>& sampleInModel) {
    return std::unique_ptr<GenoReader>(new PlinkReader(spec, sampleInModel));
  });
  registerReader("pgen", [](const Rcpp::List& spec,
                            const std::vector<std::string>& sampleInModel) {
    return std::unique_ptr<GenoReader>(new PgenReader(spec, sampleInModel));
  });
}

} // namespace

void registerReader(const std::string& genoType, ReaderFactory factory) {
  registry()[genoType] = std::move(factory);
}

std::unique_ptr<GenoReader> makeReader(const std::string& genoType,
                                       const Rcpp::List& spec,
                                       const std::vector<std::string>& sampleInModel) {
  registerBuiltinReaders();
  auto it = registry().find(genoType);
  if (it == registry().end()) {
    Rcpp::stop("Unknown genoType '%s': no genotype reader registered for it", genoType);
  }
  auto reader = it->second(spec, sampleInModel);
  if (reader->formatName() != genoType) {
    Rcpp::stop("Reader registered for '%s' reports formatName '%s'",
               genoType, reader->formatName());
  }
  return reader;
}

} // namespace Geno
