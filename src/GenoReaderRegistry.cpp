// Reader registry + the PLINK implementation (wrapping PLINK::PlinkClass by
// composition — the decode logic itself is untouched; see GenoReader.hpp).

#include "GenoReader.hpp"
#include "PLINK.hpp"
#include <map>
#include <stdexcept>

namespace Geno {

namespace {

std::map<std::string, ReaderFactory>& registry() {
  static std::map<std::string, ReaderFactory> r;
  return r;
}

class PlinkReader : public GenoReader {
public:
  PlinkReader(const Rcpp::List& spec, const std::vector<std::string>& sampleInModel)
      : m_impl(Rcpp::as<std::string>(spec["bimFile"]),
               Rcpp::as<std::string>(spec["famFile"]),
               Rcpp::as<std::string>(spec["bedFile"]),
               Rcpp::as<std::string>(spec["AlleleOrder"])) {
    std::vector<std::string> ids = sampleInModel;
    m_impl.setPosSampleInPlink(ids);
  }

  std::string formatName() const override { return "plink"; }
  uint32_t getN0() const override { return const_cast<PLINK::PlinkClass&>(m_impl).getN0(); }
  uint32_t getN() const override { return const_cast<PLINK::PlinkClass&>(m_impl).getN(); }

  bool getOneMarker(uint64_t t_gIndex, MarkerMeta& t_meta,
                    bool t_isOutputIndexForMissing, std::vector<uint>& t_indexForMissing,
                    bool t_isOnlyOutputNonZero, std::vector<uint>& t_indexForNonZero,
                    arma::vec& t_GVec) override {
    // prev=0 makes PlinkClass take the absolute-seek branch for every
    // t_gIndex > 0 — position-identical to the historical relative-skip mix —
    // and preserves the LEDGER.md #16 no-seek-at-row-0 quirk byte-for-byte.
    // The 3e fix replaces this with an unconditional absolute seek.
    uint64_t prev = 0;
    bool isTrueGenotype = true;  // PLINK-private; always true on the step-2 path
    m_impl.getOneMarker(prev, t_gIndex, t_meta.ref, t_meta.alt, t_meta.id,
                        t_meta.pos, t_meta.chr,
                        t_meta.altFreq, t_meta.altCounts, t_meta.missingRate,
                        t_meta.imputeInfo,
                        t_isOutputIndexForMissing, t_indexForMissing,
                        t_isOnlyOutputNonZero, t_indexForNonZero,
                        isTrueGenotype, t_GVec);
    return true;
  }

  void close() override { m_impl.closegenofile(); }

  std::unique_ptr<GenoReader> clone() const override {
    // Deferred to the Phase-4 shared-metadata implementation.
    throw std::logic_error("PlinkReader::clone() not implemented yet");
  }

private:
  PLINK::PlinkClass m_impl;
};

void registerBuiltinReaders() {
  static bool done = false;
  if (done) return;
  done = true;
  registerReader("plink", [](const Rcpp::List& spec,
                             const std::vector<std::string>& sampleInModel) {
    return std::unique_ptr<GenoReader>(new PlinkReader(spec, sampleInModel));
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
