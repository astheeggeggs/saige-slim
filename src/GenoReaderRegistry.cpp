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

  bool getOneMarker(uint64_t& t_gIndex_prev, uint64_t& t_gIndex,
                    std::string& t_ref, std::string& t_alt, std::string& t_marker,
                    uint32_t& t_pd, std::string& t_chr,
                    double& t_altFreq, double& t_altCounts, double& t_missingRate,
                    double& t_imputeInfo,
                    bool t_isOutputIndexForMissing, std::vector<uint>& t_indexForMissing,
                    bool t_isOnlyOutputNonZero, std::vector<uint>& t_indexForNonZero,
                    arma::vec& t_GVec) override {
    bool isTrueGenotype = true;  // PLINK-private; always true on the step-2 path
    m_impl.getOneMarker(t_gIndex_prev, t_gIndex, t_ref, t_alt, t_marker, t_pd, t_chr,
                        t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
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
