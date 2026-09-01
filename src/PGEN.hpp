
// All C++ code related to PGEN (plink2 native format) file manipulation.
// Ported from main's pre-GenoReader-refactor PgenClass (see LEDGER.md /
// CLAUDE.md Phase 7) -- mechanics unchanged, relocated to sit alongside the
// GenoReader interface. Not a drop-in: GenoReaderRegistry.cpp wraps this by
// composition (PgenReader : GenoReader), analogous to PlinkReader/PlinkClass.

#ifndef PGEN_HPP
#define PGEN_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <unordered_map>
#include <memory>
#include <pgenlib_read.h>
#include <pvar_ffi_support.h>
#include <pgenlib_ffi_support.h>

namespace PGEN {

typedef unsigned int uint;

class PgenClass {
private:
    // information from pvar file
    std::vector<std::string> chr;              // Chromosome code
    std::vector<std::string> variantId;        // Variant identifier
    std::vector<uint32_t> position;            // Base-pair coordinate (1-based)
    std::vector<std::string> ref;              // REF allele
    std::vector<std::string> alt;              // ALT allele
    size_t M;

    // information from psam file
    std::unordered_map<std::string, size_t> sampleInPgen;
    uint32_t N, N0;

    // PGEN file
    std::unique_ptr<plink2::PgenFileInfo> _info_ptr;
    std::unique_ptr<plink2::PgenReader> _state_ptr;
    std::string pgenFile, psamFile, pvarFile;
    std::vector<size_t> posSampleInPgen;

    plink2::RefcountedWptr* _allele_idx_offsetsp;
    plink2::RefcountedWptr* _nonref_flagsp;
    uintptr_t* _subset_include_vec;
    uintptr_t* _subset_include_interleaved_vec;
    uint32_t* _subset_cumulative_popcounts;
    plink2::PgrSampleSubsetIndex _subset_index;
    uint32_t _subset_size;
    plink2::PgenVariant _pgv;
    uintptr_t* _raregeno_buf;
    uint32_t* _difflist_sample_ids_buf;
    plink2::VecW* _transpose_batch_buf;
    uintptr_t* _multivar_vmaj_geno_buf;
    uintptr_t* _multivar_vmaj_phasepresent_buf;
    uintptr_t* _multivar_vmaj_phaseinfo_buf;
    uintptr_t* _multivar_smaj_geno_batch_buf;
    uintptr_t* _multivar_smaj_phaseinfo_batch_buf;
    uintptr_t* _multivar_smaj_phasepresent_batch_buf;
    std::vector<double> readBuffer;

    // internal functions
    void readPvarFile();
    void readPsamFile();
    void loadPgen(uint32_t raw_sample_ct);
public:

    PgenClass(std::string pgenFile,
        std::string psamFile,
        std::string pvarFile,
        std::vector<std::string>& sampleInModel);

    void setPosSampleInPgen(std::vector<std::string>& SampleInModel);
    void Read(std::vector<double>& buf, int variant_idx, int allele_idx);

    void getOneMarker(
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
        bool& t_isOutputIndexForMissing,
        std::vector<uint>& t_indexForMissing,
        bool& t_isOnlyOutputNonZero,
        std::vector<uint>& t_indexForNonZero,
        arma::vec& OneMarkerG1);

    uint32_t getN0() { return N0; }
    uint32_t getN() { return N; }
    uint32_t getM() { return M; }

    void closegenofile();

};

}

#endif
