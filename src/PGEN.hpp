// All C++ code related to PGEN (plink2 native format) file manipulation.
// Ported from main's pre-GenoReader-refactor PgenClass (see LEDGER.md /
// CLAUDE.md Phase 7) -- mechanics unchanged, relocated to sit alongside the
// GenoReader interface. Not a drop-in: GenoReaderRegistry.cpp wraps this by
// composition (PgenReader : GenoReader), analogous to PlinkReader/PlinkClass.
//
// S8 SPLIT THE CLASS IN TWO, so PGEN can be read from several threads at once:
//
//   PgenShared  everything that is read-only once construction finishes -- the
//               pvar/psam metadata, the analysis sample subset, and the
//               plink2::PgenFileInfo with its allocations. One instance per
//               opened .pgen, shared by every reader onto it.
//   PgenClass   one DECODE handle: its own plink2::PgenReader (its own FILE*,
//               its own fread buffer, its own LD cache and workspaces) plus the
//               per-read output buffers. Not shared, ever.
//
// That split is pgenlib's own contract, not an invention here: PgrInit()'s
// "mode 2" (per-variant fread, fname non-null) documents that many PgenReaders
// may be initialised against ONE PgenFileInfo, the first one stealing the
// header's FILE* and each later one opening its own. PgrInit copies the
// PgenFileInfo by value into the reader, and no read path writes through the
// pointers that copy shares (var_fpos/vrtypes/allele_idx_offsets/nonref_flags
// are written only by PgfiInitPhase2), so concurrent readers cannot interfere.
//
// This is deliberately NOT the shape PlinkReader uses. There, clones share one
// PlinkClass because getOneMarker is const and the .bed is an mmap; here a
// read mutates the reader's own decode state, so sharing one handle would be a
// data race and each clone gets its own.

#ifndef PGEN_HPP
#define PGEN_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <unordered_map>
#include <memory>
#include <mutex>
#include <string>
#include <pgenlib_read.h>
#include <pvar_ffi_support.h>
#include <pgenlib_ffi_support.h>

namespace PGEN {

typedef unsigned int uint;

// Immutable after construction; shared by a PgenClass and all its clones.
// Everything a read needs but never writes lives here.
class PgenShared {
public:
    PgenShared(const std::string& t_pgenFile,
               const std::string& t_psamFile,
               const std::string& t_pvarFile,
               const std::vector<std::string>& t_sampleInModel);

    // ALL file samples, in file order: posSampleInPgen is the identity and
    // N == N0. Steps 0/1 want this -- genoClass does its own subsetting, in
    // genotype-file row order, via ptrsubSampleInGeno and
    // indicatorGenoSamplesWithPheno_in (both built in R from the psam), so it
    // needs every psam row and never consults posSampleInPgen. Using this
    // overload also keeps steps 0/1 clear of readPsamFile's sample-ID map,
    // which keys on psam COLUMN 0 -- FID, not IID, whenever the psam carries a
    // "#FID IID ..." header (LEDGER #52: masked today only because both
    // committed fixtures have FID == IID).
    PgenShared(const std::string& t_pgenFile,
               const std::string& t_psamFile,
               const std::string& t_pvarFile);

    ~PgenShared();

    PgenShared(const PgenShared&) = delete;
    PgenShared& operator=(const PgenShared&) = delete;

    // information from the pvar file, one entry per variant
    std::vector<std::string> chr;              // Chromosome code
    std::vector<std::string> variantId;        // Variant identifier
    std::vector<uint32_t> position;            // Base-pair coordinate (1-based)
    std::vector<std::string> ref;              // REF allele
    std::vector<std::string> alt;              // ALT allele
    size_t M;

    // information from the psam file
    std::unordered_map<std::string, size_t> sampleInPgen;
    std::vector<size_t> posSampleInPgen;       // model sample -> row in the pgen
    // 2026-09-16: the inverse, row in the pgen -> model sample (-1 when the
    // row is not in the model), so a difference list's raw sample ids map to
    // analysis indices in O(1) each (getOneMarkerSparse).
    std::vector<int32_t> analysisIndexOfFileSample;
    uint32_t N, N0;

    std::string pgenFile, psamFile, pvarFile;

    // What PgrInit() needs from PgfiInitPhase2, kept so each reader can be
    // built without re-reading the header.
    plink2::PgenFileInfo* info() const { return _info_ptr.get(); }
    uint32_t maxVrecWidth() const { return _max_vrec_width; }
    uintptr_t pgrAllocCachelineCt() const { return _pgr_alloc_cacheline_ct; }
    uint32_t fileSampleCt() const { return _file_sample_ct; }
    bool isNotPlink1Bed() const { return _is_not_plink1_bed; }
    bool hasDosage() const { return _has_dosage; }
    bool isMultiallelic() const { return _is_multiallelic; }

    // PgrInit() moves shared_ff out of the PgenFileInfo for the FIRST reader,
    // so two readers must not be initialised concurrently. Reader construction
    // takes this; nothing on the read path does.
    std::mutex initMutex;

private:
    void readPvarFile();
    void readPsamFile();
    size_t countPsamSamples() const;
    void setPosSampleInPgen(const std::vector<std::string>& t_sampleInModel);
    void loadPgfi();

    std::unique_ptr<plink2::PgenFileInfo> _info_ptr;
    unsigned char* _pgfi_alloc;                // the block var_fpos/vrtypes live in
    plink2::RefcountedWptr* _allele_idx_offsetsp;
    plink2::RefcountedWptr* _nonref_flagsp;
    uint32_t _max_vrec_width;
    uintptr_t _pgr_alloc_cacheline_ct;
    uint32_t _file_sample_ct;
    bool _is_not_plink1_bed;
    bool _has_dosage = false;        // a dosage track is present (header gflags)
    bool _is_multiallelic = false;   // max_allele_ct != 2
};

// One decode handle onto a PgenShared. Copies are not made; additional handles
// come from cloneHandle(), which builds a second reader against the same
// shared state.
class PgenClass {
private:
    // Declared FIRST so it is destroyed LAST: this reader's pgenlib state is
    // released before the shared file info it was initialised against.
    std::shared_ptr<PgenShared> m_shared;

    // PGEN reader state -- per handle, mutated by every read.
    std::unique_ptr<plink2::PgenReader> _state_ptr;
    unsigned char* _pgr_alloc;
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

    // 2026-09-16, the carrier-space read: the longest difference list the
    // reader accepts (2 * (N/8), the raregeno/sample-id buffers' capacity and
    // plink2 --glm's own cap; 0 = never), and the SAIGE_PGEN_SPARSE=0 kill
    // switch, read once at construction.
    uint32_t _max_difflist_len = 0;
    bool _sparseEnabled = true;

    void initReader();
    void fillMeta(uint64_t t_gIndex, std::string& t_ref, std::string& t_alt,
                  std::string& t_marker, uint32_t& t_pd, std::string& t_chr) const;
    void expandFromReadBuffer(double& t_altFreq, double& t_altCounts,
                              double& t_missingRate, double& t_imputeInfo,
                              bool t_isOutputIndexForMissing, std::vector<uint>& t_indexForMissing,
                              bool t_isOnlyOutputNonZero, std::vector<uint>& t_indexForNonZero,
                              arma::vec& OneMarkerG1);
    bool sparseSupported() const;
public:

    PgenClass(std::string pgenFile,
        std::string psamFile,
        std::string pvarFile,
        std::vector<std::string>& sampleInModel);

    // All file samples, in file order (see the PgenShared overload above).
    PgenClass(std::string pgenFile,
        std::string psamFile,
        std::string pvarFile);

    // A second handle onto the same open pgen: same metadata, its own
    // plink2::PgenReader. Safe to call while other handles are reading.
    explicit PgenClass(std::shared_ptr<PgenShared> t_shared);
    ~PgenClass();

    PgenClass(const PgenClass&) = delete;
    PgenClass& operator=(const PgenClass&) = delete;

    std::shared_ptr<PgenShared> shared() const { return m_shared; }

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

    // The carrier-space read (2026-09-16; GenoReader.hpp SparseMarker states
    // the contract). Returns 0 when this file has no difference lists to
    // offer (plink1-compatible mode, a dosage track, multiallelic sites, the
    // kill switch), 1 when the record was stored densely -- or with a common
    // code other than 0/2 -- and has been expanded into the getOneMarker
    // outputs (t_GVec, t_indexForMissing, the counts) by the SAME loop
    // getOneMarker runs, and 2 when t_commonCode / t_idx / t_code hold the
    // record and the counts describe the analysis samples.
    int getOneMarkerSparse(
        uint64_t t_gIndex,
        std::string& t_ref,
        std::string& t_alt,
        std::string& t_marker,
        uint32_t& t_pd,
        std::string& t_chr,
        double& t_altFreq,
        double& t_altCounts,
        double& t_missingRate,
        double& t_imputeInfo,
        uint32_t& t_commonCode,
        std::vector<uint32_t>& t_idx,
        std::vector<uint8_t>& t_code,
        std::vector<uint>& t_indexForMissing,
        arma::vec& t_GVec);

    uint32_t getN0() const { return m_shared->N0; }
    uint32_t getN() const { return m_shared->N; }
    uint32_t getM() const { return m_shared->M; }

    void closegenofile();

};

}

#endif
