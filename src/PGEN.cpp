#include <boost/algorithm/string.hpp>
#include "PGEN.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <numeric>

namespace PGEN {

static const double kGenoRDoublePairs[32] ALIGNV16 = PAIR_TABLE16(0.0, 1.0, 2.0, NA_REAL);

// ---- PgenShared: metadata + PgenFileInfo, read-only once built -------------

PgenShared::PgenShared(const std::string& t_pgenFile,
    const std::string& t_psamFile,
    const std::string& t_pvarFile,
    const std::vector<std::string>& t_sampleInModel)
    : M(0), N(0), N0(0),
      pgenFile(t_pgenFile), psamFile(t_psamFile), pvarFile(t_pvarFile),
      _pgfi_alloc(nullptr), _allele_idx_offsetsp(nullptr), _nonref_flagsp(nullptr),
      _max_vrec_width(0), _pgr_alloc_cacheline_ct(0), _file_sample_ct(0),
      _is_not_plink1_bed(false)
{
    readPvarFile();
    readPsamFile();
    setPosSampleInPgen(t_sampleInModel);

    N0 = sampleInPgen.size();
    N = posSampleInPgen.size();
    M = variantId.size();

    loadPgfi();

    // 2026-09-16: row in the pgen -> model sample, for the carrier-space read.
    analysisIndexOfFileSample.assign(_file_sample_ct, -1);
    for(size_t i = 0; i < posSampleInPgen.size(); i++) {
        const size_t row = posSampleInPgen[i];
        if(row >= _file_sample_ct) {
            Rcpp::stop("psam row " + std::to_string(row) + " is beyond the pgen's sample count (" +
                       std::to_string(_file_sample_ct) + ").");
        }
        analysisIndexOfFileSample[row] = static_cast<int32_t>(i);
    }
}

// All file samples, in file order. Deliberately does NOT call
// setPosSampleInPgen: no sample-ID matching happens, so this path cannot be
// affected by readPsamFile's FID/IID column-0 keying.
PgenShared::PgenShared(const std::string& t_pgenFile,
    const std::string& t_psamFile,
    const std::string& t_pvarFile)
    : M(0), N(0), N0(0),
      pgenFile(t_pgenFile), psamFile(t_psamFile), pvarFile(t_pvarFile),
      _pgfi_alloc(nullptr), _allele_idx_offsetsp(nullptr), _nonref_flagsp(nullptr),
      _max_vrec_width(0), _pgr_alloc_cacheline_ct(0), _file_sample_ct(0),
      _is_not_plink1_bed(false)
{
    readPvarFile();

    // Counted, NOT taken from sampleInPgen.size(): that map keys on psam
    // column 0, so any two samples sharing an FID would collapse into one
    // entry and undercount N0. Nothing on this path needs the map, so it is
    // never built.
    N0 = countPsamSamples();
    posSampleInPgen.reserve(N0);
    for(size_t i = 0; i < N0; i++) {
        posSampleInPgen.push_back(i);
    }
    N = N0;
    M = variantId.size();

    loadPgfi();

    analysisIndexOfFileSample.resize(_file_sample_ct);
    std::iota(analysisIndexOfFileSample.begin(), analysisIndexOfFileSample.end(), 0);
}

// Data rows in the psam: every line that is neither a '#' header nor blank.
size_t PgenShared::countPsamSamples() const
{
    std::ifstream psam(psamFile);
    if(!psam.is_open()) {
        Rcpp::stop("Cannot open psam file (" + psamFile + ").");
    }
    std::string line;
    size_t n = 0;
    while(getline(psam, line)) {
        if(line.empty() || line[0] == '#') continue;
        if(line.size() == 1 && line[0] == '\r') continue;
        n++;
    }
    return n;
}

void PgenShared::loadPgfi()
{
    // similar to pgenlibr/src/pgenlibr.cpp; the second half (PgrInit and the
    // buffers carved out of pgr_alloc) is per-reader and lives in
    // PgenClass::initReader.
    const uint32_t raw_sample_ct = N0;
    plink2::PgenHeaderCtrl header_ctrl;
    _info_ptr = std::unique_ptr<plink2::PgenFileInfo>(new plink2::PgenFileInfo);
    uintptr_t pgfi_alloc_cacheline_ct;
    char errstr_buf[plink2::kPglErrstrBufBlen];

    plink2::PreinitPgfi(_info_ptr.get());
    uint32_t cur_sample_ct = raw_sample_ct;
    uint32_t cur_variant_ct = UINT32_MAX;
    if(plink2::PgfiInitPhase1(pgenFile.c_str(), nullptr, cur_variant_ct, cur_sample_ct, &header_ctrl, _info_ptr.get(), &pgfi_alloc_cacheline_ct, errstr_buf) != plink2::kPglRetSuccess) {
        Rcpp::stop(&(errstr_buf[7]));
    }
    const uint32_t raw_variant_ct = _info_ptr->raw_variant_ct;
    // LEDGER #87 FIXED. The pgen header carries its own variant count and it was
    // never compared with the pvar's row count, so a .pvar missing an interior
    // row relabelled every marker after it at exit 0 (measured: rs5 deleted, 68
    // of 70 rows carrying the previous marker's statistics). Too LONG a .pvar
    // was caught by the per-read bound; too short by nothing. Now both.
    if(raw_variant_ct != M) {
        Rcpp::stop("The .pgen file (" + pgenFile + ") holds " + std::to_string(raw_variant_ct) +
                   " variants but the .pvar (" + pvarFile + ") describes " + std::to_string(M) +
                   ". The two files do not describe the same dataset -- a .pvar row missing or "
                   "added anywhere would relabel every marker after it.");
    }
    if(header_ctrl & 0x30) {
        _allele_idx_offsetsp = plink2::CreateRefcountedWptr(cur_variant_ct + 1);
        _info_ptr->allele_idx_offsets = _allele_idx_offsetsp->p;
    }
    _info_ptr->max_allele_ct = 2;
    if((header_ctrl & 0xc0) == 0xc0) {
        const uintptr_t raw_variant_ctl = plink2::DivUp(cur_variant_ct, plink2::kBitsPerWord);
        _nonref_flagsp = plink2::CreateRefcountedWptr(raw_variant_ctl + 1);
        _info_ptr->nonref_flags = _nonref_flagsp->p;
    }
    _file_sample_ct = _info_ptr->raw_sample_ct;
    unsigned char* pgfi_alloc = nullptr;
    if(plink2::cachealigned_malloc(pgfi_alloc_cacheline_ct * plink2::kCacheline, &pgfi_alloc)) {
        Rcpp::stop("Out of memory");
    }
    uint32_t max_vrec_width;
    uintptr_t pgr_alloc_cacheline_ct;
    if(plink2::PgfiInitPhase2(header_ctrl, 1, 0, 0, 0, raw_variant_ct, &max_vrec_width, _info_ptr.get(), pgfi_alloc, &pgr_alloc_cacheline_ct, errstr_buf)) {
        if(pgfi_alloc && (!_info_ptr->vrtypes)) {
            plink2::aligned_free(pgfi_alloc);
        }
        Rcpp::stop(&(errstr_buf[7]));
    }
    if((!_allele_idx_offsetsp) && (_info_ptr->gflags & 4)) {
        Rcpp::stop("Multiallelic variants and phase/dosage info simultaneously present; not supported by SAIGE");
    }
    // var_fpos/vrtypes point INTO this block, and every reader's copy of the
    // PgenFileInfo shares those pointers -- so it must outlive them all, which
    // is what holding it here (rather than in a reader) buys.
    _pgfi_alloc = pgfi_alloc;
    _max_vrec_width = max_vrec_width;
    _pgr_alloc_cacheline_ct = pgr_alloc_cacheline_ct;
    _is_not_plink1_bed = (_info_ptr->vrtypes != nullptr);
    _has_dosage = (_info_ptr->gflags & plink2::kfPgenGlobalDosagePresent) != 0;
    _is_multiallelic = (_info_ptr->max_allele_ct != 2);
}

PgenShared::~PgenShared()
{
    if(_info_ptr) {
        CondReleaseRefcountedWptr(&_allele_idx_offsetsp);
        CondReleaseRefcountedWptr(&_nonref_flagsp);
        plink2::PglErr reterr = plink2::kPglRetSuccess;
        plink2::CleanupPgfi(_info_ptr.get(), &reterr);
        _info_ptr.reset();
    }
    if(_pgfi_alloc) {
        plink2::aligned_free(_pgfi_alloc);
        _pgfi_alloc = nullptr;
    }
}

static void to_upper_inplace(std::string& str) {
    std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c) { return std::toupper(c); });
}

// LEDGER #89: a malformed numeric field names the file, line and column instead
// of aborting through Rcpp's generic path as "stoi: no conversion".
static uint32_t parsePosField(const std::string& s, const std::string& file, size_t lineNo, const char* col)
{
    try {
        size_t used = 0;
        long long v = std::stoll(s, &used);
        if(used != s.size() || v < 0 || v > 0xFFFFFFFFLL) throw std::invalid_argument("range");
        return static_cast<uint32_t>(v);
    } catch(const std::exception&) {
        Rcpp::stop("Cannot parse the " + std::string(col) + " field '" + s + "' on line " +
                   std::to_string(lineNo) + " of " + file + " as a position.");
    }
}

void PgenShared::readPvarFile()
{
    // LEDGER #88 FIXED: the same hardening readPsamFile has -- refuse an
    // unopenable file, strip \r, skip blank lines, refuse a row shorter than
    // the header needs -- where a blank line used to index a split vector out
    // of range before failing downstream in std::stoi.
    std::ifstream bim(pvarFile);
    if(!bim.is_open()) {
        Rcpp::stop("Cannot open pvar file (" + pvarFile + ").");
    }
    std::string line;
    size_t lineNo = 0;

    int col_pos = -1;
    int col_id = -1;
    int col_ref = -1;
    int col_alt = -1;

    while(getline(bim, line)) {
        lineNo++;
        boost::replace_all(line, "\r", "");
        if(line.empty()) continue;
        std::vector<std::string> line_elements;
        boost::split(line_elements, line, boost::is_any_of("\t "));

        if(line[0] == '#') {
            if(line.rfind("#CHROM", 0) == 0) {
                for(size_t i = 0; i < line_elements.size(); i++) {
                    if(line_elements[i] == "POS") col_pos = i;
                    if(line_elements[i] == "ID") col_id = i;
                    if(line_elements[i] == "REF") col_ref = i;
                    if(line_elements[i] == "ALT") col_alt = i;
                }
                if(col_pos == -1 || col_id == -1 || col_ref == -1 || col_alt == -1) {
                    Rcpp::stop("pvar file does not have pos/id/ref/alt columns (" + pvarFile + ").");
                }
            }
            continue;
        }
        if(col_pos == -1) {
            Rcpp::stop("pvar file has a data row before any '#CHROM' header line (" + pvarFile + ", line " +
                       std::to_string(lineNo) + ").");
        }
        const size_t need = static_cast<size_t>(std::max(std::max(col_pos, col_id), std::max(col_ref, col_alt))) + 1;
        if(line_elements.size() < need) {
            Rcpp::stop("pvar file line " + std::to_string(lineNo) + " has " + std::to_string(line_elements.size()) +
                       " fields but the header needs at least " + std::to_string(need) + " (" + pvarFile + ").");
        }

        chr.push_back(line_elements.at(0));
        variantId.push_back(line_elements.at(col_id));
        position.push_back(parsePosField(line_elements.at(col_pos), pvarFile, lineNo, "POS"));
        to_upper_inplace(line_elements[col_ref]);
        to_upper_inplace(line_elements[col_alt]);
        ref.push_back(line_elements[col_ref]);
        alt.push_back(line_elements[col_alt]);
    }
}

// LEDGER #52: this used to key on line_elements[0] with the comment "// IID".
// Column 0 of a "#FID IID SEX" psam is FID, so step 2 failed outright on any
// psam with FID != IID, and two samples sharing an FID collapsed into one map
// entry -- which also undercounts N0, since that is sampleInPgen.size(). Both
// committed PGEN fixtures have FID == IID, which is the only reason the goldens
// ever passed. The IID column is now found BY NAME from the header, mirroring
// R/Geno.R::readGenoSampleIDs, which has always done it correctly.
void PgenShared::readPsamFile()
{
    std::ifstream psam(psamFile);
    if(!psam.is_open()) {
        Rcpp::stop("Cannot open psam file (" + psamFile + ").");
    }

    std::string line;
    int col_iid = -1;
    size_t idx = 0;
    while(getline(psam, line)) {
        boost::replace_all(line, "\r", "");
        if(line.empty()) continue;

        std::vector<std::string> line_elements;
        boost::split(line_elements, line, boost::is_any_of("\t "));

        if(line[0] == '#') {
            // '##' lines are metadata; the single-'#' line is the column header.
            if(col_iid == -1 && line.rfind("##", 0) != 0) {
                line_elements[0] = line_elements[0].substr(1);
                for(size_t i = 0; i < line_elements.size(); i++) {
                    if(line_elements[i] == "IID") col_iid = i;
                }
                if(col_iid == -1) {
                    Rcpp::stop("psam file does not have an IID column (" + psamFile + ").");
                }
            }
            continue;
        }

        if(col_iid == -1) {
            Rcpp::stop("psam file has no '#'-prefixed column header (" + psamFile + ").");
        }
        if(static_cast<size_t>(col_iid) >= line_elements.size()) {
            Rcpp::stop("psam file has a row with fewer columns than the header (" +
                       psamFile + ").");
        }
        const std::string& iid = line_elements[col_iid];
        // Silently overwriting here is what turned a duplicate into an undercount.
        if(sampleInPgen.count(iid)) {
            Rcpp::stop("psam file has a duplicate IID (" + iid + ") in " + psamFile + ".");
        }
        sampleInPgen[iid] = idx;
        idx++;
    }
}

void PgenShared::setPosSampleInPgen(const std::vector<std::string>& sampleInModel)
{
    auto m_N = sampleInModel.size();
    posSampleInPgen.reserve(m_N);

    for(auto& sample : sampleInModel) {
        auto res = sampleInPgen.find(sample);
        if(res == sampleInPgen.end()) {
            Rcpp::stop("At least one subject requested is not in psam file ("+sample+").");
        }
        posSampleInPgen.push_back(res->second);
    }
}

// ---- PgenClass: one decode handle ------------------------------------------

PgenClass::PgenClass(std::string pgenFile,
    std::string psamFile,
    std::string pvarFile,
    std::vector<std::string>& sampleInModel)
    : m_shared(std::make_shared<PgenShared>(pgenFile, psamFile, pvarFile, sampleInModel)),
      _pgr_alloc(nullptr), _subset_size(0)
{
    initReader();
}

PgenClass::PgenClass(std::string pgenFile,
    std::string psamFile,
    std::string pvarFile)
    : m_shared(std::make_shared<PgenShared>(pgenFile, psamFile, pvarFile)),
      _pgr_alloc(nullptr), _subset_size(0)
{
    initReader();
}

PgenClass::PgenClass(std::shared_ptr<PgenShared> t_shared)
    : m_shared(std::move(t_shared)), _pgr_alloc(nullptr), _subset_size(0)
{
    if(!m_shared) {
        Rcpp::stop("PgenClass: no shared pgen state");
    }
    initReader();
}

PgenClass::~PgenClass()
{
    closegenofile();
}

void PgenClass::initReader()
{
    // similar to pgenlibr/src/pgenlibr.cpp. Everything allocated here is this
    // handle's own: the plink2::PgenReader (its FILE*, fread buffer, LD cache
    // and workspaces) and every buffer carved out of pgr_alloc below.
    //
    // The lock is for PgrInit alone: the first reader onto a file MOVES
    // shared_ff out of the PgenFileInfo, and a second thread arriving at that
    // moment would either lose the handle or fopen a second one. Reads take
    // nothing.
    std::lock_guard<std::mutex> lock(m_shared->initMutex);

    plink2::PgenFileInfo* info_ptr = m_shared->info();
    const uint32_t file_sample_ct = m_shared->fileSampleCt();
    const uint32_t max_vrec_width = m_shared->maxVrecWidth();
    const uintptr_t pgr_alloc_cacheline_ct = m_shared->pgrAllocCachelineCt();
    char errstr_buf[plink2::kPglErrstrBufBlen];

    _state_ptr = std::unique_ptr<plink2::PgenReader>(new plink2::PgenReader);

    if(! _state_ptr) {
        Rcpp::stop("Out of memory");
    }
    plink2::PreinitPgr(_state_ptr.get());
    plink2::PgrSetFreadBuf(nullptr, _state_ptr.get());
    const uintptr_t pgr_alloc_main_byte_ct = pgr_alloc_cacheline_ct * plink2::kCacheline;
    const uintptr_t sample_subset_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerVec) * plink2::kBytesPerVec;
    const uintptr_t cumulative_popcounts_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerWord * plink2::kInt32PerVec) * plink2::kBytesPerVec;
    const uintptr_t genovec_byte_ct = plink2::DivUp(file_sample_ct, plink2::kNypsPerVec) * plink2::kBytesPerVec;
    const uint32_t is_not_plink1_bed = m_shared->isNotPlink1Bed();
    uintptr_t raregeno_byte_ct = 0;
    uintptr_t difflist_sample_ids_byte_ct = 0;
    if(is_not_plink1_bed) {
        const uint32_t max_stored_single_difflist_len = file_sample_ct / plink2::kPglMaxDifflistLenDivisor;
        raregeno_byte_ct = plink2::DivUp(2 * max_stored_single_difflist_len, plink2::kNypsPerVec) * plink2::kBytesPerVec;
        difflist_sample_ids_byte_ct = plink2::RoundUpPow2(3 * max_stored_single_difflist_len * sizeof(int32_t), plink2::kBytesPerVec);
    }
    const uintptr_t ac_byte_ct = plink2::RoundUpPow2(file_sample_ct * sizeof(plink2::AlleleCode), plink2::kBytesPerVec);
    const uintptr_t ac2_byte_ct = plink2::RoundUpPow2(file_sample_ct * 2 * sizeof(plink2::AlleleCode), plink2::kBytesPerVec);
    uintptr_t multiallelic_hc_byte_ct = 0;
    if(info_ptr->max_allele_ct != 2) {
        multiallelic_hc_byte_ct = 2 * sample_subset_byte_ct + ac_byte_ct + ac2_byte_ct;
    }
    const uintptr_t dosage_main_byte_ct = plink2::DivUp(file_sample_ct, (2 * plink2::kInt32PerVec)) * plink2::kBytesPerVec;
    unsigned char* pgr_alloc;
    if(plink2::cachealigned_malloc(pgr_alloc_main_byte_ct + (2 * plink2::kPglNypTransposeBatch + 5) * sample_subset_byte_ct + cumulative_popcounts_byte_ct + (1 + plink2::kPglNypTransposeBatch) * genovec_byte_ct + raregeno_byte_ct + difflist_sample_ids_byte_ct + multiallelic_hc_byte_ct + dosage_main_byte_ct + plink2::kPglBitTransposeBufbytes + 4 * (plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8), &pgr_alloc)) {
        Rcpp::stop("Out of memory");
    }
    // The filename must be the PGEN, not the pvar. It goes unread by the FIRST
    // reader onto a file -- that one moves the header's already-open shared_ff
    // across -- which is why passing pvarFile here was invisible for as long as
    // there was only ever one reader. A clone is the second reader: it fopen()s
    // this name, so a wrong one would decode a text file as genotypes.
    plink2::PglErr reterr = plink2::PgrInit(m_shared->pgenFile.c_str(), max_vrec_width, info_ptr, _state_ptr.get(), pgr_alloc);
    if(reterr != plink2::kPglRetSuccess) {
        if(!plink2::PgrGetFreadBuf(_state_ptr.get())) {
            plink2::aligned_free(pgr_alloc);
        }
        snprintf(errstr_buf, plink2::kPglErrstrBufBlen, "PgrInit() error %d", static_cast<int>(reterr));
        Rcpp::stop(errstr_buf);
    }
    _pgr_alloc = pgr_alloc;
    unsigned char* pgr_alloc_iter = &(pgr_alloc[pgr_alloc_main_byte_ct]);
    _subset_include_vec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _subset_include_interleaved_vec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);

#ifdef USE_AVX2
    _subset_include_interleaved_vec[-3] = 0;
    _subset_include_interleaved_vec[-2] = 0;
#endif
    _subset_include_interleaved_vec[-1] = 0;

    _subset_cumulative_popcounts = reinterpret_cast<uint32_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[cumulative_popcounts_byte_ct]);
    _pgv.genovec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[genovec_byte_ct]);
    if(is_not_plink1_bed) {
        _raregeno_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[raregeno_byte_ct]);
        _difflist_sample_ids_buf = reinterpret_cast<uint32_t*>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[difflist_sample_ids_byte_ct]);
    } else {
        _raregeno_buf = nullptr;
        _difflist_sample_ids_buf = nullptr;
    }
    if(multiallelic_hc_byte_ct) {
        _pgv.patch_01_set = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
        _pgv.patch_01_vals = reinterpret_cast<plink2::AlleleCode*>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[ac_byte_ct]);
        _pgv.patch_10_set = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
        _pgv.patch_10_vals = reinterpret_cast<plink2::AlleleCode*>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[ac2_byte_ct]);
    } else {
        _pgv.patch_01_set = nullptr;
        _pgv.patch_01_vals = nullptr;
        _pgv.patch_10_set = nullptr;
        _pgv.patch_10_vals = nullptr;
    }
    _pgv.phasepresent = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _pgv.phaseinfo = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _pgv.dosage_present = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _pgv.dosage_main = reinterpret_cast<uint16_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[dosage_main_byte_ct]);
    _subset_size = file_sample_ct;

    _transpose_batch_buf = reinterpret_cast<plink2::VecW*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglBitTransposeBufbytes]);
    _multivar_vmaj_geno_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * genovec_byte_ct]);
    _multivar_vmaj_phasepresent_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * sample_subset_byte_ct]);
    _multivar_vmaj_phaseinfo_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * sample_subset_byte_ct]);
    _multivar_smaj_geno_batch_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 4]);
    _multivar_smaj_phaseinfo_batch_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8]);
    _multivar_smaj_phasepresent_batch_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    // pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8]);

    readBuffer.resize(_subset_size);

    // 2026-09-16: the carrier-space read's cap. 2 * (N/8) is exactly what the
    // raregeno buffer above was sized for (pgenlibr's sizing, kept), and what
    // plink2 --glm accepts; pgenlib requires it below raw_sample_ct.
    _max_difflist_len = (file_sample_ct >= 8) ? 2 * (file_sample_ct / 8) : 0;
    const char* e = std::getenv("SAIGE_PGEN_SPARSE");
    _sparseEnabled = !((e != nullptr) && (e[0] != '\0') && (std::strcmp(e, "0") == 0));
}

void PgenClass::Read(std::vector<double>& buf, int variant_idx, int allele_idx) {
    // similar to pgenlibr/src/pgenlibr.cpp
    if(!m_shared || !m_shared->info()) {
        Rcpp::stop("pgen is closed");
    }
    if(static_cast<uint32_t>(variant_idx) >= m_shared->info()->raw_variant_ct) {
        Rcpp::stop("variant_num out of range ("+std::to_string(variant_idx + 1)+"; must be 1.."+std::to_string(m_shared->info()->raw_variant_ct)+")");
    }
    if(buf.size() != _subset_size) {
        Rcpp::stop("buf has wrong length ("+std::to_string(buf.size())+"; "+std::to_string(_subset_size)+" expected)");
    }
    uint32_t dosage_ct;
    plink2::PglErr reterr = plink2::PgrGet1D(_subset_include_vec, _subset_index, _subset_size, variant_idx, allele_idx, _state_ptr.get(), _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, &dosage_ct);
    if(reterr != plink2::kPglRetSuccess) {
        Rcpp::stop("PgrGet1D() error " + std::to_string(static_cast<int>(reterr)));
    }
    plink2::Dosage16ToDoubles(kGenoRDoublePairs, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct, &buf[0]);
}

// Releases THIS handle. The shared file info goes with the LAST handle onto it,
// so closing one reader can never invalidate another (the same ownership
// contract PlinkReader states, reached a different way -- there the clones
// share one decoder, here each has its own).
void PgenClass::closegenofile() {
    if(_state_ptr) {
        plink2::PglErr reterr = plink2::kPglRetSuccess;
        plink2::CleanupPgr(_state_ptr.get(), &reterr);
        if(PgrGetFreadBuf(_state_ptr.get())) {
            plink2::aligned_free(PgrGetFreadBuf(_state_ptr.get()));
            plink2::PgrSetFreadBuf(nullptr, _state_ptr.get());
            _pgr_alloc = nullptr;   // fread_buf IS the base of pgr_alloc
        }
        _state_ptr.reset();
    }
    if(_pgr_alloc) {
        plink2::aligned_free(_pgr_alloc);
        _pgr_alloc = nullptr;
    }
    _subset_size = 0;
    readBuffer.clear();
    m_shared.reset();
}

void PgenClass::fillMeta(uint64_t t_gIndex, std::string& t_ref, std::string& t_alt,
                         std::string& t_marker, uint32_t& t_pd, std::string& t_chr) const
{
    const PgenShared& sh = *m_shared;
    t_chr = sh.chr[t_gIndex];
    t_pd = sh.position[t_gIndex];
    t_ref = sh.ref[t_gIndex];
    t_alt = sh.alt[t_gIndex];
    t_marker = sh.variantId[t_gIndex];
}

// The per-sample expansion of readBuffer (one double per FILE sample, NaN when
// missing) into the analysis-ordered outputs. This is getOneMarker's loop,
// moved out unchanged on 2026-09-16 so that the dense fallback of
// getOneMarkerSparse runs the identical code.
void PgenClass::expandFromReadBuffer(double& t_altFreq, double& t_altCounts,
                                     double& t_missingRate, double& t_imputeInfo,
                                     bool t_isOutputIndexForMissing, std::vector<uint>& t_indexForMissing,
                                     bool t_isOnlyOutputNonZero, std::vector<uint>& t_indexForNonZero,
                                     arma::vec& OneMarkerG1)
{
    const PgenShared& sh = *m_shared;
    uint32_t numMissing = 0;

    t_indexForMissing.clear();
    t_indexForNonZero.clear();

    uint j = 0;
    t_altCounts = 0;
    const uint32_t N = sh.N;
    for(uint32_t i = 0; i < N; i++)
    {
        auto sampleIdx = sh.posSampleInPgen[i];
        auto numAltReads = readBuffer[sampleIdx];
        if(std::isnan(numAltReads)) {
            numMissing++;
            if(t_isOutputIndexForMissing) {
                t_indexForMissing.push_back(i);
            }
        } else {
            t_altCounts += numAltReads;
        }

        if(numAltReads > 0) {
            t_indexForNonZero.push_back(i);
        }
        if(t_isOnlyOutputNonZero){
            if(numAltReads > 0) {
                OneMarkerG1[j] = numAltReads;
                j = j + 1;
            }
        }else{
            OneMarkerG1[i] = numAltReads;
        }
    }
    int count = N - numMissing;
    t_missingRate = (double)numMissing / (double)N;
    t_imputeInfo = 1;

    if(count > 0){
        t_altFreq = t_altCounts / (double)count / 2;
    }else{
        t_altFreq = 0;
    }

    if(t_isOnlyOutputNonZero){
        OneMarkerG1.resize(j);
    }
}

void PgenClass::getOneMarker(
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
    arma::vec& OneMarkerG1)
{
    if(!m_shared) {
        Rcpp::stop("PgenClass::getOneMarker called after the pgen file was closed.");
    }
    fillMeta(t_gIndex, t_ref, t_alt, t_marker, t_pd, t_chr);

    // pvar is always REF/ALT-ordered (no ref-first/alt-first toggle like
    // bim), and allele_idx = 1 reads ALT-allele dosage directly -- so, unlike
    // PLINK's ref-first branch (LEDGER #28/#29), there is no round trip here:
    // t_altCounts/t_altFreq are accumulated straight from the decoded counts.
    int alleleIdx = 1;
    Read(readBuffer, t_gIndex, alleleIdx);

    expandFromReadBuffer(t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
                         t_isOutputIndexForMissing, t_indexForMissing,
                         t_isOnlyOutputNonZero, t_indexForNonZero, OneMarkerG1);
}

// Whether this file can offer difference lists at all. A plink1-compatible
// .pgen stores every record densely (and initReader allocated no raregeno
// buffer for it); a dosage track is not carried by the sparse form and
// PgrGet1D's per-allele counting is not reproduced for a multiallelic site,
// so both fall back to the dense read for every marker.
bool PgenClass::sparseSupported() const
{
    if(!m_shared || !_sparseEnabled) return false;
    if(_raregeno_buf == nullptr || _difflist_sample_ids_buf == nullptr) return false;
    if(_max_difflist_len == 0 || _max_difflist_len >= _subset_size) return false;
    if(m_shared->hasDosage() || m_shared->isMultiallelic()) return false;
    return true;
}

// The carrier-space read. PgrGetDifflistOrGenovec hands back the stored
// difference list when the record is one and it is no longer than the cap,
// else the dense genovec -- pgenlib never converts a dense record to a list.
// With sample_ct == raw_sample_ct (this handle decodes every file sample and
// maps to analysis order itself, as Read() does) the list's sample ids are
// pgen rows, mapped through analysisIndexOfFileSample; a row outside the
// model is dropped, and the survivors are sorted into analysis order (the
// order the dense loop visits samples in), which is a no-op when the model's
// samples are in file order. The counts are the dense loop's counts: an
// integer sum of 0/1/2 over the analysis samples, so accumulating them from
// the four genotype tallies is the same double.
int PgenClass::getOneMarkerSparse(
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
    arma::vec& t_GVec)
{
    if(!m_shared) {
        Rcpp::stop("PgenClass::getOneMarkerSparse called after the pgen file was closed.");
    }
    if(!sparseSupported()) {
        return 0;
    }
    const PgenShared& sh = *m_shared;
    if(t_gIndex >= sh.info()->raw_variant_ct) {
        Rcpp::stop("variant_num out of range ("+std::to_string(t_gIndex + 1)+"; must be 1.."+std::to_string(sh.info()->raw_variant_ct)+")");
    }
    fillMeta(t_gIndex, t_ref, t_alt, t_marker, t_pd, t_chr);

    uint32_t common = UINT32_MAX;
    uint32_t len = 0;
    plink2::PglErr reterr = plink2::PgrGetDifflistOrGenovec(
        _subset_include_vec, _subset_index, _subset_size, _max_difflist_len,
        static_cast<uint32_t>(t_gIndex), _state_ptr.get(), _pgv.genovec,
        &common, _raregeno_buf, _difflist_sample_ids_buf, &len);
    if(reterr != plink2::kPglRetSuccess) {
        Rcpp::stop("PgrGetDifflistOrGenovec() error " + std::to_string(static_cast<int>(reterr)));
    }

    if(common == UINT32_MAX || (common != 0 && common != 2)) {
        // Dense: expand exactly as Read() + getOneMarker do for a hard-call
        // record (no dosage track, so dosage_ct = 0 and Dosage16ToDoubles is
        // the plain 2-bit lookup), then the same per-sample loop.
        if(common != UINT32_MAX) {
            plink2::PgrDifflistToGenovecUnsafe(_raregeno_buf, _difflist_sample_ids_buf, common,
                                               _subset_size, len, _pgv.genovec);
        }
        plink2::Dosage16ToDoubles(kGenoRDoublePairs, _pgv.genovec, _pgv.dosage_present,
                                  _pgv.dosage_main, _subset_size, 0, &readBuffer[0]);
        std::vector<uint> dummyNonZero;
        expandFromReadBuffer(t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
                             true, t_indexForMissing, false, dummyNonZero, t_GVec);
        return 1;
    }

    t_idx.clear();
    t_code.clear();
    const std::vector<int32_t>& inv = sh.analysisIndexOfFileSample;
    uint32_t cnt[4] = {0, 0, 0, 0};
    bool sorted = true;
    uint32_t last = 0;
    for(uint32_t e = 0; e < len; e++) {
        const uint32_t row = _difflist_sample_ids_buf[e];
        const int32_t i = inv[row];
        if(i < 0) continue;
        const uint8_t c = static_cast<uint8_t>(
            (_raregeno_buf[e / plink2::kBitsPerWordD2] >> (2 * (e % plink2::kBitsPerWordD2))) & 3);
        if(!t_idx.empty() && static_cast<uint32_t>(i) < last) sorted = false;
        last = static_cast<uint32_t>(i);
        t_idx.push_back(static_cast<uint32_t>(i));
        t_code.push_back(c);
        cnt[c]++;
    }
    if(!sorted) {
        std::vector<uint32_t> order(t_idx.size());
        std::iota(order.begin(), order.end(), 0u);
        std::sort(order.begin(), order.end(), [&](uint32_t x, uint32_t y){ return t_idx[x] < t_idx[y]; });
        std::vector<uint32_t> idx2(t_idx.size());
        std::vector<uint8_t> code2(t_code.size());
        for(size_t o = 0; o < order.size(); o++) { idx2[o] = t_idx[order[o]]; code2[o] = t_code[order[o]]; }
        t_idx.swap(idx2);
        t_code.swap(code2);
    }

    const uint32_t N = sh.N;
    const uint32_t listed = static_cast<uint32_t>(t_idx.size());
    cnt[common] += N - listed;
    const uint32_t numMissing = cnt[3];
    const int count = static_cast<int>(N) - static_cast<int>(numMissing);
    t_altCounts = static_cast<double>(cnt[1] + 2u * cnt[2]);
    t_missingRate = (double)numMissing / (double)N;
    t_imputeInfo = 1;
    if(count > 0){
        t_altFreq = t_altCounts / (double)count / 2;
    }else{
        t_altFreq = 0;
    }
    t_commonCode = common;
    return 2;
}

}
