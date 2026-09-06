// This includes the main codes to connect C++ and R

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <vector>
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
// std::this_thread::sleep_for (std::chrono::seconds(1));
#include <cstdio>         // std::remove
#include <fstream>
#include <string.h>
// Currently, omp does not work well, will check it later
// error: SET_VECTOR_ELT() can only be applied to a 'list', not a 'character'
// remove all Rcpp::List to check if it works
// #include <omp.h>
// // [[Rcpp::plugins(openmp)]]]

#include <RcppParallel.h>   // S6: parallelFor for the marker scan
#include <mutex>
#include "Main.hpp"
#include <cstring>
#include <cstdlib>
#include <functional>
#include "Binary_global.hpp"   // S6: SKATMainThread, the ER handback
#include "Binary_HyperGeo.hpp" // lchoose_threadsafe, for its regression test
#include "GenoReader.hpp"
#include "SAIGE_test.hpp"
#include "UTIL.hpp"
#include "CCT.hpp"

#include <Rcpp.h>
#include "getMem.hpp"

#include <boost/math/distributions/beta.hpp>

// the genotype reader for the current analysis (see GenoReader.hpp)

// M1: the marker scan can carry several null models over one genotype read.
// g_models owns them in MANIFEST order; g_readers owns ONE reader per
// sample-set group (models whose sampleID vectors are identical() share a
// reader, a decode and a QC pass); g_modelGroup maps a model to its group.
//
// g_model and g_reader are NON-OWNING aliases of element 0. Everything outside
// the marker loop -- the region path, the conditioning setup,
// set_flagSparseGRM_cur_SAIGE_org, the S5 test hooks -- keeps using them and
// keeps meaning exactly what it meant when there could only be one model, which
// is why P > 1 is refused for those paths in R rather than half-supported here.
static std::vector<std::unique_ptr<Geno::GenoReader>> g_readers;
static Geno::GenoReader* g_reader = nullptr;
// global objects for different analysis methods
static std::vector<std::unique_ptr<SAIGE::SAIGEClass>> g_models;
static SAIGE::SAIGEClass* g_model = nullptr;
static std::vector<int> g_modelGroup;

// Install the FIRST reader (a new run) or an additional one (a second
// sample-set group). The owned objects never move, so g_reader survives the
// vector reallocating.
static void resetReaders(std::unique_ptr<Geno::GenoReader> t_first)
{
  g_readers.clear();
  g_readers.push_back(std::move(t_first));
  g_reader = g_readers[0].get();
}

static void addReader(std::unique_ptr<Geno::GenoReader> t_next)
{
  if(g_readers.empty()){
    Rcpp::stop("A genotype reader for the first sample-set group must be set up first.");
  }
  g_readers.push_back(std::move(t_next));
  g_reader = g_readers[0].get();
}

// The fastTest ESCALATION signal for region tests, set from R between calls to
// mainRegionInCPP (set_flagSparseGRM_cur_SAIGE / _org): false on the cheap
// first pass, true on the accurate re-run whose answer is kept.
//
// S12: this used to be SAIGEClass::m_flagSparseGRM_cur, where it sat next to the
// per-marker variance-ratio choice and shared a name with it -- so the region
// loops, which select a variance ratio per marker from MAC, overwrote the
// caller's intent before the write decisions could read it. That collision is
// LEDGER #27 and its twin #50. It is step-2 loop state rather than model state,
// it is written only between region calls, and mainRegionInCPP reads it exactly
// once, at entry.
static bool g_regionEscalatedPass = false;
//single, SAIGE
//Region, SAIGE-GENE+

// global variables for analysis
struct Step2Config {
  std::string impute_method;      // "mean", "minor", or "drop", //drop is not allowed
  double missingRate_cutoff = 0;
  //unsigned int g_omp_num_threads;
  double marker_minMAF_cutoff = 0;
  double marker_minMAC_cutoff = 0;
  double region_minMAC_cutoff = 0;    // for Rare Variants (RVs) whose MAC < this value, we aggregate these variants like SAIGE-GENE+
  double marker_minINFO_cutoff = 0;
  arma::vec region_maxMAF_cutoff;
  double min_gourpmac_for_burdenonly = 0;
  double maxMAFLimit = 0;
  unsigned int region_maxMarkers_cutoff = 0;   // maximal number of markers in one chunk, only used for region-based analysis to reduce memory usage
  bool isOutputMoreDetails = false;
  int marker_chunksize = 0;


  double dosage_zerod_MAC_cutoff = 0;
  double dosage_zerod_cutoff = 0;
  arma::vec weights_beta = arma::vec(2);

  double MACCutoffforER = 0;
  bool is_rewrite_XnonPAR_forMales = false;


  arma::uvec indexInModel_male;

  arma::umat X_PARregion_mat;
};
static Step2Config g_cfg;


// One single-variant output file per null model. unique_ptr because an open
// std::ofstream must never be moved, and a vector reallocating on push_back
// would do exactly that.
struct SingleWriter {
  std::ofstream stream;
  std::string prefix;
};

struct Step2Output {
  std::ofstream group, singleInGroup, singleInGroupTemp;
  std::string prefixGroup, prefixSingleInGroup, prefixSingleInGroupTemp;
  std::vector<std::unique_ptr<SingleWriter>> singles;
};
static Step2Output g_out;





// [[Rcpp::export]]
void setAssocTest_GlobalVarsInCPP(std::string t_impute_method,
		double t_missing_cutoff,
                               double t_min_maf_marker,
                               double t_min_mac_marker,
                               double t_min_info_marker,
			       double t_dosage_zerod_cutoff,
                               double t_dosage_zerod_MAC_cutoff,
			       arma::vec & t_weights_beta, 
			       std::string t_outputFilePrefix,
			       double t_MACCutoffforER)
{
  g_cfg.impute_method = t_impute_method;
  g_cfg.missingRate_cutoff = t_missing_cutoff;
  g_cfg.marker_minMAF_cutoff = t_min_maf_marker;
  g_cfg.marker_minMAC_cutoff = t_min_mac_marker;
  g_cfg.marker_minINFO_cutoff = t_min_info_marker;
  g_cfg.dosage_zerod_cutoff = t_dosage_zerod_cutoff;
  g_cfg.dosage_zerod_MAC_cutoff = t_dosage_zerod_MAC_cutoff;
  g_cfg.weights_beta = t_weights_beta;
  g_out.prefixGroup = t_outputFilePrefix;
  g_out.prefixSingleInGroup = t_outputFilePrefix + ".singleAssoc.txt";
  g_out.prefixSingleInGroupTemp = t_outputFilePrefix + ".singleAssoc.txt_temp";
  // One writer by default: the single-model run, whose file is the run prefix.
  // setSingleOutputPrefixesInCPP resizes this when a manifest asks for more.
  g_out.singles.clear();
  g_out.singles.push_back(std::unique_ptr<SingleWriter>(new SingleWriter()));
  g_out.singles[0]->prefix = t_outputFilePrefix;
  g_cfg.MACCutoffforER = t_MACCutoffforER;
}

// [[Rcpp::export]]
void setAssocTest_GlobalVarsInCPP_indexInModel_male(arma::uvec & t_indexInModel_male){
  g_cfg.indexInModel_male = t_indexInModel_male;
  g_cfg.is_rewrite_XnonPAR_forMales = true;
}

// [[Rcpp::export]]
void setAssocTest_GlobalVarsInCPP_X_PARregion_mat(arma::umat & t_X_PARregion_mat){
  g_cfg.X_PARregion_mat = t_X_PARregion_mat;
}



// [[Rcpp::export]]
void setMarker_GlobalVarsInCPP(
			       bool t_isOutputMoreDetails,
			       int t_marker_chunksize
			       )

{
  g_cfg.isOutputMoreDetails = t_isOutputMoreDetails;
  g_cfg.marker_chunksize = t_marker_chunksize;
}


//double t_DosageCutoff_for_UltraRarePresence,
			       //std::string t_method_to_CollapseUltraRare,

// [[Rcpp::export]]
void setRegion_GlobalVarsInCPP(
                               arma::vec t_max_maf_region,
                               unsigned int t_max_markers_region,
			       double t_MACCutoff_to_CollapseUltraRare, 
			       double t_min_gourpmac_for_burdenonly)
{
  g_cfg.region_maxMAF_cutoff = t_max_maf_region;
  g_cfg.maxMAFLimit = g_cfg.region_maxMAF_cutoff.max();
  g_cfg.region_maxMarkers_cutoff = t_max_markers_region;
  g_cfg.region_minMAC_cutoff = t_MACCutoff_to_CollapseUltraRare;
  g_cfg.min_gourpmac_for_burdenonly = t_min_gourpmac_for_burdenonly;
}



static void applyMaleXnonPAR(arma::vec& t_GVec, uint32_t& t_pd, double& t_altCounts, double& t_altFreq, const std::string& t_chr);

// A printed p-value back to a double, tolerating what the writer can print
// and std::stod cannot read: a subnormal (out_of_range) reads as 0, and a
// string that is not a number at all reads as 0 with a note (LEDGER #72).
static double parsePrintedPval(const std::string& s){
  try {
    return std::stod(s);
  } catch (const std::out_of_range&) {
    return 0.0;
  } catch (const std::invalid_argument&) {
    std::cerr << "Cannot parse the printed p-value '" << s << "'; taking 0\n";
    return 0.0;
  }
}

// S6: the read, addressed to an EXPLICIT reader rather than to the g_reader
// singleton. Unified_getOneMarker (the Rcpp export, further down) is now a thin
// wrapper that passes *g_reader; a worker thread passes its own clone. Nothing
// else changes -- same dispatch, same male-X rewrite, same return.
static bool getOneMarkerThrough(Geno::GenoReader& t_reader,
                                const std::string& t_genoType,
                                uint64_t t_gIndex,
                                std::string& t_ref, std::string& t_alt,
                                std::string& t_marker, uint32_t& t_pd,
                                std::string& t_chr, double& t_altFreq,
                                double& t_altCounts, double& t_missingRate,
                                double& t_imputeInfo,
                                bool t_isOutputIndexForMissing,
                                std::vector<uint>& t_indexForMissing,
                                bool t_isOnlyOutputNonZero,
                                std::vector<uint>& t_indexForNonZero,
                                arma::vec& t_GVec,
                                bool t_isImputation)
{
  bool isBoolRead = true;
  if(t_genoType == "plink" || t_genoType == "pgen"){
   Geno::MarkerMeta meta;
   isBoolRead = t_reader.getOneMarker(t_gIndex, meta,
                                      t_isOutputIndexForMissing, t_indexForMissing,
                                      t_isOnlyOutputNonZero, t_indexForNonZero,
                                      t_GVec);
   t_ref = meta.ref; t_alt = meta.alt; t_marker = meta.id;
   t_pd = meta.pos; t_chr = meta.chr;
   t_altFreq = meta.altFreq; t_altCounts = meta.altCounts;
   t_missingRate = meta.missingRate; t_imputeInfo = meta.imputeInfo;
  }

  applyMaleXnonPAR(t_GVec, t_pd, t_altCounts, t_altFreq, t_chr);

  return isBoolRead;
}


// Phase 0 instrumentation: an OFF-BY-DEFAULT stopwatch over the three phases of
// testOneMarker, so the shared/per-model split that multi-model scanning turns
// on can be predicted from a measurement rather than from arithmetic.
//
// Gated on the environment variable SAIGE_PROFILE_MARKER (precedent:
// SAIGE_PLINK_NO_MMAP, src/PLINK.cpp). When unset the only cost is a null
// pointer check per phase and nothing is printed, so no golden can move. It is
// reported for the SERIAL arm only -- the accumulators are plain statics, and
// summing wall time across workers would answer a different question anyway.
namespace {

struct MarkerProfile {
  bool enabled = false;
  double prepare = 0.0;   // read + QC gate 1 + impute/flip + QC gate 2
  double test = 0.0;      // the score/SPA/ER calls and the fastTest re-test
  double summary = 0.0;   // per-model row fields, incl. the case/ctrl gather
  long long nPrepare = 0;
  long long nTest = 0;
  long long nSummary = 0;

  void reset() {
    prepare = test = summary = 0.0;
    nPrepare = nTest = nSummary = 0;
  }
};

static MarkerProfile g_markerProfile;

static bool markerProfileRequested()
{
  const char* e = std::getenv("SAIGE_PROFILE_MARKER");
  return (e != nullptr) && (e[0] != '\0') && (std::strcmp(e, "0") != 0);
}

// The region-path twin, gated on SAIGE_PROFILE_REGION. mainRegionInCPP is
// called once per region, so the buckets accumulate across calls and one
// cumulative line is printed per call; the last such line in a log is the
// run's total. Buckets follow the function's passes: the per-marker loop
// (prepare = read + impute/flip + QC + weights; rareArm = everything on the
// non-ultra-rare branch INCLUDING test; test = the getMarkerPval calls inside
// it; urCollapse = the ultra-rare branch; flush = chunk-file writes), the
// ultra-rare pseudo-marker pass (urPass, with its own getMarkerPval calls in
// urTest), the VarMat assembly (varmat), and the emitters (emit, with the
// BURDEN getMarkerPval calls in emitTest). total is the whole call.
struct RegionProfile {
  bool enabled = false;
  long long calls = 0;
  double prepare = 0.0, rareArm = 0.0, test = 0.0, urCollapse = 0.0, flush = 0.0;
  double urPass = 0.0, urTest = 0.0, varmat = 0.0, emit = 0.0, emitTest = 0.0, total = 0.0;
  long long nPrepare = 0, nRareArm = 0, nTest = 0, nUrCollapse = 0, nFlush = 0;
  long long nUrPass = 0, nUrTest = 0, nVarmat = 0, nEmit = 0, nEmitTest = 0, nTotal = 0;
};

static RegionProfile g_regionProfile;

static bool regionProfileRequested()
{
  const char* e = std::getenv("SAIGE_PROFILE_REGION");
  return (e != nullptr) && (e[0] != '\0') && (std::strcmp(e, "0") != 0);
}

// Stops on destruction so the early returns inside the prepare phase are timed
// too; stop() makes it explicit where a phase ends mid-function.
class ProfileScope {
public:
  ProfileScope(double* t_acc, long long* t_cnt)
    : m_acc(t_acc), m_cnt(t_cnt)
  {
    if(m_acc != nullptr){
      m_t0 = std::chrono::steady_clock::now();
    }
  }
  ~ProfileScope() { stop(); }
  void stop()
  {
    if(m_acc == nullptr){
      return;
    }
    *m_acc += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - m_t0).count();
    ++(*m_cnt);
    m_acc = nullptr;
  }
  ProfileScope(const ProfileScope&) = delete;
  ProfileScope& operator=(const ProfileScope&) = delete;
private:
  double* m_acc;
  long long* m_cnt;
  std::chrono::steady_clock::time_point m_t0;
};

}  // namespace

// S4: one marker's analysis, lifted out of mainMarkerInCPP's loop body. Every
// value it produces goes into t_res (one slot per marker, never shared); every
// buffer it reuses comes from t_scratch (shared across markers by design, so
// one per worker once this is threaded). Nothing that outlives the call is
// written: the model is const on this path since S3, and the reader is addressed
// by absolute marker index rather than by a stream position.
//
// Returns false only on a short read, which is what the loop's `break` used to
// respond to -- so the caller stops the scan and the remaining slots keep their
// untested defaults, exactly as the q-length vectors did.
// S4 split it out of the loop; M1 splits it again, into the half that depends
// only on the GENOTYPE and the half that depends on a null MODEL. That is the
// whole of multi-model scanning: prepare once per sample-set group, test once
// per model.
//
// The model-independent half: the read, the two QC gates, the mean imputation
// and flip, and the two index-vector copies. Everything a model needs from it
// travels in MarkerPrep; everything the WRITER needs is written into every row
// of the group. Marker fields are duplicated per row rather than shared,
// because two models in DIFFERENT groups legitimately disagree about
// AC/AF/missingRate/N -- they are looking at different samples.
struct MarkerPrep {
  int n = 0;
  double altFreq = 0.0;
  double altCounts = 0.0;
  double MAF = 0.0;
  double MAC = 0.0;
  bool flip = false;
};

// Returns -1 on a short read (the caller stops the scan, which is what the
// loop's `break` used to do), 0 when a QC gate dropped the marker (the rows keep
// their untested defaults, so pval stays "NA"), 1 when the marker is ready.
static int prepareOneMarker(Geno::GenoReader & t_reader,
                            const std::string & t_genoType,
                            const std::string & t_genoIndexStr,
                            bool t_isImputation,
                            int t_n,
                            Marker::MarkerScratch & t_scratch,
                            const std::vector<Marker::MarkerResult*> & t_rows,
                            MarkerPrep & t_prep)
{
  int n = t_n;

  ProfileScope profPrepare(g_markerProfile.enabled ? &g_markerProfile.prepare : nullptr,
                           &g_markerProfile.nPrepare);

  // information of marker
  double altFreq, altCounts, missingRate, imputeInfo;
  std::string chr, ref, alt, marker;
  uint32_t pd;

  bool flip = false;
  char* end;
  uint64_t gIndex = std::strtoull( t_genoIndexStr.c_str(), &end,10 );

  bool isOutputIndexForMissing = true;
  bool isOnlyOutputNonZero = false;

  //clear vectors
  t_scratch.indexZeroVec.clear();
  t_scratch.indexNonZeroVec.clear();
  t_scratch.indexForMissing.clear();

  // Item 3 (2026-09-05): QC gate 1 from the code tallies, before any genotype is
  // expanded. The reader hands back exactly the quantities the full read would
  // (same bytes, same integer counts, same arithmetic), so the decision below
  // is the decision the gate after the read would have made; a marker that
  // passes is then read in full as before, and a marker that fails leaves the
  // rows as the full path left them (basic info set, pval "NA", never written).
  // Skipped when the male X non-PAR rewrite is on: it runs after the read and
  // recomputes altCounts/altFreq from the rewritten genotypes.
  if(!g_cfg.is_rewrite_XnonPAR_forMales){
    Geno::MarkerMeta pre;
    uint32_t numMissingPre = 0;
    if(t_reader.countOneMarker(gIndex, pre, numMissingPre)){
      double MAFpre = std::min(pre.altFreq, 1 - pre.altFreq);
      double nNonMissingPre = (double)(n - (int)numMissingPre);
      double MACpre = std::min(pre.altCounts, 2*nNonMissingPre - pre.altCounts);
      if((pre.missingRate > g_cfg.missingRate_cutoff) || (MAFpre < g_cfg.marker_minMAF_cutoff) || (MACpre < g_cfg.marker_minMAC_cutoff || pre.imputeInfo < g_cfg.marker_minINFO_cutoff)){
        std::string pdsPre = std::to_string(pre.pos);
        for(std::size_t r = 0; r < t_rows.size(); r++){
          Marker::MarkerResult & t_res = *t_rows[r];
          t_res.chr = pre.chr;
          t_res.pos = pdsPre;
          t_res.ref = pre.ref;
          t_res.alt = pre.alt;
          t_res.marker = pre.id;
          t_res.altFreq = pre.altFreq;
          t_res.missingRate = pre.missingRate;
          t_res.imputationInfo = pre.imputeInfo;
        }
        return 0;
      }
    }
  }

  bool isReadMarker = getOneMarkerThrough(t_reader, t_genoType, gIndex, ref, alt, marker, pd, chr,
                                          altFreq, altCounts, missingRate, imputeInfo,
                                          isOutputIndexForMissing, // bool t_isOutputIndexForMissing,
                                          t_scratch.indexForMissing,
                                          isOnlyOutputNonZero, // bool t_isOnlyOutputNonZero,
                                          t_scratch.indexNonZeroVec, t_scratch.GVec, t_isImputation);

  if(!isReadMarker){
    return -1;
  }

  std::string pds = std::to_string(pd);
  // The CHR:POS:REF:ALT `info` string the loop used to build here was written
  // into an infoVec that the writer was never handed; it is dropped with the
  // slot vectors rather than carried into MarkerResult.

  // record basic information for the marker
  for(std::size_t r = 0; r < t_rows.size(); r++){
    Marker::MarkerResult & t_res = *t_rows[r];
    t_res.chr = chr;
    t_res.pos = pds;
    t_res.ref = ref;
    t_res.alt = alt;
    t_res.marker = marker;                    // marker IDs
    t_res.altFreq = altFreq;                  // allele frequencies of ALT allele, this is not always < 0.5.
    t_res.missingRate = missingRate;
    t_res.imputationInfo = imputeInfo;
  }

  // MAF and MAC are for Quality Control (QC)
  double MAF = std::min(altFreq, 1 - altFreq);
  // LEDGER #29 (second of two lossy steps): MAC used to be derived from the
  // FREQUENCY as MAF * n * (1-missingRate) * 2. Both factors are float round
  // trips -- MAF via 1-f when ALT is major, and n*(1-missingRate) via
  // missingRate = numMissing/n -- so a marker sitting exactly on an integer
  // --minMAC cutoff was dropped by the strictly-less-than gate below whenever
  // the error happened to round low. Gate 2 already recomputes this exactly
  // from the counts, but only AFTER gate 1 has discarded the marker.
  // Derive it from counts here too. indexForMissing is populated on this path
  // (isOutputIndexForMissing is set true above), so the non-missing count is
  // exact -- unlike n*(1-missingRate).
  double nNonMissing = (double)(n - (int)t_scratch.indexForMissing.size());
  double MAC = std::min(altCounts, 2*nNonMissing - altCounts);

  // Quality Control (QC) based on missing rate, MAF, and MAC
  if((missingRate > g_cfg.missingRate_cutoff) || (MAF < g_cfg.marker_minMAF_cutoff) || (MAC < g_cfg.marker_minMAC_cutoff || imputeInfo < g_cfg.marker_minINFO_cutoff)){
    return 0;
  }

  t_scratch.indexZeroVec.clear();
  t_scratch.indexNonZeroVec.clear();

  flip = imputeGenoAndFlip(t_scratch.GVec, altFreq, altCounts, t_scratch.indexForMissing, g_cfg.impute_method, g_cfg.dosage_zerod_cutoff, g_cfg.dosage_zerod_MAC_cutoff, MAC, t_scratch.indexZeroVec, t_scratch.indexNonZeroVec);
  MAC = std::min(altCounts, 2*n-altCounts);
  MAF = std::min(altFreq, 1 - altFreq);

  if((MAF < g_cfg.marker_minMAF_cutoff) || (MAC < g_cfg.marker_minMAC_cutoff)){
    return 0;
  }

  for(std::size_t r = 0; r < t_rows.size(); r++){
    Marker::MarkerResult & t_res = *t_rows[r];
    t_res.altFreq = altFreq;                  // allele frequencies of ALT allele, this is not always < 0.5.
    t_res.altCounts = altCounts;
  }

  t_scratch.indexZeroVec_arma = arma::conv_to<arma::uvec>::from(t_scratch.indexZeroVec);
  t_scratch.indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(t_scratch.indexNonZeroVec);

  t_scratch.indexZeroVec.clear();
  t_scratch.indexNonZeroVec.clear();
  // P2Vec / G1tilde_P_G2tilde_Vec are NOT cleared here. Both are fully assigned
  // before any read (SAIGE_test.cpp:166/174/930/934 and :786 respectively), so
  // clearing them only freed the scratch buffers and made the next marker
  // re-allocate n doubles. Values are unaffected.

  t_prep.n = n;
  t_prep.altFreq = altFreq;
  t_prep.altCounts = altCounts;
  t_prep.MAF = MAF;
  t_prep.MAC = MAC;
  t_prep.flip = flip;
  return 1;
}


// The model-dependent half. Body unchanged from S4's testOneMarker apart from
// g_model-> becoming t_model. and the marker quantities coming from t_prep.
//
// getMarkerPval is const and never assigns t_GVec or the index vectors (that is
// what S3/S12 bought), so every model in a group tests THE SAME prepared buffer
// without copying it. The Unified_getMarkerPval wrapper is not used here: all it
// added was a stop on t_isOnlyOutputNonZero, which is a literal false at all
// three call sites.
static void testOneModel(const SAIGE::SAIGEClass & t_model,
                         bool t_isMoreOutput,
                         bool t_isSingleVarianceRatio,
                         const MarkerPrep & t_prep,
                         Marker::MarkerScratch & t_scratch,
                         Marker::MarkerResult & t_res)
{
  const std::string & t_traitType = t_model.m_traitType;
  const bool t_isCondition = t_model.m_isCondition;
  const int n = t_prep.n;
  const double altFreq = t_prep.altFreq;
  const double MAC = t_prep.MAC;
  const bool flip = t_prep.flip;

  // analysis results for single-marker
  double Beta, seBeta, Tstat, varT, gy;
  double Beta_c, seBeta_c, Tstat_c, varT_c;
  std::string pval, pval_noSPA, pval_c, pval_noSPA_c;
  bool isSPAConverge, is_gtilde, is_Firth, is_FirthConverge;
  double AF_case, AF_ctrl;
  uint32_t N_case, N_ctrl;

  // S3: the marker path owns its test context instead of parking it on the
  // model between markers. Same three values, same order of derivation --
  // fastTest forces the non-sparse first pass, then the variance ratio is
  // selected for that choice -- but nothing is written to the model, so no
  // marker can read what another one left behind (LEDGER #32) and two markers
  // could be tested at once without racing.
  bool hasVarRatio = true;
  SAIGE::TestContext ctx;
  ctx.useSparseGRM = t_model.m_isFastTest ? false : t_model.m_flagSparseGRM;
  ctx.noadjCov = t_model.m_isnoadjCov;
  if(t_isSingleVarianceRatio){
      ctx.varRatioVal = t_model.selectSingleVarianceRatio(ctx.useSparseGRM, ctx.noadjCov);
      hasVarRatio = true;
  }else{
      ctx.varRatioVal = t_model.selectVarianceRatio(MAC, ctx.useSparseGRM, ctx.noadjCov, hasVarRatio);
  }

  bool is_region = false;

  // Step B: the blanket ER deferral is GONE, and with it the last reason the
  // in-process path was slower than serial on a rare-variant binary scan.
  //
  // S6 excluded the whole ER branch because Rf_lchoose is called on EVERY ER
  // marker and aborts off the main thread. That was never about the arithmetic:
  // R's C-stack check measures depth against the MAIN thread's stack base, so a
  // worker reads as impossibly deep. lchoose is now computed from Rf_lbeta
  // bit-identically (HyperGeo::lchoose_threadsafe, pinned over 2,009,004
  // argument pairs by tests/testthat/test-unit-lchoose.R), so the
  // exact-enumeration arm runs on a worker like any other marker.
  //
  // What remains is SL_setseed, on the RESAMPLING arms alone, and it is guarded
  // where it lives (Binary_ComputeExact.cpp) rather than predicted from here.
  // That is the point S6 already made and this change finally leans on: the
  // THROW is the safety property; a predicate is only ever an optimisation. A
  // drawing marker now does its enumeration work, throws NeedsMainThread, and is
  // re-run whole on the main thread by post-pass 3. That wastes the work on a
  // marker rare enough to be 1 of 128,868 even at max_MAC_for_ER = 200, which is
  // a far better trade than serialising every ER marker to avoid it.
  //
  // Gated on harness/tsan.sh, both modes, WITH THIS DEFERRAL REMOVED -- which is
  // the first time a worker has ever run the ER machinery, and so the first time
  // SL_calloc and the Binary_* internals have been under a race detector.

  ProfileScope profTest(g_markerProfile.enabled ? &g_markerProfile.test : nullptr,
                        &g_markerProfile.nTest);

  if(MAC <= g_cfg.MACCutoffforER && t_traitType == "binary"){
    t_model.getMarkerPval(
                  t_scratch.GVec,
                        t_scratch.indexNonZeroVec_arma, t_scratch.indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA,
                        altFreq, Tstat, gy, varT, isSPAConverge, t_scratch.gtildeVec, is_gtilde, is_region, t_scratch.P2Vec, t_isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, t_scratch.G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, true, flip, ctx, &t_scratch.ws);
  }else{
    t_model.getMarkerPval(
                  t_scratch.GVec,
                        t_scratch.indexNonZeroVec_arma, t_scratch.indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA,
                        altFreq, Tstat, gy, varT, isSPAConverge, t_scratch.gtildeVec, is_gtilde, is_region, t_scratch.P2Vec, t_isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, t_scratch.G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, false, flip, ctx, &t_scratch.ws);
  }

  const double pval_num = parsePrintedPval(pval); // the one guarded parse (LEDGER #72)

  if((t_traitType == "binary" && MAC > g_cfg.MACCutoffforER) || t_traitType != "binary"){
    if(t_model.m_isFastTest && pval_num < (t_model.m_pval_cutoff_for_fastTest)){
     // The re-test escalates: a marker whose fast p-value cleared the cutoff is
     // re-tested with the covariate adjustment on, and on the sparse GRM unless
     // its MAC is above the last category bound. Same derivation as before, into
     // a second context rather than back onto the model.
     SAIGE::TestContext ctxRetest;
     ctxRetest.useSparseGRM = (MAC > t_model.m_cateVarRatioMinMACVecExclude.back()) ? false : t_model.m_flagSparseGRM;
     ctxRetest.noadjCov = false;
     if(!t_isSingleVarianceRatio){
        ctxRetest.varRatioVal = t_model.selectVarianceRatio(MAC, ctxRetest.useSparseGRM, ctxRetest.noadjCov, hasVarRatio);
     }else{
        ctxRetest.varRatioVal = t_model.selectSingleVarianceRatio(ctxRetest.useSparseGRM, ctxRetest.noadjCov);
        hasVarRatio = true;
     }

      t_model.getMarkerPval(
                    t_scratch.GVec,
                          t_scratch.indexNonZeroVec_arma, t_scratch.indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA,
                          altFreq, Tstat, gy, varT, isSPAConverge, t_scratch.gtildeVec, is_gtilde, is_region, t_scratch.P2Vec, t_isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, t_scratch.G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, false, flip, ctxRetest, &t_scratch.ws);
     }
  }//if((t_traitType == "binary" && MAC > g_cfg.MACCutoffforER) || t_traitType != "binary"){

  profTest.stop();

  ProfileScope profSummary(g_markerProfile.enabled ? &g_markerProfile.summary : nullptr,
                           &g_markerProfile.nSummary);

  // The two Firth counters were the loop's only cross-marker accumulators;
  // recorded per marker here and summed by the caller after the scan. Reading
  // is_FirthConverge stays nested inside is_Firth, which is what keeps it off
  // the uninitialised path (it is only written by the Firth fit itself).
  if(t_traitType == "binary"){
    if(is_Firth){
      t_res.isFirth = true;
      if(is_FirthConverge){
        t_res.isFirthConverge = true;
      }
    }
  }

  t_res.Beta = Beta * (1 - 2*flip);  // Beta if flip = false, -1*Beta is flip = true
  t_res.seBeta = seBeta;
  t_res.pval = pval;
  t_res.pvalNA = pval_noSPA;
  t_res.Tstat = Tstat * (1 - 2*flip);
  t_res.varT = varT;

  if(t_isCondition){
    t_res.Beta_c = Beta_c * (1 - 2*flip);  // Beta if flip = false, -1*Beta is flip = true
    t_res.seBeta_c = seBeta_c;
    t_res.pval_c = pval_c;
    t_res.pvalNA_c = pval_noSPA_c;
    t_res.Tstat_c = Tstat_c * (1 - 2*flip);
    t_res.varT_c = varT_c;
  }

  if(t_traitType == "binary"){
    arma::vec dosage_case = t_scratch.GVec.elem(t_model.m_case_indices);
    arma::vec dosage_ctrl = t_scratch.GVec.elem(t_model.m_ctrl_indices);
    AF_case = arma::mean(dosage_case) /2;
    AF_ctrl = arma::mean(dosage_ctrl) /2;
    N_case = dosage_case.n_elem;
    N_ctrl = dosage_ctrl.n_elem;
    if(flip){
       AF_case = 1-AF_case;
       AF_ctrl = 1-AF_ctrl;
    }
    t_res.isSPAConverge = isSPAConverge;
    t_res.AF_case = AF_case;
    t_res.AF_ctrl = AF_ctrl;

    t_res.N_case = N_case;
    t_res.N_ctrl = N_ctrl;

    arma::uvec N_case_ctrl_het_hom0;
    if(t_isMoreOutput){
      N_case_ctrl_het_hom0 = arma::find(dosage_case <= 2 && dosage_case >=1.5);
      t_res.N_case_hom = N_case_ctrl_het_hom0.n_elem;
      N_case_ctrl_het_hom0 = arma::find(dosage_case < 1.5 && dosage_case >= 0.5);
      t_res.N_case_het = N_case_ctrl_het_hom0.n_elem;
      N_case_ctrl_het_hom0 = arma::find(dosage_ctrl <= 2 && dosage_ctrl >=1.5);
      t_res.N_ctrl_hom = N_case_ctrl_het_hom0.n_elem;
      N_case_ctrl_het_hom0 = arma::find(dosage_ctrl < 1.5 && dosage_ctrl >= 0.5);
      t_res.N_ctrl_het = N_case_ctrl_het_hom0.n_elem;
      if(flip){
        t_res.N_case_hom = N_case - t_res.N_case_het - t_res.N_case_hom;
        t_res.N_ctrl_hom = N_ctrl - t_res.N_ctrl_het - t_res.N_ctrl_hom;
      }
    }
  }else if(t_traitType == "quantitative"){
    t_res.N = n;
  }
}


// One marker, one sample-set group, every model listed in t_modelIdxs.
//
// Returns false only on a short read, which is what the loop's `break` used to
// respond to -- so the caller stops the scan and the remaining slots keep their
// untested defaults, exactly as the q-length vectors did.
static bool testOneMarker(Geno::GenoReader & t_reader,
                          const std::string & t_genoType,
                          const std::string & t_genoIndexStr,
                          bool t_isMoreOutput,
                          bool t_isImputation,
                          const std::vector<unsigned char> & t_singleVR,
                          const std::vector<int> & t_modelIdxs,
                          Marker::MarkerScratch & t_scratch,
                          std::vector<Marker::MarkerResult*> & t_rows)
{
  MarkerPrep prep;
  // n is the group's sample count. Every model in a group has the same sampleID
  // vector by construction (that is what a group IS), so m_n of the first is
  // m_n of all of them -- and it is the same value the single-model loop read
  // from g_model->m_n.
  int n = g_models.at(t_modelIdxs.at(0))->m_n;
  int status = prepareOneMarker(t_reader, t_genoType, t_genoIndexStr, t_isImputation,
                                n, t_scratch, t_rows, prep);
  if(status < 0){
    return false;
  }
  if(status == 0){
    return true;   // QC dropped it; every row keeps pval == "NA"
  }

  for(std::size_t k = 0; k < t_modelIdxs.size(); k++){
    testOneModel(*g_models.at(t_modelIdxs[k]), t_isMoreOutput,
                 t_singleVR.at(t_modelIdxs[k]) != 0, prep, t_scratch, *t_rows[k]);
  }

  return true;
}


//////// ---------- Main function for marker-level analysis --------- ////////////
                           //std::vector<uint32_t> & t_genoIndex,

// S6: one worker's private state. The scratch is MarkerScratch's whole reason
// for existing (S4); the reader is an S5 clone, which is a shared_ptr copy onto
// the same open .bed rather than a second open, so making one is nearly free.
namespace {

// S6/S7: the happens-before edge that parallelFor's join really provides, made
// visible to ThreadSanitizer. TSAN cannot see it on its own: the joining code
// lives in libtbb, which is not instrumented, so every main-thread read of data
// a worker wrote (the writer's pass over MarkerResult, the pool's destructor)
// is reported as a race against a write that in fact strictly precedes it.
// Measured on the 128,868-marker fixture: 13-35 such reports per run, EVERY one
// of them main-thread-vs-worker and none worker-vs-worker, with the worker side
// holding the pool mutex -- and all of them gone once this edge is declared,
// while the deliberately-shared-scratch defect still produces 205-1130 races,
// worker-vs-worker, and a wrong answer, in every threaded cell. That pair is
// what makes a clean TSAN run mean something; see harness/tsan.sh, which runs
// both directions.
// Compiled out of every normal build -- the macro is defined only by
// harness/tsan.sh -- so this cannot affect a golden.
#ifdef SAIGE_TSAN_ANNOTATIONS
extern "C" void AnnotateHappensBefore(const char *f, int l, void *addr);
extern "C" void AnnotateHappensAfter(const char *f, int l, void *addr);
#define SAIGE_TSAN_HAPPENS_BEFORE(p) AnnotateHappensBefore(__FILE__, __LINE__, (void*)(p))
#define SAIGE_TSAN_HAPPENS_AFTER(p)  AnnotateHappensAfter(__FILE__, __LINE__, (void*)(p))
#else
#define SAIGE_TSAN_HAPPENS_BEFORE(p) ((void)0)
#define SAIGE_TSAN_HAPPENS_AFTER(p)  ((void)0)
#endif

// One worker's private state: a scratch buffer and a reader clone PER
// SAMPLE-SET GROUP. Per group, not per model: models in a group share a decode,
// so they share the buffer that decode lands in; models in different groups do
// not, so they cannot.
struct MarkerWorkerState {
  std::vector<std::unique_ptr<Marker::MarkerScratch>> scratch;
  std::vector<std::unique_ptr<Geno::GenoReader>> readers;
};

// A checkout pool rather than one state per range. RcppParallel's parallelFor
// hands out ranges by an auto partitioner, so the number of operator() calls is
// not the number of threads and a state cannot be indexed by range. Checking
// out costs one mutex acquisition per RANGE (a few hundred times per chunk),
// against allocating ~n doubles per range if each made its own.
class MarkerStatePool {
public:
  MarkerStatePool(const std::vector<std::vector<int>> & t_groups, int t_nStates)
    : m_groups(t_groups) {
    for(int i = 0; i < t_nStates; i++){
      m_free.push_back(make());
    }
  }

  std::unique_ptr<MarkerWorkerState> checkout() {
    {
      std::lock_guard<std::mutex> lock(m_mutex);
      if(!m_free.empty()){
        std::unique_ptr<MarkerWorkerState> s = std::move(m_free.back());
        m_free.pop_back();
        return s;
      }
    }
    // More concurrent ranges than states: should not happen when the pool is
    // sized to numThreads, but making one is cheap and correct.
    return make();
  }

  void checkin(std::unique_ptr<MarkerWorkerState> t_state) {
    std::lock_guard<std::mutex> lock(m_mutex);
    m_free.push_back(std::move(t_state));
  }

private:
  std::unique_ptr<MarkerWorkerState> make() const {
    std::unique_ptr<MarkerWorkerState> s(new MarkerWorkerState());
    for(std::size_t g = 0; g < m_groups.size(); g++){
      const SAIGE::SAIGEClass & first = *g_models.at(m_groups[g].at(0));
      s->scratch.push_back(std::unique_ptr<Marker::MarkerScratch>(
          new Marker::MarkerScratch(first.m_n, first.m_numMarker_cond)));
      s->readers.push_back(g_readers.at(g)->clone());
    }
    return s;
  }

  const std::vector<std::vector<int>> & m_groups;
  std::mutex m_mutex;
  std::vector<std::unique_ptr<MarkerWorkerState>> m_free;
};

// Sets the thread-local main-thread-only refusal for the duration of a worker
// range, so it cannot survive an early return or a throw.
struct RefuseMainThreadGuard {
  RefuseMainThreadGuard(){ SKATMainThread::setRefuse(true); }
  ~RefuseMainThreadGuard(){ SKATMainThread::setRefuse(false); }
};

struct MarkerScanWorker : public RcppParallel::Worker {
  MarkerScanWorker(MarkerStatePool& t_pool,
                   const std::string& t_genoType,
                   const std::vector<std::string>& t_genoIndex,
                   bool t_isMoreOutput, bool t_isImputation,
                   const std::vector<unsigned char>& t_singleVR,
                   const std::vector<std::vector<int>>& t_groups,
                   std::vector<std::vector<Marker::MarkerResult>>& t_results)
    : m_pool(t_pool), m_genoType(t_genoType),
      m_genoIndex(t_genoIndex), m_isMoreOutput(t_isMoreOutput),
      m_isImputation(t_isImputation), m_singleVR(t_singleVR),
      m_groups(t_groups), m_results(t_results) {}

  void operator()(std::size_t t_begin, std::size_t t_end) override {
    std::unique_ptr<MarkerWorkerState> state = m_pool.checkout();
    RefuseMainThreadGuard noRAPI;
    std::vector<Marker::MarkerResult*> rows;
    for(std::size_t i = t_begin; i < t_end; i++){
      // Each marker's failure is recorded in ITS OWN slot and the scan
      // continues, so which marker fails is not a function of which thread got
      // there first. The main thread then acts on the lowest index, which is
      // what a serial run would have stopped at. With several models, "its own
      // slot" means the row of EVERY model in the group that failed: a read is
      // per group, so the outcome of a read is too.
      for(std::size_t g = 0; g < m_groups.size(); g++){
        const std::vector<int> & members = m_groups[g];
        rows.clear();
        for(std::size_t k = 0; k < members.size(); k++){
          rows.push_back(&m_results.at(members[k]).at(i));
        }
        try {
          bool isReadMarker = testOneMarker(*state->readers[g], m_genoType,
                                            m_genoIndex.at(i), m_isMoreOutput,
                                            m_isImputation, m_singleVR, members,
                                            *state->scratch[g], rows);
          if(!isReadMarker){
            for(std::size_t k = 0; k < members.size(); k++){
              *rows[k] = Marker::MarkerResult();
              rows[k]->readFailed = true;
            }
          }
        } catch (const SKATMainThread::NeedsMainThread&) {
          // Only a MODEL's test can raise this (the ER resampling arm), never
          // the shared prepare -- but which model is not recorded, so the whole
          // group's rows are re-run on the main thread. Re-preparing for a
          // subset is byte-neutral because prepare is deterministic.
          for(std::size_t k = 0; k < members.size(); k++){
            *rows[k] = Marker::MarkerResult();
            rows[k]->deferredToMainThread = true;
          }
        } catch (const std::exception& e) {
          for(std::size_t k = 0; k < members.size(); k++){
            *rows[k] = Marker::MarkerResult();
            rows[k]->hasError = true;
            rows[k]->errMsg = e.what();
          }
        }
      }
    }
    m_pool.checkin(std::move(state));
    SAIGE_TSAN_HAPPENS_BEFORE(&m_pool);
  }

  MarkerStatePool& m_pool;
  const std::string& m_genoType;
  const std::vector<std::string>& m_genoIndex;
  bool m_isMoreOutput, m_isImputation;
  const std::vector<unsigned char>& m_singleVR;
  const std::vector<std::vector<int>>& m_groups;
  std::vector<std::vector<Marker::MarkerResult>>& m_results;
};

}  // namespace

// [[Rcpp::export]]
void mainMarkerInCPP(
                           std::string & t_genoType,     // "PLINK"
			   std::vector<std::string> & t_genoIndex,
			   bool & t_isMoreOutput,
			   bool & t_isImputation,
			   bool & t_isFirth,
			   int t_nThreadsMarker)
{

  int q = t_genoIndex.size();  // number of markers

  // --- the models, and the sample-set groups they fall into --------------
  const int P = (int)g_models.size();
  if(P < 1){
    Rcpp::stop("No null model has been set up.");
  }
  if((int)g_out.singles.size() != P){
    Rcpp::stop("There are " + std::to_string(P) + " null models but " +
               std::to_string(g_out.singles.size()) + " output files.");
  }
  if((int)g_modelGroup.size() != P){
    Rcpp::stop("There are " + std::to_string(P) + " null models but " +
               std::to_string(g_modelGroup.size()) + " sample-set group assignments.");
  }
  const int G = (int)g_readers.size();
  if(G < 1){
    Rcpp::stop("No genotype reader has been set up.");
  }
  std::vector<std::vector<int>> groups(G);
  for(int k = 0; k < P; k++){
    int g = g_modelGroup[k];
    if(g < 0 || g >= G){
      Rcpp::stop("Null model " + std::to_string(k + 1) + " names sample-set group " +
                 std::to_string(g + 1) + ", but only " + std::to_string(G) +
                 " genotype reader(s) were set up.");
    }
    groups[g].push_back(k);
  }
  for(int g = 0; g < G; g++){
    if(groups[g].empty()){
      Rcpp::stop("Sample-set group " + std::to_string(g + 1) +
                 " has a genotype reader but no null model.");
    }
  }

  // One output slot per (model, marker) (S4 for the marker axis, M1 for the
  // model axis). Model-major, so each model's rows are contiguous and the
  // writer streams one model's vector at a time exactly as it did when there
  // was one.
  std::vector<std::vector<Marker::MarkerResult>> results(
      P, std::vector<Marker::MarkerResult>(q));

  // S3 moved the marker path's variance ratio into the per-marker TestContext;
  // the assignSingleVarianceRatio call that used to stand here wrote a model
  // member nothing on this path reads any more, and it is gone with the member
  // (S12). All that survives is the single-vs-categorical decision itself --
  // per model, because two models can have been fitted differently.
  std::vector<unsigned char> singleVR(P);
  for(int k = 0; k < P; k++){
    singleVR[k] = ((g_models[k]->m_varRatio_null).n_elem == 1) ? 1 : 0;
  }

  int nThreads = (t_nThreadsMarker < 1) ? 1 : t_nThreadsMarker;

  // Phase 0: the profile is serial-arm only, so the decision is taken after
  // nThreads has been resolved (a reader that refuses to clone degrades to 1
  // below, which is still a serial scan and is still worth profiling).
  g_markerProfile.reset();
  g_markerProfile.enabled = (nThreads == 1) && markerProfileRequested();

  if(nThreads > 1){
    // A reader that cannot clone cannot be threaded. Since S8 that is one case
    // only: the PLINK pread fallback, which refuses (S5). PGEN clones -- each
    // worker gets its own pgenlib reader over shared metadata. Degrade to one
    // thread with a notice rather than failing the run -- the same posture as a
    // group test, and NOT the silent no-op of LEDGER #30: the reason is printed.
    // EVERY group's reader is probed: one uncloneable reader is enough to make
    // the scan serial, and finding that out mid-`parallelFor` is not an option.
    try {
      for(int g = 0; g < G; g++){
        std::unique_ptr<Geno::GenoReader> probe = g_readers[g]->clone();
      }
    } catch (const std::exception& e) {
      std::cout << "WARNING: the genotype reader cannot be cloned ("
                << e.what() << "); the marker scan will run on one thread."
                << std::endl;
      nThreads = 1;
    }
  }

  // Say what this scan is, but only when it is not the single-model scan every
  // existing golden was produced by -- P == 1 keeps its log verbatim. The
  // standing test asserts this line for the same reason the threading one does
  // (LEDGER #30): a manifest that failed to arrive would otherwise pass as a
  // green multi-model run.
  if(P > 1){
    std::cout << "Marker scan will test " << P << " models in " << G
              << " sample-set group(s) per marker." << std::endl;
  }

  // How far the scan got. Serially this is where the `break` fired; in the
  // threaded path the main thread derives it after the fact.
  int nScanned = q;

  if(nThreads == 1){
    // S6: the serial path is LITERALLY serial -- no worker, no pool, no TBB
    // touched at all -- so the default configuration cannot be affected by any
    // of the threading machinery, including its initialisation.
    std::vector<std::unique_ptr<Marker::MarkerScratch>> scratch;
    for(int g = 0; g < G; g++){
      const SAIGE::SAIGEClass & first = *g_models.at(groups[g].at(0));
      scratch.push_back(std::unique_ptr<Marker::MarkerScratch>(
          new Marker::MarkerScratch(first.m_n, first.m_numMarker_cond)));
    }
    std::vector<Marker::MarkerResult*> rows;
    for(int i = 0; i < q; i++){
      if((i+1) % g_cfg.marker_chunksize == 0){
        std::cout << "Completed " << (i+1) << "/" << q << " markers in the chunk." << std::endl;
      }

      bool isReadMarker = true;
      for(int g = 0; g < G && isReadMarker; g++){
        const std::vector<int> & members = groups[g];
        rows.clear();
        for(std::size_t k = 0; k < members.size(); k++){
          rows.push_back(&results[members[k]][i]);
        }
        isReadMarker = testOneMarker(*g_readers[g],
                                     t_genoType,
                                     t_genoIndex.at(i),
                                     t_isMoreOutput,
                                     t_isImputation,
                                     singleVR,
                                     members,
                                     *scratch[g],
                                     rows);
      }
      if(!isReadMarker){
        nScanned = i;
        break;
      }
    }
  }else{
    // Every marker's arithmetic depends only on its own genotype, the model
    // (const on this path since S3) and the config; its outputs land in its own
    // slot (S4); and its reader is its own handle onto shared, read-only file
    // state (S5 for PLINK, S8 for PGEN).
    // So the only thing order can change is which THREAD computes a row, not
    // what the row contains -- which is why byte-identity with the serial run
    // is the right gate rather than a tolerance.
    // Say so. Without this the threaded path is invisible from outside, and a
    // test that meant to exercise it but never got the knob through would pass
    // while gating nothing -- which is what --nThreads itself did for a year
    // (LEDGER #30). The standing test asserts this line, so a silent fallback
    // to serial cannot masquerade as a green threaded run.
    std::cout << "Marker scan will use " << nThreads
              << " in-process threads." << std::endl;
    MarkerStatePool pool(groups, nThreads);
    MarkerScanWorker worker(pool, t_genoType, t_genoIndex,
                            t_isMoreOutput, t_isImputation,
                            singleVR, groups, results);
    // Contiguous grain so each worker streams a run of the .bed rather than
    // stepping through it. numThreads is passed EXPLICITLY: RcppParallel's
    // resolveValue returns an explicit value whenever it is > 0 and not the
    // default, so it overrides RCPP_PARALLEL_NUM_THREADS -- which harness/env.sh
    // pins to 1, and which would otherwise silently make every threaded test
    // run serial and prove nothing.
    std::size_t grain = 64;
    RcppParallel::parallelFor(0, (std::size_t)q, worker, grain, nThreads);
    SAIGE_TSAN_HAPPENS_AFTER(&pool);

    // --- post-pass 1: report the error a serial run would have reported ---
    // Serial order is (marker, group, model), so that is the order to scan.
    for(int i = 0; i < q; i++){
      for(int g = 0; g < G; g++){
        for(std::size_t k = 0; k < groups[g].size(); k++){
          if(results[groups[g][k]][i].hasError){
            throw std::runtime_error(results[groups[g][k]][i].errMsg);
          }
        }
      }
    }

    // --- post-pass 2: truncate at the earliest short read ---
    // Serially a failed read breaks the loop, so every LATER marker is untested
    // FOR EVERY MODEL -- the loop is one loop however many models it carries.
    for(int i = 0; i < q && nScanned == q; i++){
      for(int k = 0; k < P; k++){
        if(results[k][i].readFailed){
          nScanned = i;
          break;
        }
      }
    }
    if(nScanned < q){
      for(int k = 0; k < P; k++){
        for(int j = nScanned; j < q; j++){
          results[k][j] = Marker::MarkerResult();
        }
      }
    }

    // --- post-pass 3: the ER markers, on the main thread ---
    // Ascending, though it does not matter for values: ComputeExact::Run
    // re-seeds with set.seed(1) at entry, so each drawing marker's sequence is
    // independent of how many drew before it (that is what makes the deferral
    // byte-neutral rather than merely safe).
    std::vector<std::unique_ptr<Marker::MarkerScratch>> scratch;
    for(int g = 0; g < G; g++){
      const SAIGE::SAIGEClass & first = *g_models.at(groups[g].at(0));
      scratch.push_back(std::unique_ptr<Marker::MarkerScratch>(
          new Marker::MarkerScratch(first.m_n, first.m_numMarker_cond)));
    }
    std::vector<Marker::MarkerResult*> rows;
    for(int i = 0; i < nScanned; i++){
      for(int g = 0; g < G; g++){
        const std::vector<int> & members = groups[g];
        std::vector<int> deferred;
        for(std::size_t k = 0; k < members.size(); k++){
          if(results[members[k]][i].deferredToMainThread){
            deferred.push_back(members[k]);
          }
        }
        if(deferred.empty()){
          continue;
        }
        rows.clear();
        for(std::size_t k = 0; k < deferred.size(); k++){
          results[deferred[k]][i] = Marker::MarkerResult();
          rows.push_back(&results[deferred[k]][i]);
        }
        bool isReadMarker = testOneMarker(*g_readers[g],
                                          t_genoType,
                                          t_genoIndex.at(i),
                                          t_isMoreOutput,
                                          t_isImputation,
                                          singleVR,
                                          deferred,
                                          *scratch[g],
                                          rows);
        if(!isReadMarker){
          nScanned = i;
          for(int k2 = 0; k2 < P; k2++){
            for(int j = i; j < q; j++){
              results[k2][j] = Marker::MarkerResult();
            }
          }
          break;
        }
      }
      if(nScanned <= i){
        break;
      }
    }

    // --- post-pass 4: the progress lines, in serial order ---
    // They go to stdout, not to any golden, but a threaded run should still
    // read like a serial one. Emitted here rather than from a worker, where the
    // interleaving would be nondeterministic.
    for(int i = 0; i < nScanned; i++){
      if((i+1) % g_cfg.marker_chunksize == 0){
        std::cout << "Completed " << (i+1) << "/" << q << " markers in the chunk." << std::endl;
      }
    }
  }

  if(g_markerProfile.enabled){
    double total = g_markerProfile.prepare + g_markerProfile.test + g_markerProfile.summary;
    std::cout << "MARKER PROFILE (serial, SAIGE_PROFILE_MARKER set):" << std::endl;
    std::cout << "  markers in chunk: " << q << ", scanned: " << nScanned << std::endl;
    std::cout << "  prepare (read+QC+impute): " << g_markerProfile.prepare
              << " s over " << g_markerProfile.nPrepare << " markers" << std::endl;
    std::cout << "  test    (score/SPA/ER)  : " << g_markerProfile.test
              << " s over " << g_markerProfile.nTest << " markers" << std::endl;
    std::cout << "  summary (row fields)   : " << g_markerProfile.summary
              << " s over " << g_markerProfile.nSummary << " markers" << std::endl;
    std::cout << "  accounted total        : " << total << " s" << std::endl;
    g_markerProfile.enabled = false;
  }

  //output, one file per model
  for(int k = 0; k < P; k++){
    int mFirth = 0;
    int mFirthConverge = 0;
    for(int i = 0; i < q; i++){
      if(results[k][i].isFirth){
        mFirth = mFirth + 1;
        if(results[k][i].isFirthConverge){
          mFirthConverge = mFirthConverge + 1;
        }
      }
    }

    writeOutfile_single(t_isMoreOutput,
        t_isImputation,
        g_models[k]->m_isCondition,
        t_isFirth,
        mFirth,
        mFirthConverge,
        g_models[k]->m_traitType,
        results[k],
        g_out.singles[k]->stream);
  }

}



// LEDGER #69 FIXED: the chromosome is consulted. processMale_XnonPAR decides
// from the POSITION alone, so with the flag on, every autosomal marker outside
// the two PAR windows was doubled for males too; the docstring told the user
// to supply an X-only file instead of the code checking. Only chromosome X
// qualifies now, under the spellings a variant file uses for it.
static bool isChromosomeX(const std::string& chr)
{
  std::string c = chr;
  if(c.size() > 3 && (c.compare(0, 3, "chr") == 0 || c.compare(0, 3, "Chr") == 0 || c.compare(0, 3, "CHR") == 0)) c = c.substr(3);
  return c == "X" || c == "x" || c == "23";
}

static void applyMaleXnonPAR(arma::vec& t_GVec, uint32_t& t_pd, double& t_altCounts, double& t_altFreq, const std::string& t_chr)
{
  if(g_cfg.is_rewrite_XnonPAR_forMales && isChromosomeX(t_chr)){
  	processMale_XnonPAR(t_GVec, t_pd, g_cfg.X_PARregion_mat);
	t_altCounts = arma::sum(t_GVec);
	t_altFreq = arma::mean(t_GVec)/2.0;	
  }
}


// a unified function to get single marker from genotype file


// [[Rcpp::export]]
bool Unified_getOneMarker(std::string & t_genoType,   // "PLINK"
                               uint64_t & t_gIndex,        // 0-based row of the variant metadata file
                               std::string& t_ref,       // REF allele
                               std::string& t_alt,       // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                               std::string& t_marker,    // marker ID extracted from genotype file
                               uint32_t& t_pd,           // base position
                               std::string& t_chr,       // chromosome
                               double& t_altFreq,        // frequency of ALT allele
                               double& t_altCounts,      // counts of ALT allele
                               double& t_missingRate,    // missing rate
                               double& t_imputeInfo,     // imputation information score, i.e., R2 (all 1 for PLINK)
                               bool & t_isOutputIndexForMissing,               // if true, output index of missing genotype data
                               std::vector<uint>& t_indexForMissing,     // index of missing genotype data
                               bool & t_isOnlyOutputNonZero,                   // if true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                               std::vector<uint>& t_indexForNonZero, //
			       arma::vec & t_GVec,
			       bool t_isImputation 
			       )     // the index of non-zero genotype in the all subjects. Only valid if t_isOnlyOutputNonZero == true.
{
  if(!g_reader){
    Rcpp::stop("No genotype reader has been set up.");
  }
  return getOneMarkerThrough(*g_reader, t_genoType, t_gIndex, t_ref, t_alt, t_marker,
                             t_pd, t_chr, t_altFreq, t_altCounts, t_missingRate,
                             t_imputeInfo, t_isOutputIndexForMissing, t_indexForMissing,
                             t_isOnlyOutputNonZero, t_indexForNonZero, t_GVec,
                             t_isImputation);
}


// [[Rcpp::export]]
void processMale_XnonPAR(arma::vec & t_GVec,  uint32_t& t_pd , arma::umat & t_XPARregion){
	bool isPAR = false;
	for(unsigned int j = 0; j < t_XPARregion.n_rows; j++){
		if(t_pd <= t_XPARregion(j,1) && t_pd >= t_XPARregion(j,0)){	
			isPAR = true;
		}
	}
	if(!isPAR){
		t_GVec.elem(g_cfg.indexInModel_male) = t_GVec.elem(g_cfg.indexInModel_male) * 2;
	}
}


// [[Rcpp::export]]
uint32_t Unified_getSampleSizeinGeno(std::string & t_genoType){
    uint32_t N0;
    if(t_genoType == "plink" || t_genoType == "pgen"){
	N0 = g_reader->getN0();
    }
    return(N0);
}	

// [[Rcpp::export]]
uint32_t Unified_getSampleSizeinAnalysis(std::string & t_genoType){
    uint32_t N;
    if(t_genoType == "plink" || t_genoType == "pgen"){
        N = g_reader->getN();
    }
    return(N);
}




// a unified function to get marker-level p-value
void Unified_getMarkerPval(
                           arma::vec & t_GVec,
                           bool t_isOnlyOutputNonZero,
			   arma::uvec & t_indexForNonZero_vec,
			   arma::uvec & t_indexForZero_vec,
                           double& t_Beta, 
                           double& t_seBeta, 
                           std::string& t_pval,
			   std::string& t_pval_noSPA,
                           double& t_Tstat,
			   double& t_gy,
                           double& t_varT,
                           double t_altFreq, 
			   bool & t_isSPAConverge,
			   arma::vec & t_gtilde,
			   bool & is_gtilde,
			   bool  is_region,
                           arma::vec & t_P2Vec,
			   bool t_isCondition,
                           double& t_Beta_c, 
                           double& t_seBeta_c, 
                           std::string& t_pval_c,
			   std::string& t_pval_noSPA_c,
                           double& t_Tstat_c,
                           double& t_varT_c,
			   arma::rowvec & t_G1tilde_P_G2tilde_Vec, 
			   bool & t_isFirth,
			   bool & t_isFirthConverge, 
			   bool t_isER,
			   bool t_flip,
			   const SAIGE::TestContext & t_ctx )
{
    if(t_isOnlyOutputNonZero == true)
      Rcpp::stop("When using SAIGE method to calculate marker-level p-values, 't_isOnlyOutputNonZero' should be false.");


    g_model->getMarkerPval(t_GVec, t_indexForNonZero_vec, t_indexForZero_vec, t_Beta, t_seBeta, t_pval, t_pval_noSPA, t_altFreq, t_Tstat, t_gy, t_varT, t_isSPAConverge, t_gtilde, is_gtilde, is_region, t_P2Vec, t_isCondition, t_Beta_c, t_seBeta_c, t_pval_c, t_pval_noSPA_c, t_Tstat_c, t_varT_c, t_G1tilde_P_G2tilde_Vec, t_isFirth, t_isFirthConverge, t_isER, t_flip, t_ctx); //SAIGE_new.cpp
    
    //t_indexForNonZero_vec.clear();
  
}

//////// ---------- Main functions to set objects for different genotype format --------- ////////////

// [[Rcpp::export]]
void setPLINKobjInCPP(std::string t_bimFile,
                      std::string t_famFile,
                      std::string t_bedFile,
                      std::vector<std::string> & t_SampleInModel,
                      std::string t_AlleleOrder)
{
  Rcpp::List spec = Rcpp::List::create(Rcpp::Named("bimFile") = t_bimFile,
                                       Rcpp::Named("famFile") = t_famFile,
                                       Rcpp::Named("bedFile") = t_bedFile,
                                       Rcpp::Named("AlleleOrder") = t_AlleleOrder);
  resetReaders(Geno::makeReader("plink", spec, t_SampleInModel));
}

// [[Rcpp::export]]
void setPGENobjInCPP(std::string t_pgenFile,
                     std::string t_pvarFile,
                     std::string t_psamFile,
                     std::vector<std::string> & t_SampleInModel)
{
  Rcpp::List spec = Rcpp::List::create(Rcpp::Named("pgenFile") = t_pgenFile,
                                       Rcpp::Named("pvarFile") = t_pvarFile,
                                       Rcpp::Named("psamFile") = t_psamFile);
  resetReaders(Geno::makeReader("pgen", spec, t_SampleInModel));
}

// A second (or later) SAMPLE-SET GROUP's reader. The first group's reader comes
// from setPLINKobjInCPP/setPGENobjInCPP as it always did; this only appends, so
// a single-model run never reaches it. The spec is the same Rcpp::List the two
// format-specific setters build, so nothing here knows about bim/pvar columns.
// [[Rcpp::export]]
void addGenoReaderInCPP(std::string t_genoType,
                        Rcpp::List t_spec,
                        std::vector<std::string> & t_SampleInModel)
{
  addReader(Geno::makeReader(t_genoType, t_spec, t_SampleInModel));
}

//////// ---------- Main functions to set objects for different analysis methods --------- ////////////

// Has setSAIGEobjInCPP been called in this process? addSAIGEobjInCPP appends to
// a list that setSAIGEobjInCPP is responsible for clearing, so calling them out
// of order would silently carry a previous run's models into this one.
static bool g_modelSetCalled = false;

// forward declaration: setSAIGEobjInCPP below forwards to it
void addSAIGEobjInCPP(arma::mat & t_XVX,
        arma::mat & t_XXVX_inv,
        arma::mat & t_XV,
        arma::mat & t_XVX_inv_XV,
	arma::mat & t_Sigma_iXXSigma_iX,
        arma::mat & t_X,
        arma::vec &  t_S_a,
        arma::vec & t_res,
        arma::vec & t_mu2,
        arma::vec & t_mu,
        arma::vec & t_varRatio_sparse,
        arma::vec & t_varRatio_null,
        arma::vec & t_varRatio_null_noXadj,
	arma::vec & t_cateVarRatioMinMACVecExclude,
        arma::vec & t_cateVarRatioMaxMACVecInclude,
        double t_SPA_Cutoff,
        arma::vec & t_tauvec,
        std::string t_traitType,
        arma::vec & t_y,
	std::string t_impute_method,
	bool t_flagSparseGRM,
	bool t_isFastTest,
	bool t_isnoadjCov,
	double t_pval_cutoff_for_fastTest,
        arma::umat & t_locationMat,
        arma::vec & t_valueVec,
        int t_dimNum,
	bool t_isCondition,
	std::vector<uint32_t> & t_condition_genoIndex,
	bool t_is_Firth_beta,
	double t_pCutoffforFirth,
	std::string t_firthSolver,
	arma::vec & t_offset,
	arma::vec & t_resout,
	int t_groupIndex);

// The FIRST null model of a run. Clears whatever a previous SPAGMMATtest call
// in this process left behind and installs model 0 in sample-set group 0;
// addSAIGEobjInCPP appends the rest. Splitting it this way is what keeps a
// single-model run's call sequence exactly what it was before M1.
// [[Rcpp::export]]
void setSAIGEobjInCPP(arma::mat & t_XVX,
        arma::mat & t_XXVX_inv,
        arma::mat & t_XV,
        arma::mat & t_XVX_inv_XV,
	arma::mat & t_Sigma_iXXSigma_iX,
        arma::mat & t_X,
        arma::vec &  t_S_a,
        arma::vec & t_res,
        arma::vec & t_mu2,
        arma::vec & t_mu,
        arma::vec & t_varRatio_sparse,
        arma::vec & t_varRatio_null,
        arma::vec & t_varRatio_null_noXadj,
	arma::vec & t_cateVarRatioMinMACVecExclude,
        arma::vec & t_cateVarRatioMaxMACVecInclude,
        double t_SPA_Cutoff,
        arma::vec & t_tauvec,
        std::string t_traitType,
        arma::vec & t_y,
	std::string t_impute_method,
	bool t_flagSparseGRM,
	bool t_isFastTest,
	bool t_isnoadjCov,
	double t_pval_cutoff_for_fastTest,
        arma::umat & t_locationMat,
        arma::vec & t_valueVec,
        int t_dimNum,
	bool t_isCondition,
	std::vector<uint32_t> & t_condition_genoIndex,
	bool t_is_Firth_beta,
	double t_pCutoffforFirth,
	std::string t_firthSolver,
	arma::vec & t_offset,
	arma::vec & t_resout)
{
  g_models.clear();
  g_modelGroup.clear();
  g_model = nullptr;
  g_modelSetCalled = true;
  addSAIGEobjInCPP(
        t_XVX,
        t_XXVX_inv,
        t_XV,
        t_XVX_inv_XV,
        t_Sigma_iXXSigma_iX,
        t_X,
        t_S_a,
        t_res,
        t_mu2,
        t_mu,
        t_varRatio_sparse,
        t_varRatio_null,
        t_varRatio_null_noXadj,
        t_cateVarRatioMinMACVecExclude,
        t_cateVarRatioMaxMACVecInclude,
        t_SPA_Cutoff,
        t_tauvec,
        t_traitType,
        t_y,
        t_impute_method,
        t_flagSparseGRM,
        t_isFastTest,
        t_isnoadjCov,
        t_pval_cutoff_for_fastTest,
        t_locationMat,
        t_valueVec,
        t_dimNum,
        t_isCondition,
        t_condition_genoIndex,
        t_is_Firth_beta,
        t_pCutoffforFirth,
        t_firthSolver,
        t_offset,
        t_resout,
        0);
}

// The SECOND and later null models, each with the index of the sample-set group
// whose reader it shares. A single-model run never calls this.
// [[Rcpp::export]]
void addSAIGEobjInCPP(arma::mat & t_XVX,
        arma::mat & t_XXVX_inv,
        arma::mat & t_XV,
        arma::mat & t_XVX_inv_XV,
	arma::mat & t_Sigma_iXXSigma_iX,
        arma::mat & t_X,
        arma::vec &  t_S_a,
        arma::vec & t_res,
        arma::vec & t_mu2,
        arma::vec & t_mu,
        arma::vec & t_varRatio_sparse,
        arma::vec & t_varRatio_null,
        arma::vec & t_varRatio_null_noXadj,
	arma::vec & t_cateVarRatioMinMACVecExclude,
        arma::vec & t_cateVarRatioMaxMACVecInclude,
        double t_SPA_Cutoff,
        arma::vec & t_tauvec,
        std::string t_traitType,
        arma::vec & t_y,
	std::string t_impute_method,
	bool t_flagSparseGRM,
	bool t_isFastTest,
	bool t_isnoadjCov,
	double t_pval_cutoff_for_fastTest,
        arma::umat & t_locationMat,
        arma::vec & t_valueVec,
        int t_dimNum,
	bool t_isCondition,
	std::vector<uint32_t> & t_condition_genoIndex,
	bool t_is_Firth_beta,
	double t_pCutoffforFirth,
	std::string t_firthSolver,
	arma::vec & t_offset,
	arma::vec & t_resout,
	int t_groupIndex)
{
  if(!g_modelSetCalled){
    Rcpp::stop("setSAIGEobjInCPP must be called for the first null model before addSAIGEobjInCPP.");
  }
  if(t_groupIndex < 0){
    Rcpp::stop("A null model's sample-set group index must be non-negative.");
  }
  // check SAIGE.cpp
  g_models.push_back(std::unique_ptr<SAIGE::SAIGEClass>(new SAIGE::SAIGEClass(
	t_XVX,
        t_XXVX_inv,
        t_XV,
        t_XVX_inv_XV,
	t_Sigma_iXXSigma_iX,
        t_X,
        t_S_a,
        t_res,
        t_mu2,
        t_mu,
	t_varRatio_sparse,
        t_varRatio_null,
	t_varRatio_null_noXadj,
	t_cateVarRatioMinMACVecExclude,
	t_cateVarRatioMaxMACVecInclude,
        t_SPA_Cutoff,
        t_tauvec,
        t_traitType,
        t_y,
	t_impute_method,
	t_flagSparseGRM,
	t_isFastTest,
	t_isnoadjCov,
	t_pval_cutoff_for_fastTest,
	t_locationMat,
	t_valueVec,
	t_dimNum, 
	t_isCondition,
	t_condition_genoIndex,
	t_is_Firth_beta,
        t_pCutoffforFirth,
        t_firthSolver,
	t_offset,
	t_resout)));
  g_modelGroup.push_back(t_groupIndex);
  g_model = g_models[0].get();
}


/*
// [[Rcpp::export]]
void setSparseSigmaInCPP(int r, arma::umat & t_locationMatinR, arma::vec & t_valueVecinR)
{
  g_model->setupSparseMat(r, t_locationMatinR, t_valueVecinR);
  g_model->m_flagSparseGRM = true;
}
*/

// [[Rcpp::export]]
Rcpp::List RegionSetUpConditional_binary_InCPP(arma::vec & t_weight_cond){

	unsigned int q_cond = (g_model->m_VarInvMat_cond).n_rows;
  	boost::math::beta_distribution<> beta_dist(g_cfg.weights_beta[0], g_cfg.weights_beta[1]);
  	arma::vec w0G2Vec_cond(q_cond);
  	double w0G2_cond, MAFG2_cond;
        for(unsigned int ci = 0; ci < q_cond; ci++){
		if(!(t_weight_cond.is_zero())){
			w0G2_cond = t_weight_cond(ci);
		}else{
                	MAFG2_cond = (g_model->m_MAF_cond)[ci];
                	w0G2_cond = boost::math::pdf(beta_dist, MAFG2_cond);
		}
                w0G2Vec_cond.at(ci) = w0G2_cond;
        }
	arma::mat m_VarMat_weighted_cond = (w0G2Vec_cond * (w0G2Vec_cond.t())) % (g_model->m_VarMat_cond);

	 Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("VarMat_G2_cond") = m_VarMat_weighted_cond,
                                          Rcpp::Named("Score_G2_cond") = g_model->m_Tstat_cond,
                                          Rcpp::Named("pval_G2_cond") = g_model->m_p_cond,
                                          Rcpp::Named("gsum_G2_cond") = g_model->m_gsum_cond,
                                          Rcpp::Named("qsum_G2_cond") = g_model->m_qsum_cond
					  );
	 return(OutList);

}



// LEDGER #57: the chunked VarMat assembly writes P1Mat/P2Mat to temporary .bin
// files and reads them back. arma::Mat::save()/.load() both RETURN a bool and
// every return value here used to be discarded; Armadillo does not cover for
// that either, because save/load failures warn at arma_warn(3, ...) and
// ARMA_WARN_LEVEL defaults to 2, so a failure was silent by construction. A
// failed load left the matrix empty, the assembly skipped that block, SKAT then
// threw for every cell of the region and get_SKAT_pvalue's try(silent = TRUE)
// turned the whole gene into NA rows -- exit 0, no warning. Observed once in
// seven N = 50,000 scans, with the volume at 100%.
//
// These stop instead. That is byte-neutral wherever the I/O succeeds, which is
// every path any golden takes (the bundled fixtures never chunk).
static void saveChunkMatOrStop(const arma::mat & t_mat, const std::string & t_file,
                               const std::string & t_region)
{
  if(!t_mat.save(t_file)){
    Rcpp::stop("Failed to WRITE the temporary chunk file '" + t_file + "' while testing region '"
               + t_region + "'. The set/group test streams the variance-matrix chunks through "
               "this file; continuing would silently drop or misplace part of the region's "
               "variance matrix. Check free space and permissions on the output directory.");
  }
}

static void loadChunkMatOrStop(arma::mat & t_mat, const std::string & t_file,
                               const std::string & t_region)
{
  if(!t_mat.load(t_file)){
    Rcpp::stop("Failed to READ the temporary chunk file '" + t_file + "' while testing region '"
               + t_region + "'. It was written earlier in this run, so it has been removed, "
               "truncated or made unreadable since. Continuing would silently drop or misplace "
               "part of the region's variance matrix. Check free space and permissions on the "
               "output directory, and that no other SAIGE run shares this --SAIGEOutputFile.");
  }
}

// ---- step 2a: the chunk store (in-memory by default, spills above a budget) ----
// Replaces the per-region round trip through <prefix>_P{1,2}Mat_Chunk_<i>.bin.
// The chunk STRUCTURE is unchanged -- one (P1 rows, P2 cols) pair per chunk,
// the ultra-rare block last -- so VarMat is assembled from exactly the same
// sub-block products with exactly the same operands, which is what keeps every
// group golden byte-identical: the .bin files held raw doubles, so a chunk read
// back from disk and a chunk kept in memory are the same matrix.
//
// Budget: SAIGE_REGION_CHUNK_MEMORY_MB (default 8192). Chunks are kept in
// memory while the running total of their bytes stays within it; a chunk that
// would exceed it goes to disk under the old file name and is loaded on
// demand, exactly as before. The tests force the spill with a budget of 0.
namespace {

class RegionChunkStore {
public:
  RegionChunkStore(const std::string & t_prefix, const std::string & t_region)
    : m_prefix(t_prefix), m_region(t_region)
  {
    const char* e = std::getenv("SAIGE_REGION_CHUNK_MEMORY_MB");
    double mb = 8192.0;
    if(e != nullptr && e[0] != '\0'){
      char* end = nullptr;
      double v = std::strtod(e, &end);
      if(end != e && v >= 0.0){ mb = v; }
    }
    m_budgetBytes = mb * 1024.0 * 1024.0;
  }
  ~RegionChunkStore(){ clear(); }

  // Store chunk `ichunk` (called once per chunk, in order).
  // By value: a `P1Mat.rows(0, k)` subview argument is materialised exactly
  // once, into the matrix that is then moved into the slot.
  void add(unsigned int t_ichunk, arma::mat t_P1, arma::mat t_P2)
  {
    double bytes = 8.0 * (double(t_P1.n_elem) + double(t_P2.n_elem));
    Slot s;
    if(m_usedBytes + bytes <= m_budgetBytes){
      s.p1 = std::move(t_P1); s.p2 = std::move(t_P2); s.inMemory = true;
      m_usedBytes += bytes;
    }else{
      saveChunkMatOrStop(t_P1, p1File(t_ichunk), m_region);
      saveChunkMatOrStop(t_P2, p2File(t_ichunk), m_region);
      s.inMemory = false;
      m_spilled = true;
    }
    if(m_slots.size() <= t_ichunk){ m_slots.resize(t_ichunk + 1); }
    m_slots[t_ichunk] = std::move(s);
  }
  // The chunk's P1 rows / P2 cols; a spilled chunk is (re)loaded into scratch.
  const arma::mat & P1(unsigned int i, arma::mat & t_scratch)
  {
    if(m_slots.at(i).inMemory){ return m_slots[i].p1; }
    loadChunkMatOrStop(t_scratch, p1File(i), m_region);
    return t_scratch;
  }
  const arma::mat & P2(unsigned int i, arma::mat & t_scratch)
  {
    if(m_slots.at(i).inMemory){ return m_slots[i].p2; }
    loadChunkMatOrStop(t_scratch, p2File(i), m_region);
    return t_scratch;
  }
  bool spilled() const { return m_spilled; }
  void clear()
  {
    for(unsigned int i = 0; i < m_slots.size(); i++){
      if(!m_slots[i].inMemory){
        std::remove(p1File(i).c_str());
        std::remove(p2File(i).c_str());
      }
    }
    m_slots.clear();
    m_usedBytes = 0.0;
  }
private:
  struct Slot { arma::mat p1, p2; bool inMemory = true; };
  std::string p1File(unsigned int i) const { return m_prefix + "_P1Mat_Chunk_" + std::to_string(i) + ".bin"; }
  std::string p2File(unsigned int i) const { return m_prefix + "_P2Mat_Chunk_" + std::to_string(i) + ".bin"; }
  std::string m_prefix, m_region;
  double m_budgetBytes = 0.0, m_usedBytes = 0.0;
  bool m_spilled = false;
  std::vector<Slot> m_slots;
};

}  // namespace

//////// ---------- Main function for region-level analysis --------- ////////////
// [[Rcpp::export]]
Rcpp::List mainRegionInCPP(
                           std::string t_genoType,     // "PLINK"
                           std::vector<std::string> & t_genoIndex,
			   arma::mat & annoIndicatorMat,
			   arma::vec & maxMAFVec, 
                           std::string t_outputFile,
			   std::string t_traitType,
                           unsigned int t_n,           // sample size  
                           arma::mat & t_P1Mat_placeholder,  // 1 x 1 since 2026-09-06; the scratch is s_regionP1 below
                           arma::mat & t_P2Mat_placeholder, 
			   std::string t_regionTestType, 
			   bool t_isImputation,
			   arma::mat & t_weight,       // q0 x Kc custom (group-file) weight sets; 0-matrix sentinel when none
			   arma::mat & t_beta_params,  // Kb x 2 beta(a,b) parameter pairs; 0-matrix sentinel when none
			   bool t_includeNoWeights,    // append an all-ones "No_Weights" set
			   std::vector<std::string> & t_weightNames, // per-set labels; only consumed when the Weight column is emitted
			   bool t_forceWeightColumn, // emit the Weight column in the BURDEN writer even at q_weight == 1 (is_output_weight_column)
			   arma::vec & t_weight_cond,
			   bool t_isSingleinGroupTest,
			   bool t_isOutputMarkerList, 
			   std::vector<std::string> & annoStringVec,
			   std::string regionName, 
			   bool t_isFastTest, 
			   bool t_isMoreOutput,
			   int t_nThreads)
{


  //create the output list
  Rcpp::List OutList = Rcpp::List::create();
  //arma::vec timeoutput1 = getTime();
  unsigned int q0 = t_genoIndex.size();                 // number of markers (before QC) in one region

  // The caller's intent for this pass: false on the cheap fastTest pass, true on
  // the accurate re-run. Both write decisions below read THIS and nothing else.
  //
  // It used to be a snapshot of g_model->m_flagSparseGRM_cur, taken here because
  // the per-marker variance-ratio selection in the loops below overwrote that
  // member from MAC and destroyed the signal (LEDGER #27, and #50 in the
  // non-BURDEN branch, which was still reading the mutated member). Since S12
  // the two jobs are separate variables: this one is written only by R between
  // region calls, and each marker's variance-ratio choice lives in a local
  // TestContext that outlives nothing.
  bool isEscalatedPass = g_regionEscalatedPass;

  // Serial arm only, like the marker profile: the buckets are plain doubles
  // that phase-1 workers would race on.
  g_regionProfile.enabled = regionProfileRequested() && (t_nThreads <= 1);
  RegionProfile & rp = g_regionProfile;
  ProfileScope rpTotal(rp.enabled ? &rp.total : nullptr, &rp.nTotal);

  // ---- weight-set layout (jmr extension) ----------------------------------
  // Weight sets are ordered: custom (group-file) sets first, then beta(a,b)
  // sets, then the optional all-ones "No_Weights" set. When the caller
  // supplies none of these (t_weight and t_beta_params are 0-matrix
  // sentinels, t_includeNoWeights false), a single IMPLICIT beta set
  // parameterized by g_cfg.weights_beta is used -- exactly the legacy
  // single-weight path.
  // Per-marker weight by set kind:
  //   custom set  : w0 = t_weight(i, col);  UR collapse max(w0*G)
  //   beta set    : w0 = pdf(beta, MAF);    UR collapse max(G), then scaled
  //                                         by pdf(beta, MAF_collapsed)
  //   No_Weights  : w0 = 1;                 UR collapse max(G)
  // NORMATIVE flat index over (anno j, maxMAF m, weight set r):
  //   jmr = j*q_maf*q_weight + m*q_weight + r      (weight fastest)
  // which reduces to the legacy jm = j*q_maf + m when q_weight == 1. The
  // *_noweights accumulators stay jm-indexed (weight-independent). The same
  // formula is documented at the R call site (run_mainRegionInCPP).
  boost::math::beta_distribution<> beta_dist(g_cfg.weights_beta[0], g_cfg.weights_beta[1]);
  bool hasCustomWeights = (!(t_weight.is_zero()) && t_weight.n_rows == q0);
  unsigned int q_weight_custom = hasCustomWeights ? t_weight.n_cols : 0;
  bool hasExplicitBeta = !(t_beta_params.is_zero());
  if(hasExplicitBeta && !(t_beta_params.n_cols == 2 && t_beta_params.n_rows >= 1)){
    // fail fast on a malformed (e.g. transposed Kb x 2 -> 2 x Kb) matrix
    // rather than silently degrading to "no explicit beta sets"; a
    // transposed 2 x 2 is shape-indistinguishable and cannot be caught here
    Rcpp::stop("mainRegionInCPP: t_beta_params must be the 0-matrix sentinel or a Kb x 2 matrix of beta(a,b) parameter pairs; got %d x %d", (int)(t_beta_params.n_rows), (int)(t_beta_params.n_cols));
  }
  unsigned int q_weight_beta = hasExplicitBeta ? t_beta_params.n_rows : 0;
  bool useImplicitBeta = (q_weight_custom + q_weight_beta + (t_includeNoWeights ? 1u : 0u)) == 0;
  if(useImplicitBeta){
    q_weight_beta = 1;
  }
  unsigned int q_weight = q_weight_custom + q_weight_beta + (t_includeNoWeights ? 1u : 0u);
  // Weight-column emission in the BURDEN writer: on for a real multi-weight
  // run, or forced at a single set by is_output_weight_column (schema only;
  // no arithmetic depends on this flag)
  bool emitWeightColumn = (q_weight > 1) || t_forceWeightColumn;
  const int SET_CUSTOM = 0, SET_BETA = 1, SET_NOWEIGHT = 2;
  std::vector<int> setKindVec(q_weight);
  std::vector<unsigned int> setColVec(q_weight, 0);  // custom column / beta row per set
  std::vector<boost::math::beta_distribution<>> betaDistVec;
  // per-set ER exact-test gate: the legacy (!isWeightCustomized ||
  // isEqualWeights) predicate evaluated per set (true for beta and
  // No_Weights sets, all-ones test for custom sets)
  std::vector<bool> erAllowedVec(q_weight, true);
  {
    unsigned int rset = 0;
    for(unsigned int c = 0; c < q_weight_custom; c++){
      setKindVec[rset] = SET_CUSTOM;
      setColVec[rset] = c;
      erAllowedVec[rset] = arma::all(t_weight.col(c) == 1.0);
      rset++;
    }
    for(unsigned int b = 0; b < q_weight_beta; b++){
      setKindVec[rset] = SET_BETA;
      setColVec[rset] = b;
      if(useImplicitBeta){
        betaDistVec.push_back(beta_dist);
      }else{
        betaDistVec.push_back(boost::math::beta_distribution<>(t_beta_params(b,0), t_beta_params(b,1)));
      }
      rset++;
    }
    if(t_includeNoWeights){
      setKindVec[rset] = SET_NOWEIGHT;
    }
  }
  arma::vec w0_vec(q_weight);   // per-set weight of the marker being processed
  // -------------------------------------------------------------------------

  unsigned int q_anno = annoIndicatorMat.n_cols;
  unsigned int q_maf = maxMAFVec.n_elem;
  unsigned int q_anno_maf = q_anno*q_maf;                    // weight-independent (jm) slots
  unsigned int q_anno_maf_weight = q_anno_maf*q_weight;      // per-set (jmr) slots
  arma::mat genoURMat(t_n, q_anno_maf_weight, arma::fill::zeros);
  arma::mat genoURMat_noweights(t_n, q_anno_maf, arma::fill::zeros);  // jm-indexed (weight-independent)
  unsigned int q = q0 + q_anno_maf_weight;
  arma::imat annoMAFIndicatorMat(q, q_anno_maf_weight, arma::fill::zeros);
  arma::ivec annoMAFIndicatorVec(q_anno_maf_weight);
  annoMAFIndicatorVec.zeros();
  arma::vec maxMAFperAnno(q_anno, arma::fill::zeros);
  arma::vec MAFIndicatorVec(maxMAFVec.n_elem);
  MAFIndicatorVec.zeros();


  //setting up conditional markers
  unsigned int q_cond = (g_model->m_VarInvMat_cond).n_rows;
  arma::rowvec G1tilde_P_G2tilde_Vec(q_cond);

  bool isCondition = g_model->m_isCondition;
  arma::vec w0G2Vec_cond(q_cond);
  double w0G2_cond, MAFG2_cond;
  if(isCondition){
	for(unsigned int ci = 0; ci < q_cond; ci++){
		if(!(t_weight_cond.is_zero())){
			w0G2_cond = t_weight_cond(ci);
		}else{	
			MAFG2_cond = (g_model->m_MAF_cond)[ci];	
  			w0G2_cond = boost::math::pdf(beta_dist, MAFG2_cond);
		}	
		w0G2Vec_cond.at(ci) = w0G2_cond;
	}
  }
  arma::mat w0G2Mat_cond(q_cond, q_cond);
  w0G2Mat_cond = w0G2Vec_cond * (w0G2Vec_cond.t());
  arma::mat genoSumMat(t_n, q_anno_maf_weight, arma::fill::zeros); //for Phi_cc for binary traits and BURDEN test
  arma::vec genoSumcount_noweight(q_anno_maf, arma::fill::zeros);  // jm-indexed (weight-independent)
  //arma::sp_mat genoSumMat_sp(t_n, q_anno_maf); //for Phi_cc for binary traits and BURDEN test
  std::vector<double> Beta_cVec(q, arma::datum::nan);         // beta value for ALT allele
  std::vector<double> seBeta_cVec(q, arma::datum::nan);
  //std::vector<double> pval_cVec(q, arma::datum::nan);
  std::vector<std::string> pval_cVec(q, "NA");
  std::vector<double> Tstat_cVec(q, arma::datum::nan);
  std::vector<double> varT_cVec(q, arma::datum::nan);
  //std::vector<double> pvalNA_cVec(q, arma::datum::nan);
  std::vector<std::string> pvalNA_cVec(q, "NA");
  arma::mat G1tilde_P_G2tilde_Weighted_Mat(q, q_cond);
  //group test output
  //arma::vec MAC_GroupVec = arma::zeros<vec>(q_anno_maf);
  arma::vec MAC_GroupVec(q_anno_maf_weight);
  MAC_GroupVec.zeros();
  arma::vec MACCase_GroupVec(q_anno_maf_weight);
  MACCase_GroupVec.zeros();
  arma::vec MACControl_GroupVec(q_anno_maf_weight);
  MACControl_GroupVec.zeros();
  arma::vec NumRare_GroupVec(q_anno_maf_weight);
  NumRare_GroupVec.zeros();
  arma::vec NumUltraRare_GroupVec(q_anno_maf_weight);
  NumUltraRare_GroupVec.zeros();
  arma::vec gtildeVec;
  double MACgroup, MACcasegroup, MACcontrolgroup, AF_case, AF_ctrl, N_case_hom, N_ctrl_het, N_case_het, N_ctrl_hom;
  uint32_t N_case, N_ctrl, N; 

  //single-variant assoc output
  arma::uvec indicatorVec(q, arma::fill::zeros);       // 0: does not pass QC, 1: non-URV, 2: URV
  std::vector<std::string> markerVec(q);
  std::vector<std::string> chrVec(q);  // marker IDs
  std::vector<std::string> posVec(q);  // marker IDs
  std::vector<std::string> refVec(q);  // marker IDs
  std::vector<std::string> altVec(q);  // marker IDs
  std::vector<std::string> infoVec(q);    // marker information: CHR:POS:REF:ALT
  std::vector<double> altFreqVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> MACVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> MAFVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> AF_caseVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> AF_ctrlVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<uint32_t> N_caseVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<uint32_t> N_ctrlVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double>  N_case_homVec(q, arma::datum::nan);
  std::vector<double>  N_ctrl_hetVec(q, arma::datum::nan);
  std::vector<double>  N_case_hetVec(q, arma::datum::nan);
  std::vector<double>  N_ctrl_homVec(q, arma::datum::nan);
  std::vector<uint32_t> N_Vec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> altCountsVec(q, arma::datum::nan);    // allele counts of ALT allele.
  std::vector<double> imputationInfoVec(q, arma::datum::nan);    // imputation info of ALT allele.
  std::vector<double> missingRateVec(q, arma::datum::nan);
  std::vector<double> BetaVec(q, arma::datum::nan);         // beta value for ALT allele
  std::vector<double> seBetaVec(q, arma::datum::nan);
  std::vector<std::string> pvalVec(q, "NA");
  //std::vector<double> pvalVec_val(q, arma::datum::nan);
  std::vector<double> TstatVec(q, arma::datum::nan);
  std::vector<double> TstatVec_flip(q, arma::datum::nan);
  std::vector<double> gyVec(q, arma::datum::nan);
  std::vector<double> varTVec(q, arma::datum::nan);
  std::vector<std::string> pvalNAVec(q, "NA");  
  std::vector<unsigned char> isSPAConvergeVec(q, 0);   // not vector<bool>: bit-packed writes race across markers


  unsigned int m1 = g_cfg.region_maxMarkers_cutoff;     // number of markers in all chunks expect for the last chunk

  // 2026-09-06: the P1/P2 working scratch persists across region calls and is
  // sized once, for the largest chunk a call can hold: m1 tested-marker rows,
  // or the q_anno_maf_weight rows of the ultra-rare pass. It is NOT zeroed:
  // every row/column is written in full before a chunk is stored, and a chunk
  // holds only the rows/columns written for it. The matrices used to arrive as
  // m1 x n copies from R on every call, be truncated, resized (zero-filled)
  // back for the ultra-rare pass, truncated again and streamed through .bin
  // files: 42% of the region C++ at N = 50,000 for a 1.6% matrix product
  // (harness/bench/NOTES.md). Completed chunks now go to chunkStore, in memory
  // unless SAIGE_REGION_CHUNK_MEMORY_MB says otherwise. R passes 1 x 1
  // placeholders (R/SAIGE_SPATest_Region.R).
  (void)t_P1Mat_placeholder; (void)t_P2Mat_placeholder;
  static arma::mat s_regionP1, s_regionP2;
  arma::mat & P1Mat = s_regionP1;
  arma::mat & P2Mat = s_regionP2;
  const unsigned int m1max = std::max(m1, q_anno_maf_weight);
  RegionChunkStore chunkStore(t_outputFile, regionName);
  if(t_regionTestType != "BURDEN"){
    if(P1Mat.n_rows != m1max || P1Mat.n_cols != t_n){ P1Mat.set_size(m1max, t_n); }
    if(P2Mat.n_rows != t_n || P2Mat.n_cols != m1max){ P2Mat.set_size(t_n, m1max); }
  }
  // P2Mat should be of dimension: t_n * m1 
  // Suppose that 
  // n is the sample size in analysis 
  // m (<q) is the number of markers that pass the marker-level QC (e.g., g_cfg.missingRate_cutoff and g_cfg.region_maxMAF_cutoff)
  // VarMat (m x m) is the variance matrix of these m markers
  // VarMat = P1Mat %*% P2Mat, where P1Mat is of (m x n) and P2Mat is of (n x m)
  // Added on 09-17-2021: we collapse all ultra-rare variants (URV) to get one "fake" marker. 
  // That part has been moved to function mainRegionURVInCPP()
  
  std::vector<unsigned int> mPassCVVec;
    
  // conduct marker-level analysis
  std::string pval, pval_noSPA,  pval_c, pval_noSPA_c;
  double Beta, seBeta, Tstat, varT, gy;
  double Beta_c, seBeta_c, Tstat_c, varT_c;
  bool isSPAConverge, is_gtilde, is_Firth, is_FirthConverge;
  arma::vec P1Vec(t_n), P2Vec(t_n);
  arma::vec GVec(t_n);
  arma::vec GZeroVec(t_n);
  std::vector<uint> indexZeroVec;
  std::vector<uint> indexNonZeroVec;
  // Hoisted out of the marker and ultra-rare loops below: conv_to fully
  // overwrites them at each use, so reuse is value-neutral.
  arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;

  // variance ratio. Every site below selects its own into a local TestContext,
  // so all that is decided here is single vs categorical; the
  // assignSingleVarianceRatio call that used to sit in the else branch wrote a
  // model member that no read in this function could reach -- each of the three
  // sites reassigns before its own read -- and it went with the member (S12).
  bool hasVarRatio = true;
  bool isSingleVarianceRatio = ((g_model->m_varRatio_null).n_elem <= 1);
  // initiate chunk information
  unsigned int nchunks = 0; //number of chunks
  unsigned int ichunk = 0; //ith chunk
  unsigned int i1InChunk = 0; //i1th marker in ith chunk
  unsigned int i1 = 0;    // index of Markers (non-URV)
  unsigned int i2 = 0;    // index of Markers (Ultra-Rare Variants, URV)
  unsigned int jm;    // legacy weight-independent index j*q_maf + m
  unsigned int jmr;   // per-set index jm*q_weight + r (weight fastest)

  double cctpval;
  double cctpval_cond;
  // ---- pass A, restructured 2026-09-06 (steps 2-3 of the region plan) -------
  // Each marker's work is split the way the marker path's was in S4:
  //   prepareAndTest(i, k, scratch, reader)  read + impute/flip + QC + weights,
  //       then (rare markers) the test context and getMarkerPval. Everything
  //       it produces lands in preps[k] (one record per marker in the block)
  //       and in column k of blockP1/blockP2. It reads the model and never
  //       writes it, touches no cross-marker state and no per-index vector, so
  //       one marker per worker is safe by construction.
  //   accumulate(i, k)  the serial part, in marker order: the slot writes into
  //       the per-index vectors, the running chunk index, the P1/P2 row copy,
  //       the float `+=` into genoSumMat (order-dependent), the max-collapse
  //       into genoURMat, the group counters, the anno x MAF indicator row and
  //       the chunk flush.
  // Markers are processed in blocks of m1 (= markers_per_chunk): phase 1 over
  // the block -- serially at t_nThreads == 1, on a parallelFor otherwise --
  // then phase 2 over the same block in order. Every expression in phase 2 is
  // the pre-restructure one over the values phase 1 stored, so the output is
  // byte-identical whatever the thread count.
  struct RegionPrep {
    bool readFailed = false;
    bool passQC = false;
    bool isRare = false;
    bool tested = false;      // test outputs below and blockP1/P2 column k are valid
    bool deferredToMainThread = false;   // the ER arm, main-thread only: re-run there
    bool hasError = false;
    std::string errMsg;
    bool flip = false;
    // marker fields (written for every marker that was read, QC or not)
    std::string chr, pos, ref, alt, marker, info;
    double altFreq = 0.0, altCounts = 0.0, missingRate = 0.0, imputeInfo = 0.0;
    double MAF = 0.0;
    double MAC = 0.0;
    double MACcasegroup = 0.0;
    double MACcontrolgroup = 0.0;
    arma::uvec nz;            // nonzero sample indices after impute/flip
    arma::vec gnz;            // their genotype values (all phase 2 needs of GVec)
    arma::vec w0;             // per weight set
    // test outputs
    double Beta = 0.0, seBeta = 0.0, Tstat = 0.0, varT = 0.0, gy = 0.0;
    std::string pval, pval_noSPA;
    bool isSPAConverge = false;
    double Beta_c = 0.0, seBeta_c = 0.0, Tstat_c = 0.0, varT_c = 0.0;
    std::string pval_c, pval_noSPA_c;
    arma::rowvec G1row;       // G1tilde_P_G2tilde_Weighted_Mat.row(i), already weighted
    // binary summary
    double AF_case = 0.0, AF_ctrl = 0.0;
    uint32_t N_case = 0, N_ctrl = 0;
    double N_case_hom = 0.0, N_case_het = 0.0, N_ctrl_hom = 0.0, N_ctrl_het = 0.0;
  };
  const bool doSingleTests = (t_regionTestType != "BURDEN" || t_isSingleinGroupTest);
  const bool keepP1P2 = (t_regionTestType != "BURDEN");
  const unsigned int blockSize = std::max(1u, m1);
  std::vector<RegionPrep> preps(blockSize);
  arma::mat blockP1, blockP2;   // n x blockSize: column k = sqrt(vr) * gtilde / sqrt(vr) * P2Vec of marker k
  if(keepP1P2){
    blockP1.set_size(t_n, blockSize);
    blockP2.set_size(t_n, blockSize);
  }
  Marker::MarkerScratch scratchA(t_n, q_cond);

  auto prepareAndTest = [&](unsigned int i, unsigned int k, Marker::MarkerScratch & sc, Geno::GenoReader & reader)
  {
    RegionPrep & pr = preps[k];
    pr = RegionPrep();
    double altFreq, altCounts, missingRate, imputeInfo;
    std::string chr, ref, alt, marker;
    uint32_t pd;
    bool flip = false;

    ProfileScope rpPrep(rp.enabled ? &rp.prepare : nullptr, &rp.nPrepare);
    sc.GVec.resize(t_n);
    sc.GVec.zeros();

    std::string t_genoIndex_str = t_genoIndex.at(i);
    char* end;
    uint64_t gIndex = std::strtoull( t_genoIndex_str.c_str(), &end,10 );

    bool isReadMarker = getOneMarkerThrough(reader, t_genoType, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          true,   // t_isOutputIndexForMissing
                                          sc.indexForMissing,
                                          false,  // t_isOnlyOutputNonZero
                                          sc.indexNonZeroVec,
                                          sc.GVec,
                                          t_isImputation);
    if(!isReadMarker){
      pr.readFailed = true;
      return;
    }
    pr.chr = chr; pr.pos = std::to_string(pd); pr.ref = ref; pr.alt = alt; pr.marker = marker;
    pr.info = chr+":"+std::to_string(pd)+":"+ref+":"+alt;

    double MAF = std::min(altFreq, 1 - altFreq);
    double MAC = MAF * 2 * t_n * (1 - missingRate);   // checked on 08-10-2021
    flip = imputeGenoAndFlip(sc.GVec, altFreq, altCounts, sc.indexForMissing, g_cfg.impute_method, g_cfg.dosage_zerod_cutoff, g_cfg.dosage_zerod_MAC_cutoff, MAC, sc.indexZeroVec, sc.indexNonZeroVec);
    MAF = std::min(altFreq, 1 - altFreq);
    MAC = std::min(altCounts, t_n *2 - altCounts);
    pr.altFreq = altFreq; pr.missingRate = missingRate; pr.altCounts = altCounts;
    pr.imputeInfo = imputeInfo; pr.MAC = MAC; pr.MAF = MAF; pr.flip = flip;

    if((missingRate > g_cfg.missingRate_cutoff) || (MAF > g_cfg.maxMAFLimit) || (MAF < g_cfg.marker_minMAF_cutoff) || (MAC < g_cfg.marker_minMAC_cutoff) || (imputeInfo < g_cfg.marker_minINFO_cutoff)){
      return;
    }
    pr.passQC = true;

    pr.w0.set_size(q_weight);
    for(unsigned int r = 0; r < q_weight; r++){
      if(setKindVec[r] == SET_CUSTOM){
        pr.w0(r) = t_weight(i, setColVec[r]);
      }else if(setKindVec[r] == SET_BETA){
        pr.w0(r) = boost::math::pdf(betaDistVec[setColVec[r]], MAF);
      }else{
        pr.w0(r) = 1.0;
      }
    }

    sc.indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(sc.indexNonZeroVec);
    pr.nz = sc.indexNonZeroVec_arma;
    pr.gnz = sc.GVec.elem(pr.nz);
    rpPrep.stop();

    if(!(MAC > g_cfg.region_minMAC_cutoff)){
      return;   // ultra-rare: collapsed in accumulate()
    }
    pr.isRare = true;
    ProfileScope rpRare(rp.enabled ? &rp.rareArm : nullptr, &rp.nRareArm);

    // S12: this marker's test context, derived and consumed here. A common
    // marker uses the non-sparse variance ratio, a rare one follows the model;
    // nothing is parked on g_model, so two markers can be tested at once.
    // hasVarRatio is per call, as in testOneModel: selectVarianceRatio writes
    // it, nothing here reads it, and a frame-shared flag written from every
    // worker is a data race (TSAN, 2026-09-06) even when the value never changes.
    bool hasVarRatio = true;
    SAIGE::TestContext ctx;
    ctx.useSparseGRM = (MAC > g_model->m_cateVarRatioMinMACVecExclude.back())
                         ? false : g_model->m_flagSparseGRM;
    ctx.noadjCov = false;
    if(!isSingleVarianceRatio){
      ctx.varRatioVal = g_model->selectVarianceRatio(MAC, ctx.useSparseGRM, false, hasVarRatio);
    }else{
      ctx.varRatioVal = g_model->selectSingleVarianceRatio(ctx.useSparseGRM, false);
    }

    if(doSingleTests){ //perform single-variant assoc tests
      bool is_gtilde, is_Firth, is_FirthConverge;
      sc.indexZeroVec_arma = arma::conv_to<arma::uvec>::from(sc.indexZeroVec);
      ProfileScope rpTest(rp.enabled ? &rp.test : nullptr, &rp.nTest);
      if(MAC <=  g_cfg.MACCutoffforER && t_traitType == "binary"){
        Unified_getMarkerPval(
                  sc.GVec,
                  false, // bool t_isOnlyOutputNonZero,
        sc.indexNonZeroVec_arma, sc.indexZeroVec_arma, pr.Beta, pr.seBeta, pr.pval, pr.pval_noSPA, pr.Tstat, pr.gy, pr.varT, altFreq, pr.isSPAConverge, sc.gtildeVec, is_gtilde, true, sc.P2Vec, isCondition, pr.Beta_c, pr.seBeta_c, pr.pval_c, pr.pval_noSPA_c, pr.Tstat_c, pr.varT_c, sc.G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, true, flip, ctx);
      }else{
        Unified_getMarkerPval(
                  sc.GVec,
                  false, // bool t_isOnlyOutputNonZero,
        sc.indexNonZeroVec_arma, sc.indexZeroVec_arma, pr.Beta, pr.seBeta, pr.pval, pr.pval_noSPA, pr.Tstat, pr.gy, pr.varT, altFreq, pr.isSPAConverge, sc.gtildeVec, is_gtilde, true, sc.P2Vec, isCondition, pr.Beta_c, pr.seBeta_c, pr.pval_c, pr.pval_noSPA_c, pr.Tstat_c, pr.varT_c, sc.G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, false, flip, ctx);
      }
      rpTest.stop();

      if(isCondition){
        if(q_weight == 1){
          // the legacy baked row: G2-side conditioning weights times this
          // marker's single-set weight w0
          pr.G1row = sc.G1tilde_P_G2tilde_Vec % w0G2Vec_cond.t() * pr.w0(0);
        }else{
          // design D7 (C11): at K>1 the exported matrix is G1-raw -- still
          // scaled by the (set-independent) G2 conditioning weights, but with
          // NO per-row w0 baked in. region_postprocess applies set r's weight
          // per (anno, max_MAF, set) cell.
          pr.G1row = sc.G1tilde_P_G2tilde_Vec % w0G2Vec_cond.t();
        }
      }

      if(keepP1P2){
        blockP1.col(k) = sqrt(ctx.varRatioVal)*sc.gtildeVec;
        blockP2.col(k) = sqrt(ctx.varRatioVal)*sc.P2Vec;
      }
      pr.tested = true;
    }//if(doSingleTests)

    arma::vec dosage_case, dosage_ctrl;
    if(t_traitType == "binary"){
      dosage_case = sc.GVec.elem(g_model->m_case_indices);
      dosage_ctrl = sc.GVec.elem(g_model->m_ctrl_indices);
      pr.MACcasegroup = arma::accu(dosage_case);
      pr.MACcontrolgroup = arma::accu(dosage_ctrl);
    }

    if(doSingleTests){
      if(t_traitType == "binary"){
        double AF_case = arma::mean(dosage_case) /2;
        double AF_ctrl = arma::mean(dosage_ctrl) /2;
        if(flip){
          AF_case = 1-AF_case;
          AF_ctrl = 1-AF_ctrl;
        }
        pr.AF_case = AF_case;
        pr.AF_ctrl = AF_ctrl;
        pr.N_case = dosage_case.n_elem;
        pr.N_ctrl = dosage_ctrl.n_elem;

        arma::uvec N_case_ctrl_het_hom0;
        if(t_isMoreOutput){
          N_case_ctrl_het_hom0 = arma::find(dosage_case <= 2 && dosage_case >=1.5);
          pr.N_case_hom = N_case_ctrl_het_hom0.n_elem;
          N_case_ctrl_het_hom0 = arma::find(dosage_case < 1.5 && dosage_case >= 0.5);
          pr.N_case_het = N_case_ctrl_het_hom0.n_elem;
          N_case_ctrl_het_hom0 = arma::find(dosage_ctrl <= 2 && dosage_ctrl >=1.5);
          pr.N_ctrl_hom = N_case_ctrl_het_hom0.n_elem;
          N_case_ctrl_het_hom0 = arma::find(dosage_ctrl < 1.5 && dosage_ctrl >= 0.5);
          pr.N_ctrl_het = N_case_ctrl_het_hom0.n_elem;
          if(flip){
                pr.N_case_hom = pr.N_case - pr.N_case_het -  pr.N_case_hom;
                pr.N_ctrl_hom = pr.N_ctrl - pr.N_ctrl_het - pr.N_ctrl_hom;
          }
        }
      }
    } //if(doSingleTests)
  };

  auto accumulate = [&](unsigned int i, unsigned int k)
  {
    const RegionPrep & pr = preps[k];
    // the marker fields, for every marker that was read (QC or not)
    chrVec.at(i) = pr.chr;
    posVec.at(i) = pr.pos;
    refVec.at(i) = pr.ref;
    altVec.at(i) = pr.alt;
    markerVec.at(i) = pr.marker;          // marker IDs
    infoVec.at(i) = pr.info;              // marker information: CHR:POS:REF:ALT
    altFreqVec.at(i) = pr.altFreq;        // allele frequencies of ALT allele, this is not always < 0.5.
    missingRateVec.at(i) = pr.missingRate;
    altCountsVec.at(i) = pr.altCounts;
    MACVec.at(i) = pr.MAC;
    MAFVec.at(i) = pr.MAF;
    imputationInfoVec.at(i) = pr.imputeInfo;
    if(!pr.passQC){
      return;
    }
    const uint nNonZero = pr.nz.n_elem;
    const bool flip = pr.flip;
    if(pr.isRare){  // not Ultra-Rare Variants
      ProfileScope rpRare(rp.enabled ? &rp.rareArm : nullptr, &rp.nRareArm);
      indicatorVec.at(i) = 1;
      if(i1InChunk == 0){
        std::cout << "Start analyzing chunk " << ichunk << "....." << std::endl;
      }
      if(pr.tested){
        BetaVec.at(i) = pr.Beta * (1 - 2*flip);  // Beta if flip = false, -1 * Beta is flip = true
        seBetaVec.at(i) = pr.seBeta;
        pvalVec.at(i) = pr.pval;
        pvalNAVec.at(i) = pr.pval_noSPA;
        TstatVec.at(i) = pr.Tstat * (1 - 2*flip);
        TstatVec_flip.at(i) = pr.Tstat;
        gyVec.at(i) = pr.gy;
        varTVec.at(i) = pr.varT;
        isSPAConvergeVec.at(i) = pr.isSPAConverge;
        if(isCondition){
          Beta_cVec.at(i) = pr.Beta_c * (1 - 2*flip);  // Beta if flip = false, -1 * Beta is flip = true
          seBeta_cVec.at(i) = pr.seBeta_c;
          pval_cVec.at(i) = pr.pval_c;
          pvalNA_cVec.at(i) = pr.pval_noSPA_c;
          Tstat_cVec.at(i) = pr.Tstat_c * (1 - 2*flip);
          varT_cVec.at(i) = pr.varT_c;
          G1tilde_P_G2tilde_Weighted_Mat.row(i) = pr.G1row;
        }
        if(keepP1P2){
          P1Mat.row(i1InChunk) = blockP1.col(k).t();
          P2Mat.col(i1InChunk) = blockP2.col(k);
        }
      }
      i1 += 1;
      i1InChunk += 1;

      MAFIndicatorVec.zeros();
      MAFIndicatorVec.elem( find(maxMAFVec >= pr.MAF) ).ones();
      annoMAFIndicatorVec.zeros();
      for(unsigned int j = 0; j < q_anno; j++){
        if(annoIndicatorMat(i,j) == 1){
          maxMAFperAnno(j) = std::max(maxMAFperAnno(j), pr.MAF);
          for(unsigned int m = 0; m < q_maf; m++){
            if(MAFIndicatorVec(m) == 1){
              jm = j*q_maf + m;
              for(unsigned int r = 0; r < q_weight; r++){
                jmr = jm*q_weight + r;
                annoMAFIndicatorVec(jmr) = 1;
                MAC_GroupVec(jmr) = MAC_GroupVec(jmr) + pr.MAC;
                if(t_traitType == "binary"){
                  MACCase_GroupVec(jmr) = MACCase_GroupVec(jmr) + pr.MACcasegroup;
                  MACControl_GroupVec(jmr) = MACControl_GroupVec(jmr) + pr.MACcontrolgroup;
                }
                for(unsigned int k2 = 0; k2 < nNonZero; k2++){
                  genoSumMat(pr.nz(k2), jmr) = genoSumMat(pr.nz(k2), jmr) + pr.w0(r)*pr.gnz(k2);
                  // jm-indexed, weight-independent: accumulate once per (j,m)
                  if(r == 0){
                    genoSumcount_noweight(jm) = genoSumcount_noweight(jm) + pr.gnz(k2);
                  }
                }
                NumRare_GroupVec(jmr) = NumRare_GroupVec(jmr) + 1;
              }
            }
          }
        }
      }
      annoMAFIndicatorMat.row(i) = annoMAFIndicatorVec.t();
      if(doSingleTests){
        if(t_traitType == "binary"){
          AF_caseVec.at(i) = pr.AF_case;
          AF_ctrlVec.at(i) = pr.AF_ctrl;
          N_caseVec.at(i) = pr.N_case;
          N_ctrlVec.at(i) = pr.N_ctrl;
          if(t_isMoreOutput){
            N_case_homVec.at(i) = pr.N_case_hom;
            N_case_hetVec.at(i) = pr.N_case_het;
            N_ctrl_homVec.at(i) = pr.N_ctrl_hom;
            N_ctrl_hetVec.at(i) = pr.N_ctrl_het;
          }
        }else if(t_traitType == "quantitative"){
          N_Vec.at(i) = t_n;
        }
      }
    }else{   // Ultra-Rare Variants (URV)
      ProfileScope rpUrc(rp.enabled ? &rp.urCollapse : nullptr, &rp.nUrCollapse);
      indicatorVec.at(i) = 2;
      arma::vec MAFIndicatorVecUR(maxMAFVec.n_elem);
      MAFIndicatorVecUR.zeros();
      MAFIndicatorVecUR.elem( find(maxMAFVec >= pr.MAF) ).ones();
      annoMAFIndicatorVec.zeros();
      for(unsigned int j = 0; j < q_anno; j++){
        if(annoIndicatorMat(i,j) == 1){
          maxMAFperAnno(j) = std::max(maxMAFperAnno(j), pr.MAF);
          for(unsigned int m = 0; m < q_maf; m++){
            if(MAFIndicatorVecUR(m) == 1){
              jm = j*q_maf + m;
              for(unsigned int r = 0; r < q_weight; r++){
                jmr = jm*q_weight + r;
                annoMAFIndicatorVec(jmr) = 2;
                if(setKindVec[r] != SET_CUSTOM){
                  // beta / No_Weights sets: unweighted collapse (per-kind table above)
                  for(unsigned int k2 = 0; k2 < nNonZero; k2++){
                    genoURMat(pr.nz(k2), jmr) = std::max(genoURMat(pr.nz(k2), jmr), pr.gnz(k2));
                  }
                  if(q_weight > 1){
                    // keep the jm-indexed unweighted collapse authoritative for
                    // every set composition (beta-only included): the once-per-jm
                    // genoSumcount_noweight accumulation reads its support. At
                    // q_weight == 1 this never runs -- legacy leaves the matrix
                    // untouched on the beta path (max is idempotent, so repeats
                    // across sets are exact).
                    for(unsigned int k2 = 0; k2 < nNonZero; k2++){
                      genoURMat_noweights(pr.nz(k2), jm) = std::max(genoURMat_noweights(pr.nz(k2), jm) , (pr.gnz(k2)));
                    }
                  }
                }else{
                  // custom set: collapse the w_r-weighted genotypes; the jm-indexed
                  // unweighted collapse is kept for the count accumulators (max is
                  // idempotent, so repeating it per custom set is exact)
                  for(unsigned int k2 = 0; k2 < nNonZero; k2++){
                    genoURMat(pr.nz(k2), jmr) = std::max(genoURMat(pr.nz(k2), jmr) , t_weight(i, setColVec[r]) * (pr.gnz(k2)));
                    genoURMat_noweights(pr.nz(k2), jm) = std::max(genoURMat_noweights(pr.nz(k2), jm) , (pr.gnz(k2)));
                  }
                }
                NumUltraRare_GroupVec(jmr) = NumUltraRare_GroupVec(jmr) + 1;
              }
            }
          }
        }
      }
      annoMAFIndicatorMat.row(i) = annoMAFIndicatorVec.t();
      i2 += 1;
    }

    if(i1InChunk == m1 ){
      ProfileScope rpFlush(rp.enabled ? &rp.flush : nullptr, &rp.nFlush);
      std::cout << "In chunks 0-" << ichunk << ", " << i2 << " markers are ultra-rare and " << i1 << " markers are not ultra-rare." << std::endl;
      if(t_regionTestType != "BURDEN"){
        chunkStore.add(ichunk, P1Mat.rows(0, i1InChunk - 1), P2Mat.cols(0, i1InChunk - 1));
      }
      mPassCVVec.push_back(m1);
      ichunk += 1;
      i1InChunk = 0;
      nchunks = nchunks + 1;
    }
  };

  // Phase 1 on a parallelFor. Each range checks out a worker state (one
  // scratch + one reader clone) from the marker path's pool, refuses the R
  // API for its duration, and records a marker's failure IN ITS OWN prep
  // record: the main thread then acts on the lowest index, which is what a
  // serial run would have stopped at (the marker path's post-passes, S6).
  struct RegionBlockWorker : public RcppParallel::Worker {
    RegionBlockWorker(MarkerStatePool & t_pool, unsigned int t_iStart,
                      std::vector<RegionPrep> & t_preps,
                      const std::function<void(unsigned int, unsigned int, Marker::MarkerScratch&, Geno::GenoReader&)> & t_fn)
      : m_pool(t_pool), m_iStart(t_iStart), m_preps(t_preps), m_fn(t_fn) {}
    void operator()(std::size_t t_begin, std::size_t t_end) override {
      std::unique_ptr<MarkerWorkerState> state = m_pool.checkout();
      RefuseMainThreadGuard noRAPI;
      for(std::size_t k = t_begin; k < t_end; k++){
        try {
          m_fn(m_iStart + k, k, *state->scratch[0], *state->readers[0]);
        } catch (const SKATMainThread::NeedsMainThread&) {
          m_preps[k] = RegionPrep();
          m_preps[k].deferredToMainThread = true;
        } catch (const std::exception& e) {
          m_preps[k] = RegionPrep();
          m_preps[k].hasError = true;
          m_preps[k].errMsg = e.what();
        }
      }
      m_pool.checkin(std::move(state));
      SAIGE_TSAN_HAPPENS_BEFORE(&m_pool);
    }
    MarkerStatePool & m_pool;
    unsigned int m_iStart;
    std::vector<RegionPrep> & m_preps;
    const std::function<void(unsigned int, unsigned int, Marker::MarkerScratch&, Geno::GenoReader&)> & m_fn;
  };
  const int nThreadsRegion = std::max(1, t_nThreads);
  std::unique_ptr<MarkerStatePool> regionPool;
  std::vector<std::vector<int>> regionGroups;
  if(nThreadsRegion > 1){
    regionGroups.push_back(std::vector<int>(1, 0));   // model 0, reader 0
    regionPool.reset(new MarkerStatePool(regionGroups, nThreadsRegion));
  }
  const std::function<void(unsigned int, unsigned int, Marker::MarkerScratch&, Geno::GenoReader&)> prepareFn = prepareAndTest;

  // The block driver. A failed read stops the region at that marker exactly
  // where the serial loop's `break` did: nothing after it is accumulated, and
  // the rest of the call proceeds on what it has.
  {
    bool shortRead = false;
    unsigned int iStart = 0;
    while(iStart < q0 && !shortRead){
      const unsigned int iEnd = std::min(q0, iStart + blockSize);
      const unsigned int nBlock = iEnd - iStart;
      unsigned int nPrepared = nBlock;
      if(nThreadsRegion == 1){
        // literally serial: no pool, no guard, ER runs inline
        nPrepared = 0;
        for(unsigned int i = iStart; i < iEnd; i++){
          prepareAndTest(i, i - iStart, scratchA, *g_reader);
          if(preps[i - iStart].readFailed){
            break;
          }
          nPrepared++;
        }
      }else{
        RegionBlockWorker worker(*regionPool, iStart, preps, prepareFn);
        std::size_t grain = 8;
        RcppParallel::parallelFor(0, (std::size_t)nBlock, worker, grain, nThreadsRegion);
        SAIGE_TSAN_HAPPENS_AFTER(regionPool.get());
        // post-pass 1: the error a serial run would have reported (lowest index)
        for(unsigned int k = 0; k < nBlock; k++){
          if(preps[k].hasError){
            throw std::runtime_error(preps[k].errMsg);
          }
        }
        // post-pass 2: truncate at the earliest short read
        for(unsigned int k = 0; k < nBlock; k++){
          if(preps[k].readFailed){
            nPrepared = k;
            break;
          }
        }
        // post-pass 3: the ER markers, on the main thread. Re-preparing is
        // byte-neutral because prepare is deterministic; ComputeExact re-seeds
        // at entry, so the order the deferred markers run in does not matter.
        for(unsigned int k = 0; k < nPrepared; k++){
          if(preps[k].deferredToMainThread){
            prepareAndTest(iStart + k, k, scratchA, *g_reader);
            if(preps[k].readFailed){
              nPrepared = k;
              break;
            }
          }
        }
      }
      for(unsigned int i = iStart; i < iStart + nPrepared; i++){
        accumulate(i, i - iStart);
      }
      if(nPrepared < nBlock){
        std::cout << "ERROR: Reading " << (iStart + nPrepared) << "th marker failed." << std::endl;
        shortRead = true;
      }
      Rcpp::checkUserInterrupt();
      iStart = iEnd;
    }
  }



//the second last chunk
  if(i1InChunk != 0){
    ProfileScope rpFlush2(rp.enabled ? &rp.flush : nullptr, &rp.nFlush);
    std::cout << "In chunks 0-" << ichunk << ", " << i2 << " markers are ultra-rare and " << i1 << " markers are not ultra-rare." << std::endl;
    if(t_regionTestType != "BURDEN"){	  
      chunkStore.add(ichunk, P1Mat.rows(0, i1InChunk - 1), P2Mat.cols(0, i1InChunk - 1));
    }
    ichunk = ichunk + 1;
    //}
    mPassCVVec.push_back(i1InChunk);
    nchunks = nchunks + 1; 
    i1InChunk = 0;
  }


//for all UR variants
if(i2 > 0){
  ProfileScope rpUr(rp.enabled ? &rp.urPass : nullptr, &rp.nUrPass);
  int m1new = std::max(m1, q_anno_maf_weight);
  //if(isWeightCustomized){
    //weightURMat_cnt.replace(0, 1); 	  
    //genoURMat = genoURMat / weightURMat_cnt;
    //weightURMat.ones();    
  //}

  // (P1Mat/P2Mat already hold m1new = m1max rows/cols: sized once per call above)
  (void)m1new;

  arma::mat XV, XXVX_inv;
  g_model->extract_XV_XXVX_inv(XV, XXVX_inv);
  //the last chunk for UR
  unsigned int i;
  for(unsigned int j = 0; j < q_anno; j++){
     for(unsigned int m = 0; m < q_maf; m++){
	jm = j*q_maf+m;
	if(q_weight > 1){
	  // Fix (C5 review #1): the jm-indexed weight-independent count must
	  // not depend on any one set's weighted support (a zero weight in
	  // custom set 0 would undercount the altCounts every set's BURDEN
	  // reads). Accumulate once per (j,m) over the unweighted collapse's
	  // own support, replicating the legacy per-sample additions. At
	  // q_weight == 1 this is skipped: the accumulation happens inside the
	  // single set's emission loop below, textually the legacy code.
	  arma::vec genoURVec_noweights_jm = genoURMat_noweights.col(jm);
	  arma::uvec indexForNonZero_nw = arma::find(genoURVec_noweights_jm != 0);
	  for(unsigned int k = 0; k < indexForNonZero_nw.n_elem; k++){
	    genoSumcount_noweight(jm) = genoSumcount_noweight(jm) + genoURVec_noweights_jm(indexForNonZero_nw(k));
	  }
	}
	for(unsigned int r = 0; r < q_weight; r++){
	jmr = jm*q_weight + r;
	arma::vec genoURVec = genoURMat.col(jmr);
	arma::vec genoURVec_noweights = genoURMat_noweights.col(jm);
	int n = genoURVec.size();
	arma::uvec indexForNonZero = arma::find(genoURVec != 0);
	i = q0 + jmr;
	markerVec.at(i) = "UR";             // marker IDs
	if(indexForNonZero.n_elem > 0){
	//URindVec.push_back(jm+1);
	  double altFreq = arma::mean(genoURVec)/2;
	  double altCounts = arma::accu(genoURVec);
	  double missingRate = 0;
	  double imputeInfo = 1;
    	  std::string chr, ref, alt, marker;
    	  bool flip = false;
	  std::string info = "UR";	
    	  double MAF = std::min(altFreq, 1 - altFreq);
	  double w0;
	  double MAC = MAF * 2 * t_n * (1 - missingRate);
	  std::vector<uint32_t> indexForMissing;
    	  flip = imputeGenoAndFlip(genoURVec, altFreq, altCounts, indexForMissing, g_cfg.impute_method, g_cfg.dosage_zerod_cutoff, g_cfg.dosage_zerod_MAC_cutoff, MAC, indexZeroVec, indexNonZeroVec);
	  // Count/emission split (C5 review #1): at q_weight == 1 each branch's
	  // first arm is textually the legacy loop -- genoSumcount_noweight
	  // gated on THIS set's collapse support (for a weighted file that is
	  // the WEIGHTED support, exactly as upstream). At q_weight > 1 the
	  // count was already accumulated once per (j,m) before the r loop, so
	  // the emission loops touch genoSumMat only.
	  if(setKindVec[r] == SET_CUSTOM){
	    // LEDGER #17 FIXED. w0 used to be left UNASSIGNED here, and the
	    // conditional export below multiplies by it -- so a K=1 custom-weighted
	    // group file plus --condition read an uninitialised double.
	    //
	    // 1.0 is the value, and it is not a placeholder: a custom set's UR
	    // collapse is WEIGHTED (`max(w_r * G)` at the genoURMat line above), so
	    // the weight is already inside genoURVec and hence inside
	    // G1tilde_P_G2tilde_Vec. Multiplying by w_r again would double-count it.
	    // That is exactly what the K>1 path does -- it exports the G1-raw row
	    // and region_postprocess applies 1 for custom sets (design D7) -- so
	    // this makes the two agree instead of merely making the read defined.
	    // The emission loops below deliberately do NOT scale by w0, for the
	    // same reason.
	    w0 = 1.0;
	    //genoSumMat.col(jm) = genoSumMat.col(jm) + genoURVec % (weightURMat.col(jm));
	    //genoSumMat.col(jm) = genoSumMat.col(jm) + genoURVec;
	    if(q_weight == 1){
	      for(unsigned int k = 0; k < indexForNonZero.n_elem; k++){
                genoSumMat(indexForNonZero(k), jmr) = genoSumMat(indexForNonZero(k), jmr) + genoURVec(indexForNonZero(k));
	 	genoSumcount_noweight(jm) = genoSumcount_noweight(jm) + genoURVec_noweights(indexForNonZero(k));
	      }
	    }else{
	      for(unsigned int k = 0; k < indexForNonZero.n_elem; k++){
                genoSumMat(indexForNonZero(k), jmr) = genoSumMat(indexForNonZero(k), jmr) + genoURVec(indexForNonZero(k));
	      }
	    }
    	  }else if(setKindVec[r] == SET_BETA){
                w0 = boost::math::pdf(betaDistVec[setColVec[r]], MAF);
	        if(q_weight == 1){
	        for(unsigned int k = 0; k < indexForNonZero.n_elem; k++){
                genoSumMat(indexForNonZero(k), jmr) = genoSumMat(indexForNonZero(k), jmr) + genoURVec(indexForNonZero(k)) * w0;
	 	genoSumcount_noweight(jm) = genoSumcount_noweight(jm) + genoURVec(indexForNonZero(k));
	        }
	        }else{
	        for(unsigned int k = 0; k < indexForNonZero.n_elem; k++){
                genoSumMat(indexForNonZero(k), jmr) = genoSumMat(indexForNonZero(k), jmr) + genoURVec(indexForNonZero(k)) * w0;
	        }
	        }
    	  }else{
	        // No_Weights: unweighted collapse contributes with weight 1.
	        // w0 = 1.0 is assigned (C5 review #3): the conditional export
	        // below reads w0, and this branch is reachable at q_weight == 1
	        // via t_includeNoWeights (unlike the pre-existing custom-branch
	        // non-assignment above, which stays as ledgered legacy behavior).
	        w0 = 1.0;
	        if(q_weight == 1){
	        for(unsigned int k = 0; k < indexForNonZero.n_elem; k++){
                genoSumMat(indexForNonZero(k), jmr) = genoSumMat(indexForNonZero(k), jmr) + genoURVec(indexForNonZero(k));
	 	genoSumcount_noweight(jm) = genoSumcount_noweight(jm) + genoURVec(indexForNonZero(k));
	        }
	        }else{
	        for(unsigned int k = 0; k < indexForNonZero.n_elem; k++){
                genoSumMat(indexForNonZero(k), jmr) = genoSumMat(indexForNonZero(k), jmr) + genoURVec(indexForNonZero(k));
	        }
	        }
    	  }

	 //genoSumMat.col(1).print("genoSumMat.col(1)");


  	 if(t_regionTestType != "BURDEN"){
	    arma::vec genoSumMatvec1 = genoSumMat.col(jmr);
	    arma::vec genoSumMatvec2 = XV * genoSumMatvec1;
	    arma::vec genoSumMatvec3 = genoSumMatvec1 - XXVX_inv * genoSumMatvec2;
	    genoSumMat.col(jmr) = genoSumMatvec3;
          }

	  MAC = MAF * 2 * t_n * (1 - missingRate);   // checked on 08-10-2021
if(t_regionTestType != "BURDEN" || t_isSingleinGroupTest){

      // S12: this marker's test context, derived and consumed here. Same two
      // steps as before -- a common marker uses the non-sparse variance ratio,
      // a rare one follows the model -- but nothing is parked on g_model, so no
      // later read can pick up a marker's choice by accident (LEDGER #27/#50)
      // and two markers could be tested at once without racing.
      SAIGE::TestContext ctx;
      ctx.useSparseGRM = (MAC > g_model->m_cateVarRatioMinMACVecExclude.back())
                           ? false : g_model->m_flagSparseGRM;
      ctx.noadjCov = false;
      if(!isSingleVarianceRatio){
        ctx.varRatioVal = g_model->selectVarianceRatio(MAC, ctx.useSparseGRM, false, hasVarRatio);
      }else{
        ctx.varRatioVal = g_model->selectSingleVarianceRatio(ctx.useSparseGRM, false);
      }




	    annoMAFIndicatorVec.zeros();
	    annoMAFIndicatorVec(jmr) = 1;
	    annoMAFIndicatorMat.row(i) = annoMAFIndicatorVec.t();
    	    infoVec.at(i) = info;                 // marker information: CHR:POS:REF:ALT
    	    altFreqVec.at(i) = altFreq;	// allele frequencies of ALT allele, this is not always < 0.5.
	    altCountsVec.at(i) = altCounts;
	    missingRateVec.at(i) = missingRate;
    	    MACVec.at(i) = MAC;
    	    MAFVec.at(i) = MAF;

            indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);
            indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);


	    ProfileScope rpUrTest(rp.enabled ? &rp.urTest : nullptr, &rp.nUrTest);
	    if(MAC <= g_cfg.MACCutoffforER && t_traitType == "binary" && erAllowedVec[r]){
              g_model->getMarkerPval(genoURVec, indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, altFreq, Tstat, gy, varT, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, true, flip, ctx);
	    }else{
              g_model->getMarkerPval(genoURVec, indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, altFreq, Tstat, gy, varT, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, false, flip, ctx);
	    }
	    rpUrTest.stop();

            BetaVec.at(i) = Beta* (1 - 2*flip);
            seBetaVec.at(i) = seBeta;
            pvalVec.at(i) = pval;
	    //pvalVec_val.at(i) = std::stod(pval);
            pvalNAVec.at(i) = pval_noSPA;
            TstatVec.at(i) = Tstat * (1 - 2*flip);
            TstatVec_flip.at(i) = Tstat;
            gyVec.at(i) = gy;
            varTVec.at(i) = varT;
            isSPAConvergeVec.at(i) = isSPAConverge;
            if(isCondition){
              Beta_cVec.at(i) = Beta_c * (1 - 2*flip);  // Beta if flip = false, -1 * Beta is flip = true
              seBeta_cVec.at(i) = seBeta_c;
              pval_cVec.at(i) = pval_c;
              pvalNA_cVec.at(i) = pval_noSPA_c;
              Tstat_cVec.at(i) = Tstat_c * (1 - 2*flip);
              varT_cVec.at(i) = varT_c;
              if(q_weight == 1){
                // w0 is assigned on every branch above (LEDGER #17 was that the
                // custom branch did not assign it and this line read garbage;
                // fixed, with 1.0 -- see the comment there for why 1 and not
                // w_r). Reachable only at q_weight == 1: a K>1 run takes the
                // raw-export branch below, which reads no w0.
                G1tilde_P_G2tilde_Weighted_Mat.row(i) = G1tilde_P_G2tilde_Vec % w0G2Vec_cond.t() * w0;
              }else{
                // design D7 (C11): G1-raw row at K>1, per-set weight applied
                // in region_postprocess (1 for custom sets -- their collapse
                // already baked w_r into genoURVec -- and pdf(beta_r, MAF)
                // for beta sets, matching what the baked path would use)
                G1tilde_P_G2tilde_Weighted_Mat.row(i) = G1tilde_P_G2tilde_Vec % w0G2Vec_cond.t();
              }
            }

    	    chrVec.at(i) = "UR";
    	    posVec.at(i) = "UR";
    	    refVec.at(i) = "UR";
    	    altVec.at(i) = "UR";


	    std::string str = std::to_string(maxMAFVec.at(m));
	    str.erase ( str.find_last_not_of('0') + 1, std::string::npos );
    	    markerVec.at(i) = regionName + ":" + annoStringVec.at(j) + ":" + str ;
	    if(q_weight > 1){
	      // per-set UR pseudo-markers get a weight-set suffix; q_weight == 1
	      // keeps the legacy name (and legacy output bytes)
	      std::string wname = (r < t_weightNames.size()) ? t_weightNames.at(r) : ("W" + std::to_string(r + 1));
	      markerVec.at(i) = markerVec.at(i) + ":" + wname;
	    }
            arma::vec dosage_case, dosage_ctrl;
            MAC_GroupVec(jmr) = MAC_GroupVec(jmr) + MAC;
            if(t_traitType == "binary"){
                        dosage_case = genoURVec.elem(g_model->m_case_indices);
                        dosage_ctrl = genoURVec.elem(g_model->m_ctrl_indices);
                        MACcasegroup = arma::accu(dosage_case);
                        MACcontrolgroup = arma::accu(dosage_ctrl);
                        MACCase_GroupVec(jmr) = MACCase_GroupVec(jmr) + MACcasegroup;
                        MACControl_GroupVec(jmr) = MACControl_GroupVec(jmr) + MACcontrolgroup;
            }

          if(t_traitType == "binary"){
            AF_case = arma::mean(dosage_case) /2;
            AF_ctrl = arma::mean(dosage_ctrl) /2;
            if(flip){
              AF_case = 1-AF_case;
              AF_ctrl = 1-AF_ctrl;
            }
          AF_caseVec.at(i) = AF_case;
          AF_ctrlVec.at(i) = AF_ctrl;
          N_caseVec.at(i) = dosage_case.n_elem;
          N_ctrlVec.at(i) = dosage_ctrl.n_elem;

         arma::uvec N_case_ctrl_het_hom0;
        if(t_isMoreOutput){
          N_case_ctrl_het_hom0 = arma::find(dosage_case <= 2 && dosage_case >=1.5);
          N_case_homVec.at(i)  = N_case_ctrl_het_hom0.n_elem;
          N_case_ctrl_het_hom0 = arma::find(dosage_case < 1.5 && dosage_case >= 0.5);
          N_case_hetVec.at(i) = N_case_ctrl_het_hom0.n_elem;
          N_case_ctrl_het_hom0 = arma::find(dosage_ctrl <= 2 && dosage_ctrl >=1.5);
          N_ctrl_homVec.at(i) = N_case_ctrl_het_hom0.n_elem;
          N_case_ctrl_het_hom0 = arma::find(dosage_ctrl < 1.5 && dosage_ctrl >= 0.5);
          N_ctrl_hetVec.at(i) = N_case_ctrl_het_hom0.n_elem;
          if(flip){
                N_case_homVec.at(i) = N_case - N_case_hetVec.at(i) -  N_case_homVec.at(i);
                N_ctrl_homVec.at(i) = N_ctrl - N_ctrl_hetVec.at(i) - N_ctrl_homVec.at(i);
          }
        }

        }else if(t_traitType == "quantitative"){
          N_Vec.at(i) = n;
        }
      if(t_regionTestType != "BURDEN"){	
        P1Mat.row(i1InChunk) = sqrt(ctx.varRatioVal)*gtildeVec.t();
	//gtildeVec.print("gtildeVec");
        P2Mat.col(i1InChunk) = sqrt(ctx.varRatioVal)*P2Vec;
	//P2Vec.print("P2Vec");

      } //if(t_regionTestType != "BURDEN"){
    }else{//if(t_regionTestType != "BURDEN" || t_isSingleinGroupTest){
            MAC_GroupVec(jmr) = MAC_GroupVec(jmr) + MAC;
	                arma::vec dosage_case, dosage_ctrl;
            if(t_traitType == "binary"){
                        dosage_case = genoURVec.elem(g_model->m_case_indices);
                        dosage_ctrl = genoURVec.elem(g_model->m_ctrl_indices);
                        MACcasegroup = arma::accu(dosage_case);
                        MACcontrolgroup = arma::accu(dosage_ctrl);
                        MACCase_GroupVec(jmr) = MACCase_GroupVec(jmr) + MACcasegroup;
                        MACControl_GroupVec(jmr) = MACControl_GroupVec(jmr) + MACcontrolgroup;
            }
    }

    i1InChunk = i1InChunk + 1;
    i1 = i1 + 1;
    }
    }// for r (weight sets)
   }
  }


  if(i1InChunk != 0){
    nchunks = nchunks + 1;
    if(t_regionTestType != "BURDEN"){
      chunkStore.add(ichunk, P1Mat.rows(0, i1InChunk - 1), P2Mat.cols(0, i1InChunk - 1));
    }
      ichunk = ichunk + 1;
//    }
    mPassCVVec.push_back(i1InChunk);
  }
   //std::cout << "P1Mat.n_rows ok2 " << P1Mat.n_rows << std::endl; 

  }// if(i2 > 0)    
  int mPassCVVecsize = mPassCVVec.size();
  nchunks = mPassCVVecsize;


arma::mat VarMat;
//(i1, i1);

if(t_regionTestType != "BURDEN"){
  ProfileScope rpVm(rp.enabled ? &rp.varmat : nullptr, &rp.nVarmat);
  VarMat.resize(i1, i1);

  // VarMat is assembled from the SAME per-chunk sub-block products the disk
  // path issued -- same shapes, same operands -- whether the chunks sit in
  // memory or were spilled, which is what keeps every group golden byte-
  // identical (the .bin files held raw doubles). A single-chunk region is the
  // diagonal product alone, i.e. the old `VarMat = P1Mat * P2Mat` on the
  // truncated scratch. The "Analyzing chunks" lines are kept for nchunks > 1
  // (stdout only; LEDGER #57's tell was their count).
  {
    int first_row = 0, first_col = 0, last_row = 0, last_col = 0;
    arma::mat chunkP1, chunkP2;   // load targets for spilled chunks only

    for(unsigned int index1 = 0; index1 < nchunks; index1++)
    {
      last_row = first_row + mPassCVVec.at(index1) - 1;

      const arma::mat & P1c = chunkStore.P1(index1, chunkP1);

      // LEDGER #57: an empty chunk used to `continue` WITHOUT advancing
      // first_row, so every later chunk's submatrix landed in the wrong rows of
      // VarMat -- a wrong answer rather than a missing one. The cursor is
      // advanced on the skip now. Unreachable in practice (mPassCVVec never
      // holds 0), but it is the kind of thing that must not be left armed.
      if(P1c.n_cols == 0){
        first_row = last_row + 1;
        first_col = 0;
        continue;
      }

      // off-diagonal sub-matrices
      for(unsigned int index2 = 0; index2 < index1; index2++)
      {
        std::cout << "Analyzing chunks (" << index1 << "/" << nchunks - 1 << ", " << index2 << "/" << nchunks - 1 << ")........" << std::endl;

        const arma::mat & P2c = chunkStore.P2(index2, chunkP2);

        // Computed BEFORE the skip below: the column cursor has to advance past
        // this chunk either way, for the same reason as the row cursor above.
        last_col = first_col + mPassCVVec.at(index2) - 1;

        if(P2c.n_cols == 0){
          first_col = last_col + 1;
          continue;
        }
        arma::mat offVarMat = P1c * P2c;

        VarMat.submat(first_row, first_col, last_row, last_col) = offVarMat;
        VarMat.submat(first_col, first_row, last_col, last_row) = offVarMat.t();
        first_col = last_col + 1;
      }

      // diagonal sub-matrix
      last_col = first_col + mPassCVVec.at(index1) - 1;
      if(nchunks > 1){
        std::cout << "Analyzing chunks (" << index1 << "/" << nchunks - 1 << ", " << index1 << "/" << nchunks - 1 << ")........" << std::endl;
      }
      const arma::mat & P2d = chunkStore.P2(index1, chunkP2);

      arma::mat diagVarMat = P1c * P2d;

      VarMat.submat(first_row, first_col, last_row, last_col) = diagVarMat;

      first_row = last_row + 1;
      first_col = 0;
      Rcpp::checkUserInterrupt();
    }
    chunkStore.clear();   // removes any spilled files; a no-op in memory
  }
}  //if(t_regionTestType != "BURDEN"){

//read and test single markers done
//arma::vec timeoutput2 = getTime();
//printTime(timeoutput1, timeoutput2, "read and test single markers done");
//check the max MAF of all markers
//if there are sets containing the same markers
ProfileScope rpEmit(rp.enabled ? &rp.emit : nullptr, &rp.nEmit);
arma::uvec q_maf_for_anno(q_anno);
for(unsigned int j = 0; j < q_anno; j++){
	arma::uvec jtemp = find(maxMAFVec >= maxMAFperAnno(j));
	q_maf_for_anno(j) = jtemp.min();
}

//If only conduct Burden test
//arma::vec BURDEN_pval_Vec(q_anno_maf);
std::vector<std::string> BURDEN_pval_Vec(q_anno_maf_weight, "NA");
//BURDEN_pval_Vec.fill(-1.0);
//arma::vec BURDEN_pval_cVec(q_anno_maf);
std::vector<std::string> BURDEN_pval_cVec(q_anno_maf_weight, "NA");
//BURDEN_pval_cVec.fill(-1.0);
std::vector<std::string> BURDEN_AnnoName_Vec(q_anno_maf_weight);
std::vector<std::string> BURDEN_maxMAFName_Vec(q_anno_maf_weight);
std::vector<double> BURDEN_Beta_Vec(q_anno_maf_weight);
std::vector<double> BURDEN_seBeta_Vec(q_anno_maf_weight);
//std::vector<double> BURDEN_pval_cVec(q_anno_maf);
std::vector<double> BURDEN_Beta_cVec(q_anno_maf_weight);
std::vector<double> BURDEN_seBeta_cVec(q_anno_maf_weight);



bool iswriteOutput = false;
bool isregion = true;
if(!g_model->m_flagSparseGRM){
	isregion = false;
}
//Rcpp::DataFrame OUT_BURDEN = Rcpp::DataFrame::create();
unsigned int i= 0;
unsigned int q_maf_m;
std::string AnnoName;
double maxMAFName;
if(t_regionTestType == "BURDEN"){
     // Carry-forward cache for maxMAF cutoffs beyond q_maf_for_anno(j): those
     // cutoffs add no markers, so the last tested (j,*,r) result is re-emitted.
     // With the weight dimension the cache and the isPolyMarker flag are per
     // set r; at q_weight == 1 they reproduce the legacy scalar flow exactly.
     std::vector<bool> isPolyMarkerVec(q_weight, true);
     std::vector<std::string> carry_pval(q_weight), carry_pval_c(q_weight);
     std::vector<double> carry_Beta(q_weight), carry_seBeta(q_weight);
     std::vector<double> carry_Beta_c(q_weight), carry_seBeta_c(q_weight);
     for(unsigned int j = 0; j < q_anno; j++){
       q_maf_m = q_maf_for_anno(j);
       AnnoName = annoStringVec[j];
       for(unsigned int r = 0; r < q_weight; r++){
         isPolyMarkerVec[r] = true;
       }
       for(unsigned int m = 0; m < q_maf; m++){
	maxMAFName = maxMAFVec(m);
	jm = j*q_maf+m;
	for(unsigned int r = 0; r < q_weight; r++){
	jmr = jm*q_weight + r;
	i = jmr;
      if(m <= q_maf_m){
        arma::vec genoSumVec = genoSumMat.col(jmr);
        int n = genoSumVec.size();
        arma::uvec indexNonZeroVec_arma = arma::find(genoSumVec != 0);
	arma::uvec indexZeroVec_arma = arma::find(genoSumVec == 0);
        //double altFreq = arma::mean(genoSumVec)/2;
        //double altCounts = arma::accu(genoSumVec);
	double altCounts =  genoSumcount_noweight(jm);   // weight-independent count
	double altFreq = altCounts/(2*t_n);
	double MAC, MAF;
	if(altFreq > 1){
		MAF = 1;
		MAC = t_n;
	}else{
		MAF = std::min(altFreq, 1 - altFreq);
		MAC = MAF*2*t_n;
	}
	//std::cout << "altCounts " << altCounts << std::endl;
	//std::cout << "altFreq " << altCounts << std::endl;
        double missingRate = 0;
        double imputeInfo = 1;
        std::string chr, ref, alt, marker;
        bool flip = false;
        std::string info = "UR";
        //double MAF = std::min(altFreq, 1 - altFreq);
        double w0;
        //double MAC = MAF * 2 * t_n * (1 - missingRate);
        if(indexNonZeroVec_arma.n_elem > 0 && MAC >= g_cfg.min_gourpmac_for_burdenonly){
	  //if(MAC >= g_cfg.min_gourpmac_for_burdenonly){
          isPolyMarkerVec[r] = true;
	  std::vector<uint32_t> indexForMissing;

      // S12: this marker's test context, derived and consumed here. Same two
      // steps as before -- a common marker uses the non-sparse variance ratio,
      // a rare one follows the model -- but nothing is parked on g_model, so no
      // later read can pick up a marker's choice by accident (LEDGER #27/#50)
      // and two markers could be tested at once without racing.
      SAIGE::TestContext ctx;
      ctx.useSparseGRM = (MAC > g_model->m_cateVarRatioMinMACVecExclude.back())
                           ? false : g_model->m_flagSparseGRM;
      ctx.noadjCov = false;
      if(!isSingleVarianceRatio){
        ctx.varRatioVal = g_model->selectVarianceRatio(MAC, ctx.useSparseGRM, false, hasVarRatio);
      }else{
        ctx.varRatioVal = g_model->selectSingleVarianceRatio(ctx.useSparseGRM, false);
      }
          //if(!isSingleVarianceRatio){
           // (the duplicate variance-ratio selection this used to repeat is above)
          //}


	  //arma::vec timeoutput_getp = getTime();
	  ProfileScope rpEmitTest(rp.enabled ? &rp.emitTest : nullptr, &rp.nEmitTest);
	  if(MAC <= g_cfg.MACCutoffforER && t_traitType == "binary" && erAllowedVec[r] ){
          g_model->getMarkerPval(genoSumVec, indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, altFreq, Tstat, gy, varT, isSPAConverge, gtildeVec, is_gtilde, isregion, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, true, flip, ctx);
	  }else{
          g_model->getMarkerPval(genoSumVec, indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, altFreq, Tstat, gy, varT, isSPAConverge, gtildeVec, is_gtilde, isregion, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, false, flip, ctx);
	  }	  
	  rpEmitTest.stop();
	  //arma::vec timeoutput_getp2 = getTime();
	  //printTime(timeoutput_getp, timeoutput_getp2, "get p  done");
	  
	  if(isCondition){
	    BURDEN_pval_cVec.at(i) = pval_c;
	    BURDEN_Beta_cVec.at(i) = Beta_c;
	    BURDEN_seBeta_cVec.at(i) = seBeta_c;
          }

	   BURDEN_AnnoName_Vec.at(i) = AnnoName;
	   std::string str = std::to_string(maxMAFName);
	   str.erase ( str.find_last_not_of('0') + 1, std::string::npos );
	   BURDEN_maxMAFName_Vec.at(i) = str;
	   BURDEN_pval_Vec.at(i) = pval;
	   BURDEN_Beta_Vec.at(i) = Beta;
	   BURDEN_seBeta_Vec.at(i) = seBeta;

	   // cache this set's result for the m > q_maf_for_anno(j) carry-forward
	   carry_pval[r] = pval;
	   carry_Beta[r] = Beta;
	   carry_seBeta[r] = seBeta;
	   if(isCondition){
	     carry_pval_c[r] = pval_c;
	     carry_Beta_c[r] = Beta_c;
	     carry_seBeta_c[r] = seBeta_c;
	   }

         }else{
	   isPolyMarkerVec[r] = false;
	 }
     }else{
	if(isPolyMarkerVec[r]){
	  if(isCondition){
            BURDEN_pval_cVec.at(i) = carry_pval_c[r];
            BURDEN_Beta_cVec.at(i) = carry_Beta_c[r];
            BURDEN_seBeta_cVec.at(i) = carry_seBeta_c[r];
          }

           BURDEN_AnnoName_Vec.at(i) = AnnoName;
           BURDEN_maxMAFName_Vec.at(i) = std::to_string(maxMAFName);
           BURDEN_pval_Vec.at(i) = carry_pval[r];
           BURDEN_Beta_Vec.at(i) = carry_Beta[r];
           BURDEN_seBeta_Vec.at(i) = carry_seBeta[r];
	 }
     }
	} // for r (weight sets)

   }
 }
	  std::vector<double> nonMissingPvalVec_std, nonMissingPvalVec_c_std;
	  double burden_p, burden_p_cond;
	  for(unsigned int i = 0; i < BURDEN_pval_Vec.size(); i++){
		if(BURDEN_pval_Vec.at(i) != "NA"){
	  		// LEDGER #72 FIXED: the same guarded parse the marker path has had
	  		// (Main.cpp, the fastTest re-test): a subnormal printed p (below
	  		// DBL_MIN, which std::stod refuses as out_of_range) used to abort the
	  		// whole run from inside mainRegionInCPP. It is taken as 0, the value
	  		// the marker path uses and the limit the Cauchy combination wants.
	  		burden_p = parsePrintedPval(BURDEN_pval_Vec.at(i));
			nonMissingPvalVec_std.push_back(burden_p);
			if(isCondition){
				burden_p_cond = parsePrintedPval(BURDEN_pval_cVec.at(i));
				nonMissingPvalVec_c_std.push_back(burden_p_cond);
			}
		}	
	  } 

	  arma::vec nonMissingPvalVec = arma::conv_to< arma::vec >::from(nonMissingPvalVec_std);
	  //arma::uvec nonMissingPvalVecInd = arma::find(BURDEN_pval_Vec >= 0);
	  //arma::vec nonMissingPvalVec = BURDEN_pval_Vec.elem(nonMissingPvalVecInd);
	   cctpval = CCT_cpp(nonMissingPvalVec);
           if(isCondition){
	        //arma::vec nonMissingPvalVec_cond = BURDEN_pval_cVec.elem(nonMissingPvalVecInd);
		arma::vec nonMissingPvalVec_cond  = arma::conv_to< arma::vec >::from(nonMissingPvalVec_c_std);
	   	cctpval_cond = CCT_cpp(nonMissingPvalVec_cond);	   
           }

	
	   if(!t_isFastTest){
		iswriteOutput = true;
	   }else{
		//LEDGER #27: the caller's intent at entry, not the mutated flag
		if(isEscalatedPass){
			iswriteOutput = true;
		}else{
		  if(cctpval >= 0.1){
			iswriteOutput = true;
		  } 	
		}
	   }
	
	   if(iswriteOutput){
		writeOutfile_BURDEN(regionName,
			BURDEN_AnnoName_Vec,
			BURDEN_maxMAFName_Vec,
			BURDEN_pval_Vec,
			BURDEN_Beta_Vec,
			BURDEN_seBeta_Vec,
			BURDEN_pval_cVec,
			BURDEN_Beta_cVec,
			BURDEN_seBeta_cVec,
			MAC_GroupVec,
			MACCase_GroupVec,
			MACControl_GroupVec,
			NumRare_GroupVec,
			NumUltraRare_GroupVec,
			cctpval,
			cctpval_cond,
			q_anno,
			q_maf,
			q_weight,
			t_weightNames,
			emitWeightColumn,
			isCondition,
			t_traitType);
	   }else{
	     OutList.push_back(iswriteOutput, "iswriteOutput");
	   }

 //arma::vec timeoutput3 = getTime();
 //printTime(timeoutput2, timeoutput3, "burden test done");
 }else{
  q_maf_for_anno = q_maf_for_anno + 1;
  OutList.push_back(MAC_GroupVec, "MAC_GroupVec");
  OutList.push_back(q_maf_for_anno, "q_maf_for_annoVec");
  if(t_traitType == "binary"){
   if(t_traitType == "binary"){
    OutList.push_back(MACCase_GroupVec, "MACCase_GroupVec");
    OutList.push_back(MACControl_GroupVec, "MACCtrl_GroupVec");
   }

    OutList.push_back(genoSumMat, "genoSumMat");
    OutList.push_back(gyVec, "gyVec");
  }

    OutList.push_back(VarMat, "VarMat");	
    OutList.push_back(MAFVec, "MAFVec");	
    OutList.push_back(TstatVec_flip, "TstatVec_flip");	
  //arma::mat scaled_m_VarInvMat_cond;
    if(isCondition){
  //std::cout << "okk5" << std::endl;
     if(q_weight == 1){
      arma::mat AdjCondMat = G1tilde_P_G2tilde_Weighted_Mat * (g_model->m_VarInvMat_cond / (w0G2Mat_cond));
      arma::mat VarMatAdjCond = AdjCondMat * (G1tilde_P_G2tilde_Weighted_Mat.t());
      arma::vec TstatAdjCond = AdjCondMat * (g_model->m_Tstat_cond % w0G2Vec_cond );
      OutList.push_back(G1tilde_P_G2tilde_Weighted_Mat, "G1tilde_P_G2tilde_Weighted_Mat");
      OutList.push_back(g_model->m_scalefactor_G2_cond, "scalefactor_G2_cond");
      OutList.push_back(g_model->m_VarInvMat_cond_scaled_weighted, "VarInvMat_G2_cond_scaled");
      OutList.push_back(g_model->m_Tstat_cond, "Tstat_G2_cond"); //m_Tstat_cond is weighted
      OutList.push_back(g_model->m_G2_Weight_cond, "G2_Weight_cond");
      OutList.push_back(TstatAdjCond, "TstatAdjCond");
      OutList.push_back(VarMatAdjCond, "VarMatAdjCond");
     }else{
      // design D7 (C11): at K>1 the baked TstatAdjCond/VarMatAdjCond would
      // carry set 0's weights only, so they are not computed here.
      // G1tilde_P_G2tilde_Weighted_Mat is the G1-raw matrix (see the fill
      // sites); region_postprocess forms, per set r,
      //   G1w_r = raw[cell rows] * AnnoWeights_r
      //   AdjCondMat_r  = G1w_r * VarInvMat_G2_cond_w0
      //   VarMatAdjCond_r = AdjCondMat_r * G1w_r'
      //   TstatAdjCond_r  = AdjCondMat_r * Tstat_G2_cond_w0
      // -- the literal q_weight == 1 arithmetic above with w0 moved to R.
      arma::mat VarInvMat_G2_cond_w0 = g_model->m_VarInvMat_cond / (w0G2Mat_cond);
      arma::vec Tstat_G2_cond_w0 = g_model->m_Tstat_cond % w0G2Vec_cond;
      OutList.push_back(G1tilde_P_G2tilde_Weighted_Mat, "G1tilde_P_G2tilde_Weighted_Mat");
      OutList.push_back(g_model->m_scalefactor_G2_cond, "scalefactor_G2_cond");
      OutList.push_back(g_model->m_VarInvMat_cond_scaled_weighted, "VarInvMat_G2_cond_scaled");
      OutList.push_back(g_model->m_Tstat_cond, "Tstat_G2_cond");
      OutList.push_back(g_model->m_G2_Weight_cond, "G2_Weight_cond");
      OutList.push_back(VarInvMat_G2_cond_w0, "VarInvMat_G2_cond_w0");
      OutList.push_back(Tstat_G2_cond_w0, "Tstat_G2_cond_w0");
     }
    }

    if(!t_isFastTest){
            iswriteOutput = true;
     }else{
        // LEDGER #50 fix, the non-BURDEN twin of #27 at :2240. This read used to
        // be g_model->m_flagSparseGRM_cur, which the marker and ultra-rare loops
        // above overwrite from each marker's MAC -- so by the time the write
        // decision runs, R's "this is the escalated pass, DO write" signal has
        // been replaced by the last pseudo-marker's variance-ratio choice. On the
        // fast pass R takes the rows from the temp bin itself (copy_singleInGroup),
        // so the escalated pass is the ONLY route by which an escalated region's
        // single-variant rows can be written, and they were silently dropped.
        if(isEscalatedPass){
            iswriteOutput = true;
        }
     }
}


 int numofUR = q_anno_maf_weight;
 int numofUR0;
 int mFirth = 0;
if(t_isSingleinGroupTest){
 // OutList.push_back(pvalVec_val, "pvalVec");
 OutList.push_back(pvalVec, "pvalVec");
if(iswriteOutput){
  numofUR0 = writeOutfile_singleInGroup(t_isMoreOutput,
      t_isImputation,
      isCondition,
      is_Firth,
      mFirth,
      is_FirthConverge,
      t_traitType,
      chrVec,
      posVec,
      markerVec,
      refVec,
      altVec,
      altCountsVec,
      altFreqVec,
      imputationInfoVec,
      missingRateVec,
      BetaVec,
      seBetaVec,
      TstatVec,
      varTVec,
      pvalVec,
      pvalNAVec,
      isSPAConvergeVec,
      Beta_cVec,
      seBeta_cVec,
      Tstat_cVec,
      varT_cVec,
      pval_cVec,
      pvalNA_cVec,
      AF_caseVec,
      AF_ctrlVec,
      N_caseVec,
      N_ctrlVec,
      N_case_homVec,
      N_ctrl_hetVec,
      N_case_hetVec,
      N_ctrl_homVec,
      N_Vec,
      g_out.singleInGroup);
 }else{
  g_out.singleInGroupTemp.open(g_out.prefixSingleInGroupTemp.c_str(), std::ofstream::out);
  numofUR0 = writeOutfile_singleInGroup(t_isMoreOutput,
      t_isImputation,
      isCondition,
      is_Firth,
      mFirth,
      is_FirthConverge,
      t_traitType,
      chrVec,
      posVec,
      markerVec,
      refVec,
      altVec,
      altCountsVec,
      altFreqVec,
      imputationInfoVec,
      missingRateVec,
      BetaVec,
      seBetaVec,
      TstatVec,
      varTVec,
      pvalVec,
      pvalNAVec,
      isSPAConvergeVec,
      Beta_cVec,
      seBeta_cVec,
      Tstat_cVec,
      varT_cVec,
      pval_cVec,
      pvalNA_cVec,
      AF_caseVec,
      AF_ctrlVec,
      N_caseVec,
      N_ctrlVec,
      N_case_homVec,
      N_ctrl_hetVec,
      N_case_hetVec,
      N_ctrl_homVec,
      N_Vec,
      g_out.singleInGroupTemp);
  g_out.singleInGroupTemp.close();

 }//iswriteOutput
 OutList.push_back(numofUR, "numofUR");
}

OutList.push_back(iswriteOutput, "iswriteOutput");

 if(t_isOutputMarkerList){
	OutList.push_back(indicatorVec, "markerIndcatorVec");
 }

  

  OutList.push_back(NumRare_GroupVec, "NumRare_GroupVec");
  OutList.push_back(NumUltraRare_GroupVec, "NumUltraRare_GroupVec");


 if(t_regionTestType != "BURDEN" || t_isOutputMarkerList){
  OutList.push_back(annoMAFIndicatorMat, "annoMAFIndicatorMat");
 }

  if(rp.enabled){
    rpEmit.stop();
    rpTotal.stop();
    rp.calls += 1;
    // one cumulative line per call; the last one in a log is the run's total
    std::cout << "REGION PROFILE cum (SAIGE_PROFILE_REGION set):"
              << " calls=" << rp.calls
              << " markers=" << rp.nPrepare
              << " tested=" << rp.nTest
              << " urCollapsed=" << rp.nUrCollapse
              << " urTested=" << rp.nUrTest
              << " burdenTested=" << rp.nEmitTest
              << " | prepare=" << rp.prepare
              << " rareArm=" << rp.rareArm
              << " test=" << rp.test
              << " rareOther=" << (rp.rareArm - rp.test)
              << " urCollapse=" << rp.urCollapse
              << " flush=" << rp.flush
              << " urPass=" << rp.urPass
              << " urTest=" << rp.urTest
              << " varmat=" << rp.varmat
              << " emit=" << rp.emit
              << " emitTest=" << rp.emitTest
              << " total=" << rp.total
              << std::endl;
  }

  return OutList;
}








// [[Rcpp::export]]
void assign_conditionMarkers_factors(
                           std::string t_genoType,     // "plink"
                           std::vector<std::string> & t_genoIndex,
                           unsigned int t_n, 
			   arma::vec & t_weight_cond,
			   int t_modelIndex
			   )           // sample size
{
  if(t_modelIndex < 0 || (std::size_t)t_modelIndex >= g_models.size()){
    Rcpp::stop("assign_conditionMarkers_factors: no null model at index " +
               std::to_string(t_modelIndex));
  }
  // Read through THIS model's sample-set group, and write the conditioning
  // factors onto THIS model. A single-model run passes 0 and lands exactly
  // where it did before M1.
  SAIGE::SAIGEClass & model = *g_models[t_modelIndex];
  Geno::GenoReader & reader = *g_readers.at(g_modelGroup.at(t_modelIndex));
  bool isImpute = false;	
  unsigned int q = t_genoIndex.size();
  arma::mat P1Mat(q, t_n);
  arma::mat P2Mat(t_n, q);
  arma::mat G2tildeMat(t_n, q); // LEDGER #84: the adjusted genotypes themselves, unscaled
  arma::mat VarInvMat(q, q);
  arma::vec TstatVec(q);
  //arma::vec pVec(q);
    std::vector<std::string> pVec(q, "NA");
  arma::vec MAFVec(q);
  arma::vec gyVec(q);
  arma::vec w0G2_cond_Vec(q);
  arma::vec gsumVec(t_n, arma::fill::zeros);
  //double beta1 = g_cfg.weights_beta[0];
  //double beta2 = g_cfg.weights_beta[1];  
  boost::math::beta_distribution<> beta_dist(g_cfg.weights_beta[0], g_cfg.weights_beta[1]);
  //boost::math::beta_distribution<> beta_dist(beta1, beta2);
  //g_cfg.weights_beta.print();
  //boost::math::beta_distribution<> beta_dist(1, 25);
  //std::vector<double> GVec0(t_n);
  arma::vec GVec(t_n);
  std::string pval, pval_noSPA;
  double Beta, seBeta, Tstat, varT, gy, w0G2_cond;
  bool isSPAConverge, is_gtilde, is_Firth, is_FirthConverge;
  arma::vec P2Vec(t_n);

  //std::vector<uint> indexZeroVec;
  //std::vector<uint> indexNonZeroVec;
  std::string pval_c, pval_noSPA_c;
  double Beta_c, seBeta_c, Tstat_c, varT_c;
  arma::rowvec G1tilde_P_G2tilde_Vec;
  bool isCondition = false;
  for(unsigned int i = 0; i < q; i++)
  {
    // marker-level information
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing;
    std::vector<uint> indexZeroVec;
    std::vector<uint> indexNonZeroVec;
    std::string chr, ref, alt, marker;
    uint32_t pd;
    bool flip = false;

    bool isOutputIndexForMissing = true;
    bool isOnlyOutputNonZero = false; 




    std::string t_genoIndex_str = t_genoIndex.at(i);
    char* end;
    uint64_t gIndex = std::strtoull( t_genoIndex_str.c_str(), &end,10 );


    bool isReadMarker = getOneMarkerThrough(reader, t_genoType, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          isOutputIndexForMissing, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          isOnlyOutputNonZero, // bool t_isOnlyOutputNonZero,
                                          indexNonZeroVec, GVec, isImpute);
     //arma::vec GVec(GVec0);
     //GVec0.clear();
    if(!isReadMarker){
      break;
    }

    std::string info = chr+":"+std::to_string(pd)+"_"+ref+"/"+alt;

  double MAF = std::min(altFreq, 1 - altFreq);
  double MAC = MAF * 2 * t_n * (1 - missingRate);

  // S12: the conditioning markers get their own context too. This function used
  // to set m_flagSparseGRM_cur at entry and leave it set for whatever ran next;
  // now the model's own m_flagSparseGRM is read directly and the choice does not
  // outlive the marker.
  bool hasVarRatio;
  SAIGE::TestContext ctx;
  ctx.useSparseGRM = model.m_flagSparseGRM;
  ctx.noadjCov = false;
  if((model.m_varRatio_null).n_elem == 1){
	ctx.varRatioVal = model.selectSingleVarianceRatio(ctx.useSparseGRM, false);
  }else{
	ctx.varRatioVal = model.selectVarianceRatio(MAC, ctx.useSparseGRM, false, hasVarRatio);
  }	  
 
  flip = imputeGenoAndFlip(GVec, altFreq, altCounts, indexForMissing, g_cfg.impute_method, g_cfg.dosage_zerod_cutoff, g_cfg.dosage_zerod_MAC_cutoff, MAC, indexZeroVec, indexNonZeroVec);


 arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;
       indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);
       indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);


  MAF = std::min(altFreq, 1 - altFreq);
  MAC = std::min(altCounts, 2*t_n-altCounts);

  // LEDGER #73 FIXED. A conditioning marker is named by the user and bypassed
  // every check, so a marker MONOMORPHIC in the analysed samples -- rs1 of the
  // bundled fixture -- was conditioned on anyway and the conditional columns
  // came out IDENTICAL to the marginal ones on 70/70 rows at exit 0, with only
  // a std::cerr line for MAF == 0 and then on. That case is a silent no-op
  // labelled as an adjustment, and it STOPS now, by name.
  //
  // A marker that merely fails the SCAN's thresholds is a different matter:
  // golden P4 conditions on rs13 (MAC 2) under --minMAC=20, and conditioning on
  // a marker rarer than the markers being scanned is a legitimate design (the
  // pre-push gate is what caught a first draft that refused it). So the scan's
  // --minMAC/--minMAF/--maxMissing/--minInfo are reported as an R warning, and
  // the marker is used as before.
  if(MAC == 0){
    Rcpp::stop("conditioning marker " + marker + " (" + info + ") is monomorphic in the analysed samples"
               " (MAC 0): conditioning on it would change nothing while the _c columns claim otherwise."
               " It cannot be conditioned on.\n");
  }
  if((missingRate > g_cfg.missingRate_cutoff) || (MAF < g_cfg.marker_minMAF_cutoff) ||
     (MAC < g_cfg.marker_minMAC_cutoff) || (imputeInfo < g_cfg.marker_minINFO_cutoff)){
    Rcpp::warning("conditioning marker " + marker + " (" + info + ") would be dropped by the scan's own QC"
                  " -- MAC " + std::to_string(MAC) + ", MAF " + std::to_string(MAF) +
                  ", missing rate " + std::to_string(missingRate) + ", imputation info " + std::to_string(imputeInfo) +
                  " against --minMAC " + std::to_string(g_cfg.marker_minMAC_cutoff) +
                  ", --minMAF " + std::to_string(g_cfg.marker_minMAF_cutoff) +
                  ", --maxMissing " + std::to_string(g_cfg.missingRate_cutoff) +
                  ", --minInfo " + std::to_string(g_cfg.marker_minINFO_cutoff) +
                  " -- and is conditioned on anyway.");
  }
  arma::vec gtildeVec;

  if(MAC <=  g_cfg.MACCutoffforER && model.m_traitType == "binary"){
     // model.getMarkerPval rather than the Unified_ wrapper: the wrapper
     // dispatches to g_models[0], which is the wrong model whenever this is
     // called for the second or later entry of a manifest. Argument order
     // follows SAIGEClass::getMarkerPval (altFreq before Tstat), which is the
     // one place it differs from the wrapper's.
     model.getMarkerPval(
                    GVec,
                    indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, altFreq, Tstat, gy, varT, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, true, flip, ctx);
  }else{
     // model.getMarkerPval rather than the Unified_ wrapper: the wrapper
     // dispatches to g_models[0], which is the wrong model whenever this is
     // called for the second or later entry of a manifest. Argument order
     // follows SAIGEClass::getMarkerPval (altFreq before Tstat), which is the
     // one place it differs from the wrapper's.
     model.getMarkerPval(
                    GVec,
                    indexNonZeroVec_arma, indexZeroVec_arma, Beta, seBeta, pval, pval_noSPA, altFreq, Tstat, gy, varT, isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c, G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge, false, flip, ctx);

  }
      P1Mat.row(i) = sqrt(ctx.varRatioVal)*gtildeVec.t();
      P2Mat.col(i) = sqrt(ctx.varRatioVal)*P2Vec;
      G2tildeMat.col(i) = gtildeVec;
      //P1Mat.row(i) = gtildeVec.t();
      //P2Mat.col(i) = P2Vec;
      MAFVec(i) = MAF;
      //w0G2_cond = boost::math::pdf(beta_dist, MAF);
      //t_weight_cond.print();
      // (the former std::cerr "Conditioning marker is monomorphic" notice is
      // subsumed by the LEDGER #73 stop above)
     if(!t_weight_cond.is_zero()){
	 w0G2_cond = t_weight_cond(i);
    }else{
	 w0G2_cond = boost::math::pdf(beta_dist, MAF);
    }
     w0G2_cond_Vec(i) = w0G2_cond;
     gyVec(i) = gy * w0G2_cond;
     gsumVec = gsumVec + GVec * w0G2_cond;
     TstatVec(i) = Tstat;
     //pVec(i) = std::stod(pval);
     pVec.at(i) = pval;
  }
  arma::mat VarMat = P1Mat * P2Mat;

  VarInvMat =  arma::pinv(VarMat);
  //VarInvMat = VarMat.i();
  double qsum = arma::accu(gyVec);
  arma::vec gsumtildeVec;


  model.getadjG(gsumVec, gsumtildeVec);
  model.assignConditionFactors(
		   			P2Mat,
					G2tildeMat,
					VarInvMat,
					VarMat,
					TstatVec,
				        w0G2_cond_Vec,	
					MAFVec,
					qsum,
					gsumtildeVec,
					pVec);

}

// [[Rcpp::export]]
void assign_conditionMarkers_factors_binary_region(
			   arma::vec & scalefactor_G2_cond){
	//std::cout << "assign_conditionMarkers_factors_binary_region" << std::endl;
	g_model->assignConditionFactors_scalefactor(scalefactor_G2_cond);
}


// LEDGER #20 FIXED (2026-09-05): the logistf-shaped parameters this port never
// read -- col_fit (column subsetting), maxhs and lconv (step-halving) -- are
// gone from the signature rather than accepted and ignored. Nothing in the
// package calls this export; it is public API via exportPattern, characterised
// in tests/testthat/test-unit-firth.R. The loop bound is `iter < maxit` so
// that maxit = 0 returns `init` untouched instead of taking one Newton step;
// for maxit >= 1 the iteration count is unchanged (the post-increment
// `iter == maxit` break fires at the same step).
// [[Rcpp::export]]
arma::vec fast_logistf_fit(arma::mat & x,
		arma::vec & y,
		arma::vec & weight,
		arma::vec & offset,
		bool firth,
    	arma::vec init, 
	int maxit, 
	int maxstep, 
	double gconv, 
	double xconv, 
	bool & isfirthconverge){
  int n = x.n_rows;
  int k = x.n_cols;
  arma::vec beta = init;
  int iter = 0;
  arma::vec pi_0 = -x * beta - offset; 
  pi_0 = arma::exp(pi_0) + 1;
  arma::vec pi = 1/pi_0;
  int evals = 1;
  arma::vec beta_old;
  arma::mat oneVec(k, 1 , arma::fill::ones);
  arma::mat XX_covs(k, k, arma::fill::zeros);
  isfirthconverge = false;
  while(iter < maxit){
	beta_old = beta;
	arma::vec wpi = weight % pi % (1 - pi);
	arma::vec wpi_sqrt = arma::sqrt(wpi);
	arma::vec W2 = weight % wpi_sqrt;
	arma::mat XW2(n, k, arma::fill::zeros);
	for(int j = 0; j < k; j++){
                XW2.col(j) = x.col(j) % W2;
        }       

	arma::mat Q;
	arma::mat R;
	arma::qr_econ(Q, R, XW2);
	arma::vec h = Q % Q * oneVec;
	arma::vec U_star(2, arma::fill::zeros);
	arma::vec ypih;
	if(firth){
		ypih = (weight % (y - pi)) + (h % (0.5 - pi));
	}else{
		ypih = (weight % (y - pi));
	}
	//ypih.print();
	arma::vec xcol(n, arma::fill::zeros);
	U_star = x.t() * ypih;
	
	arma::mat XX_XW2(n, k, arma::fill::zeros);
	for(int j = 0; j < k; j++){
		xcol = x.col(j);
		XX_XW2.col(j) = xcol % wpi_sqrt; 	
	}
	arma::mat XX_Fisher = XX_XW2.t() * (XX_XW2);
	bool isinv = arma::inv_sympd (XX_covs, XX_Fisher); 
	if(!isinv){
		break;
	}	
	//}
	arma::vec delta = XX_covs * U_star;
	delta.replace(arma::datum::nan, 0);	

	double mx = arma::max(arma::abs(delta))/maxstep;
	if(mx > 1){
		delta = delta/mx;
	}
	evals = evals + 1;
	iter = iter + 1;
	beta = beta + delta;
	pi_0 = -x * beta - offset;
  	pi_0 = arma::exp(pi_0) + 1;
  	pi = 1/pi_0;
	if((iter == maxit) || ( (arma::max(arma::abs(delta)) <= xconv) & (abs(U_star).is_zero(gconv)))){
		isfirthconverge = true;
		break;
	}
  }
	arma::mat var;
	if(XX_covs.has_nan()){
		var = XX_covs;
		beta = arma::datum::nan;
	}
	return beta;
}

// [[Rcpp::export]]
void closeGenoFile(std::string & t_genoType)
{
  if(t_genoType == "plink" || t_genoType == "pgen"){
    // Every sample-set group's reader, not just the first: with a manifest
    // there is one per group and closing only g_readers[0] would leak the rest.
    for(std::size_t g = 0; g < g_readers.size(); g++){
      g_readers[g]->close();
    }
  }
}

// [[Rcpp::export]]
bool openOutfile(std::string t_traitType, bool isappend, unsigned int t_nWeightSets, bool t_forceWeightColumn){
	bool isopen;
	if(!isappend){
	g_out.group.open(g_out.prefixGroup.c_str());
	isopen = g_out.group.is_open();
	if(isopen){
		g_out.group << "Region\tGroup\tmax_MAF\t";
		// additive multi-weight schema (design D3): the Weight label column
		// sits right after max_MAF and exists only when more than one weight
		// set is in play (or the user forces it at a single set with
		// is_output_weight_column); otherwise the emitted bytes are exactly
		// the legacy header.
		if(t_nWeightSets > 1 || t_forceWeightColumn){
			g_out.group << "Weight\t";
		}
		g_out.group << "Pvalue_Burden\tBETA_Burden\tSE_Burden\t";
		if(g_model->m_isCondition){
			g_out.group << "Pvalue_Burden_c\tBeta_Burden_c\tseBeta_Burden_c\t";
		}
		g_out.group << "MAC\t";
		if(t_traitType == "binary"){	
			g_out.group << "MAC_case\tMAC_control\t";
		}

		g_out.group << "Number_rare\tNumber_ultra_rare\n";
	}
     }else{
	g_out.group.open(g_out.prefixGroup.c_str(), std::ofstream::out | std::ofstream::app) ;				
	isopen = g_out.group.is_open();
     }
     return(isopen);
}

// [[Rcpp::export]]
bool openOutfile_singleinGroup(std::string t_traitType, bool t_isImputation, bool isappend, bool t_isMoreOutput){
        bool isopen;
     if(!isappend){
        g_out.singleInGroup.open(g_out.prefixSingleInGroup.c_str());
	isopen = g_out.singleInGroup.is_open();
        if(isopen){
                g_out.singleInGroup << "CHR\tPOS\tMarkerID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\t";
		if(t_isImputation){
			g_out.singleInGroup << "imputationInfo\t";
		}else{
			
			g_out.singleInGroup << "MissingRate\t";
		}
		g_out.singleInGroup << "BETA\tSE\tTstat\tvar\tp.value\t";
	        if(t_traitType == "binary"){
                        g_out.singleInGroup << "p.value.NA\tIs.SPA\t";
                }	

                if(g_model->m_isCondition){
                        g_out.singleInGroup << "BETA_c\tSE_c\tTstat_c\tvar_c\tp.value_c\t";
			if(t_traitType == "binary"){
				g_out.singleInGroup << "p.value.NA_c\t";		
			}
                }
		
	        if(t_traitType == "binary"){
                        g_out.singleInGroup << "AF_case\tAF_ctrl\tN_case\tN_ctrl";
 			if(t_isMoreOutput){
                                g_out.singleInGroup << "\tN_case_hom\tN_case_het\tN_ctrl_hom\tN_ctrl_het";
                        }
			g_out.singleInGroup << "\n";
                }else if(t_traitType == "quantitative"){
			g_out.singleInGroup << "N\n";	
			
		}

        }
    }else{
       g_out.singleInGroup.open(g_out.prefixSingleInGroup.c_str(), std::ofstream::out | std::ofstream::app);
       isopen = g_out.singleInGroup.is_open();
    }

    //if(t_testType != "BURDEN" && t_isfastTest){
    //  g_out.singleInGroupTemp.open(g_out.prefixSingleInGroupTemp.c_str(), std::ofstream::out | std::ofstream::app);
    //}
      return(isopen);
}


// Open (or reopen for append) ONE model's single-variant output file. The
// header is derived from THAT model rather than passed in, because with several
// models in one scan the trait type and the conditional columns differ per file.
// [[Rcpp::export]]
bool openOutfile_single(int t_modelIndex, bool t_isImputation, bool isappend, bool t_isMoreOutput){
      if(t_modelIndex < 0 || (std::size_t)t_modelIndex >= g_out.singles.size()){
        Rcpp::stop("openOutfile_single: no output file for model index " + std::to_string(t_modelIndex));
      }
      if((std::size_t)t_modelIndex >= g_models.size()){
        Rcpp::stop("openOutfile_single: no null model at index " + std::to_string(t_modelIndex));
      }
      SingleWriter & w = *g_out.singles[t_modelIndex];
      const SAIGE::SAIGEClass & model = *g_models[t_modelIndex];
      const std::string t_traitType = model.m_traitType;
      bool isopen;
      if(!isappend){
        w.stream.open(w.prefix.c_str());
        isopen = w.stream.is_open();
        if(isopen){
                w.stream << "CHR\tPOS\tMarkerID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\t";
                if(t_isImputation){
                        w.stream << "imputationInfo\t";
                }else{

                        w.stream << "MissingRate\t";
                }
                w.stream << "BETA\tSE\tTstat\tvar\tp.value\t";
                if(t_traitType == "binary"){
                        w.stream << "p.value.NA\tIs.SPA\t";
                }

                if(model.m_isCondition){
                        w.stream << "BETA_c\tSE_c\tTstat_c\tvar_c\tp.value_c\t";
			if(t_traitType == "binary"){
				w.stream << "p.value.NA_c\t";
			}



                }
		

                if(t_traitType == "binary"){
                        w.stream << "AF_case\tAF_ctrl\tN_case\tN_ctrl";

			if(t_isMoreOutput){	
				w.stream << "\tN_case_hom\tN_case_het\tN_ctrl_hom\tN_ctrl_het";
			}
			w.stream << "\n";

                }else if(t_traitType == "quantitative"){
                        w.stream << "N\n";

                }

        }
     }else{
        w.stream.open(w.prefix.c_str(), std::ofstream::out | std::ofstream::app);
        isopen = w.stream.is_open();
     }

     return(isopen);
}

// Flush and close the single-variant output stream.
//
// g_out.single is a file-scope std::ofstream. In a serial run it is never
// closed explicitly and survives only because static destruction at normal
// process exit flushes it. That is not available to a forked worker:
// parallel:::mcexit ends a child with _exit(), which runs no static
// destructors, so the child would silently lose everything still sitting in
// the stream buffer (measured: the tail past the 192 KiB buffer boundary,
// with the child still reporting success). Serially this call is
// byte-neutral -- it only moves the same flush earlier.
// [[Rcpp::export]]
void closeOutfile_single(){
  for(std::size_t k = 0; k < g_out.singles.size(); k++){
    if(g_out.singles[k]->stream.is_open()){
      g_out.singles[k]->stream.close();
    }
  }
}

// Flush the single-variant output stream, so that everything written so far is
// on disk rather than sitting in the stream's buffer.
//
// LEDGER #47: the `.index` progress file and this stream flush at INDEPENDENT
// points -- the index once per chunk from R, the stream whenever its ~192 KiB
// buffer fills -- and those boundaries do not align. A killed run therefore
// leaves rows on disk belonging to chunks LATER than the last indexed one, and
// `--is_overwrite_output=FALSE` resumes at the indexed chunk and appends them a
// second time (measured: 3,826 duplicated markers). The chunk loop calls this
// at each boundary so R can record a byte count that means exactly "everything
// on disk belongs to chunks <= i"; a restart truncates the file to it.
//
// Returns false if the stream is not open, so the caller can tell "nothing to
// flush" from "flushed zero bytes".
// [[Rcpp::export]]
bool flushOutfile_singleInCPP(){
  if(g_out.singles.empty()){
    return false;
  }
  bool allOpen = true;
  for(std::size_t k = 0; k < g_out.singles.size(); k++){
    if(!g_out.singles[k]->stream.is_open()){
      allOpen = false;
      continue;
    }
    g_out.singles[k]->stream.flush();
  }
  return allOpen;
}

// Re-point the single-variant output files. One prefix per null model, in
// manifest order; the first is the one setAssocTest_GlobalVarsInCPP already
// installed, and a single-model run never calls this at all (so its file naming
// is byte- and name-identical to before M1). This REPLACES the fork-era
// setSingleOutputPrefixInCPP, which had no R caller.
// [[Rcpp::export]]
void setSingleOutputPrefixesInCPP(std::vector<std::string> t_outputFilePrefixes){
  if(t_outputFilePrefixes.empty()){
    Rcpp::stop("At least one single-variant output prefix is required.");
  }
  g_out.singles.clear();
  for(std::size_t k = 0; k < t_outputFilePrefixes.size(); k++){
    g_out.singles.push_back(std::unique_ptr<SingleWriter>(new SingleWriter()));
    g_out.singles[k]->prefix = t_outputFilePrefixes[k];
  }
}

// Whether the active genotype reader may be used from forked children. Nothing
// in R calls this any more -- S7 retired the forked path -- but it remains the
// reader's answer to a real question and is exercised by
// tests/testthat/test-unit-genoreader.R.
// [[Rcpp::export]]
bool isGenoReaderForkSafeInCPP(){
  if(!g_reader){
    Rcpp::stop("No genotype reader has been set up.");
  }
  return g_reader->isForkSafe();
}

// Whether the active reader can hand out per-worker clones, i.e. whether the
// marker scan can actually thread. Both built-in readers do since S8; the one
// refusal left is the PLINK pread fallback (S5). The C++ marker loop degrades
// to one thread on its own if this is false; R asks FIRST only so the user gets
// a real R warning() rather than a line of stdout that scrolls past in a
// cluster log.
// [[Rcpp::export]]
bool isGenoReaderCloneableInCPP(){
  if(g_readers.empty()){
    Rcpp::stop("No genotype reader has been set up.");
  }
  // EVERY sample-set group's reader. One uncloneable reader makes the whole
  // scan serial, so "can this scan thread?" is the conjunction, not the answer
  // for group 1.
  try {
    for(std::size_t g = 0; g < g_readers.size(); g++){
      std::unique_ptr<Geno::GenoReader> probe = g_readers[g]->clone();
    }
  } catch (const std::exception&) {
    return false;
  }
  return true;
}


// ---- S5 test hooks -------------------------------------------------------
// Not used by the pipeline. They exist because the reader's ownership and
// cloning behaviour is invisible from R otherwise: no golden can see whether a
// clone decodes the same bytes as its primary, or whether the map survives the
// primary's close. Same reasoning as the S2 seed counter -- a change the
// goldens cannot see carries its own measurement.

namespace {

// Decode a list of bim rows through one reader. Returns the genotype vectors
// alongside the metadata so a caller can compare decodes, not just summaries.
Rcpp::List readMarkersThrough(Geno::GenoReader& t_reader,
                              const std::vector<int>& t_gIndexVec)
{
  uint32_t n = t_reader.getN();
  size_t q = t_gIndexVec.size();
  Rcpp::List gvecs(q);
  Rcpp::CharacterVector ids(q);
  Rcpp::NumericVector altFreq(q), altCounts(q), missingRate(q);
  for(size_t k = 0; k < q; k++){
    Geno::MarkerMeta meta;
    std::vector<uint> indexForMissing, indexForNonZero;
    arma::vec GVec(n);
    t_reader.getOneMarker((uint64_t)t_gIndexVec[k], meta,
                          true, indexForMissing, false, indexForNonZero, GVec);
    gvecs[k] = Rcpp::wrap(GVec);
    ids[k] = meta.id;
    altFreq[k] = meta.altFreq;
    altCounts[k] = meta.altCounts;
    missingRate[k] = meta.missingRate;
  }
  return Rcpp::List::create(Rcpp::Named("g") = gvecs,
                            Rcpp::Named("id") = ids,
                            Rcpp::Named("altFreq") = altFreq,
                            Rcpp::Named("altCounts") = altCounts,
                            Rcpp::Named("missingRate") = missingRate);
}

}  // namespace

// Read the given bim rows through the reader currently set up.
// [[Rcpp::export]]
Rcpp::List testReadMarkersInCPP(std::vector<int> t_gIndexVec){
  if(!g_reader){
    Rcpp::stop("No genotype reader has been set up.");
  }
  return readMarkersThrough(*g_reader, t_gIndexVec);
}

// Read the same rows three ways: through the primary, through a clone, and
// through the clone AFTER the primary has been closed. The third is the one
// that matters -- it is the claim that the map dies with the last owner rather
// than with whoever calls close() first.
// NOTE: this closes the primary, so the reader must be set up again before any
// further use in the same process.
// [[Rcpp::export]]
Rcpp::List testGenoReaderCloneInCPP(std::vector<int> t_gIndexVec){
  if(!g_reader){
    Rcpp::stop("No genotype reader has been set up.");
  }
  std::unique_ptr<Geno::GenoReader> cloned = g_reader->clone();
  Rcpp::List fromPrimary = readMarkersThrough(*g_reader, t_gIndexVec);
  Rcpp::List fromClone = readMarkersThrough(*cloned, t_gIndexVec);
  g_reader->close();
  Rcpp::List fromCloneAfterClose = readMarkersThrough(*cloned, t_gIndexVec);
  return Rcpp::List::create(Rcpp::Named("primary") = fromPrimary,
                            Rcpp::Named("clone") = fromClone,
                            Rcpp::Named("cloneAfterPrimaryClose") = fromCloneAfterClose);
}


// ---- LEDGER #42 test hook --------------------------------------------------
// The two-sided quantile the SE back-calculation uses for a log-scale p. No
// golden can reach it: the log format is only chosen when the score p-value
// UNDERFLOWS to exactly 0 in double, which needs |z| ~ 40, and the bundled
// 1,000-sample fixtures top out around |z| = 2.6. So, as with the S2 seed
// counter and the S5 reader hooks, the change carries its own measurement.
// [[Rcpp::export]]
double testQnormTwoSidedLogPInCPP(double t_logp){
  return SAIGE::qnormTwoSidedLogP(t_logp);
}

// LEDGER #44: an effect size the score test declined to compute is printed as
// the literal "NA", not as a NaN. NaN spelling is platform-dependent -- glibc
// writes a sign-bit-set NaN as "-nan" where macOS writes "nan" -- which made
// exactly these rows non-byte-portable and cost two red CI runs before the
// underlying bug was fixed. R_IsNA distinguishes the deliberate NA_REAL the
// score tests now emit from an ordinary NaN arriving from anywhere else, so
// this changes those rows and nothing else.
static inline void writeNumOrNA(std::ofstream & t_os, double t_x){
  if(R_IsNA(t_x)){
    t_os << "NA";
  }else{
    t_os << t_x;
  }
}

// S4: the writer now streams MarkerResult slots instead of ~30 parallel
// q-length vectors. Same rows, same order, same skip-on-"NA" rule, same
// numtest; it just reads one object per marker rather than indexing thirty.
void writeOutfile_single(bool t_isMoreOutput,
			bool t_isImputation,
			bool t_isCondition,
			bool t_isFirth,
                         int mFirth,
                         int mFirthConverge,
                        std::string t_traitType,
                        const std::vector<Marker::MarkerResult> & t_resultVec,
                        std::ofstream & t_stream

){
  int numtest = 0;
  for(unsigned int k = 0; k < t_resultVec.size(); k++){
        const Marker::MarkerResult & res = t_resultVec.at(k);
        if(res.pval != "NA"){
                numtest = numtest + 1;
                t_stream << res.chr;
                t_stream << "\t";
                t_stream << res.pos;
                t_stream << "\t";
                t_stream << res.marker;
                t_stream << "\t";
                t_stream << res.ref;
                t_stream << "\t";
                t_stream << res.alt;
                t_stream << "\t";
                t_stream << res.altCounts;
                t_stream << "\t";
                t_stream << res.altFreq;
                t_stream << "\t";

                if(t_isImputation){
                        t_stream << res.imputationInfo;
                        t_stream << "\t";
                }else{
                        t_stream << res.missingRate;
                        t_stream << "\t";

                }
                writeNumOrNA(t_stream, res.Beta);
                t_stream << "\t";
                writeNumOrNA(t_stream, res.seBeta);
                t_stream << "\t";
                t_stream << res.Tstat;
                t_stream << "\t";
                t_stream << res.varT;
                t_stream << "\t";
                t_stream << res.pval;
                t_stream << "\t";

                if(t_traitType == "binary"){
                        t_stream << res.pvalNA;
                        t_stream << "\t";
                        t_stream << std::boolalpha << res.isSPAConverge;
                        t_stream << "\t";
                }
                if(t_isCondition){
                        writeNumOrNA(t_stream, res.Beta_c);
                        t_stream << "\t";
                        writeNumOrNA(t_stream, res.seBeta_c);
                        t_stream << "\t";
                        t_stream << res.Tstat_c;
                        t_stream << "\t";
                        t_stream << res.varT_c;
                        t_stream << "\t";
                        t_stream << res.pval_c;
                        t_stream << "\t";
			if(t_traitType == "binary"){	
                        	t_stream << res.pvalNA_c;
                        	t_stream << "\t";
			}	
                }
                if(t_traitType == "binary"){
                        t_stream << res.AF_case;
                        t_stream << "\t";
                        t_stream << res.AF_ctrl;
                        t_stream << "\t";
                        t_stream << res.N_case;
                        t_stream << "\t";
                        t_stream << res.N_ctrl;

                        if(t_isMoreOutput){
                                t_stream << "\t";
                                t_stream << res.N_case_hom;
                                t_stream << "\t";
                                t_stream << res.N_case_het;
                                t_stream << "\t";
                                t_stream << res.N_ctrl_hom;
                                t_stream << "\t";
                                t_stream << res.N_ctrl_het;
                        }
                        t_stream << "\n";
                }else if(t_traitType == "quantitative"){
                        t_stream << res.N;
                        t_stream << "\n";

                }
        }
  }
  std::cout << numtest << " markers were tested." << std::endl;
  if(t_traitType == "binary"){
      if(t_isFirth){
        std::cout << "Firth approx was applied to " << mFirth << " markers. " << mFirthConverge << " sucessfully converged." <<std::endl;
      }
   }
}


// The two R entry points that raise and lower the escalation signal. The names
// are R's ABI and are kept; what they write is g_regionEscalatedPass, not a
// model member (S12).
// [[Rcpp::export]]
void set_flagSparseGRM_cur_SAIGE(bool t_flagSparseGRM_cur){
	g_regionEscalatedPass = t_flagSparseGRM_cur;
}

// [[Rcpp::export]]
void set_flagSparseGRM_cur_SAIGE_org(){
	g_regionEscalatedPass = g_model->m_flagSparseGRM;
}



void writeOutfile_BURDEN(std::string regionName,
			std::vector<std::string>  & BURDEN_AnnoName_Vec,
			std::vector<std::string> & BURDEN_maxMAFName_Vec,
			std::vector<std::string> & BURDEN_pval_Vec,
			std::vector<double> & BURDEN_Beta_Vec,
			std::vector<double> & BURDEN_seBeta_Vec,
			std::vector<std::string> & BURDEN_pval_cVec,
			std::vector<double> & BURDEN_Beta_cVec,
			std::vector<double> & BURDEN_seBeta_cVec,
			arma::vec & MAC_GroupVec,
			arma::vec & MACCase_GroupVec,
			arma::vec & MACControl_GroupVec,
			arma::vec & NumRare_GroupVec,
			arma::vec & NumUltraRare_GroupVec,
			double cctpval,
			double cctpval_cond,
			unsigned int q_anno,
			unsigned int q_maf,
			unsigned int q_weight,
			std::vector<std::string> & t_weightNames,
			bool t_emitWeightColumn,
			bool isCondition,
			std::string t_traitType){
     unsigned int i;
     for(unsigned int j = 0; j < q_anno; j++){
       for(unsigned int m = 0; m < q_maf; m++){
         for(unsigned int r = 0; r < q_weight; r++){
           // jmr index, weight fastest -> rows come out weight-innermost;
           // q_weight == 1 emits the legacy rows byte for byte
           i = (j*q_maf+m)*q_weight + r;
	   if(BURDEN_pval_Vec.at(i) != "NA"){
           g_out.group << regionName;
           g_out.group << "\t";
           g_out.group << BURDEN_AnnoName_Vec.at(i);
           g_out.group << "\t";
           g_out.group << BURDEN_maxMAFName_Vec.at(i);
           g_out.group << "\t";
           if(t_emitWeightColumn){
             // per-set Weight label right after max_MAF (design D3); same
             // name fallback as the UR pseudo-marker suffix in mainRegionInCPP
             g_out.group << ((r < t_weightNames.size()) ? t_weightNames.at(r) : ("W" + std::to_string(r + 1)));
             g_out.group << "\t";
           }
           g_out.group << BURDEN_pval_Vec.at(i);
           g_out.group << "\t";
           // LEDGER #66: same reason as writeOutfile_singleInGroup above.
           writeNumOrNA(g_out.group, BURDEN_Beta_Vec.at(i));
           g_out.group << "\t";
           writeNumOrNA(g_out.group, BURDEN_seBeta_Vec.at(i));
           g_out.group << "\t";
           if(isCondition){
               g_out.group << BURDEN_pval_cVec.at(i);
               g_out.group << "\t";
               writeNumOrNA(g_out.group, BURDEN_Beta_cVec.at(i));
               g_out.group << "\t";
               writeNumOrNA(g_out.group, BURDEN_seBeta_cVec.at(i));
               g_out.group << "\t";
	   }
	   g_out.group << MAC_GroupVec(i);
           g_out.group << "\t";
           if(t_traitType == "binary"){
               g_out.group << MACCase_GroupVec(i);
               g_out.group << "\t";
               g_out.group << MACControl_GroupVec(i);
               g_out.group << "\t";
           }
           g_out.group << NumRare_GroupVec(i);
           g_out.group << "\t";
           g_out.group << NumUltraRare_GroupVec(i);
           g_out.group << "\n";
	}
         }
     }
   }
     g_out.group << regionName;
     g_out.group << "\tCauchy\tNA\t";
     if(t_emitWeightColumn){
       // the single pooled Cauchy row is not tied to one weight set
       g_out.group << "NA\t";
     }
     g_out.group << cctpval;
     g_out.group << "\tNA\tNA\t";	
     if(isCondition){
	g_out.group << cctpval_cond;
	g_out.group << "\tNA\tNA\t";
     }
     g_out.group << "NA\t";
     if(t_traitType == "binary"){
        g_out.group << "NA\t";
        g_out.group << "NA\t";
     }
     g_out.group << "NA\t";
     g_out.group << "NA\n";	
}


int writeOutfile_singleInGroup(bool t_isMoreOutput,
                        bool t_isImputation,
                        bool t_isCondition,
                        bool t_isFirth,
                         int mFirth,
                         int mFirthConverge,
                        std::string t_traitType,
                        std::vector<std::string> & chrVec,
                        std::vector<std::string> & posVec,
                        std::vector<std::string> & markerVec,
                        std::vector<std::string> & refVec,
                        std::vector<std::string> & altVec,
                        std::vector<double> & altCountsVec,
                        std::vector<double> & altFreqVec,
                        std::vector<double> & imputationInfoVec,
                        std::vector<double> & missingRateVec,
                        std::vector<double> & BetaVec,
                        std::vector<double> & seBetaVec,
                        std::vector<double> & TstatVec,
                        std::vector<double> & varTVec,
                        std::vector<std::string> & pvalVec,
                        std::vector<std::string> & pvalNAVec,
                        std::vector<unsigned char> & isSPAConvergeVec,
                        std::vector<double> & Beta_cVec,
                        std::vector<double> & seBeta_cVec,
                        std::vector<double> & Tstat_cVec,
                        std::vector<double> & varT_cVec,
                        std::vector<std::string> & pval_cVec,
                        std::vector<std::string> & pvalNA_cVec,
                        std::vector<double> & AF_caseVec,
                        std::vector<double> & AF_ctrlVec,
                        std::vector<uint32_t> & N_caseVec,
                        std::vector<uint32_t> & N_ctrlVec,
                        std::vector<double>  & N_case_homVec,
                        std::vector<double>  & N_ctrl_hetVec,
                        std::vector<double>  & N_case_hetVec,
                        std::vector<double>  & N_ctrl_homVec,
                        std::vector<uint32_t> & N_Vec,
                        std::ofstream & t_OutFile_singleInGroup){
  int numofUR = 0;
  for(unsigned int k = 0; k < pvalVec.size(); k++){
        //if(std::isfinite(pvalVec.at(k))){
        //if(!std::isnan(pvalVec.at(k))){
        if(pvalVec.at(k) != "NA"){
		if(chrVec.at(k) == "UR"){
                        numofUR = numofUR + 1;
                }
                t_OutFile_singleInGroup << chrVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << posVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << markerVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << refVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << altVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << altCountsVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << altFreqVec.at(k);
                t_OutFile_singleInGroup << "\t";

                if(t_isImputation){
                        t_OutFile_singleInGroup << imputationInfoVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                }else{
                        t_OutFile_singleInGroup << missingRateVec.at(k);
                        t_OutFile_singleInGroup << "\t";

                }
                // LEDGER #66: writeNumOrNA, not a raw <<. The score tests
                // (LEDGER #44) and now the conditional block emit NA_REAL for
                // an effect size they declined to compute, and this writer --
                // unlike writeOutfile_single -- printed it as a bare NaN,
                // whose spelling is platform-dependent. Byte-neutral on every
                // committed golden: none carries a non-finite effect size.
                writeNumOrNA(t_OutFile_singleInGroup, BetaVec.at(k));
                t_OutFile_singleInGroup << "\t";
                writeNumOrNA(t_OutFile_singleInGroup, seBetaVec.at(k));
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << TstatVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << varTVec.at(k);
                t_OutFile_singleInGroup << "\t";
                t_OutFile_singleInGroup << pvalVec.at(k);
                t_OutFile_singleInGroup << "\t";

                if(t_traitType == "binary"){
                        t_OutFile_singleInGroup << pvalNAVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << std::boolalpha << (isSPAConvergeVec.at(k) != 0);
                        t_OutFile_singleInGroup << "\t";
                }
                if(t_isCondition){
                        writeNumOrNA(t_OutFile_singleInGroup, Beta_cVec.at(k));
                        t_OutFile_singleInGroup << "\t";
                        writeNumOrNA(t_OutFile_singleInGroup, seBeta_cVec.at(k));
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << Tstat_cVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << varT_cVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << pval_cVec.at(k);
                        t_OutFile_singleInGroup << "\t";
			if(t_traitType == "binary"){
                        	t_OutFile_singleInGroup << pvalNA_cVec.at(k);
                        	t_OutFile_singleInGroup << "\t";
			}
                }
                if(t_traitType == "binary"){
                        t_OutFile_singleInGroup << AF_caseVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << AF_ctrlVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << N_caseVec.at(k);
                        t_OutFile_singleInGroup << "\t";
                        t_OutFile_singleInGroup << N_ctrlVec.at(k);

                        if(t_isMoreOutput){
                                t_OutFile_singleInGroup << "\t";
                                t_OutFile_singleInGroup << N_case_homVec.at(k);
                                t_OutFile_singleInGroup << "\t";
                                t_OutFile_singleInGroup << N_case_hetVec.at(k);
                                t_OutFile_singleInGroup << "\t";
                                t_OutFile_singleInGroup << N_ctrl_homVec.at(k);
                                t_OutFile_singleInGroup << "\t";
                                t_OutFile_singleInGroup << N_ctrl_hetVec.at(k);
                        }
                        t_OutFile_singleInGroup << "\n";
                }else if(t_traitType == "quantitative"){
                        t_OutFile_singleInGroup << N_Vec.at(k);
                        t_OutFile_singleInGroup << "\n";

                }
        }
  }
  return(numofUR);

}

// [[Rcpp::export]]
void copy_singleInGroup(){
  std::ifstream ini_file;
  ini_file.open(g_out.prefixSingleInGroupTemp.c_str());
  if (!ini_file)
  {
    std::cout << "Error in Opening the temp file!" << std::endl;
  }
  std::string str;
  while (getline(ini_file, str))
  {
    g_out.singleInGroup << str << "\n"; 	
  }
  ini_file.close();
}

// Exhaustively compare lchoose_threadsafe against Rf_lchoose over the integer
// domain HyperGeo::lCombinations can reach, and report the number of pairs that
// differ AS BITS (memcmp, so NaN/Inf/-0.0 are held to the same standard rather
// than compared with ==, which would call every NaN pair a difference and every
// -0.0/0.0 pair a match).
//
// This exists because the goldens cannot do the job. They DO execute the
// function -- replacing its main branch with a hard error fails the standard
// tier -- but perturbing that branch by +1e-9 leaves the entire tier green, so
// they bound its numerics only loosely. The swap's whole justification is
// bit-identity, so bit-identity is what gets tested.
// [[Rcpp::export]]
int lchooseThreadsafeDiffCountInCPP(int t_nmax)
{
  int ndiff = 0;
  for(int n = 0; n <= t_nmax; n++){
    for(int k = -3; k <= n; k++){
      double a = Rf_lchoose((double)n, (double)k);
      double b = lchoose_threadsafe((double)n, (double)k);
      if(std::memcmp(&a, &b, sizeof(double)) != 0) ndiff++;
    }
  }
  return ndiff;
}
