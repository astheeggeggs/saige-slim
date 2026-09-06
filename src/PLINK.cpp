
// All C++ codes related to PLINK file manipulation

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <boost/algorithm/string.hpp>
#include <array>
#include <algorithm>
#include "PLINK.hpp"
#include "UTIL.hpp"
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <RcppArmadilloExtensions/sample.h> // sample
#include <string>

// make a global variable for future usage
// static PLINK::PlinkClass* ptr_gPLINKobj = NULL;

// ptr_gPLINKobj = new PLINK::PlinkClass(t_bimFile,
//                                       t_famFile,
//                                       t_bedFile,
//                                       t_SampleInModel);

namespace PLINK {

PlinkClass::PlinkClass(std::string t_bimFile,
                       std::string t_famFile,
                       std::string t_bedFile,
                       //std::vector<std::string> t_SampleInModel,
                       std::string t_AlleleOrder)
{
  setPlinkobj(t_bimFile, t_famFile, t_bedFile);
  //setPosSampleInPlink(t_SampleInModel);
  m_AlleleOrder = t_AlleleOrder;
}



// set PlinkClass by reading PLINK files
void PlinkClass::setPlinkobj(std::string t_bimFile,
                             std::string t_famFile,
                             std::string t_bedFile)
{
  m_bimFile = t_bimFile;
  m_famFile = t_famFile;
  m_bedFile = t_bedFile;
  
  // setChrMaps();
  readBimFile();
  readFamFile();
  
  // LEDGER #58: the .bed must be at least as long as the .bim and .fam imply.
  // Checked ONCE, here, rather than per read: it fails up front and can name
  // both files and the shortfall. Without it the default mmap path indexes
  // `m_map + posSeek` with no bound (m_mapLen is set only to be handed back to
  // munmap), so a .bim claiming more markers than the .bed holds read as zeros
  // in the last page's slack -- decoding to genotype 2 for every sample, hence
  // AF 1.0, MAC 0, and a silent drop by the --minMAC gate -- and SIGSEGV'd
  // beyond it. Measured on genotype_100markers: one phantom .bim row gave
  // output BYTE-IDENTICAL to the intact run at exit 0; three gave exit 139.
  // The pread fallback caught it per read, so the loud path was the one nothing
  // executes and the silent one was the default.
  //
  // LEDGER #87 FIXED: the check is now an EQUALITY. #58 allowed a longer .bed
  // on the argument that every marker the .bim names is then present and
  // correctly placed; that holds only when the .bim is a PREFIX of the true
  // marker set. A .bim missing an INTERIOR row gives exactly the same length
  // and relabels every marker after the gap -- measured: rs5 deleted, 70 rows
  // at exit 0, 68 of them carrying the statistics of the marker before them.
  // plink itself refuses a size mismatch in either direction; so do we.
  {
    struct stat stCheck;
    const uint64_t needed = (uint64_t)3 +
      (uint64_t)m_numBytesofEachMarker0 * (uint64_t)m_M0;
    if(stat(t_bedFile.c_str(), &stCheck) != 0){
      Rcpp::stop("Cannot stat the .bed file (" + t_bedFile + ").");
    }
    if((uint64_t)stCheck.st_size > needed){
      Rcpp::stop("The .bed file (" + t_bedFile + ") is LONGER than the .bim and .fam "
                 "that accompany it imply: " + std::to_string(m_M0) + " markers x " +
                 std::to_string(m_N0) + " samples need " + std::to_string(needed) +
                 " bytes, but the file is " + std::to_string((uint64_t)stCheck.st_size) +
                 ". The .bim describes fewer markers than the .bed holds -- a row missing "
                 "from the .bim would relabel every marker after it -- so the three files "
                 "do not describe the same dataset.");
    }
    if((uint64_t)stCheck.st_size < needed){
      Rcpp::stop("The .bed file (" + t_bedFile + ") is too short for the .bim and .fam "
                 "that accompany it: " + std::to_string(m_M0) + " markers x " +
                 std::to_string(m_N0) + " samples need " + std::to_string(needed) +
                 " bytes, but the file is " + std::to_string((uint64_t)stCheck.st_size) +
                 ". The three files do not describe the same dataset -- check that the "
                 ".bed is complete and that the .bim was not regenerated from a "
                 "different marker set.");
    }
  }

  m_fin = fopen(t_bedFile.c_str(), "rb");
  // LEDGER #54. The mode byte IS checked below and always was; the two MAGIC
  // bytes were read into a 64 KiB stack buffer and never looked at, so a file
  // that is not a PLINK .bed at all was accepted as one -- measured on a copy of
  // genotype_100markers with 0x6c 0x1b overwritten by zeros, which produced
  // output BYTE-IDENTICAL to the intact run at exit 0, while plink2 refuses the
  // same file ("first two bytes don't match the magic number"). The fread return
  // values were discarded too, which would leave magicNumber3 uninitialised on a
  // short file; the #58 size check above already refuses one, but the read is
  // checked rather than left resting on that.
  unsigned char magicNumber_2[2];
  size_t nHdr = fread(magicNumber_2, 1, 2, m_fin);
  char magicNumber3;
  nHdr += fread(&magicNumber3, 1, 1, m_fin);
  fseek(m_fin, 3, SEEK_SET);
  if(nHdr != 3)
    Rcpp::stop("The .bed file (" + t_bedFile + ") is shorter than its 3-byte header.");
  if(magicNumber_2[0] != 0x6c || magicNumber_2[1] != 0x1b)
    Rcpp::stop("The .bed file (" + t_bedFile + ") does not start with the PLINK magic number "
               "0x6c 0x1b, so it is not a PLINK .bed file. Check that the --bedFile argument "
               "names the right file.");

  // Phase 4: map the whole .bed read-only. Falls back to pread on failure.
  // S5: SAIGE_PLINK_NO_MMAP forces that fallback. It exists because the
  // fallback was reachable only by an mmap failure -- which does not happen for
  // a small local file -- so NOTHING had ever exercised it, on any platform, in
  // any golden. The same species of gap as F8/S0: untested is not the same as
  // absent, and the two paths must decode identically.
  if(getenv("SAIGE_PLINK_NO_MMAP") == NULL){
    int fd = fileno(m_fin);
    struct stat st;
    if(fstat(fd, &st) == 0 && st.st_size > 3){
      void* p = mmap(NULL, (size_t)st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
      if(p != MAP_FAILED){
        m_map = (const unsigned char*)p;
        m_mapLen = (size_t)st.st_size;
      }
    }
  }

/* 
  m_ibedFile.open(m_bedFile.c_str(), std::ios::binary);
  m_ibedFile.seekg(2);
  char magicNumber3;
  m_ibedFile.read(&magicNumber3, 1);
*/

  if(magicNumber3 != 1)
    Rcpp::stop("The third magic number of the plink bed file is not 00000001. Please use SNP-major plink (plink version >= 1.9) files.");
}

static uint32_t parseBimPosField(const std::string& s, const std::string& file, size_t lineNo)
{
  try {
    size_t used = 0;
    long long v = std::stoll(s, &used);
    if(used != s.size() || v < 0 || v > 0xFFFFFFFFLL) throw std::invalid_argument("range");
    return static_cast<uint32_t>(v);
  } catch(const std::exception&) {
    Rcpp::stop("Cannot parse the position field '" + s + "' on line " + std::to_string(lineNo) +
               " of " + file + " (column 4 of a .bim).");
  }
}

void PlinkClass::readBimFile()
{
  std::cout << "Reading bim file...." << std::endl;
  std::ifstream bim(m_bimFile);
  m_M0 = 0;
  std::string line;
  
  while(getline(bim, line)){
    m_M0++;
    std::vector<std::string> line_elements;
    boost::split(line_elements, line, boost::is_any_of("\t "));
    boost::replace_all(line_elements[line_elements.size() - 1], "\r", "");
    // LEDGER #89 FIXED. A .bim row has six fields; a shorter one used to index
    // the split vector out of range. Column 3 (centimorgans) was parsed with a
    // bare std::stof into m_gd, which nothing in the tree ever read -- so a
    // '.' there refused a usable dataset -- and is no longer parsed at all.
    // Column 4 (position) is parsed with a message naming file, line, column.
    if(line_elements.size() < 6){
      Rcpp::stop("bim file line " + std::to_string(m_M0) + " has " + std::to_string(line_elements.size()) +
                 " fields; a .bim row needs 6 (" + m_bimFile + ").");
    }
    m_chr.push_back(line_elements[0]);
    m_MarkerInPlink.push_back(line_elements[1]);
    m_pd.push_back(parseBimPosField(line_elements[3], m_bimFile, m_M0));
    std::transform(line_elements[4].begin(), line_elements[4].end(), line_elements[4].begin(), toupper);
    std::transform(line_elements[5].begin(), line_elements[5].end(), line_elements[5].begin(), toupper);
    m_alt.push_back(line_elements[4]);  // allele 1, usually minor allele, alt-first
    m_ref.push_back(line_elements[5]);  // allele 2, usually major allele
  }
  m_M = m_M0;
}

void PlinkClass::readFamFile()
{
  std::cout << "Reading fam file...." << std::endl;
  std::ifstream fam(m_famFile);
  m_N0 = 0;
  std::string line;
  while(getline(fam, line)){
    m_N0 ++;
    std::vector<std::string> line_elements;
    boost::split(line_elements, line, boost::is_any_of("\t "));
    m_SampleInPlink.push_back(line_elements[1]);  // put IID to m_SampleInPlink
  }
  m_N = m_N0;
  m_numBytesofEachMarker0 = (m_N0 + 3) / 4;
} 


void PlinkClass::setPosSampleInPlink(std::vector<std::string> & t_SampleInModel)
{
  std::cout << "Setting position of samples in PLINK files...." << std::endl;
  m_N = t_SampleInModel.size();
  m_numBytesofEachMarker = (m_N + 3) / 4;
  
  // convert from std::vector<std::string> to Rcpp::CharacterVector
  Rcpp::CharacterVector SampleInModel(m_N);
  for(uint32_t i = 0; i < m_N; i++)
    SampleInModel(i) = t_SampleInModel.at(i);
  
  // convert from std::vector<std::string> to Rcpp::CharacterVector
  uint32_t N_plink = m_SampleInPlink.size();
  Rcpp::CharacterVector SampleInPlink(N_plink);
  for(uint32_t i = 0; i < N_plink; i++)
    SampleInPlink(i) = m_SampleInPlink.at(i);
    
  // Rcpp::match is much faster than loop
  Rcpp::IntegerVector posSampleInPlink = Rcpp::match(SampleInModel, SampleInPlink);
  m_posSampleInPlink.resize(m_N);
  for(uint32_t i = 0; i < m_N; i++){
    if(Rcpp::IntegerVector::is_na(posSampleInPlink.at(i)))
      Rcpp::stop("At least one subject requested is not in Plink file.");
    m_posSampleInPlink.at(i) = posSampleInPlink.at(i) - 1;   // convert "starting from 1" to "starting from 0"
    //std::cout << "m_posSampleInPlink.at(i) " << i << " " << m_posSampleInPlink.at(i) << std::endl;
     }
  m_isIdentityPos = (m_N == m_N0);
  for(uint32_t i = 0; i < m_N && m_isIdentityPos; i++){
    if(m_posSampleInPlink[i] != i) m_isIdentityPos = false;
  }


}

// std::vector<uint32_t> PlinkClass::getPosMarkerInPlink(std::vector<std::string> t_MarkerReqstd)
// {
//   int M = t_MarkerReqstd.size();
//   std::vector<uint32_t> posMarkerInPlink;
//   for(int i = 0; i < M; i++){
//     std::string marker = t_MarkerReqstd.at(i);
//     auto pos = std::find(m_MarkerInPlink.begin(), m_MarkerInPlink.end(), marker);
//     if(pos != m_MarkerInPlink.end()){
//       posMarkerInPlink.push_back(pos - m_MarkerInPlink.begin());
//     }else{
//       Rcpp::warning("Marker %s is not found in plink file.", marker);
//     }
//   }
//   return posMarkerInPlink;
// }

// S5: idempotent. The map was already guarded; the fclose was not, so a second
// close was undefined behaviour -- which mattered the moment a clone could
// outlive the owner that closed. Both handles are now cleared as they go.
void PlinkClass::closegenofile(){
        if(m_map){ munmap((void*)m_map, m_mapLen); m_map = NULL; m_mapLen = 0; }
        if(m_fin){ fclose(m_fin); m_fin = NULL; }
	//m_ibedFile.close();
}

PlinkClass::~PlinkClass(){
        closegenofile();
}


const unsigned char* PlinkClass::markerBytes(uint64_t t_gIndex, std::vector<unsigned char>& t_buf) const
{
  // LEDGER #16 fix: seek unconditionally to the requested record. The old
  // gate (if t_gIndex > 0) skipped the seek for bim row 0, so a mid-stream
  // request for row 0 decoded whatever record the cursor sat on and reported
  // it under row 0's marker ID (demonstrated via conditional analysis).
  uint64_t posSeek = 3 + m_numBytesofEachMarker0 * t_gIndex;
  // S5: the fallback was fseek + fread into the shared member buffer, which put
  // TWO pieces of mutable state on every read -- the FILE*'s offset and the
  // buffer. pread carries the offset as an argument and never touches the
  // file description's own, so the fd can be shared; the buffer is a local. The
  // bytes are identical (same offset, same length), and the mmap branch above
  // is untouched.
  if(m_map){
    return m_map + posSeek;               // no syscalls, no copy
  }
  if(!m_fin)
    Rcpp::stop("PlinkClass::getOneMarker called after the .bed file was closed.");
  t_buf.resize(m_numBytesofEachMarker0);
  ssize_t nread = pread(fileno(m_fin), (void*)t_buf.data(),
                        (size_t)m_numBytesofEachMarker0, (off_t)posSeek);
  if(nread != (ssize_t)m_numBytesofEachMarker0)
    Rcpp::stop("Short read from the .bed file at marker index %s.",
               std::to_string(t_gIndex));
  return t_buf.data();
}

void PlinkClass::summariseCounts(const int t_counts[4], double& t_altFreq, double& t_altCounts,
                                 double& t_missingRate, double& t_imputeInfo) const
{
  int numMissing = t_counts[MISSING];
  int count = m_N - numMissing;
  t_missingRate = (double)numMissing / (double)m_N;
  t_imputeInfo = 1;
  t_altCounts = (double)(t_counts[HET] + 2*t_counts[HOM_ALT]);

  if(count > 0){
        t_altFreq = t_altCounts/ (double)count / 2;
  }else{
        t_altFreq = 0;
  }

  // updated on 03/14/2021
  if(m_AlleleOrder == "ref-first"){
    // LEDGER #29 (first of two lossy steps): the ALT count used to be
    // reconstructed FROM the frequency as 2*count*(1-f), a float round trip
    // landing up to ~500 ulps off the true integer -- so under ref-first
    // t_altCounts reached the QC gates already lossy. counts[] are integers and
    // count is exact, so take the complement by exact integer arithmetic and
    // derive the frequency from it, not the other way round.
    if(count > 0){
      t_altCounts = 2.0 * (double)count - t_altCounts;
      t_altFreq = t_altCounts / (double)count / 2;
    }else{
      // Degenerate all-missing marker: preserved verbatim (altFreq 1-0 = 1,
      // altCounts 0). Dropped by the missingRate gate long before this matters.
      t_altFreq = 1 - t_altFreq;
      t_altCounts = 2 * (double)count * t_altFreq;
    }
  }
}

void PlinkClass::countOneMarker(uint64_t t_gIndex,
                                std::string& t_ref, std::string& t_alt, std::string& t_marker,
                                uint32_t& t_pd, std::string& t_chr,
                                double& t_altFreq, double& t_altCounts, double& t_missingRate,
                                double& t_imputeInfo, uint32_t& t_numMissing) const
{
  // Item 3 (2026-09-05). On the N = 50,000 fixture 210,515 markers were decoded
  // to doubles and 70,709 passed --minMAC=20; the decode was 15% of the binary
  // scan's CPU and 60% of the quantitative one. QC gate 1 needs only the four
  // code tallies, so tally first and let the caller decide whether to expand.
  std::vector<unsigned char> oneMarkerG4;
  const unsigned char* src = markerBytes(t_gIndex, oneMarkerG4);

  t_marker = m_MarkerInPlink[t_gIndex];
  t_pd = m_pd[t_gIndex];
  t_chr = m_chr[t_gIndex];
  if(m_AlleleOrder == "alt-first"){
    t_ref = m_ref[t_gIndex];
    t_alt = m_alt[t_gIndex];
  }
  if(m_AlleleOrder == "ref-first"){
    t_ref = m_alt[t_gIndex];
    t_alt = m_ref[t_gIndex];
  }

  int counts[] = {0, 0, 0, 0};
  if(m_isIdentityPos){
    // Whole file in file order: one table lookup and one add per byte. Each
    // table word packs the byte's four code tallies in four 8-bit lanes, so a
    // lane holds at most 4 per byte and 63 bytes can be summed before any lane
    // could overflow; the block sum is unpacked into the 32-bit counts. The
    // first version of this loop did four masked increments per byte and was
    // no faster than the decode it replaces (quant scan 61 s vs 62 s).
    // The last, partial byte is walked per sample so its padding bits (zero,
    // which is HOM_ALT's code) are never counted.
    static const std::array<uint32_t, 256> lanes = [](){
      std::array<uint32_t, 256> t{};
      for(int v = 0; v < 256; v++){
        uint32_t w = 0;
        for(int k = 0; k < 4; k++){
          w += 1u << (8 * ((v >> (2 * k)) & 0x3));
        }
        t[v] = w;
      }
      return t;
    }();
    const uint32_t nFull = m_N / 4;
    uint32_t b = 0;
    while(b < nFull){
      const uint32_t end = std::min(nFull, b + 63);
      uint32_t acc = 0;
      for(; b < end; b++){
        acc += lanes[src[b]];
      }
      counts[0] += (int)(acc & 0xFF);
      counts[1] += (int)((acc >> 8) & 0xFF);
      counts[2] += (int)((acc >> 16) & 0xFF);
      counts[3] += (int)((acc >> 24) & 0xFF);
    }
    for(uint32_t i = nFull * 4; i < m_N; i++){
      size_t g;
      getGenotype(src[i/4], i%4, g);
      counts[g]++;
    }
  }else{
    for(uint32_t i = 0; i < m_N; i++){
      auto ind = m_posSampleInPlink[i];
      size_t g;
      getGenotype(src[ind/4], ind%4, g);
      counts[g]++;
    }
  }
  t_numMissing = counts[MISSING];
  summariseCounts(counts, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo);
}

void PlinkClass::getOneMarker(				uint64_t t_gIndex,        // different meanings for different genoType
                                   std::string& t_ref,       // REF allele
                                   std::string& t_alt,       // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                                   std::string& t_marker,    // marker ID extracted from genotype file
                                   uint32_t& t_pd,           // base position
                                   std::string& t_chr,       // chromosome
                                   double& t_altFreq,        // frequency of ALT allele
                                   double& t_altCounts,      // counts of ALT allele
                                   double& t_missingRate,    // missing rate
                                   double& t_imputeInfo,     // imputation information score, i.e., R2 (all 1 for PLINK)
                                   bool t_isOutputIndexForMissing,               // if true, output index of missing genotype data
                                   std::vector<uint>& t_indexForMissing,     // index of missing genotype data
                                   bool t_isOnlyOutputNonZero,                   // if true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                                   std::vector<uint>& t_indexForNonZero,     // only used when t_isOnlyOutputNonZero = TRUE
                                   bool t_isTrueGenotype,    // only used in PLINK. check m_genoMaps for details about the genotype mapping in PLINK.
				   arma::vec & OneMarkerG1) const
				   //std::vector<double>& OneMarkerG1)
{		
  
  //std::vector<double> OneMarkerG1;
 //arma::vec timeoutput1 = getTime(); 
  // t_isTrueGenotype = FALSE is used only when calculating GRM
  if(!t_isTrueGenotype){
    if(t_isOutputIndexForMissing)
      Rcpp::stop("Check PlinkClass::getOneMarker, if t_isTrueGenotype = FALSE, then t_isOutputIndexForMissing should be FALSE.");
    if(t_isOnlyOutputNonZero)
      Rcpp::stop("Check PlinkClass::getOneMarker, if t_isTrueGenotype = FALSE, then t_isOnlyOutputNonZero should be FALSE.");
  }

  // The record bytes come from markerBytes (LEDGER #16 seek, S5 pread fallback).
  std::vector<unsigned char> oneMarkerG4;
  const unsigned char* src = markerBytes(t_gIndex, oneMarkerG4);

  t_indexForMissing.clear();
  t_indexForNonZero.clear();
  
  t_marker = m_MarkerInPlink[t_gIndex];
  t_pd = m_pd[t_gIndex];
  t_chr = m_chr[t_gIndex];
  
  //std::map<int8_t, int8_t> genoMaps;
  const int8_t* genoMaps = m_genoMaps_alt_first.data();  // Phase 4: no per-marker vector copy
  
  if(m_AlleleOrder == "alt-first"){
    t_ref = m_ref[t_gIndex];
    t_alt = m_alt[t_gIndex];
    genoMaps = m_genoMaps_alt_first.data();
  }
    
  if(m_AlleleOrder == "ref-first"){
    t_ref = m_alt[t_gIndex];
    t_alt = m_ref[t_gIndex];
    genoMaps = m_genoMaps_ref_first.data();
  }
 //arma::vec timeoutput2 = getTime();
  uint j = 0;
  int counts[] = {0, 0, 0, 0};
  for(uint32_t i = 0; i < m_N; i++)
  {
    auto ind = m_posSampleInPlink[i];             // C++ start from 0
    unsigned char bufferG4 = src[ind/4];    // unsigned char: 1 byte for 4 genotypes (4 samples)
    size_t bufferG1;                                     // int: 1 genotype (1 sample)
    
    // https://www.cog-genomics.org/plink/1.9/formats#bed
    getGenotype(bufferG4, ind%4, bufferG1);          // bufferG4 -> bufferG1
    
    counts[bufferG1]++;

    if(bufferG1 == MISSING && t_isOutputIndexForMissing){
      t_indexForMissing.push_back(i);
    }
    
    
    if(t_isTrueGenotype)
      bufferG1 = genoMaps[bufferG1];


    // Only built when the caller actually asked for the packed non-zero form.
    // Both step-2 call sites pass t_isOnlyOutputNonZero = false and then have
    // imputeGenoAndFlip rebuild this vector from scratch (it clears it first,
    // UTIL.cpp:70-71), so the per-sample push_back here was pure waste.
    if(t_isOnlyOutputNonZero && bufferG1 > 0){
      t_indexForNonZero.push_back(i);
    }
    if(t_isOnlyOutputNonZero){
      if(bufferG1 > 0){
        OneMarkerG1[j] = bufferG1;
	      j = j + 1;
      }
    }else{
      OneMarkerG1[i] = bufferG1;
    }
  }

  summariseCounts(counts, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo);
 
  if(t_isOnlyOutputNonZero){
    OneMarkerG1.resize(j);
  }
 //arma::vec timeoutput4 = getTime();
  //return OneMarkerG1;
   //  printTime(timeoutput1, timeoutput2, "Unified_getOneMarker 1");
   //  printTime(timeoutput2, timeoutput3, "Unified_getOneMarker 2");
   //  printTime(timeoutput3, timeoutput4, "Unified_getOneMarker 3");
}

// C++ version of which(). Note: start from 0, not 1 
std::vector<unsigned int> whichCPP(std::vector<std::string> & strVec, 
                                   std::string strValue)
{
  std::vector<unsigned int> indexVec;
  for(unsigned int i = 0; i < strVec.size(); i++){
    if(strVec.at(i)==strValue)
      indexVec.push_back(i);
  }
  return(indexVec);
}


}

// make a global variable for future usage
// static PLINK::PlinkClass* ptr_gPLINKobj = NULL;

// ptr_gPLINKobj = new PLINK::PlinkClass(t_bimFile,
//                                       t_famFile,
//                                       t_bedFile,
//                                       t_SampleInModel);


// // [[Rcpp::export]]
// void setPLINKobjInR(std::string t_bimFile,
//                     std::string t_famFile,
//                     std::string t_bedFile,
//                     std::vector<std::string> t_SampleInModel)
// {
//   ptr_gPLINKobj = new PLINK::PlinkClass(t_bimFile,
//                                         t_famFile,
//                                         t_bedFile,
//                                         t_SampleInModel);
//   
//   int n = ptr_gPLINKobj->getN();
//   std::cout << "n:\t" << n << std::endl;
// }



