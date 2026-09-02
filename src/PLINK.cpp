
// All C++ codes related to PLINK file manipulation

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <boost/algorithm/string.hpp>
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
  
  m_fin = fopen(t_bedFile.c_str(), "rb");
  char magicNumber_2[65536]; 
  fread(magicNumber_2, 2, 1, m_fin);
  char magicNumber3;
  fread(&magicNumber3, 1, 1, m_fin);
  fseek(m_fin, 3, SEEK_SET);    

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
    
    m_chr.push_back(line_elements[0]);
    m_MarkerInPlink.push_back(line_elements[1]);
    m_gd.push_back(std::stof(line_elements[2]));
    m_pd.push_back(std::stoi(line_elements[3]));
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

  // LEDGER #16 fix: seek unconditionally to the requested record. The old
  // gate (if t_gIndex > 0) skipped the seek for bim row 0, so a mid-stream
  // request for row 0 decoded whatever record the cursor sat on and reported
  // it under row 0's marker ID (demonstrated via conditional analysis).
  uint64_t posSeek = 3 + m_numBytesofEachMarker0 * t_gIndex;
  const unsigned char* src;
  // S5: the fallback was fseek + fread into the shared member buffer, which put
  // TWO pieces of mutable state on every read -- the FILE*'s offset and the
  // buffer. pread carries the offset as an argument and never touches the
  // file description's own, so the fd can be shared; the buffer is a local. The
  // bytes are identical (same offset, same length), and the mmap branch above
  // is untouched.
  std::vector<unsigned char> oneMarkerG4;
  if(m_map){
    src = m_map + posSeek;               // no syscalls, no copy
  }else{
    if(!m_fin)
      Rcpp::stop("PlinkClass::getOneMarker called after the .bed file was closed.");
    oneMarkerG4.resize(m_numBytesofEachMarker0);
    ssize_t nread = pread(fileno(m_fin), (void*)oneMarkerG4.data(),
                          (size_t)m_numBytesofEachMarker0, (off_t)posSeek);
    if(nread != (ssize_t)m_numBytesofEachMarker0)
      Rcpp::stop("Short read from the .bed file at marker index %s.",
                 std::to_string(t_gIndex));
    src = oneMarkerG4.data();
  }

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

  int numMissing = counts[MISSING];
  
 //arma::vec timeoutput3 = getTime();
  int count = m_N - numMissing;
  t_missingRate = (double)numMissing / (double)m_N;
  t_imputeInfo = 1;
  t_altCounts = (double)(counts[HET] + 2*counts[HOM_ALT]);

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



