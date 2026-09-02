
// All C++ codes related to PLINK file manipulation

#ifndef PLINK_HPP
#define PLINK_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

namespace PLINK {

class PlinkClass{
private:
  
  // added on 03/14/2021
  std::string m_AlleleOrder;           // "ref-first" (default for BGEN) or "alt-first" (default for PLINK)
  
  // information from bim file
  uint32_t m_M0, m_M;
  std::vector<std::string> m_chr;               // Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
  std::vector<std::string> m_MarkerInPlink;     // Variant identifier
  std::vector<float> m_gd;                      // Position in morgans or centimorgans (safe to use dummy value of '0')
  std::vector<uint32_t> m_pd;                   // Base-pair coordinate (1-based; limited to 2^31-2)
  std::vector<std::string> m_alt;               // Allele 1 (corresponding to clear bits in .bed; usually minor)
  std::vector<std::string> m_ref;               // Allele 2 (corresponding to set bits in .bed; usually major)
  
  // information from fam file
  std::vector<std::string> m_SampleInPlink;
  uint32_t m_N0, m_N;
  unsigned long long int m_numBytesofEachMarker0, m_numBytesofEachMarker;
  
  // input file stream of .bed file
  //std::ifstream m_ibedFile;
  FILE *m_fin = NULL;
  // Phase 4: read-only memory map of the .bed file (NULL if mmap failed;
  // the pread path below is the fallback). Avoids a read syscall and a buffer
  // copy per marker.
  const unsigned char* m_map = NULL;
  size_t m_mapLen = 0;

  // PLINK files
  std::string m_bimFile, m_famFile, m_bedFile;
  std::vector<uint32_t> m_posSampleInPlink;
  
  // https://www.cog-genomics.org/plink/1.9/formats#bed
  // PLINK format
  // The two-bit genotype codes have the following meanings:
  // 00	Homozygous for first allele in .bim file
  // 01	Missing genotype
  // 10	Heterozygous
  // 11	Homozygous for second allele in .bim file
  
  const static unsigned char HOM_REF = 0x3;  // 0b11 ;
  const static unsigned char HET = 0x2;      // 0b10 ;
  const static unsigned char HOM_ALT = 0x0;  // 0b00 ;
  const static unsigned char MISSING = 0x1;  // 0b01 ;
  
  // or use "arma::datum::nan"
  // std::vector<int8_t> m_genoMaps = {{2, -1, 1, 0}};
  std::vector<int8_t> m_genoMaps_alt_first = {{2, -1, 1, 0}};
  std::vector<int8_t> m_genoMaps_ref_first = {{0, -1, 1, 2}};
  
  // S5: the fallback's read buffer used to be the member m_OneMarkerG4, which
  // is the whole reason getOneMarker could not be const -- one shared scratch
  // buffer that every read wrote through. It is now a local in the fallback
  // branch, which runs only when mmap failed and already pays a syscall per
  // marker, so one small allocation there is not measurable.

  // void setChrMaps();
  void readBimFile();
  void readFamFile();
  
  // extract geno (0,1,2,3) at specific pos (0,1,2,3) of address c (1 byte)  
  inline void getGenotype(const unsigned char c, const uint32_t pos, size_t& geno) const {
    geno = (c >> (pos << 1)) & 0x3;  // 0b11 = 0x3 
  }
  
public:
  
  PlinkClass(std::string t_bimFile,
             std::string t_famFile,
             std::string t_bedFile,
             std::string t_AlleleOrder);

  // S5: ownership, not convention, releases the map and the FILE*. Nothing in
  // R ever called the step-2 closeGenoFile(), so before this both survived to
  // process exit. Copying would double-close, so the class is non-copyable --
  // clones share ONE instance through a shared_ptr (see GenoReaderRegistry).
  ~PlinkClass();
  PlinkClass(const PlinkClass&) = delete;
  PlinkClass& operator=(const PlinkClass&) = delete;
  
             //std::vector<std::string> t_SampleInModel,dd
  // setup PlinkClass
  void setPlinkobj(std::string t_bimFile,
                   std::string t_famFile,
                   std::string t_bedFile);
  
  void setPosSampleInPlink(std::vector<std::string> & t_SampleInModel);
  
  // S5: const. It writes NO member -- the mmap branch is pointer arithmetic and
  // the fallback preads into a local -- which is exactly what lets several
  // threads run it on one shared PlinkClass. The compiler is the audit here, as
  // it was for the model in S3: do NOT add a mutable member to get around it.
  void getOneMarker(  			uint64_t t_gIndex,        // different meanings for different genoType
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
                         bool t_isOnlyOutputNonZero,                   // is true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                         std::vector<uint>& t_indexForNonZero,
                         bool t_isTrueGenotype, // only used in PLINK. check m_genoMaps for details about the genotype mapping in PLINK
			 arma::vec & OneMarkerG1) const;
			 //std::vector<double>& OneMarkerG1);

  // True when the .bed is read through the mmap path (no per-read seek on
  // the shared FILE*), i.e. when reads are safe from forked children.
  bool isMemMapped() const {return m_map != NULL;}

  uint32_t getN0() const {return m_N0;}
  uint32_t getN() const {return m_N;}
  uint32_t getM0() const {return m_M0;}
  uint32_t getM() const {return m_M;}
  uint32_t getnumBytesofEachMarker0() const {return m_numBytesofEachMarker0;}
  uint32_t getnumBytesofEachMarker() const {return m_numBytesofEachMarker;}

  Rcpp::StringVector getChrVec() const
  {
    Rcpp::StringVector chrVec(m_M0);
    for(unsigned int i = 0; i < m_M0; i++){
      chrVec[i] = Rcpp::String(m_chr.at(i));
    }
    return chrVec;
  }

  void closegenofile();

};

}

#endif
