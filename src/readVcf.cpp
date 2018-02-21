#include <Rcpp.h>     // Rcpp-specific part
#include <string>     // to use std::string
#include <vector>     // to use std::vector
#include <zlib.h>     // to read gzip file

using namespace Rcpp; // to avoid typing Rcpp::
using namespace std;  // to avoid typing std::

#define MAX_LINE_LENGTH 1000000 // maximum line length

// C++ function to tokenize a string
std::vector<std::string> tokenizeLine(std::string line) {
  std::vector<std::string> toks;
  int pos = 0, nextpos = 0;
  while( (nextpos = line.find_first_of("\t\n\r", pos)) != string::npos ) {
    toks.push_back( line.substr(pos,nextpos-pos) );
    pos = nextpos+1;
  }
  if ( pos < (int)line.size() ) 
    toks.push_back(line.substr(pos));
  return toks;
}

//' Read VCF compressed with gzip/bgzip
//'
//' This function reads a compressed gzip/bgzip files and returns 
//' a genotype matrix of size [# variants] x [# samples].
//' 
//' @param filename A file name of a compressed VCF file ending with .gz or .bgz
//' @return A integer genotype matrix of a size [# variants] x [# samples]. 
//'         Only GT field will be recognized. The row names will be variant IDs, 
//'         in the format of [CHROM]:[POS]:[REF]:[ALT], and the column names 
//'         will be the sample IDs
//' @examples
//' fpath = system.file("extdata",
//'   "1000G_exome_chr20_example_softFiltered.calls.vcf.gz",
//'   package="tinyVcfReader")
//' mat = readVcf(fpath)
//' dim(mat)
//' @export
// [[Rcpp::export]]
IntegerMatrix readVcf(std::string filename) {
  gzFile file = gzopen(filename.c_str(), "r");
  if (!file) stop("Cannot open VCF file");
  
  char line[MAX_LINE_LENGTH];
  
  StringVector sampleIDs;
  StringVector markerIDs;
  vector<int> genotypes;
  int nvar = 0;
  while( gzgets(file, line, MAX_LINE_LENGTH) != NULL ) {
    if ( line[0] == '#' ) { // header line
      if ( line[1] == '#' ) continue; // ignore meta line
      vector<string> toks = tokenizeLine(line);
      for(int i=9; i < toks.size(); ++i)
        sampleIDs.push_back(toks[i]);
    }
    else {
      vector<std::string> toks = tokenizeLine(line);
      if ( sampleIDs.size() + 9 != toks.size() )
        stop("The number of lines do not match");
      markerIDs.push_back( toks[0] + ":" + toks[1] + ":" + toks[3] + ":" + toks[4] );
      for(int i=9; i < (int)toks.size(); ++i) {
        if ( toks[i][0] == '.' ) genotypes.push_back(-9);
        else  // assume that allele is 0 to 9
          genotypes.push_back((toks[i][0] - '0') + (toks[i][2] - '0')); 
      }
      ++nvar;
    }
  }
  int nsamples = (int)sampleIDs.size();
  IntegerMatrix mat(nvar,nsamples);
  for(int i=0, k=0; i < nvar; ++i) {
    for(int j=0; j < nsamples; ++j, ++k) {
      if ( genotypes[k] < 0 )
        mat(i,j) = NumericVector::get_na();
      else
        mat(i,j) = genotypes[k];
    }
  }
  rownames(mat) = markerIDs;
  colnames(mat) = sampleIDs;
  return mat;
}