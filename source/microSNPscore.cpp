#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include "mRNA.h"
#include "miRNA.h"
#include "sequenceFile.h"
#include "alignment.h"
#include "SNP.h"
#include "conservationList.h"
#include "filepath.h"

using namespace microSNPscore;

int main(int argc, char * argv[])
{
   /*******************************\ 
  | Define help and usage messages: |
   \*******************************/
  const std::string usage(std::string(argv[0])+" [mRNA file] [miRNA file] [conservation file] [SNP file] [prediction file]\n");
  const std::string help(usage);
   /*****************************\ 
  | Parse command line arguments: |
   \*****************************/
  if(argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) // help requested
  {
    std::cout << help;
    return 1;
  } // help requested
  else if(argc < 6) // bad call
  {
    std::cerr << usage;
    return 0;
  } // bad call
  else // good call
  {
    std::string mRNA_file_path(argv[1]);
    std::string miRNA_file_path(argv[2]);
    std::string conservation_file_path(argv[3]);
    std::string SNP_file_path(argv[4]);
    std::string prediction_file_path(argv[5]);
     /***************************\ 
    | Read data from input files: |
     \***************************/
    sequenceFile mRNA_file(mRNA_file_path);
    mRNA_file.read();
    sequenceFile miRNA_file(miRNA_file_path);
    miRNA_file.read();
    conservationList conservation_list(conservation_file_path);
    std::map<sequenceID,mRNA> mRNAs;
    for(sequenceFile::const_iterator mRNA_it(mRNA_file.begin());mRNA_it!=mRNA_file.end();++mRNA_it)
    {
      mRNA the_mRNA(mRNA_it->get_mRNA(conservation_list));
      mRNAs.insert(std::pair<sequenceID,mRNA>(the_mRNA.get_ID(),the_mRNA));
    }
    std::map<sequenceID,miRNA> miRNAs;
    for(sequenceFile::const_iterator miRNA_it(miRNA_file.begin());miRNA_it!=miRNA_file.end();++miRNA_it)
    {
      miRNA the_miRNA(miRNA_it->get_miRNA(conservation_list));
      miRNAs.insert(std::pair<sequenceID,miRNA>(the_miRNA.get_ID(),the_miRNA));
    }
    std::map<SNPID,SNP> SNPS;
    return 0;
  } // good call
} // int main
