#include <string>
#include <iostream>
#include <fstream>
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
    return 0;
  } // good call
} // int main
