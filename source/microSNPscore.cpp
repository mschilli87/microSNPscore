#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <regex.h>
#include "mRNA.h"
#include "miRNA.h"
#include "sequenceFile.h"
#include "alignment.h"
#include "SNP.h"
#include "conservationList.h"
#include "filepath.h"

using namespace microSNPscore;

struct prediction
{
  sequenceID         miRNA;
  sequenceID         mRNA;
  chromosomePosition three_prime;
  SNPID              SNP;
};

void read_sequences(std::map<sequenceID,mRNA> & mRNA_map,filePath mRNA_path,
                    std::map<sequenceID,miRNA> & miRNA_map,filePath miRNA_path,
                    const conservationList & conservations)
{
     /****************************************************************\ 
    | Read the given files and insert the corresponing sequences into  |
    | their maps reusing the sequenceFile object to hold only one file |
    | at a time in memory:                                             |
     \****************************************************************/
    sequenceFile the_file(mRNA_path);
    the_file.read();
    for(sequenceFile::const_iterator mRNA_it(the_file.begin());mRNA_it!=the_file.end();++mRNA_it)
    {
      mRNA the_mRNA(mRNA_it->get_mRNA(conservations));
      mRNA_map.insert(std::pair<sequenceID,mRNA>(the_mRNA.get_ID(),the_mRNA));
    }
    the_file=sequenceFile(miRNA_path);
    the_file.read();
    for(sequenceFile::const_iterator miRNA_it(the_file.begin());miRNA_it!=the_file.end();++miRNA_it)
    {
      miRNA the_miRNA(miRNA_it->get_miRNA(conservations));
      miRNA_map.insert(std::pair<sequenceID,miRNA>(the_miRNA.get_ID(),the_miRNA));
    }
}

void read_SNPs(std::map<SNPID,SNP> & map, filePath path)
{
   /********************************************************\ 
  | Try to open an input file stream associated to the given |
  | file path stating an error in the case of failure:       |
   \********************************************************/
  std::ifstream file(path.c_str());
  if(file.fail())
  {
    std::cerr << "microSNPscore::read_SNPs\n";
    std::cerr << " ==> Cannot open file to read from: ";
    std::cerr << path << std::endl;
    std::cerr << "  --> no SNPs will be read from the file\n";
  }
  else
  {
     /********************************************************\ 
    | Try to initialize extended regular expression matching a |
    | valid SNP file line stating error in case of failure:    |
     \********************************************************/
    regex_t line_regex;
    char line_pattern[] = "^(.+)\t(.+)\t(.*)\t(.+)\t([-+])\t([[:digit:]]+)$";
    int error_code = regcomp(&line_regex,line_pattern,REG_EXTENDED);
    if(error_code != 0)
    {
      std::cerr << "microSNPscore::read_SNPs\n";
      std::cerr << " ==> compiling line regular expression failed:\n";
      const size_t error_len(regerror(error_code,&line_regex,NULL,0));
      char error_message[error_len];
      regerror(error_code,&line_regex,error_message,error_len);
      std::cerr << error_message << std::endl;
      std::cerr << "  --> no SNPs will be read from the file\n";
    }
    else
    {
       /**************************************************************\ 
      | Read the content of the file linewise into a string and try to |
      | match the initialized regular expression stating error in case |
      | of failure:                                                    |
       \**************************************************************/
      std::string line_string;
      size_t line_nmatch(7);
      regmatch_t line_pmatch[line_nmatch];
      while(getline(file,line_string).good())
      {
        error_code = regexec(&line_regex,line_string.c_str(),line_nmatch,line_pmatch,0);
        if(error_code != 0)
        {
              std::cerr << "microSNPscore::read_SNPs\n";
              std::cerr << " ==> no valid SNP file line:\n";
              std::cerr << line_string << std::endl;
              std::cerr << "     error message:\n";
              const size_t error_len(regerror(error_code,&line_regex,NULL,0));
              char error_message[error_len];
              regerror(error_code,&line_regex,error_message,error_len);
              std::cerr << error_message << std::endl;
              std::cerr << "  --> omitting line\n";
        }
        else
        {
           /*****************************************************************\ 
          | Extract the subsequences matching the regular expression's groups |
          | from the line assigning them to the corresponding parameters      |
          | (converting them via a stringstream) to create a SNP:             |
           \*****************************************************************/
          SNPID ID(line_string.substr(line_pmatch[1].rm_so,line_pmatch[1].rm_eo-line_pmatch[1].rm_so));
          std::string reference(line_string.substr(line_pmatch[2].rm_so,line_pmatch[2].rm_eo-line_pmatch[2].rm_so));
          std::string alternative(line_string.substr(line_pmatch[3].rm_so,line_pmatch[3].rm_eo-line_pmatch[3].rm_so));
          chromosomeType chromosome(line_string.substr(line_pmatch[4].rm_so,line_pmatch[4].rm_eo-line_pmatch[4].rm_so));
          strandType strand(line_string.substr(line_pmatch[5].rm_so,line_pmatch[5].rm_eo-line_pmatch[5].rm_so) == "+" ? Plus : Minus);
          chromosomePosition position;
          std::istringstream stream_position(line_string.substr(line_pmatch[6].rm_so,line_pmatch[6].rm_eo-line_pmatch[6].rm_so));
          stream_position >> position;
          SNP line_SNP(ID,reference,alternative,chromosome,strand,position);
           /******************************\ 
          | Insert the new SNP in the map: |
           \******************************/
          map.insert(std::pair<SNPID,SNP>(ID,line_SNP));
        } // regexec(&line_regex,line_string.c_str(),line_nmatch,line_pmatch,0) == 0
      } // getline(file,line_string).good()
    } // regcomp(&line_regex,line_pattern,REG_EXTENDED) == 0
  } // !file.fail()
}

void read_predictions(std::vector<prediction> & vector,filePath path)
{
   /********************************************************\ 
  | Try to open an input file stream associated to the given |
  | file path stating an error in the case of failure:       |
   \********************************************************/
  std::ifstream file(path.c_str());
  if(file.fail())
  {
    std::cerr << "microSNPscore::read_predictions\n";
    std::cerr << " ==> Cannot open file to read from: ";
    std::cerr << path << std::endl;
    std::cerr << "  --> no predictions will be read from the file\n";
  }
  else
  {
     /**************************************************************\ 
    | Try to initialize extended regular expression matching a valid |
    | prediction file line stating error in case of failure:         |
     \**************************************************************/
    regex_t line_regex;
    char line_pattern[] = "^([^\t]+)\t([^\t]+)\t([[:digit:]]+)\t([^\t]+)$";
    int error_code = regcomp(&line_regex,line_pattern,REG_EXTENDED);
    if(error_code != 0)
    {
      std::cerr << "microSNPscore::read_predictions\n";
      std::cerr << " ==> compiling line regular expression failed:\n";
              const size_t error_len(regerror(error_code,&line_regex,NULL,0));
              char error_message[error_len];
              regerror(error_code,&line_regex,error_message,error_len);
              std::cerr << error_message << std::endl;
      std::cerr << "  --> no predictions will be read from the file\n";
    }
    else
    {
       /**************************************************************\ 
      | Read the content of the file linewise into a string and try to |
      | match the initialized regular expression stating error in case |
      | of failure:                                                    |
       \**************************************************************/
      std::string line_string;
      size_t line_nmatch(5);
      regmatch_t line_pmatch[line_nmatch];
      while(getline(file,line_string).good())
      {
        error_code = regexec(&line_regex,line_string.c_str(),line_nmatch,line_pmatch,0);
        if(error_code != 0)
        {
              std::cerr << "microSNPscore::read_predictions\n";
              std::cerr << " ==> no valid prediction file line:\n";
              std::cerr << line_string << std::endl;
              std::cerr << "     error message:\n";
              const size_t error_len(regerror(error_code,&line_regex,NULL,0));
              char error_message[error_len];
              regerror(error_code,&line_regex,error_message,error_len);
              std::cerr << error_message << std::endl;
              std::cerr << "  --> omitting line\n";
        }
        else
        {
           /*****************************************************************\ 
          | Extract the subsequences matching the regular expression's groups |
          | from the line assigning them to the corresponding parameters      |
          | (converting them via a stringstream) to create a prediction:      |
           \*****************************************************************/
          sequenceID miRNA(line_string.substr(line_pmatch[1].rm_so,line_pmatch[1].rm_eo-line_pmatch[1].rm_so));
          sequenceID mRNA(line_string.substr(line_pmatch[2].rm_so,line_pmatch[2].rm_eo-line_pmatch[2].rm_so));
          chromosomePosition three_prime;
          std::istringstream stream_three_prime(line_string.substr(line_pmatch[3].rm_so,line_pmatch[3].rm_eo-line_pmatch[3].rm_so));
          stream_three_prime >> three_prime;
          SNPID SNP(line_string.substr(line_pmatch[4].rm_so,line_pmatch[4].rm_eo-line_pmatch[4].rm_so));
          prediction line_prediction;
          line_prediction.miRNA = miRNA;
          line_prediction.mRNA = mRNA;
          line_prediction.three_prime = three_prime;
          line_prediction.SNP = SNP;
           /****************************************\ 
          | Append the new prediction to the vector: |
           \****************************************/
          vector.push_back(line_prediction);
        } // regexec(&line_regex,line_string.c_str(),line_nmatch,line_pmatch,0) == 0
      } // getline(file,line_string).good()
    } // regcomp(&line_regex,line_pattern,REG_EXTENDED) == 0
  } // !file.fail()
}

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
    conservationList conservation_list(conservation_file_path);
    std::map<sequenceID,mRNA> mRNAs;
    std::map<sequenceID,miRNA> miRNAs;
    std::map<SNPID,SNP> SNPs;
    std::vector<prediction> predictions;
    read_sequences(mRNAs,mRNA_file_path,miRNAs,miRNA_file_path,conservation_list);
    read_SNPs(SNPs,SNP_file_path);
    read_predictions(predictions,prediction_file_path);
     /**************************************************\ 
    | Iterate over the predictions printing the original |
    | line followed by the deregulation score:           |
     \**************************************************/
    for(std::vector<prediction>::const_iterator prediction_it(predictions.begin());prediction_it!=predictions.end();++prediction_it)
    {
      std::cout << prediction_it->miRNA << '\t';
      std::cout << prediction_it->mRNA << '\t';
      std::cout << prediction_it->three_prime << '\t';
      std::cout << prediction_it->SNP << '\t';
      std::cout << SNPs[prediction_it->SNP].get_deregulation_score(miRNAs[prediction_it->miRNA],
                                                                   mRNAs[prediction_it->mRNA],
                                                                   prediction_it->three_prime) << std::endl;
    }
    return 0;
  } // good call
} // int main
