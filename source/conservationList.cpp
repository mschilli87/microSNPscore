
#include <iostream>
// for std::cerr and std::endl (error stating)
#include <algorithm>
// for std::lower_bound (binary search for ranges)
#include <fstream>
// for std::ifstream (file reading)
#include <regex.h>
// for regex_t, regmatch_t, regcomp and regexec (regular expressions)
#include <sstream>
// for std::istringstream (type conversion)

#include "conservationList.h"

namespace microSNPscore {

    /*****************************************************************//**
    * @brief constructor - Do not call without parameter values!
    *
    * This is used to create an instance of the class conservation range.
    * The default values are not intended to be used directly.
    * They are only provided to allow array allocation but you will need
    * to assign a valid object created by giving those parameters a value
    * to actually use it. This is done by containers like std::vector and
    * the reason for providing those default values is to allow using
    * containers containing objects of this class.
    *
    * @param the_chromosome (pseudo-optional) chromosomeType that
    *     represents the the chromosome the conservation range is located
    *     on - Defaults to ""
    * @param the_start (pseudo-optional) chromosomePosition that
    *     represents the position on the chromosome (the 5' end of the +
    *     strand (i.e. the 3' end of the - strand) beeing position 1)
    *     where the conservation range starts - Defaults to 0
    * @param the_score (pseudo-optional) conservationScore that represents
    *     the score of the conservation range - Defaults to 0
    *
    * @return a conservation range with the given score located at
    *     the given position on the given chromosome
    *********************************************************************/
    
    conservationRange::conservationRange(chromosomeType the_chromosome, chromosomePosition the_start, conservationScore the_score)
    :chromosome(the_chromosome),start(the_start),score(the_score) {
}

    /*****************************************************************//**
    * @brief constructor
    *
    * This is used to create an instance of the class conservationList
    * from a given file.
    * The file should contain a tab-separated table with the columns
    * chromosome name, start position and conservation score for each
    * range sorted ascending by chromosome name and start position.
    * If the file does not exist (or cannot be opened for reading) an
    * error is raised and an empty list is created.
    * If a line does not match the format or is not in order, an error is
    * raised and the line is ignored.
    *
    * @param conservation_file file path of the input file
    *
    * @return a conservationList containing the ranges given in the file
    *********************************************************************/
    conservationList::conservationList(const filePath & conservation_file) {
       /********************************************************\ 
      | Try to open an input file stream associated to the given |
      | file path stating an error in the case of failure:       |
       \********************************************************/
      std::ifstream infile(conservation_file.c_str());
      if(infile.fail())
      {
        std::cerr << "microSNPscore::conservationList::conservationList\n";
        std::cerr << " ==> Cannot open file to read from: ";
        std::cerr << conservation_file << std::endl;
        std::cerr << "  --> no conservations will be read from the file\n";
      }
      else
      {
         /**************************************************************\ 
        | Try to initialize extended regular expression matching a valid |
        | conservation range line stating error in case of failure:      |
         \**************************************************************/
        regex_t line_regex;
        char line_pattern[] = "^([[:alnum:]]+)\t([[:digit:]]+)\t([[:digit:].]+)$";
        if(regcomp(&line_regex,line_pattern,REG_EXTENDED) != 0)
        {
          std::cerr << "microSNPscore::conservationList::conservationList\n";
          std::cerr << " ==> compiling line regular expression failed\n";
          std::cerr << "  --> no conservations will be read from the file\n";
        }
        else
        {
           /**************************************************************\ 
          | Read the content of the file linewise into a string and try to |
          | match the initialized regular expression stating error in case |
          | of failure:                                                    |
           \**************************************************************/
          std::string line_string;
          size_t line_nmatch(3);
          regmatch_t line_pmatch[line_nmatch];
          while(getline(infile,line_string).good())
          {
            if(regexec(&line_regex,line_string.c_str(),line_nmatch,line_pmatch,0) != 0)
            {
                  std::cerr << "microSNPscore::conservationList::conservationList\n";
                  std::cerr << " ==> no valid conservation range:\n";
                  std::cerr << line_string << std::endl;
                  std::cerr << "  --> omitting line\n";
            }
            else
            {
               /*****************************************************************\ 
              | Extract the subsequences matching the regular expression's groups |
              | from the line assigning them to the corresponding parameters      |
              | (converting them via a stringstream) create conservation range:   |
               \*****************************************************************/
              chromosomeType line_chromosome(line_string.substr(line_pmatch[1].rm_so,line_pmatch[1].rm_eo-line_pmatch[1].rm_so));
              chromosomePosition line_start;
              std::istringstream stream_start(line_string.substr(line_pmatch[2].rm_so,line_pmatch[2].rm_eo-line_pmatch[2].rm_so));
              stream_start >> line_start;
              conservationScore line_score;
              std::istringstream stream_score(line_string.substr(line_pmatch[3].rm_so,line_pmatch[3].rm_eo-line_pmatch[3].rm_so));
              stream_score >> line_score;
              conservationRange line_range(line_chromosome,line_start,line_score);
               /*******************************************************\ 
              | If this is not the first line check whether is in order |
              | with its predecessor stating an error if not:           |
               \*******************************************************/
              if(ranges.begin() != ranges.end() && line_range <= *(ranges.end()-1))
              {
                  std::cerr << "microSNPscore::conservationList::conservationList\n";
                  std::cerr << " ==> conservation range out of order:\n";
                  std::cerr << line_string << std::endl;
                  std::cerr << "  --> omitting line\n";
              }
              else
              {
                 /***********************************\ 
                | Append the new range to the vector: |
                 \***********************************/
                ranges.push_back(line_range);
              } // ranges.begin() == ranges.end() || line_range > *(ranges.end()-1)
            } // regexec(&line_regex,line_string.c_str(),line_nmatch,line_pmatch,0) == 0
          } // getline(infile,line_string).good()
        } // regcomp(&line_regex,line_pattern,REG_EXTENDED) == 0
      } // !infile.fail()
}

    /*****************************************************************//**
    * @brief conservation list begin
    *
    * This is used to get the first conservation range in the list.
    *
    * @return const_iterator pointing to the first conservation range
    *********************************************************************/
    conservationList::const_iterator conservationList::begin() const {
      return ranges.begin();
}

    /*****************************************************************//**
    * @brief conservation list end
    *
    * This is used to get the end of the conservation list.
    *
    * @return const_iterator pointing behind the last conservation range
    *********************************************************************/
    conservationList::const_iterator conservationList::end() const {
      return ranges.end();
}

    /*****************************************************************//**
    * @brief get conservation score by chromosome and position
    *
    * This method is used to access the conservation score of a given
    * position on a given chromosome.
    * If the chromosome is unknown, an error is raised and 0 is returned.
    *
    * @param chromosome the chromosome the position of interest is
    * located on
    * @param position the position of interest on the given chromosome
    *     (the 5' end of the + strand (i.e. the 3' end of the - end)
    *     beeing position 1)
    *
    * @return the conservation score of the given position on the given
    *     chromosome
    *********************************************************************/
    
    conservationScore conservationList::get_score(const chromosomeType & chromosome, const chromosomePosition & position) const {
       /************************************************************\ 
      | Search for the first range not beginning before the searched |
      | position and return its score if it is a perfect match or    |
      | the one of its predecessor (if one exists):                  |
       \************************************************************/
      const std::vector<conservationRange>::const_iterator possible_border(std::lower_bound(ranges.begin(),ranges.end(),conservationRange(chromosome,position)));
      if(possible_border!=ranges.end())
      {
        if(possible_border->get_chromosome() == chromosome && possible_border->get_start() == position)
        {
          return possible_border->get_score();
        }
        else if(possible_border != ranges.begin() && (possible_border-1)->get_chromosome() == chromosome)
        {
          return (possible_border-1)->get_score();
        }
      } // else case for both if statements:
       /***************************************************\ 
      | If no score is defined state an error and return 0: |
       \***************************************************/
      std::cerr << "microSNPscore::conservationList::get_score\n";
      std::cerr << " ==> Unkown chromosome: " << chromosome << std::endl;
      std::cerr << "  --> assuming zero conservation\n";
      return 0;
      
}

/*****************************************************************//**
* @brief output stream conservation range insertion operator
*
* This operator is used to insert a conservation range to an output
* stream (e.g. to print it on screen).
* The conservation range will be represented by the tab separated
* values chromosome, start and score.
*
* @param the_stream output stream the nucleotide should be inserted in
* @param the_range conservation range to be inserted in the output
*     stream
*
* @return output stream with the inserted conservation range
*********************************************************************/
std::ostream & operator<<(std::ostream & the_stream, const conservationRange & the_range)
{
  the_stream << the_range.get_chromosome() << '\t';
  the_stream << the_range.get_start() << '\t';
  the_stream << the_range.get_score();
  return the_stream;
}

/*****************************************************************//**
* @brief output stream conservation list insertion operator
*
* This operator is used to insert a conservation list to an output
* stream (e.g. to print it on screen).
* The conservation list will be represented by its ranges each one on
* a single line.
*
* @param the_stream output stream the nucleotide should be inserted in
* @param the_list conservation list to be inserted in the output
*     stream
*
* @return output stream with the inserted conservation list
*********************************************************************/
std::ostream & operator<<(std::ostream & the_stream, const conservationList & the_list)
{
  for(conservationList::const_iterator range_it(the_list.begin());range_it!=the_list.end();++range_it)
  {
    the_stream << *range_it << std::endl;
  }
  return the_stream;
}

} // namespace microSNPscore
