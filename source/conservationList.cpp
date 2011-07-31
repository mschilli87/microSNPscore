
#include <iostream>
// for std::cerr and std::endl (error stating)
#include <algorithm>
// for std::lower_bound (binary search for ranges)
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
    
    conservationScore conservationList::get_score(const chromosomeType & chromosome, const chromosomePosition & position) {
       /************************************************************\ 
      | Search for the first range not beginning before the searched |
      | and return its score if it is a perfect macht or the one     |
      | position of its predecessor (if one exists):                 |
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


} // namespace microSNPscore
