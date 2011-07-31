
#include <iostream>
// for std::cerr and std::endl (error stating)
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
}

    /*****************************************************************//**
    * @brief get conservation range by chromosome and position
    *
    * This method is used to access the conservation range of a given
    * position on a given chromosome.
    * If the list is empty its end will be returned.
    * If the chromosome is unknown, the last range of the preceeding 
    * chromosome in order is returned.
    *
    * @param chromosome the chromosome the position of interest is
    * located on
    * @param position the position of interest on the given chromosome
    *     (the 5' end of the + strand (i.e. the 3' end of the - end)
    *     beeing position 1)
    *
    * @return iterator pointing to the conservation range of the given
    *     position on the given chromosome
    *********************************************************************/
    
    conservationRange conservationList::get_range(const chromosomeType & chromosome, const chromosomePosition & position) {
       /**************************\ 
      | Initialize lookup borders: |
       \**************************/
      iterator begin = ranges.begin();
      iterator end = ranges.end();
       /**********************************************\ 
      | As long as more than one position is possible: |
       \**********************************************/
      while(end-begin>1)
      {
         /***************************************\ 
        | Look at the to elements in the middle:  |
         \***************************************/
        const iterator mid = begin + ((sequenceLength)(end-begin-1)/2);
        const chromosomeType mid_chr = mid->get_chromosome();
        const chromosomePosition mid_start = mid->get_start();
        const iterator next = mid+1;
        const chromosomeType next_chr = next->get_chromosome();
        const chromosomePosition next_start = next->get_start();
         /************************************************************\ 
        | If the chromosome does not exist or the position is between  |
        | the two elements it has to be the left middle element:       |
         \************************************************************/
        if((mid_chr < chromosome && next_chr > chromosome) ||
           (mid_chr == chromosome && mid_start <= position &&
           (next_chr > chromosome || next_start > position))) // mid <= position && next > position
        {
          begin=mid;
          end=next;
        } // mid =< position && next > position
         /***********************************\ 
        | If even the right middle element is |
        | too the target cannot be before it: |
         \***********************************/
        else if(next_chr < chromosome || (next_chr == chromosome && next_start < position)) // next < position
        {
          begin=next;
        } // next < position
         /*********************************************\ 
        | If the target is not at nor right of the left |
        | middle element, it has to be in the left :    |
         \*********************************************/
        else // mid > position
        {
          end=mid;
        }// mid > position
      } //begin != end
       /****************************************\ 
      | In the end return the remaining element: |
       \****************************************/
      return begin;
}


} // namespace microSNPscore
