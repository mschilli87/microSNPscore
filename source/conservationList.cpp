
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
    * @brief add conservation range to the list
    *
    * This method is used to insert a new conservation range to the list.
    *
    * @param conservation_range conservationRange to insert in the list
    *********************************************************************/
    void conservationList::add(const conservationRange & conservation_range) {
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
    
    conservationScore conservationList::get_range(const chromosomeType & chromosome, const chromosomePosition & position) {
}


} // namespace microSNPscore
