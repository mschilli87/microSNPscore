#ifndef MICROSNPSCORE_CONSERVATIONLIST_H
#define MICROSNPSCORE_CONSERVATIONLIST_H


#include "sequence.h"
#include "nucleotide.h"
#include "filepath.h"
#include <vector>

namespace microSNPscore {

/*****************************************************************//**
* @brief conservation score type
*
* This represents a score of a conservation range
*********************************************************************/
typedef double conservationScore;
/*****************************************************************//**
* @brief conservation range
*
* This represent a range on a chromosome with the same conservation.
*********************************************************************/
class conservationRange {
  public:
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
    
    conservationRange(chromosomeType the_chromosome = "", chromosomePosition the_start = 0, conservationScore the_score = 0);

    /*****************************************************************//**
    * @brief less than operator
    *
    * This returns if the range is located before another given range.
    *
    * @param the_range range to compare with
    *
    * @return true if this range is before the other one, false otherwise
    *********************************************************************/
    inline bool operator<(const conservationRange & the_range) const;

    /*****************************************************************//**
    * @brief less than operator
    *
    * This returns if the range is located before another given range.
    *
    * @param the_range range to compare with
    *
    * @return true if this range is before the other one, false otherwise
    *********************************************************************/
    inline bool operator <=(const conservationRange & the_range) const;

    /*****************************************************************//**
    * @brief get method for chromosome attribute
    *
    * This method is used to access the chromosome the conservation range
    * is located on.
    *
    * @return the chromosome of the conservation range
    *********************************************************************/
    inline const chromosomeType get_chromosome() const;

    /*****************************************************************//**
    * @brief get method for start attribute
    *
    * This method is used to access the position on the chromosome (the 5'
    * end of the + strand (i.e. the 3' end of the - strand) beeing
    * position 1) where the conservation range starts.
    *
    * @return the start position of the conservation range
    *********************************************************************/
    inline const chromosomePosition get_start() const;

    /*****************************************************************//**
    * @brief get method for score attribute
    *
    * This method is used to access the score of the conservation range.
    *
    * @return the score of the conservation range
    *********************************************************************/
    inline const conservationScore get_score() const;


  private:
    /*****************************************************************//**
    * @brief chromosome
    *
    * This is the chromosome the conservation range is located on.
    *********************************************************************/
    chromosomeType chromosome;

    /*****************************************************************//**
    * @brief start position
    *
    * This is the position on the chromosome (the 5' end of the + strand
    * (i.e. the 3' end of the - strand) beeing position 1) where the
    * conservation range starts.
    *********************************************************************/
    chromosomePosition start;

    /*****************************************************************//**
    * @brief conservation score
    *
    * This is the score of the conservation range.
    *********************************************************************/
    conservationScore score;

};
    /*****************************************************************//**
    * @brief less than operator
    *
    * This returns if the range is located before another given range.
    *
    * @param the_range range to compare with
    *
    * @return true if this range is before the other one, false otherwise
    *********************************************************************/
    inline bool conservationRange::operator<(const conservationRange & the_range) const {
      return get_chromosome()<the_range.get_chromosome() || (get_chromosome() == the_range.get_chromosome() && get_start() < the_range.get_start());
}

    /*****************************************************************//**
    * @brief less than operator
    *
    * This returns if the range is located before another given range.
    *
    * @param the_range range to compare with
    *
    * @return true if this range is before the other one, false otherwise
    *********************************************************************/
    inline bool conservationRange::operator <=(const conservationRange & the_range) const {
      return get_chromosome()<the_range.get_chromosome() || (get_chromosome() == the_range.get_chromosome() && get_start() <= the_range.get_start());
}

    /*****************************************************************//**
    * @brief get method for chromosome attribute
    *
    * This method is used to access the chromosome the conservation range
    * is located on.
    *
    * @return the chromosome of the conservation range
    *********************************************************************/
    inline const chromosomeType conservationRange::get_chromosome() const {
      return chromosome;
    }

    /*****************************************************************//**
    * @brief get method for start attribute
    *
    * This method is used to access the position on the chromosome (the 5'
    * end of the + strand (i.e. the 3' end of the - strand) beeing
    * position 1) where the conservation range starts.
    *
    * @return the start position of the conservation range
    *********************************************************************/
    inline const chromosomePosition conservationRange::get_start() const {
      return start;
    }

    /*****************************************************************//**
    * @brief get method for score attribute
    *
    * This method is used to access the score of the conservation range.
    *
    * @return the score of the conservation range
    *********************************************************************/
    inline const conservationScore conservationRange::get_score() const {
      return score;
    }

/*****************************************************************//**
* @brief conservation list
*
* This represent a searchable list range on a chromosome with the same
* conservation.
*********************************************************************/
class conservationList {
  public:
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
    conservationList(const filePath & conservation_file);

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
    
    conservationScore get_score(const chromosomeType & chromosome, const chromosomePosition & position) const;


  private:
    /*****************************************************************//**
    * @brief conservation range list
    *
    * This is a vector containing the conservation list's entries.
    *********************************************************************/
    std::vector<conservationRange> ranges;

};

} // namespace microSNPscore
#endif
