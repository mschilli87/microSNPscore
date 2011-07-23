#ifndef MICROSNPSCORE_ALIGNMENT_H
#define MICROSNPSCORE_ALIGNMENT_H


#include "nucleotide.h"
#include <vector>

namespace microSNPscore {

/*****************************************************************//**
* @brief Alignment column class
*
* This represent a column of an alignment of mRNA and miRNA.
*********************************************************************/
class alignmentColumn {
  public:
    /*****************************************************************//**
    * @brief constructor - Do not use without parameter values!
    *
    * This is used to create an instance of the class alignmentColumn.
    * The default values for @p the_mRNA_nucleotide and
    * @p the_miRNA_nucleotide are not intended to be used directly.
    * They are only provided to allow array allocation but you will need
    * to assign a valid object created by giving those parameters a value
    * to actually use it. This is done by containers like std::vector and
    * the reason for providing those default values is to allow using
    * containers containing objects of this class.
    *
    * @param the_mRNA_nucleotide (pseudo-optional) const nucleotide
    *     reference to the nucleotide of the messenger RNA in that
    *     alignment column - Defaults to uninitialized nucleotide
    *     constructed without parameters
    * @param the_miRNA_nucleotide (pseudo-optional) const nucleotide
    *     reference to the nucleotide of the microRNA in that alignment
    *     column - Defaults to uninitialized nucleotide constructed
    *     without parameters
    * @param position (optional) matchPosition indicating whether the
    *     match occurs in the seed (Seed) or in the 3' region of the miRNA
    *     (ThreePrime) - Defaults to ThreePrime
    * @param indel_type (optional) indelType telling whether this match
    *     would be the first continuing indel (i.e. Open) if it would be
    *     an indel or if it is directly following an existing indel (i.e.
    *     Extend) - Defaults to Open
    *
    * @return alignmentColumn representing an alignment column aligning
    *     the given nucleotides
    *
    * @see nucleotide::nucleotide(nucleoBase,sequencePosition,
    *    chromosomePosition)
    *********************************************************************/
    
    alignmentColumn(const nucleotide & the_mRNA_nucleotide = nucleotide(), const nucleotide & the_miRNA_nucleotide = nucleotide(), matchPosition position = ThreePrime, IndelType indel_type = Open);

    /*****************************************************************//**
    * @brief get method for messenger RNA nucleotide attribute
    *
    * This method is used to access the nucleotide of the messenger RNA
    * that is aligned in that alignment column.
    *
    * @return the mRNA nucleotide aligned in that alignment column
    *********************************************************************/
    inline const nucleotide get_mRNA_nucleotide() const;

    /*****************************************************************//**
    * @brief get method for microRNA nucleotide attribute
    *
    * This method is used to access the nucleotide of the microRNA
    * that is aligned in that alignment column.
    *
    * @return the miRNA nucleotide aligned in that alignment column
    *********************************************************************/
    inline const nucleotide get_miRNA_nucleotide() const;

    /*****************************************************************//**
    * @brief get method for microRNA nucleotide attribute
    *
    * This method is used to access the nucleotide of the microRNA
    * that is aligned in that alignment column.
    *
    * @return the miRNA nucleotide aligned in that alignment column
    *********************************************************************/
    inline const matchType get_match() const;


  private:
    /*****************************************************************//**
    * @brief messenger RNA nucleotide
    *
    * This is the nucleotide of the messenger RNA that is aligned in that
    * alignment column.
    * It should be const but because alignmentColumns shall be used in a
    * vector and std::vector tries to assign its elements to an internal
    * array it needs a working assignment operator which has to change the
    * object's members and therefore they cannot be declared const.
    * Nevertheless this attribute is not intended to be changed in any
    * other context than assigning an initialized object to an unitialized
    * one produced by the standard constructor which is in fact not
    * designed to be used directly but only provided to allow array
    * allocation which is needed to create containers, too.
    *
    * @see alignmentColumn()
    *********************************************************************/
    nucleotide mRNA_nucleotide;

    /*****************************************************************//**
    * @brief microRNA nucleotide
    *
    * This is the nucleotide of the microRNA that is aligned in that
    * alignment column.
    * It should be const but because alignmentColumns shall be used in a
    * vector and std::vector tries to assign its elements to an internal
    * array it needs a working assignment operator which has to change the
    * object's members and therefore they cannot be declared const.
    * Nevertheless this attribute is not intended to be changed in any
    * other context than assigning an initialized object to an unitialized
    * one produced by the standard constructor which is in fact not
    * designed to be used directly but only provided to allow array
    * allocation which is needed to create containers, too.
    *
    * @see alignmentColumn()
    *********************************************************************/
    nucleotide miRNA_nucleotide;

    /*****************************************************************//**
    * @brief match state
    *
    * This is the match state (IndelOpen, IndelExtend, Mismatch, Masked,
    * Wobble, Match) of the alignment column.
    * It should be const but because alignmentColumns shall be used in a
    * vector and std::vector tries to assign its elements to an internal
    * array it needs a working assignment operator which has to change the
    * object's members and therefore they cannot be declared const.
    * Nevertheless this attribute is not intended to be changed in any
    * other context than assigning an initialized object to an unitialized
    * one produced by the standard constructor which is in fact not
    * designed to be used directly but only provided to allow array
    * allocation which is needed to create containers, too.
    *
    * @see alignmentColumn()
    *********************************************************************/
    
    matchType match;

};
    /*****************************************************************//**
    * @brief get method for messenger RNA nucleotide attribute
    *
    * This method is used to access the nucleotide of the messenger RNA
    * that is aligned in that alignment column.
    *
    * @return the mRNA nucleotide aligned in that alignment column
    *********************************************************************/
    inline const nucleotide alignmentColumn::get_mRNA_nucleotide() const {
      return mRNA_nucleotide;
    }

    /*****************************************************************//**
    * @brief get method for microRNA nucleotide attribute
    *
    * This method is used to access the nucleotide of the microRNA
    * that is aligned in that alignment column.
    *
    * @return the miRNA nucleotide aligned in that alignment column
    *********************************************************************/
    inline const nucleotide alignmentColumn::get_miRNA_nucleotide() const {
      return miRNA_nucleotide;
    }

    /*****************************************************************//**
    * @brief get method for microRNA nucleotide attribute
    *
    * This method is used to access the nucleotide of the microRNA
    * that is aligned in that alignment column.
    *
    * @return the miRNA nucleotide aligned in that alignment column
    *********************************************************************/
    inline const matchType alignmentColumn::get_match() const {
      return match;
    }

/*****************************************************************//**
* @brief seed type type
*
* This represents the seed type (sixMer, sevenMerAOne, sevenMerMEight,
* eightMer) of a mRNA:miRNA alignement.
*********************************************************************/

enum seedType {
  sixMer,
  sevenMerAOne,
  sevenMerMEight,
  eightMer

};
/*****************************************************************//**
* @brief mRNA:miRNA alignment class
*
* This represents an alignment of mRNA and miRNA.
*********************************************************************/
class alignment {
  public:
    /*****************************************************************//**
    * @brief const iterartor type
    *
    * This type is used to access the alignment's columns.
    *********************************************************************/
    typedef std::vector<alignmentColumn>::const_iterator const_iterator;

    /*****************************************************************//**
    * @brief standard constructor
    *
    * This is used to create an instance of the class alignment.
    *
    * @return an empty alignment
    *********************************************************************/
    alignment();

    /*****************************************************************//**
    * @brief alignment begin
    *
    * This is used to get the first column (beeing the 5' end of the miRNA
    * and thus the 3' end of the mRNA) of the alignment.
    *
    * @return const_iterator pointing to the first alignment column
    *********************************************************************/
    inline const_iterator begin() const;

    /*****************************************************************//**
    * @brief alignment begin
    *
    * This is used to get the end of the alignment's column vector (beeing
    * the 3' end of the miRNA and thus the 5' end of the mRNA).
    *
    * @return const_iterator pointing behind the last alignment column
    *********************************************************************/
    inline const_iterator end() const;

    /*****************************************************************//**
    * @brief seed type calculation
    *
    * This method is used to calculate the alignment's seed type (sixMer,
    * sevenMerAOne, sevenMerMEight, eightMer).
    *
    * @return the seed type of the alignment
    *********************************************************************/
    seedType get_seed_type() const;


  private:
    /*****************************************************************//**
    * @brief alignment columns
    *
    * A vector containing the alignment's columns from miRNA 5' to 3'
    *********************************************************************/
    std::vector<alignmentColumn> columns;

};
    /*****************************************************************//**
    * @brief alignment begin
    *
    * This is used to get the first column (beeing the 5' end of the miRNA
    * and thus the 3' end of the mRNA) of the alignment.
    *
    * @return const_iterator pointing to the first alignment column
    *********************************************************************/
    inline alignment::const_iterator alignment::begin() const {
      return columns.end();
}

    /*****************************************************************//**
    * @brief alignment begin
    *
    * This is used to get the end of the alignment's column vector (beeing
    * the 3' end of the miRNA and thus the 5' end of the mRNA).
    *
    * @return const_iterator pointing behind the last alignment column
    *********************************************************************/
    inline alignment::const_iterator alignment::end() const {
      return columns.begin();
}

/*****************************************************************//**
* @brief output stream seed type insertion operator
*
* This operator is used to insert a seed type to an output stream
* (e.g. to print it on screen).
* The seed type will be represented by its name.
*
* @param the_stream output stream the seed type should be inserted in
* @param seed_type seedType to be inserted in the output stream
*
* @return output stream with the inserted seed type
*********************************************************************/
std::ostream & operator<<(std::ostream & the_stream, const seedType & seed_type);

} // namespace microSNPscore
#endif
