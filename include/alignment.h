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
    * @brief constructor
    *
    * This is used to create an instance of the class alignmentColumn.
    *
    * @param the_mRNA_nucleotide const nucleotide reference to the
    *        nucleotide of the messenger RNA in that alignment column
    * @param the_miRNA_nucleotide const nucleotide reference to the
    *        nucleotide of the microRNA in that alignment column
    * @param position (optional) matchPosition indicating whether the
    *     match occurs in the seed (Seed) or in the 3' region of the miRNA
    *     (ThreePrime) - Defaults to ThreePrime
    * @param indel_type (optional) indelType telling whether this match
    *     would be the first continuing indel (i.e. Open) if it would be
    *     an indel or if it is directly following an existing indel (i.e.
    *     Extend) - Defaults to Open
    * @return alignmentColumn representing an alignment column aligning
    *         the given nucleotides
    *********************************************************************/
    alignmentColumn(const nucleotide & the_mRNA_nucleotide, const nucleotide & the_miRNA_nucleotide, matchPosition position = ThreePrime, IndelType indel_type = Open);

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

    inline const matchType get_match() const;


  private:
    /*****************************************************************//**
    * @brief messenger RNA nucleotide
    *
    * This is the nucleotide of the messenger RNA that is aligned in that
    * alignment column
    *********************************************************************/
    const nucleotide mRNA_nucleotide;

    /*****************************************************************//**
    * @brief microRNA nucleotide
    *
    * This is the nucleotide of the microRNA that is aligned in that
    * alignment column
    *********************************************************************/
    const nucleotide miRNA_nucleotide;

    /*****************************************************************//**
    * @brief match state
    *
    * This is the match state (Indel, Mismatch, Wobble, Match) of the
    * alignment column
    *********************************************************************/
    const matchType match;

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

inline const matchType alignmentColumn::get_match() const {
  return match;
}

/*****************************************************************//**
* @brief mRNA:miRNA alignment class
*
* This represents an alignment of mRNA and miRNA.
*********************************************************************/
class alignment {
  public:
    /*****************************************************************//**
    * @brief standard constructor
    *
    * This is used to create an instance of the class alignment.
    *
    * @return an empty alignment
    *********************************************************************/
    alignment();


  private:
    /*****************************************************************//**
    * @brief alignment columns
    *
    * A vector containing the alignment's columns from miRNA 5' to 3'
    *********************************************************************/
    std::vector<alignmentColumn> columns;

};

} // namespace microSNPscore
#endif
