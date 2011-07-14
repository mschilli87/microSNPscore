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
    * @brief standard constructor
    *
    * This is used to create an instance of the class alignmentColumn.
    *
    * @return an empty alignment column
    *********************************************************************/
    alignmentColumn();

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


  private:
    /*****************************************************************//**
    * @brief messenger RNA nucleotide
    *
    * This is the nucleotide of the messenger RNA that is aligned in that
    * alignment column
    *********************************************************************/
    nucleotide mRNA_nucleotide;

    /*****************************************************************//**
    * @brief microRNA nucleotide
    *
    * This is the nucleotide of the microRNA that is aligned in that
    * alignment column
    *********************************************************************/
    nucleotide miRNA_nucleotide;

    /*****************************************************************//**
    * @brief match state
    *
    * This is the match state (Indel, Mismatch, Wobble, Match) of the
    * alignment column
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
