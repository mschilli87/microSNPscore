#ifndef MICROSNPSCAN_ALIGNMENT_H
#define MICROSNPSCAN_ALIGNMENT_H


#include "sequence.h"

namespace microSNPscan {

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

    /*****************************************************************//**
    * @brief destructor
    *
    * This is used to remove an instance of the alignment class.
    *
    *********************************************************************/
    ~alignment();

    /*****************************************************************//**
    * @brief const copy constructor
    *
    * This is used to copy an instance of the alignment class.
    *
    * @param source alignment object to copy
    * @return a copy of the source alignment object
    *********************************************************************/
    alignment(const alignment & source);

    /*****************************************************************//**
    * @brief assignment operator
    *
    * This is used to copy an instance of  the alignment class by
    * assigning it to another one
    *
    * @param source the alignment object to be copied
    * @return a copy of the source alignment object
    *********************************************************************/
    alignment & operator=(const alignment & source);


  private:
    /*****************************************************************//**
    * @brief alignment columns
    *
    * A vector containing the alignment's columns from miRNA 5' to 3'
    *********************************************************************/
    std::vector<alignmentColumn> columns;

};
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
    * @brief destructor
    *
    * This is used to remove an instance of the alignmentColumn class.
    *
    *********************************************************************/
    ~alignmentColumn();

    /*****************************************************************//**
    * @brief const copy constructor
    *
    * This is used to copy an instance of the alignmentColumn class.
    *
    * @param source alignmentColumn object to copy
    * @return a copy of the source alignmentColumn object
    *********************************************************************/
    alignmentColumn(const alignmentColumn & source);

    /*****************************************************************//**
    * @brief assignment operator
    *
    * This is used to copy an instance of  the alignmentColumn class by
    * assigning it to another one
    *
    * @param source the alignmentColumn object to be copied
    * @return a copy of the source alignmentColumn object
    *********************************************************************/
    alignmentColumn & operator=(const alignmentColumn & source);

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


} // namespace microSNPscan
#endif
