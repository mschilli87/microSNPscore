
#include "sequence.h"

namespace microSNPscan {

/*****************************************************************//**
* @brief standard constructor
*
* This is used to create an instance of the class sequence.
*
* @return an empty sequence
*********************************************************************/

sequence::sequence() {
}

/*****************************************************************//**
* @brief destructor
*
* This is used to remove an instance of the sequence class.
*
*********************************************************************/

sequence::~sequence() {
}

/*****************************************************************//**
* @brief const copy constructor
*
* This is used to copy an instance of the sequence class.
*
* @param source sequence object to copy
* @return a copy of the source sequence object
*********************************************************************/

sequence::sequence(const sequence & source) {
}

/*****************************************************************//**
* @brief assignment operator
*
* This is used to copy an instance of  the sequence class by
* assigning it to another one
*
* @param source the sequece object to be copied
* @return a copy of the source sequence object
*********************************************************************/

sequence & sequence::operator=(const sequence & source) {
}

/*****************************************************************//**
* @brief standard constructor
*
* This is used to create an instance of the class exon.
*
* @return an empty exon
*********************************************************************/
exon::exon() {
}

/*****************************************************************//**
* @brief destructor
*
* This is used to remove an instance of the exon class.
*
*********************************************************************/
exon::~exon() {
}

/*****************************************************************//**
* @brief const copy constructor
*
* This is used to copy an instance of the exon class.
*
* @param source exon object to copy
* @return a copy of the source exon object
*********************************************************************/
exon::exon(const exon & source) {
}

/*****************************************************************//**
* @brief assignment operator
*
* This is used to copy an instance of  the exon class by
* assigning it to another one
*
* @param source the exon object to be copied
* @return a copy of the source exon object
*********************************************************************/
exon & exon::operator=(const exon & source) {
}

/*****************************************************************//**
* @brief standard constructor
*
* This is used to create an instance of the class nucleotide.
*
* @return an empty exon
*********************************************************************/
nucleotide::nucleotide() {
}

/*****************************************************************//**
* @brief destructor
*
* This is used to remove an instance of the nucleotide class.
*
*********************************************************************/
nucleotide::~nucleotide() {
}

/*****************************************************************//**
* @brief const copy constructor
*
* This is used to copy an instance of the nucleotide class.
*
* @param source nucleotide object to copy
* @return a copy of the source nucleotide object
*********************************************************************/
nucleotide::nucleotide(const nucleotide & source) {
}

/*****************************************************************//**
* @brief assignment operator
*
* This is used to copy an instance of  the nucleotide class by
* assigning it to another one
*
* @param source the nucleotide object to be copied
* @return a copy of the source nucleotide object
*********************************************************************/
nucleotide & nucleotide::operator=(const nucleotide & source) {
}


} // namespace microSNPscan
