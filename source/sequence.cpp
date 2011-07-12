
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
* @brief get subsequence from chromosome position
*
* This method can be used to extract a subsequence of a given length
* starting (i.e. 5' end) at a given chromosome position from the
* sequence.
* If the length is too high so that the queried subsequence would
* reach over the (3') end of the sequence, a shorter subsequence
* starting at the disired start position and ending at the (3') end of
* the sequence will be returned.
* If the given chromosome position is not part of the sequence, an
* empty sequence will be returned.
*
* @param from the start position on the chromosome of the subsequence
* @param len the (maximal) length of the subsequence
* @return the subsequence starting at the given position and ending
*         after the given length (5' to 3') or at the end of the
*         sequence
*********************************************************************/

sequence sequence::get_subsequence_chr_from(const chromosomePosition & from, unsigned short len) {
}

/*****************************************************************//**
* @brief get subsequence from sequence position
*
* This method can be used to extract a subsequence of a given length
* starting (i.e. 5' end) at a given position from the sequence.
* If the length is too high so that the queried subsequence would
* reach over the (3') end of the sequence, a shorter subsequence
* starting at the disired start position and ending at the (3') end of
* the sequence will be returned.
* If the given position is not part of the sequence, an empty sequence
* will be returned.
*
* @param from the start position in the sequence of the subsequence
* @param len the (maximal) length of the subsequence
* @return the subsequence starting at the given position and ending
*         after the given length (5' to 3') or at the end of the
*         sequence
*********************************************************************/

sequence sequence::get_subsequence_from(const sequencePosition & from, unsigned short len) {
}

/*****************************************************************//**
* @brief get subsequence to chromosome position
*
* This method can be used to extract a subsequence of a given length
* ending (i.e. 3' end) at a given chromosome position from the
* sequence.
* If the length is too high so that the queried subsequence would
* reach over the (5') end of the sequence, a shorter subsequence
* ending at the disired end position and starting at the (5') end of
* the sequence will be returned.
* If the given chromosome position is not part of the sequence, an
* empty sequence will be returned.
*
* @param to the end position on the chromosome of the subsequence
* @param len the (maximal) length of the subsequence
* @return the subsequence ending at the given position and starting
*         after the given length (3' to 5') or at the start of the
*         sequence
*********************************************************************/

sequence sequence::get_subsequence_chr_to(const chromosomePosition & to, unsigned short len) {
}

/*****************************************************************//**
* @brief get subsequence to sequence position
*
* This method can be used to extract a subsequence of a given length
* ending (i.e. 3' end) at a given position from the sequence.
* If the length is too high so that the queried subsequence would
* reach over the (5') end of the sequence, a shorter subsequence
* ending at the disired end position and starting at the (5') end of
* the sequence will be returned.
* If the given position is not part of the sequence, an empty sequence
* will be returned.
*
* @param to the end position in the sequence of the subsequence
* @param len the (maximal) length of the subsequence
* @return the subsequence ending at the given position and starting
*         after the given length (3' to 5') or at the start of the
*         sequence
*********************************************************************/

sequence sequence::get_subsequence_to(const sequencePosition & to, unsigned short len) {
}

/*****************************************************************//**
* @brief get subsequence between chromosome positions
*
* This method can be used to extract a subsequence starting (i.e. 5'
* end) end ending (i.e. 3' end) at given chromosome positions from the
* sequence.
* If at least one of the given chromosome positions is not part of the
* sequence, an empty sequence will be returned.
*
* @param from the start position on the chromosome of the subsequence
* @param to the end position on the chromosome of the subsequence
* @return the subsequence starting and ending at the given positions
*********************************************************************/

sequence sequence::get_subsequence_chr_from_to(const chromosomePosition & from, const chromosomePosition & to) {
}

/*****************************************************************//**
* @brief get subsequence between sequence positions
*
* This method can be used to extract a subsequence starting (i.e. 5'
* end) end ending (i.e. 3' end) at given positions from the sequence.
* If at least one of the given positions is not part of the sequence,
* an empty sequence will be returned.
*
* @param from the start position in the sequence of the subsequence
* @param to the end position in the sequence of the subsequence
* @return the subsequence starting and ending at the given positions
*********************************************************************/

sequence sequence::get_subsequence_from_to(const sequencePosition & from, const sequencePosition & to) {
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
