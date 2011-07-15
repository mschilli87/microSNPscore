
#include "alignment.h"

namespace microSNPscore {

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

alignmentColumn::alignmentColumn(const nucleotide & the_mRNA_nucleotide, const nucleotide & the_miRNA_nucleotide, matchPosition position, IndelType indel_type)
:mRNA_nucleotide(the_mRNA_nucleotide),miRNA_nucleotide(the_miRNA_nucleotide),match(the_mRNA_nucleotide.get_match(the_miRNA_nucleotide,position,indel_type)) {
}

/*****************************************************************//**
* @brief standard constructor
*
* This is used to create an instance of the class alignment.
*
* @return an empty alignment
*********************************************************************/
alignment::alignment() {
}


} // namespace microSNPscore
