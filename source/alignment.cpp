
#include "alignment.h"

namespace microSNPscore {

/*****************************************************************//**
* @brief constructor
*
* This is used to create an instance of the class alignmentColumn.
*
* @param the_mRNA_nucleotide const nucleotide reference to the
*        nucleotide of the messenger RNA in that alignment column
* @param the_miRNA_nucleotide const nucleotide reference to the
*        nucleotide of the microRNA in that alignment column
* @return alignmentColumn representing an alignment column aligning
*         the given nucleotides
*********************************************************************/
alignmentColumn::alignmentColumn(const nucleotide & the_mRNA_nucleotide, const nucleotide & the_miRNA_nucleotide)
:mRNA_nucleotide(the_mRNA_nucleotide),miRNA_nucleotide(the_miRNA_nucleotide),match(the_mRNA_nucleotide.get_match(the_miRNA_nucleotide)) {
  // Bouml preserved body begin 00022C13
  // Bouml preserved body end 00022C13
}

/*****************************************************************//**
* @brief standard constructor
*
* This is used to create an instance of the class alignment.
*
* @return an empty alignment
*********************************************************************/
alignment::alignment() {
  // Bouml preserved body begin 00022E13
  // Bouml preserved body end 00022E13
}


} // namespace microSNPscore
