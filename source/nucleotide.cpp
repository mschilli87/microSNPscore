
#include "nucleotide.h"

namespace microSNPscore {

matchType::matchType(matchIdentifier match_type) {
  // Bouml preserved body begin 00031F13
  // Bouml preserved body end 00031F13
}

/*****************************************************************//**
* @brief standard constructor
*
* A match type must have an identifier.
* Calling a standard constructor to create a match type is therefore
* not intended.
*
* @todo remove this constructor when calling constructors are finished
*********************************************************************/

matchType::matchType() {
  // Bouml preserved body begin 00031F93
  // Bouml preserved body end 00031F93
}

/*****************************************************************//**
* @brief constructor
*
* This is used to create an instance of the class nucleotide.
*
* @param the_base nucleoBase that represents the nucleo base of the
*        nucleotide
* @param the_sequence_position sequencePosition that represents the
*        position in sequence of the nucleotide, the 5' end beeing
*	 position 1 (gaps should be given the position of their
*        predecessor in the alignment)
* @param the_chromosome_position chromosomePosition that represents
*        the position on the chromosome of the nucleotide, the 5' end
*        of the + strand (i.e. the 3' end of the - strand) beeing
*	 position 1 (gaps should be given the position of their
*        predecessor in the alignment)
* @return a nucleotide containing the given nucleo base and located at
*         the given positions on chromosome and in sequence
*********************************************************************/

nucleotide::nucleotide(nucleoBase the_base, sequencePosition the_sequence_position, chromosomePosition the_chromosome_position)
:base(the_base),sequence_position(the_sequence_position),chromosome_position(the_chromosome_position) {
  // Bouml preserved body begin 00024B13
return;
  // Bouml preserved body end 00024B13
}

/*****************************************************************//**
* @brief match calculation
*
* This method returns the match type of the nucleotide and a given
* other one.
* This operation is commutative, meaning that n1.get_match(n2) is
* always the same as n2.get_match(n1).
*
* @param matching_nucleotide const nucleotide reference to the
*        nucleotide that is paired with this one
* @return matchType representing the match between this and the given
*         nucleotide.
*********************************************************************/

matchType nucleotide::get_match(const nucleotide & matching_nucleotide) {
  // Bouml preserved body begin 00032013
  // Bouml preserved body end 00032013
}

/*****************************************************************//**
* @brief standard constructor
*
* A nucleotide must have a base, and a position on chromosome and in
* sequence.
* Calling a standard constructor to create a nucleotide is therefore
* not intended.
*
* @todo remove this constructor when calling constructors are finished
*********************************************************************/

nucleotide::nucleotide() {
  // Bouml preserved body begin 0002EB13
  // Bouml preserved body end 0002EB13
}


} // namespace microSNPscore
