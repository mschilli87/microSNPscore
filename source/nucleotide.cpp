
#include <iostream>
//for std::cerr and std::endl (error stating)
#include "nucleotide.h"

namespace microSNPscore {

/*****************************************************************//**
* @brief constructor
*
* This is used to create an instance of the class match type.
*
* @param match_type matchIdentifier representing the match type
* (Indel, Mismatch, Masked, Wobble, Match)
*********************************************************************/
matchType::matchType(matchIdentifier match_type)
:identifier(match_type),score(calculate_score(match_type)) {
}

/*****************************************************************//**
* @brief score initialization function
*
* This is used set the score according to the match identifier.
* The scoring scheme is taken from miRanda because mirSVR was trained
* with miRanda alignments (IndelOpen: -9 IndelExtend: -4 Mismatch: -3,
* Masked: -1, Wobble: -1, Match: +5).
* Because score is const is has to be initialized before the
* constructor runs and because the object is not completed at this
* time, this function is declared static to rule out any side effects.
*
* @param the_identifier matchIdentifier representing the match type
*        (Indel, Mismatch, Masked, Wobble, Match).
* @return matchScore representing the score of a match of the given
*         type
*********************************************************************/

matchScore matchType::calculate_score(matchIdentifier the_identifier)
{
   /*****************************************************************\ 
  | Returning the score corresponding to the match type as used in    |
  | miRanda. The cases are ordered from common to uncommon to reduce  |
  | comparisms as much as possible. Of course the default case should |
  | never be reached. Because return exits the function there is no   |
  | break statement needed after the cases.                           |
   \*****************************************************************/
  switch(the_identifier)
  {
    case Match: return 5;
    case Mismatch: return -3;
    case IndelExtend: return -4;
    case IndelOpen: return -9;
    case Wobble: return -1;
    case Masked: return -1;
    default:
     std::cerr << "microSNPscore::matchType::calculateScore\n";
     std::cerr << " ==> Undefined match type identifier: ";
     std::cerr << the_identifier << std:: endl;
     std::cerr << "  --> assuming Masked --> returning -1\n";
     return -1;
  }
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
return;
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
* @param indel_type (optional) indelType telling whether this match
*        would be the first continuing indel (i.e. Open) if it would be
*        an indel or if it is directly following an existing indel (i.e.
*        Extend) - Defaults to Open
* @return matchType representing the match between this and the given
*         nucleotide.
*********************************************************************/
matchType nucleotide::get_match(const nucleotide & matching_nucleotide, IndelType indel_type) const {
   /****************************************************************\ 
  | Calling both get methods only once and deciding by switch which  |
  | match type to return because this allows short runtime by due to |
  | large code fragments beeing skipped. The cases are ordered from  |
  | common to uncommon to reduce comparisms as much as possible. Of  |
  | course the default case should never be reached. Because return  |
  | exits the function there is no break statement needed after the  |
  | cases.                                                           |
   \****************************************************************/
  nucleoBase this_base(this->get_base());
  nucleoBase match_base(matching_nucleotide.get_base());
  switch(this_base)
  {
    case Uracil:
      switch(match_base)
      {
        case Uracil: return Mismatch;
        case Adenine: return Match;
        case Guanine: return Wobble;
        case Cytosine: return Mismatch;
        case Gap: return (indel_type==Open ? IndelOpen : IndelExtend);
        case Mask: return Masked;
        default:
          std::cerr << "microSNPscore::nucleotide::get_match\n";
          std::cerr << " ==> Undefined nucleo base: ";
          std::cerr << match_base << std:: endl;
          std::cerr << "  --> assuming Masked\n";
          return Masked;
      }
    case Adenine:
      switch(match_base)
      {
        case Uracil: return Match;
        case Adenine: return Mismatch;
        case Guanine: return Mismatch;
        case Cytosine: return Mismatch;
        case Gap: return (indel_type==Open ? IndelOpen : IndelExtend);
        case Mask: return Masked;
        default:
          std::cerr << "microSNPscore::nucleotide::get_match\n";
          std::cerr << " ==> Undefined nucleo base: ";
          std::cerr << match_base << std:: endl;
          std::cerr << "  --> assuming Masked\n";
          return Masked;
      }
    case Guanine:
      switch(match_base)
      {
        case Uracil: return Wobble;
        case Adenine: return Mismatch;
        case Guanine: return Mismatch;
        case Cytosine: return Match;
        case Gap: return (indel_type==Open ? IndelOpen : IndelExtend);
        case Mask: return Masked;
        default:
          std::cerr << "microSNPscore::nucleotide::get_match\n";
          std::cerr << " ==> Undefined nucleo base: ";
          std::cerr << match_base << std:: endl;
          std::cerr << "  --> assuming Masked\n";
          return Masked;
      }
    case Cytosine:
      switch(match_base)
      {
        case Uracil: return Mismatch;
        case Adenine: return Mismatch;
        case Guanine: return Match;
        case Cytosine: return Mismatch;
        case Gap: return (indel_type==Open ? IndelOpen : IndelExtend);
        case Mask: return Masked;
        default:
          std::cerr << "microSNPscore::nucleotide::get_match\n";
          std::cerr << " ==> Undefined nucleo base: ";
          std::cerr << match_base << std:: endl;
          std::cerr << "  --> assuming Masked\n";
          return Masked;
      }
    case Gap: return (indel_type==Open ? IndelOpen : IndelExtend);
    case Mask: return Masked;
    default:
      std::cerr << "microSNPscore::nucleotide::get_match\n";
      std::cerr << " ==> Undefined nucleo base: ";
      std::cerr << match_base << std:: endl;
      std::cerr << "  --> assuming Masked\n";
      return Masked;
  }
}


} // namespace microSNPscore
