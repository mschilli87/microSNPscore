
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
* @param position matchPosition indicating whether the match occurs in
*        the seed (Seed) or in the 3' region of the miRNA (ThreePrime)
*********************************************************************/

matchType::matchType(matchIdentifier match_type, matchPosition position)
:identifier(match_type),score(calculate_score(match_type,position)) {
}

/*****************************************************************//**
* @brief score initialization function
*
* This is used set the score according to the match identifier.
* The scoring scheme is taken from miRanda because mirSVR was trained
* with miRanda alignments (IndelOpen: -9 IndelExtend: -4 Mismatch: -3,
* Masked: -1, Wobble: -1, Match: +5, Seed: *4).
* Because score is const is has to be initialized before the
* constructor runs and because the object is not completed at this
* time, this function is declared static to rule out any side effects.
*
* @param the_identifier matchIdentifier representing the match type
*        (Indel, Mismatch, Masked, Wobble, Match).
* @param position matchPosition indicating whether the match occurs in
*        the seed (Seed) or in the 3' region of the miRNA (ThreePrime)
* @return matchScore representing the score of a match of the given
*         type
*********************************************************************/

matchScore matchType::calculate_score(matchIdentifier the_identifier, matchPosition position)
{
   /*****************************************************************\ 
  | Returning the score corresponding to the match type as used in    |
  | miRanda. The cases are ordered from common to uncommon to reduce  |
  | comparisms as much as possible. Of course the default case should |
  | never be reached. Because return exits the function there is no   |
  | break statement needed after the cases.                           |
   \*****************************************************************/
  const matchScore weight(position==Seed ? 4 : 1);
  switch(the_identifier)
  {
    case Match: return (5*weight);
    case Mismatch: return (-3*weight);
    case IndelExtend: return (-4*weight);
    case IndelOpen: return (-9*weight);
    case Wobble: return (-1*weight);
    case Masked: return (-1*weight);
    default:
     std::cerr << "microSNPscore::matchType::calculateScore\n";
     std::cerr << " ==> Undefined match type identifier: ";
     std::cerr << the_identifier << std:: endl;
     std::cerr << "  --> assuming Masked --> returning -1\n";
     return (-1*weight);
  }
}

/*****************************************************************//**
* @brief constructor - Do not call without parameter values!
*
* This is used to create an instance of the class nucleotide.
* The default values are not intended to be used directly.
* They are only provided to allow array allocation but you will need
* to assign a valid object created by giving those parameters a value
* to actually use it. This is done by containers like std::vector and
* the reason for providing those default values is to allow using
* containers containing objects of this class.
*
* @param the_base (pseudo-optional) nucleoBase that represents the
*     nucleo base of the nucleotide - Defaults to Mask
* @param the_sequence_position (pseudo-optional) sequencePosition that
*     represents the position in sequence of the nucleotide, the 5'
*     end beeing position 1 (gaps should be given the position of
*     their predecessor in the alignment) - Defaults to 0
* @param the_chromosome_position (pseudo-optional) chromosomePosition
*     that represents the position on the chromosome of th
*     nucleotide, the 5' end of the + strand (i.e. the 3' end of
*     the - strand) beeing position 1 (gaps should be given the
*     position of their predecessor in the alignment) - Defaults to 0
*
* @return a nucleotide containing the given nucleo base and located at
*     the given positions on chromosome and in sequence
*********************************************************************/
nucleotide::nucleotide(nucleoBase the_base, sequencePosition the_sequence_position, chromosomePosition the_chromosome_position)
:base(the_base),sequence_position(the_sequence_position),chromosome_position(the_chromosome_position) {
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
    *     nucleotide that is paired with this one
    * @param position (optional) matchPosition indicating whether the
    *     match occurs in the seed (Seed) or in the 3' region of the miRNA
    *     (ThreePrime) - Defaults to ThreePrime
    * @param indel_type (optional) indelType telling whether this match
    *     would be the first continuing indel (i.e. Open) if it would be
    *     an indel or if it is directly following an existing indel (i.e.
    *     Extend) - Defaults to Open
    * @return matchType representing the match between this and the given
    *     nucleotide.
    *********************************************************************/
    matchType nucleotide::get_match(const nucleotide & matching_nucleotide, matchPosition position, IndelType indel_type) const {
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
            case Uracil: return matchType(Mismatch,position);
            case Adenine: return matchType(Match,position);
            case Guanine: return matchType(Wobble,position);
            case Cytosine: return matchType(Mismatch,position);
            case Gap: return matchType((indel_type==Open ? IndelOpen : IndelExtend),position);
            case Mask: return matchType(Masked,position);
            default:
              std::cerr << "microSNPscore::nucleotide::get_match\n";
              std::cerr << " ==> Undefined nucleo base: ";
              std::cerr << match_base << std:: endl;
              std::cerr << "  --> assuming Masked\n";
              return matchType(Masked,position);
          }
        case Adenine:
          switch(match_base)
          {
            case Uracil: return matchType(Match,position);
            case Adenine: return matchType(Mismatch,position);
            case Guanine: return matchType(Mismatch,position);
            case Cytosine: return matchType(Mismatch,position);
            case Gap: return matchType((indel_type==Open ? IndelOpen : IndelExtend),position);
            case Mask: return matchType(Masked,position);
            default:
              std::cerr << "microSNPscore::nucleotide::get_match\n";
              std::cerr << " ==> Undefined nucleo base: ";
              std::cerr << match_base << std:: endl;
              std::cerr << "  --> assuming Masked\n";
              return matchType(Masked,position);
          }
        case Guanine:
          switch(match_base)
          {
            case Uracil: return matchType(Wobble,position);
            case Adenine: return matchType(Mismatch,position);
            case Guanine: return matchType(Mismatch,position);
            case Cytosine: return matchType(Match,position);
            case Gap: return matchType((indel_type==Open ? IndelOpen : IndelExtend),position);
            case Mask: return matchType(Masked,position);
            default:
              std::cerr << "microSNPscore::nucleotide::get_match\n";
              std::cerr << " ==> Undefined nucleo base: ";
              std::cerr << match_base << std:: endl;
              std::cerr << "  --> assuming Masked\n";
              return matchType(Masked,position);
          }
        case Cytosine:
          switch(match_base)
          {
            case Uracil: return matchType(Mismatch,position);
            case Adenine: return matchType(Mismatch,position);
            case Guanine: return matchType(Match,position);
            case Cytosine: return matchType(Mismatch,position);
            case Gap: return matchType((indel_type==Open ? IndelOpen : IndelExtend),position);
            case Mask: return matchType(Masked,position);
            default:
              std::cerr << "microSNPscore::nucleotide::get_match\n";
              std::cerr << " ==> Undefined nucleo base: ";
              std::cerr << match_base << std:: endl;
              std::cerr << "  --> assuming Masked\n";
              return matchType(Masked,position);
          }
        case Gap: return matchType((indel_type==Open ? IndelOpen : IndelExtend),position);
        case Mask: return matchType(Masked,position);
        default:
          std::cerr << "microSNPscore::nucleotide::get_match\n";
          std::cerr << " ==> Undefined nucleo base: ";
          std::cerr << match_base << std:: endl;
          std::cerr << "  --> assuming Masked\n";
          return matchType(Masked,position);
      }
}

/*****************************************************************//**
* @brief output stream insertion operator
*
* This operator is used to insert a nucleotide to an output stream
* (e.g. to print it on screen).
* The nucleotide will be represented by its letter code (Adenine: A,
* Cytosine: C, Guanine: G, Uracil: U, Masked: X).
*
* @param the_stream output stream the nucleotide should be inserted in
* @param the_nucleotide nucleotide to be inserted in the output stream
*
* @return output stream with the inserted nucleotide
*********************************************************************/
std::ostream & operator<<(std::ostream & the_stream, const nucleotide & the_nucleotide)
{
   /******************************************************************\ 
  | Call base get method only once and append letter code depending on |
  | the nucleo base.                                                   |
  | The cases are ordered from common to uncommon to reduce comparisms |
  | as much as possible. Of course the default case should never be    |
  | reached. Because return exits the function there is no break       |
  | statement needed after the cases.                                  | 
   \******************************************************************/
  nucleoBase the_base(the_nucleotide.get_base());
  switch(the_base)
  {
    case Adenine: return the_stream << 'A';
    case Uracil: return the_stream << 'U';
    case Cytosine: return the_stream << 'C';
    case Guanine: return the_stream << 'G';
    case Gap: return the_stream << '-';
    case Mask: return the_stream << 'X';
    default:
      std::cerr << "microSNPscore::operator<<(nucleotide)\n";
      std::cerr << " ==> unkown base: ";
      std::cerr << the_base << std::endl;
      std::cerr << "  --> assuming Mask\n";
      return the_stream << 'X';
  }
}

/*****************************************************************//**
* @brief output stream nucleo base insertion operator
*
* This operator is used to insert a nucleo base to an output stream
* (e.g. to print it on screen).
* The nucle base will be represented by its full name.
*
* @param the_stream output stream the nucle base should be inserted in
* @param the_base nucleoBase to be inserted in the output stream
*
* @return output stream with the inserted nucleo base
*********************************************************************/
std::ostream & operator<<(std::ostream & the_stream, const nucleoBase & the_base)
{
   /******************************************************************\ 
  | Append name depending the nucleo base.                             |
  | The cases are ordered from common to uncommon to reduce comparisms |
  | as much as possible. Of course the default case should never be    |
  | reached. Because return exits the function there is no break       |
  | statement needed after the cases.                                  | 
   \******************************************************************/
  switch(the_base)
  {
    case Adenine: return the_stream << "Adenine";
    case Uracil: return the_stream << "Uracil";
    case Cytosine: return the_stream << "Cytosine";
    case Guanine: return the_stream << "Guanine";
    case Gap: return the_stream << "Gap";
    case Mask: return the_stream << "Mask";
    default:
      std::cerr << "microSNPscore::operator<<(nucleoBase)\n";
      std::cerr << " ==> unkown base: ";
      std::cerr << the_base << std::endl;
      std::cerr << "  --> assuming Mask\n";
      return the_stream << "Mask";
  }
}

/*****************************************************************//**
* @brief output stream match insertion operator
*
* This operator is used to insert a match to an output stream (e.g. to
* print it on screen).
* The match will be represented by the full name of its identifier.
*
* @param the_stream output stream the match should be inserted in
* @param the_match matchType to be inserted in the output stream
*
* @return output stream with the inserted match
*********************************************************************/
std::ostream & operator<<(std::ostream & the_stream, const matchType & the_match)
{
   /******************************************************************\ 
  | Call identifier get method only once and append name depending on  |
  | the identifier.                                                    |
  | The cases are ordered from common to uncommon to reduce comparisms |
  | as much as possible. Of course the default case should never be    |
  | reached. Because return exits the function there is no break       |
  | statement needed after the cases.                                  | 
   \******************************************************************/
  matchIdentifier the_identifier(the_match.get_identifier());
  switch(the_identifier)
  {
    case Match: return the_stream << "Match";
    case Mismatch: return the_stream << "Mismatch";
    case IndelExtend: return the_stream << "IndelExtend";
    case IndelOpen: return the_stream << "IndelOpen";
    case Wobble: return the_stream << "Wobble";
    case Masked: return the_stream << "Masked";
    default:
      std::cerr << "microSNPscore::operator<<(matchType)\n";
      std::cerr << " ==> unkown match type identifier: ";
      std::cerr << the_identifier << std::endl;
      std::cerr << "  --> assuming Masked\n";
      return the_stream << 'Masked';
  }
}

} // namespace microSNPscore
