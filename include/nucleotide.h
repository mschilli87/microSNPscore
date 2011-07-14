#ifndef MICROSNPSCORE_NUCLEOTIDE_H
#define MICROSNPSCORE_NUCLEOTIDE_H


namespace microSNPscore {

/*****************************************************************//**
* @brief chromosome position type
*
* This represents a position on a chromosome, the 5' end of the
* + strand (i.e. the 3' end of the - strand) beeing position 1.
*********************************************************************/

typedef unsigned int chromosomePosition;
/*****************************************************************//**
* @brief chromosome position type
*
* This represents a position on a sequence, the 5' end of the
* sequence beeing position 1.
*********************************************************************/
typedef unsigned short sequencePosition;
/*****************************************************************//**
* @brief nucleo base type
*
* This represents the nucleo bases Adenine, Cytosine, Guanine & Uracil
* as well as a gap or a mask
*********************************************************************/

enum nucleoBase {
  Adenine,
  Cytosine,
  Guanine,
  Uracil,
  Mask,
  Gap

};
/*****************************************************************//**
* @brief nucleotide class
*
* This represents a nucleotide.
*********************************************************************/

class nucleotide {
  public:
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
    
    nucleotide(nucleoBase the_base, sequencePosition the_sequence_position, chromosomePosition the_chromosome_position);

    /*****************************************************************//**
    * @brief get method for nucleo base attribute
    *
    * This method is used to access the nucleo base of the nucleotide
    *
    * @return the nucleo base of the nucleotide
    *********************************************************************/
    inline const nucleoBase get_base() const;

    /*****************************************************************//**
    * @brief get method for sequence position attribute
    *
    * This method is used to access the position of the nucleotide in the
    * sequence, the 5' end beeing position 1. Gaps will return the
    * position of their predecessor in the alignment.
    *
    * @return the position of the nucleotide in the sequence
    *********************************************************************/
    
    inline const sequencePosition get_sequence_position() const;

    /*****************************************************************//**
    * @brief get method for chromosome position attribute
    *
    * This method is used to access the position of the nucleotide on its
    * chromosome, the 5' end of the + strand (i.e. the 3' end of the -
    * strand) beeing position 1. Gaps will return the position of their
    * predecessor in the alignment.
    *
    * @return the position of the nucleotide on its chromosome
    *********************************************************************/
    
    inline const chromosomePosition get_chromosome_position() const;

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
    
    nucleotide();


  private:
    /*****************************************************************//**
    * @brief nucleo base
    *
    * This is the nucleo base of the nucleotide.
    *********************************************************************/
    
    nucleoBase base;

    /*****************************************************************//**
    * @brief position in sequence
    *
    * This is the nucleotide's position in the sequence, the 5' end
    * beeing position 1.
    * Gaps are given the position of their predecessors in the alignment.
    *********************************************************************/
    sequencePosition sequence_position;

    /*****************************************************************//**
    * @brief position on chromosome
    *
    * This is the nucleotide's position on its chromosome, the 5' end of
    * the + strand (i.e. the 3' end of the - strand) beeing position 1.
    * Gaps are given the position of their predecessors in the alignment.
    *********************************************************************/
    chromosomePosition chromosome_position;

};
/*****************************************************************//**
* @brief get method for nucleo base attribute
*
* This method is used to access the nucleo base of the nucleotide
*
* @return the nucleo base of the nucleotide
*********************************************************************/
inline const nucleoBase nucleotide::get_base() const {
  return base;
}

/*****************************************************************//**
* @brief get method for sequence position attribute
*
* This method is used to access the position of the nucleotide in the
* sequence, the 5' end beeing position 1. Gaps will return the
* position of their predecessor in the alignment.
*
* @return the position of the nucleotide in the sequence
*********************************************************************/

inline const sequencePosition nucleotide::get_sequence_position() const {
  return sequence_position;
}

/*****************************************************************//**
* @brief get method for chromosome position attribute
*
* This method is used to access the position of the nucleotide on its
* chromosome, the 5' end of the + strand (i.e. the 3' end of the -
* strand) beeing position 1. Gaps will return the position of their
* predecessor in the alignment.
*
* @return the position of the nucleotide on its chromosome
*********************************************************************/

inline const chromosomePosition nucleotide::get_chromosome_position() const {
  return chromosome_position;
}

/*****************************************************************//**
* @brief match state identifier type
*
* This represents the identifiers of the match states (Indel,
* Mismatch, Masked, Wobble, Match) of a pair of nucleotides.
*********************************************************************/

enum matchIdentifier {
  Indel,
  Mismatch,
  Masked,
  Wobble,
  Match

};
/*****************************************************************//**
* @brief match state score type
*
* This represents the score of the match states (Indel, Mismatch,
* Masked, Wobble, Match) of a pair of nucleotides.
*********************************************************************/

typedef short matchScore;
/*****************************************************************//**
* @brief match state type
*
* This represents the match states (Indel, Mismatch, Masked, Wobble,
* Match) of a pair of nucleotides.
*********************************************************************/

class matchType {
  public:
    matchType(matchIdentifier match_type);

    /*****************************************************************//**
    * @brief standard constructor
    *
    * A match type must have an identifier.
    * Calling a standard constructor to create a match type is therefore
    * not intended.
    *
    * @todo remove this constructor when calling constructors are finished
    *********************************************************************/
    
    matchType();

    /*****************************************************************//**
    * @brief get method for identifier attribute
    *
    * This method is used to access the identifier (Indel, Mismatch,
    * Masked, Wobble, Match) of the match type.
    *
    * @return the identifier (Indel, Mismatch, Masked, Wobble, Match)
    * of the match type.
    *********************************************************************/
    inline const matchIdentifier get_identifier() const;

    /*****************************************************************//**
    * @brief get method for score attribute
    *
    * This method is used to access the score of the match type.
    *
    * @return the score of the match type.
    *********************************************************************/
    inline const matchScore get_score() const;


  private:
    /*****************************************************************//**
    * @brief match state identifier
    *
    * This is the identifiers of the match state (Indel, Mismatch, Masked,
    * Wobble, Match).
    *********************************************************************/
    
    matchIdentifier identifier;

    /*****************************************************************//**
    * @brief match state score
    *
    * This is the score of the match state.
    *********************************************************************/
    matchScore score;

};
/*****************************************************************//**
* @brief get method for identifier attribute
*
* This method is used to access the identifier (Indel, Mismatch,
* Masked, Wobble, Match) of the match type.
*
* @return the identifier (Indel, Mismatch, Masked, Wobble, Match)
* of the match type.
*********************************************************************/
inline const matchIdentifier matchType::get_identifier() const {
  return identifier;
}

/*****************************************************************//**
* @brief get method for score attribute
*
* This method is used to access the score of the match type.
*
* @return the score of the match type.
*********************************************************************/
inline const matchScore matchType::get_score() const {
  return score;
}


} // namespace microSNPscore
#endif
