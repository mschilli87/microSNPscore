#ifndef MICROSNPSCORE_NUCLEOTIDE_H
#define MICROSNPSCORE_NUCLEOTIDE_H


#include <iostream>
//for std::ostream (operator<<)
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
* @brief match state identifier type
*
* This represents the identifiers of the match states (Indel,
* Mismatch, Masked, Wobble, Match) of a pair of nucleotides.
*********************************************************************/
enum matchIdentifier {
  IndelOpen,
  IndelExtend,
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
* @brief Indel type
*
* This represents the type (Open or Extend) of an Indel.
*********************************************************************/
enum IndelType {
  Open,
  Extend

};
/*****************************************************************//**
* @brief match position type
*
* This represents the positon (Seed or ThreePrime) of a match.
*********************************************************************/
enum matchPosition {
  Seed,
  ThreePrime

};
/*****************************************************************//**
* @brief match state type
*
* This represents the match states (Indel, Mismatch, Masked, Wobble,
* Match) of a pair of nucleotides.
*********************************************************************/
class matchType {
  public:
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
    
    matchType(matchIdentifier match_type, matchPosition position);

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
    * The scoring scheme is taken from miRanda because mirSVR was trained
    * with miRanda alignments (IndelOpen: -9 IndelExtend: -4 Mismatch: -3,
    * Masked: -1, Wobble: -1, Match: +5).
    *
    * @return the score of the match type.
    *********************************************************************/
    inline const matchScore get_score() const;


  private:
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
    
    static matchScore calculate_score(const matchIdentifier the_identifier, matchPosition position);

    /*****************************************************************//**
    * @brief match state identifier
    *
    * This is the identifiers of the match state (Indel, Mismatch, Masked,
    * Wobble, Match).
    * It should be const but because match types shall be used in a vector
    * and std::vector tries to assign its elements to an internal array it
    * needs a working assignment operator which has to change the object's
    * members and therefore they cannot be declared const.
    * Nevertheless this attribute is not intended to be changed in any
    * other context than assigning an initialized object to an unitialized
    * one produced by the standard constructor which is in fact not
    * designed to be used directly but only provided to allow array
    * allocation which is needed to create containers, too.
    *********************************************************************/
    matchIdentifier identifier;

    /*****************************************************************//**
    * @brief match state score
    *
    * This is the score of the match state according to the miRanda
    * scoring scheme (IndelOpen: -9 IndelExtend: -4 Mismatch: -3,
    * Masked: -1, Wobble: -1, Match: +5).
    * It should be const but because match types shall be used in a vector
    * and std::vector tries to assign its elements to an internal array it
    * needs a working assignment operator which has to change the object's
    * members and therefore they cannot be declared const.
    * Nevertheless this attribute is not intended to be changed in any
    * other context than assigning an initialized object to an unitialized
    * one produced by the standard constructor which is in fact not
    * designed to be used directly but only provided to allow array
    * allocation which is needed to create containers, too.
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
    * The scoring scheme is taken from miRanda because mirSVR was trained
    * with miRanda alignments (IndelOpen: -9 IndelExtend: -4 Mismatch: -3,
    * Masked: -1, Wobble: -1, Match: +5).
    *
    * @return the score of the match type.
    *********************************************************************/
    inline const matchScore matchType::get_score() const {
      return score;
    }

/*****************************************************************//**
* @brief nucleotide class
*
* This represents a nucleotide.
*********************************************************************/
class nucleotide {
  public:
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
    nucleotide(nucleoBase the_base = Mask, sequencePosition the_sequence_position = 0, chromosomePosition the_chromosome_position = 0);

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
    matchType get_match(const nucleotide & matching_nucleotide, matchPosition position = ThreePrime, IndelType indel_type = Open) const;


  private:
    /*****************************************************************//**
    * @brief nucleo base
    *
    * This is the nucleo base of the nucleotide.
    * It should be const but because nucleotides shall be used in a vector
    * and std::vector tries to assign its elements to an internal array it
    * needs a working assignment operator which has to change the object's
    * members and therefore they cannot be declared const.
    * Nevertheless this attribute is not intended to be changed in any
    * other context than assigning an initialized object to an unitialized
    * one produced by the standard constructor which is in fact not
    * designed to be used directly but only provided to allow array
    * allocation which is needed to create containers, too.
    *
    * @see nucleotide()
    *********************************************************************/
    
    nucleoBase base;

    /*****************************************************************//**
    * @brief position in sequence
    *
    * This is the nucleotide's position in the sequence, the 5' end
    * beeing position 1.
    * Gaps are given the position of their predecessors in the alignment.
    * It should be const but because nucleotides shall be used in a vector
    * and std::vector tries to assign its elements to an internal array it
    * needs a working assignment operator which has to change the object's
    * members and therefore they cannot be declared const.
    * Nevertheless this attribute is not intended to be changed in any
    * other context than assigning an initialized object to an unitialized
    * one produced by the standard constructor which
    * is in fact not designed to be used directly but only provided to
    * allow array allocation which is needed to create containers, too.
    *
    * @see nucleotide()
    *********************************************************************/
    sequencePosition sequence_position;

    /*****************************************************************//**
    * @brief position on chromosome
    *
    * This is the nucleotide's position on its chromosome, the 5' end of
    * the + strand (i.e. the 3' end of the - strand) beeing position 1.
    * Gaps are given the position of their predecessors in the alignment.
    * It should be const but because nucleotides shall be used in a vector
    * and std::vector tries to assign its elements to an internal array it
    * needs a working assignment operator which has to change the object's
    * members and therefore they cannot be declared const.
    * Nevertheless this attribute is not intended to be changed in any
    * other context than assigning an initialized object to an unitialized
    * one produced by the standard constructor which
    * is in fact not designed to be used directly but only provided to
    * allow array allocation which is needed to create containers, too.
    *
    * @see nucleotide()
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
* @brief output stream nucleotide insertion operator
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
std::ostream & operator<<(std::ostream & the_stream, const nucleotide & the_nucleotide);

/*****************************************************************//**
* @brief output stream nucleo base insertion operator
*
* This operator is used to insert a nucleo base to an output stream
* (e.g. to print it on screen).
* The nucle base will be represented by its full name.
*
* @param the_stream output stream the nucleo base should be inserted
* in
* @param the_base nucleoBase to be inserted in the output stream
*
* @return output stream with the inserted nucleo base
*********************************************************************/
std::ostream & operator<<(std::ostream & the_stream, const nucleoBase & the_base);

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
std::ostream & operator<<(std::ostream & the_stream, const matchType & the_match);

} // namespace microSNPscore
#endif
