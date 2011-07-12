#ifndef MICROSNPSCORE_SEQUENCE_H
#define MICROSNPSCORE_SEQUENCE_H


#include <string>
#include <vector>

namespace microSNPscore {

/*****************************************************************//**
* @brief chromosome type
*
* This represents a chromosome. It is defined as string to handle
* different notations (like "chr1" or "1" and "MIT" or "24") (but this
* is done without any consistency checking) and special 'chromosomes'
* (like "HSCHR12_3_CTG2_1" and "GL000195.1").
*********************************************************************/

typedef std::string chromosomeType;
/*****************************************************************//**
* @brief strand type
*
* This represents the strand of a sequence (Plus or Minus).
*********************************************************************/

enum strandType {
  Plus,
  Minus

};
/*****************************************************************//**
* @brief exon class
*
* This represents an exon as a pair of positions (start/stop).
*********************************************************************/

class exon {
  public:
    /*****************************************************************//**
    * @brief standard constructor
    *
    * This is used to create an instance of the class exon.
    *
    * @return an empty exon
    *********************************************************************/
    exon();

    /*****************************************************************//**
    * @brief destructor
    *
    * This is used to remove an instance of the exon class.
    *
    *********************************************************************/
    ~exon();

    /*****************************************************************//**
    * @brief const copy constructor
    *
    * This is used to copy an instance of the exon class.
    *
    * @param source exon object to copy
    * @return a copy of the source exon object
    *********************************************************************/
    exon(const exon & source);

    /*****************************************************************//**
    * @brief assignment operator
    *
    * This is used to copy an instance of  the exon class by
    * assigning it to another one
    *
    * @param source the exon object to be copied
    * @return a copy of the source exon object
    *********************************************************************/
    exon & operator=(const exon & source);

};
/*****************************************************************//**
* @brief nucleo base type
*
* This represents the nucleo bases Adenine, Cytosine, Guanine & Uracil
*********************************************************************/

enum nucleoBase {
  Adenine,
  Cytosine,
  Guanine,
  Uracil,
  Masked,
  Gap

};
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
* @brief nucleotide class
*
* This represents a nucleotide.
*********************************************************************/

class nucleotide {
  public:
    /*****************************************************************//**
    * @brief standard constructor
    *
    * This is used to create an instance of the class nucleotide.
    *
    * @return an empty exon
    *********************************************************************/
    nucleotide();

    /*****************************************************************//**
    * @brief destructor
    *
    * This is used to remove an instance of the nucleotide class.
    *
    *********************************************************************/
    ~nucleotide();

    /*****************************************************************//**
    * @brief const copy constructor
    *
    * This is used to copy an instance of the nucleotide class.
    *
    * @param source nucleotide object to copy
    * @return a copy of the source nucleotide object
    *********************************************************************/
    nucleotide(const nucleotide & source);

    /*****************************************************************//**
    * @brief assignment operator
    *
    * This is used to copy an instance of  the nucleotide class by
    * assigning it to another one
    *
    * @param source the nucleotide object to be copied
    * @return a copy of the source nucleotide object
    *********************************************************************/
    nucleotide & operator=(const nucleotide & source);

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
* @brief sequence class
*
* This shall become the representation for sequences (DNA and RNA)
* and is used as my first example to generate code with BOUML
*********************************************************************/

class sequence {
  public:
    /*****************************************************************//**
    * @brief standard constructor
    *
    * This is used to create an instance of the class sequence.
    *
    * @return an empty sequence
    *********************************************************************/
    
    sequence();

    /*****************************************************************//**
    * @brief destructor
    *
    * This is used to remove an instance of the sequence class.
    *
    *********************************************************************/
    
    ~sequence();

    /*****************************************************************//**
    * @brief const copy constructor
    *
    * This is used to copy an instance of the sequence class.
    *
    * @param source sequence object to copy
    * @return a copy of the source sequence object
    *********************************************************************/
    
    sequence(const sequence & source);

    /*****************************************************************//**
    * @brief assignment operator
    *
    * This is used to copy an instance of  the sequence class by
    * assigning it to another one
    *
    * @param source the sequece object to be copied
    * @return a copy of the source sequence object
    *********************************************************************/
    
    sequence & operator=(const sequence & source);

    /*****************************************************************//**
    * @brief get method for chromosome attribute
    *
    * This method is used to access the chromosome the sequence is on.
    *
    * @return the chromosome of the sequence
    *********************************************************************/
    
    inline const chromosomeType get_chromosome() const;

    /*****************************************************************//**
    * @brief get method for strand attribute
    *
    * This method is used to access the strand of the chromosome the
    * sequence is on.
    *
    * @return the strand of the sequence (Plus or Minus)
    *********************************************************************/
    inline const strandType get_strand() const;

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
    
    sequence get_subsequence_from(sequencePosition from, unsigned short len);

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
    
    sequence get_subsequence_to(sequencePosition to, unsigned short len);

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
    
    sequence get_subsequence_from_to(sequencePosition from, sequencePosition to);

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
    
    sequence get_subsequence_chr_from(chromosomePosition from, unsigned short len);

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
    
    sequence get_subsequence_chr_to(chromosomePosition to, unsigned short len);

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
    
    sequence get_subsequence_chr_from_to(chromosomePosition from, chromosomePosition to);


  private:
    /*****************************************************************//**
    * @brief chromosome
    *
    * This is the chromosome the sequence is on.
    *********************************************************************/
    
    chromosomeType chromosome;

    /*****************************************************************//**
    * @brief strand on chromosome
    *
    * This is the chromosome strand (Plus or Minus) the sequence is on.
    *********************************************************************/
    
    strandType strand;

    /*****************************************************************//**
    * @brief nucleotide sequence
    *
    * A vector containing the sequence's nucleotides from 5' to 3'
    *********************************************************************/
    std::vector<nucleotide> nucleotides;

};
/*****************************************************************//**
* @brief get method for chromosome attribute
*
* This method is used to access the chromosome the sequence is on.
*
* @return the chromosome of the sequence
*********************************************************************/

inline const chromosomeType sequence::get_chromosome() const {
  return chromosome;
}

/*****************************************************************//**
* @brief get method for strand attribute
*
* This method is used to access the strand of the chromosome the
* sequence is on.
*
* @return the strand of the sequence (Plus or Minus)
*********************************************************************/
inline const strandType sequence::get_strand() const {
  return strand;
}


} // namespace microSNPscore
#endif
