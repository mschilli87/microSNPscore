#ifndef MICROSNPSCAN_SEQUENCE_H
#define MICROSNPSCAN_SEQUENCE_H


#include <vector>
using namespace std;
#include "alignment.h"

namespace microSNPscan {

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


  private:
    /*****************************************************************//**
    * @brief nucleotide sequence
    *
    * A vector containing the sequence's nucleotides from 5' to 3'
    *********************************************************************/
    vector<nucleotide> nucleotides;

    /*****************************************************************//**
    * @brief strand on chromosome
    *
    * This is the chromosome strand (Plus or Minus) the sequence is on.
    *********************************************************************/
    
    strand strand;


  public:
    /*****************************************************************//**
    * @brief get method for strand attribute
    *
    * This method is used to access the strand of the chromosome the
    * sequence is on.
    *
    * @return the strand of the sequence (Plus or Minus)
    *********************************************************************/
    
    inline const strand get_strand() const;

};
/*****************************************************************//**
* @brief get method for strand attribute
*
* This method is used to access the strand of the chromosome the
* sequence is on.
*
* @return the strand of the sequence (Plus or Minus)
*********************************************************************/

inline const strand sequence::get_strand() const {
  return strand;
}

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


  private:
    /*****************************************************************//**
    * @brief position on chromosome
    *
    * This is the nucleotide's position on its chromosome, the 5' end of
    * the + strand (i.e. the 3' end of the - strand) beeing position 1.
    * Gaps are given the position of their predecessors in the alignment.
    *********************************************************************/
    
    unsigned int chromosome_position;


  public:
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
    
    inline const unsigned int get_chromosome_position() const;


  private:
    /*****************************************************************//**
    * @brief nucleo base
    *
    * This is the nucleo base of the nucleotide.
    *********************************************************************/
    
    nucleoBase base;


  public:
    /*****************************************************************//**
    * @brief get method for nucleo base attribute
    *
    * This method is used to access the nucleo base of the nucleotide
    *
    * @return the nucleo base of the nucleotide
    *********************************************************************/
    inline const nucleoBase get_base() const;


  private:
    alignmentColumn ;

    /*****************************************************************//**
    * @brief position in sequence
    *
    * This is the nucleotide's position in the sequence, the 5' end
    * beeing position 1.
    * Gaps are given the position of their predecessors in the alignment.
    *********************************************************************/
    
    unsigned int sequence_position;


  public:
    /*****************************************************************//**
    * @brief get method for sequence position attribute
    *
    * This method is used to access the position of the nucleotide in the
    * sequence, the 5' end beeing position 1. Gaps will return the
    * position of their predecessor in the alignment.
    *
    * @return the position of the nucleotide in the sequence
    *********************************************************************/
    
    inline const unsigned int get_sequence_position() const;

};
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

inline const unsigned int nucleotide::get_chromosome_position() const {
  return chromosome_position;
}

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

inline const unsigned int nucleotide::get_sequence_position() const {
  return sequence_position;
}

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
* @brief strand type
*
* This represents the strand of a sequence (Plus or Minus).
*********************************************************************/

${template}class strand${inherit} {
${members}};
${inlines}

} // namespace microSNPscan
#endif
