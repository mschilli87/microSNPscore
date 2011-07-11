#ifndef MICROSNPSCAN_SEQUENCE_H
#define MICROSNPSCAN_SEQUENCE_H


#include <vector>
using namespace std;

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

};

} // namespace microSNPscan
#endif
