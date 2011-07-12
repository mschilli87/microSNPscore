#ifndef MICROSNPSCORE_MRNA_H
#define MICROSNPSCORE_MRNA_H


#include "sequence.h"

namespace microSNPscore {

/*****************************************************************//**
* @brief messenger RNA class
*
* This represents a messenger RNA sequence.
* It is a specialisation of the sequence class.
*********************************************************************/

class mRNA : public sequence {
  public:
    /*****************************************************************//**
    * @brief standard constructor
    *
    * This is used to create an instance of the class mRNA.
    *
    * @return an empty mRNA
    *********************************************************************/
    mRNA();

    /*****************************************************************//**
    * @brief destructor
    *
    * This is used to remove an instance of the mRNA class.
    *
    *********************************************************************/
    ~mRNA();

    /*****************************************************************//**
    * @brief const copy constructor
    *
    * This is used to copy an instance of the mRNA class.
    *
    * @param source mRNA object to copy
    * @return a copy of the source mRNA object
    *********************************************************************/
    mRNA(const mRNA & source);

    /*****************************************************************//**
    * @brief assignment operator
    *
    * This is used to copy an instance of  the mRNA class by
    * assigning it to another one
    *
    * @param source the mRNA object to be copied
    * @return a copy of the source mRNA object
    *********************************************************************/
    mRNA & operator=(const mRNA & source);

};

} // namespace microSNPscore
#endif
