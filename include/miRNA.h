#ifndef MICROSNPSCORE_MIRNA_H
#define MICROSNPSCORE_MIRNA_H


#include "sequence.h"

namespace microSNPscore {

/*****************************************************************//**
* @brief microRNA class
*
* This represents a microRNA sequence.
* It is a specialisation of the sequence class.
*********************************************************************/

class miRNA : public sequence {
  public:
    /*****************************************************************//**
    * @brief standard constructor
    *
    * This is used to create an instance of the class miRNA.
    *
    * @return an empty miRNA
    *********************************************************************/
    miRNA();

    /*****************************************************************//**
    * @brief destructor
    *
    * This is used to remove an instance of the miRNA class.
    *
    *********************************************************************/
    ~miRNA();

    /*****************************************************************//**
    * @brief const copy constructor
    *
    * This is used to copy an instance of the miRNA class.
    *
    * @param source miRNA object to copy
    * @return a copy of the source miRNA object
    *********************************************************************/
    miRNA(const miRNA & source);

    /*****************************************************************//**
    * @brief assignment operator
    *
    * This is used to copy an instance of  the miRNA class by
    * assigning it to another one
    *
    * @param source the miRNA object to be copied
    * @return a copy of the source miRNA object
    *********************************************************************/
    miRNA & operator=(const miRNA & source);

};

} // namespace microSNPscore
#endif
