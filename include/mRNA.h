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

};

} // namespace microSNPscore
#endif
