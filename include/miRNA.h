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

};

} // namespace microSNPscore
#endif
