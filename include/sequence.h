#ifndef MICROSNPSCAN_SEQUENCE_H
#define MICROSNPSCAN_SEQUENCE_H


namespace microSNPscan {

/*****************************************************************//**
* @brief sequence class                                              *
*                                                                    *
* This shall become the representation for sequences (DNA and RNA)   *
* and is used as my first example to generate code with BOUML        *
*********************************************************************/

class sequence {
  public:
    sequence();

    ~sequence();

    sequence(const sequence & source);

    sequence & operator=(const sequence & source);

};

} // namespace microSNPscan
#endif
