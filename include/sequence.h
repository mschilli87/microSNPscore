#ifndef MICROSNPSCAN_SEQUENCE_H
#define MICROSNPSCAN_SEQUENCE_H


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
    * This is used to copy an instance of  the sequnce class by
    * assigning it to another one
    *
    * @param source the sequece object to be copied
    * @return a copy of the source sequence object
    *********************************************************************/
    
    sequence & operator=(const sequence & source);

};

} // namespace microSNPscan
#endif
