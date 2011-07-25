
#include "alignment.h"

namespace microSNPscore {

/*****************************************************************//**
* @brief constructor - Do not use without parameter values!
*
* This is used to create an instance of the class alignmentColumn.
* The default values for @p the_mRNA_nucleotide and
* @p the_miRNA_nucleotide are not intended to be used directly.
* They are only provided to allow array allocation but you will need
* to assign a valid object created by giving those parameters a value
* to actually use it. This is done by containers like std::vector and
* the reason for providing those default values is to allow using
* containers containing objects of this class.
*
* @param the_mRNA_nucleotide (pseudo-optional) const nucleotide
*     reference to the nucleotide of the messenger RNA in that
*     alignment column - Defaults to uninitialized nucleotide
*     constructed without parameters
* @param the_miRNA_nucleotide (pseudo-optional) const nucleotide
*     reference to the nucleotide of the microRNA in that alignment
*     column - Defaults to uninitialized nucleotide constructed
*     without parameters
* @param position (optional) matchPosition indicating whether the
*     match occurs in the seed (Seed) or in the 3' region of the miRNA
*     (ThreePrime) - Defaults to ThreePrime
* @param indel_type (optional) indelType telling whether this match
*     would be the first continuing indel (i.e. Open) if it would be
*     an indel or if it is directly following an existing indel (i.e.
*     Extend) - Defaults to Open
*
* @return alignmentColumn representing an alignment column aligning
*     the given nucleotides
*
* @see nucleotide::nucleotide(nucleoBase,sequencePosition,
*    chromosomePosition)
*********************************************************************/

alignmentColumn::alignmentColumn(const nucleotide & the_mRNA_nucleotide, const nucleotide & the_miRNA_nucleotide, matchPosition position, IndelType indel_type)
:mRNA_nucleotide(the_mRNA_nucleotide),miRNA_nucleotide(the_miRNA_nucleotide),match(the_mRNA_nucleotide.get_match(the_miRNA_nucleotide,position,indel_type)) {
}

/*****************************************************************//**
* @brief constructor
*
* This is used to create an instance of the class alignment.
* The default values are not intended to be used directly.
* They are only provided to allow array allocation but you will need
* to assign a valid object created by giving those parameters a value
* to actually use it. This is done by containers like std::vector and
* the reason for providing those default values is to allow using
* containers containing objects of this class.
*
* @param the_columns (pseudo-optional) vector containing the
*     alignment's coumns - Defaults to empty
* @param the_score (pseudo-optinal) the overall score of the alignment
*     (It is NOT checked against the given coloumn score to improve
*     performance) - Defaults to 0
* @return alignment containing the given columns with the given score
*********************************************************************/
alignment::alignment(const std::vector<alignmentColumn> & the_columns, alignmentScore the_score):
columns(the_columns),score(the_score),seed_type(sixMer) {
   /******************************************************************\ 
  | Iterating over the first eight alignment columns checking for the  |
  | first position to be Adenine in mRNA, positions 2 to 7 and 8 to be |
  | matches and update the seed type if necessary:                     |
   \******************************************************************/
  const_iterator column_it(the_columns.begin());
  if(column_it != the_columns.end())
  {
    bool AOne(column_it->get_mRNA_nucleotide().get_base()==Adenine); // A1?
    while(++column_it != the_columns.end() &&
          column_it - the_columns.begin() != 7 &&
          column_it->get_match().get_identifier() == Match) {}
    if(column_it - the_columns.begin() == 7) // m2-m7?
    {
      if(column_it != the_columns.end() && column_it->get_match().get_identifier() == Match) // m8?
      {
        if(AOne)
        {
          seed_type = eightMer;
        }
        else
        {
          seed_type = sevenMerMEight;
        }
      } // m8
      else
      {
        if(AOne)
        {
          seed_type = sevenMerAOne;
        }
      }
    } // m2-m7
  } // the_columns.begin() == the_columns.end()
}

/*****************************************************************//**
* @brief output stream seed type insertion operator
*
* This operator is used to insert a seed type to an output stream
* (e.g. to print it on screen).
* The seed type will be represented by its name.
*
* @param the_stream output stream the seed type should be inserted in
* @param seed_type seedType to be inserted in the output stream
*
* @return output stream with the inserted seed type
*********************************************************************/
std::ostream & operator<<(std::ostream & the_stream, const seedType & seed_type)
{
   /******************************************************************\ 
  | Append name depending on the seed type.                            |
  | The cases are ordered from common to uncommon to reduce comparisms |
  | as much as possible. Of course the default case should never be    |
  | reached. Because return exits the function there is no break       |
  | statement needed after the cases.                                  | 
   \******************************************************************/
  switch(seed_type)
  {

    case eightMer: return the_stream << "8mer";
    case sevenMerAOne: return the_stream << "7mer-A1";
    case sevenMerMEight: return the_stream << "7mer-m8";
    case sixMer: return the_stream << "6mer";
    default:
      std::cerr << "microSNPscore::operator<<(seedType)\n";
      std::cerr << " ==> unkown seed type: ";
      std::cerr << seed_type << std::endl;
      std::cerr << "  --> assuming sixMer\n";
      return the_stream << "6mer";
  }
}

} // namespace microSNPscore
