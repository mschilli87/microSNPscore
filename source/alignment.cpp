
#include "alignment.h"
#include "mRNA.h"
#include "miRNA.h"

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
    * @brief constructor
    *
    * This is used to create an instance of the class alignmentMatrixCell.
    * The default values are not intended to be used directly.
    * They are only provided to allow array allocation but you will need
    * to assign a valid object created by giving those parameters a value
    * to actually use it. The reason for providing those default values is
    * to allow using (two-dimensional) arrays (beeing the alignment
    * matrices) containing objects of this class.
    *
    * @param the_mRNA_nucleotide (pseudo-optional) pointer to the mRNA
    *     nucleotide corresponding to the alignment matrix row holding the
    *     cell - Defaults to NULL
    * @param the_miRNA_nucleotide (pseudo-optional) pointer to the miRNA
    *     nucleotide corresponding to the alignment matrix column holding
    *     the cell - Defaults to NULL
    * @param the_score (pseudo-optional) score of an optimal aligment up
    *     to this cell - Defaults to 0
    * @param the_predecessors (pseudo-optional) vector containing pointers
    *     to the cells that are part of an optimal alignment up to the
    *     cell - Defaults to empty
    *
    * @return an alignment matrix cell with the given attributes
    *********************************************************************/
    
    alignmentMatrixCell::alignmentMatrixCell(const nucleotide * the_mRNA_nucleotide, const nucleotide * the_miRNA_nucleotide, alignmentScore the_score, const std::vector<const alignmentMatrixCell *> & the_predecessors)
    :mRNA_nucleotide(the_mRNA_nucleotide),miRNA_nucleotide(the_miRNA_nucleotide),score(the_score),predecessors(the_predecessors) {
}

    /*****************************************************************//**
    * @brief constructor
    *
    * This is used to create an instance of the class optimalAlignmentList
    * by aligning a given mRNA to a given miRNA.
    *
    * @param the_mRNA mRNA to be aligned
    * @param the_miRNA miRNA to be aligned
    *
    * @return an optimal alignment list for the given mRNA:miRNA-duplex
    *********************************************************************/
    optimalAlignmentList::optimalAlignmentList(const mRNA & the_mRNA, const miRNA & the_miRNA)
    :alignments(std::vector<alignment>()) {
       /****************************************************************\ 
      | If both sequences contain at least one nucleotide:               |
      | Allocate memory for matrices, fill them, search the last column  |
      | for best-scoring alignment ends and backtrace them appending the |
      | backtraced alignments to the alignments attribute:               |
       \****************************************************************/
      sequenceLength mRNA_length(the_mRNA.get_length());
      sequenceLength miRNA_length(the_miRNA.get_length());
      if(mRNA_length != 0 && miRNA_length != 0)
      {
        alignmentMatrixCell mRNA_gap_matrix[mRNA_length][miRNA_length];
        alignmentMatrixCell miRNA_gap_matrix[mRNA_length][miRNA_length];
        alignmentMatrixCell score_matrix[mRNA_length][miRNA_length];
        alignmentScore max_score(fill_matrices(&mRNA_gap_matrix[0][0],&miRNA_gap_matrix[0][0],&score_matrix[0][0],the_mRNA,the_miRNA));
        for(unsigned short int row=0;row!=mRNA_length;++row)
        {
          if(score_matrix[row][miRNA_length-1].get_score()==max_score)
          {
            backtrace_alignments(&score_matrix[row][miRNA_length-1],alignments);
          }
        }
      } // mRNA_length != 0 && miRNA_length != 0
}

    /*****************************************************************//**
    * @brief fill alignment matrices
    *
    * This method is used to calculate the values of the cells of the
    * alignment matrices as well as the optimal alignment's score.
    * The given matrices are assumed to have proper dimensions, otherwise
    * the behavior is undefined.
    *
    * @param matrix_mRNA_gap pointer to the [0][0] element of the
    *     alignment matrix that should hold the values for the optimal
    *     alignments up to each cells coordinates where there is an open
    *     gap in the mRNA in the last alignment column
    * @param matrix_miRNA_gap pointer to the [0][0] element of the
    *     alignment matrix that should hold the values for the optimal
    *     alignments up to each cells coordinates where there is an open
    *     gap in the miRNA in the last alignment column
    * @param matrix_overall pointer to the [0][0] element of the alignment
    *     matrix that should hold the values for the optimal alignments up
    *     to each cells coordinates
    * @param the_mRNA const reference to the mRNA to align
    * @param the_miRNA const reference to the miRNA to align
    *
    * @return the optimal alignment score
    *********************************************************************/
    alignmentScore optimalAlignmentList::fill_matrices(alignmentMatrixCell * matrix_mRNA_gap, alignmentMatrixCell * matrix_miRNA_gap, alignmentMatrixCell * matrix_overall, const mRNA & the_mRNA, const miRNA & the_miRNA)
    {
}

    /*****************************************************************//**
    * @brief recursive alignment calculation
    *
    * This method is used to calculate the optimal alignments going
    * through a given alignment matrix cell ending with a given postfix
    * and a given overall score and appending them to a given alignment
    * vector. The score is NOT checked against the alignment column scores
    * to improve performance.
    * The parameters @p first_call, @p the_score and @p postfix should be
    * left to the standard values. They only need to be changed for the
    * (intern) recursive calls.
    *
    * @param cell pointer to the alingment matrix cell up to which the
    *     backtrace should be performed
    * @param alignment_vector reference to the vector the backtraced
    *     alignments should be added to
    * @param (optional/intern) first_call boolean telling whether this is
    *     the first call (i.e. top of the recursion) - Defaults to True
    * @param (optional/intern) the_score alignment score of the backtraced
    *     alignments - Defaults to 0
    * @param (optional/intern) postfix pointer to a vector containing the
    *     alignment columns following the backtraced part - Defaults to
    *     empty
    *********************************************************************/
    
    void optimalAlignmentList::backtrace_alignments(alignmentMatrixCell * cell, std::vector<alignment> & alignment_vector, bool first_call, alignmentScore the_score, std::vector<alignmentColumn> * postfix)
    {
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
