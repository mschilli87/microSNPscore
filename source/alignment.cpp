
#include <algorithm>
// for std::max (alignment score calculation) and std::back_copy (alignment traceback)
#include <iostream>
// for std::cerr and std::endl (error stating)
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
    * This is used to create an instance of the class
    * alignmentMatrixCellEntry.
    * The default values are not intended to be used directly.
    * They are only provided to allow array allocation but you will need
    * to assign a valid object created by giving those parameters a value
    * to actually use it. This is done by containers like std::vector and
    * the reason for providing those default values is to allow using
    * containers containing objects of this class.
    *
    * @param the_column (pseudo-optional) alignmentColumn corresponding to
    *     the entry - Defaults to uninitialized
    * @param the_predecessor (pseudo-optional) pointer to the alignment
    *     matrix cell preceeding the entry in an optimal alignment up to
    *     its cell - Defaults to NULL
    *
    * @return an alignment matrix cell entry with the given attributes
    *********************************************************************/
    
    alignmentMatrixCellEntry::alignmentMatrixCellEntry(const alignmentColumn & the_column, const alignmentMatrixCell * the_predecessor)
    :column(the_column),predecessor(the_predecessor) {
}

    /*****************************************************************//**
    * @brief constructor
    *
    * This is used to create an instance of the class alignmentMatrixCell.
    * Because it is not intended to bo used directly, it is declared
    * protected. It is only provided to have a standard constructor the
    * inherited classes can use. The reason for providing those
    * constructors is to allow using (two-dimensional) arrays (beeing the
    * alignment matrices) containing objects of those classes.
    *
    * @param the_entries (optional) a vector containing the entrys of the
    *      alignment cell - Defaults to empty
    * @param the_score (optional) score the alignment cell - Defaults to 0
    *
    * @return an alignment matrix cell with the given parameters.
    *********************************************************************/
    
    alignmentMatrixCell::alignmentMatrixCell(const std::vector<alignmentMatrixCellEntry> & the_entries, alignmentScore the_score)
    :entries(the_entries),score(the_score) {
}

    /*****************************************************************//**
    * @brief constructor
    *
    * This is used to construct a cell in a open-gap-matrix (except the
    * first row and column as well as the second row of the open-miRNA-
    * gap-matrix and the the second column of the open-mRNA-gap-matrix)
    * (i.e. where there could be a preceding gap in the corresponding
    * sequence).
    * It will have entries for all steps leading to a best scoring
    * alignment up to the cell ending with a gap in the corresponding
    * sequence.
    *
    * @param preceding_open_gap_cell const reference to open-gap-matrix
    *     cell that comes before the new inserted gap (one row up for
    *     open-miRNA-gap-matrix or one column left for open-mRNA-gap-
    *     matrix)
    * @param preceding_overall_cell const reference to overall-matrix cell
    *     that comes before the new inserted gap (one row up for open-
    *     miRNA-gap-matrix or one column left for open-mRNA-gap-matrix)
    * @param miRNA_nucleotide const reference to the miRNA nucleotide
    *     corresponing to the column of the open-gap-matrix cell that
    *     should be created (beeing a gap for open-miRNA-gap-matrix)
    * @param mRNA_nucleotide const reference to the mRNA nucleotide
    *     corresponing to the row of the open-gap-matrix cell that schould
    *     be created (beeing a gap for open-mRNA-gap-matrix)
    * @param match_position matchPosition indicating whether the score
    *    of the overall-matrix cell should be weighted (Seed) or not
    *    (ThreePrime)
    *
    * @return an open-gap-matrix cell with one entry for every step that
    *     is part of an optimal alignment up to the cell ending with a gap
    *     in the corresponding sequence
    *********************************************************************/
    
    openGapMatrixCell::openGapMatrixCell(const openGapMatrixCell & preceding_open_gap_cell, const overallMatrixCell & preceding_overall_cell, const nucleotide & miRNA_nucleotide, const nucleotide & mRNA_nucleotide, matchPosition match_position) {
       /*****************************************************************\ 
      | Calculate the best score of an alignment up to the cell opening a |
      | new gap and that of one extending an open gap:                    |
       \*****************************************************************/
      alignmentColumn gap_open_column(mRNA_nucleotide,miRNA_nucleotide,match_position,Open);
      alignmentColumn gap_extend_column(mRNA_nucleotide,miRNA_nucleotide,match_position,Extend);
      alignmentScore gap_open_score(gap_open_column.get_match().get_score() + preceding_overall_cell.get_score());
      alignmentScore gap_extend_score(gap_extend_column.get_match().get_score() + preceding_open_gap_cell.get_score());
       /***************************************************************\ 
      | Calculate the best score of an alingnment up to the cell ending |
      | with a gap and add entries for those steps that score that good |
      | before updating the score:                                      |
       \***************************************************************/
      alignmentScore best_score(std::max(gap_open_score,gap_extend_score));
      if(gap_open_score == best_score)
      {
        entries.push_back(alignmentMatrixCellEntry(gap_open_column,&preceding_overall_cell));
      }
      if(gap_extend_score == best_score)
      {
        entries.push_back(alignmentMatrixCellEntry(gap_extend_column,&preceding_open_gap_cell));
      }
      score = best_score;
}

    /*****************************************************************//**
    * @brief first row and column constructor
    *
    * This is used to construct a cell in the first row of the open-
    * mRNA-gap-matrix (except the first two columns) or the first column
    * of the open-mRNA-gap-matrix (except the first two rows) (i.e. the
    * first definded row/column where there must be a preceding gap in the
    * corresponding sequence).
    * It will have one entry extending the existing gap in the
    * corresponding sequence with the given gap.
    *
    * @param preceding_open_gap_cell const reference to open-gap-matrix
    *     cell that comes before the new inserted gap (one row up for
    *     open-miRNA-gap-matrix or one column left for open-mRNA-gap-
    *     matrix)
    * @param miRNA_nucleotide const reference to the miRNA nucleotide
    *     corresponing to the column of the open-gap-matrix cell that
    *     should be created (beeing a gap for open-miRNA-gap-matrix)
    * @param mRNA_nucleotide const reference to the mRNA nucleotide
    *     corresponing to the row of the open-gap-matrix cell that schould
    *     be created (beeing a gap for open-mRNA-gap-matrix)
    * @param match_position matchPosition indicating whether the score
    *    of the overall-matrix cell should be weighted (Seed) or not
    *    (ThreePrime)
    *
    * @return an open-gap-matrix cell with one entry extending an open gap
    *********************************************************************/
    openGapMatrixCell::openGapMatrixCell(const openGapMatrixCell & preceding_open_gap_cell, const nucleotide & miRNA_nucleotide, const nucleotide & mRNA_nucleotide, matchPosition match_position) {
       /*****************************************************************\ 
      | Calculate the best score of an alignment up to the cell extending |
      | an open gap and create a new entry for it before updating the     |
      | score:                                                            |
       \*****************************************************************/
      alignmentColumn gap_extend_column(mRNA_nucleotide,miRNA_nucleotide,match_position,Extend);
      alignmentScore gap_extend_score(gap_extend_column.get_match().get_score() + preceding_open_gap_cell.get_score());
      entries.push_back(alignmentMatrixCellEntry(gap_extend_column,&preceding_open_gap_cell));
      score = gap_extend_score;
}

    /*****************************************************************//**
    * @brief second row and column constructor
    *
    * This is used to construct a cell in the second row of the open-
    * miRNA-gap-matrix or the second column of the open-mRNA-gap-matrix
    * (i.e. the first definded row/column where there cannot be a
    * preceding gap in the corresponding sequence).
    * It will have one entry inserting the given gap in the corresponding
    * sequence.
    *
    * @param preceding_overall_cell const reference to overall-matrix cell
    *     that comes before the new inserted gap (one row up for open-
    *     miRNA-gap-matrix or one column left for open-mRNA-gap-matrix)
    * @param miRNA_nucleotide const reference to the miRNA nucleotide
    *     corresponing to the column of the open-gap-matrix cell that
    *     should be created (beeing a gap for open-miRNA-gap-matrix)
    * @param mRNA_nucleotide const reference to the mRNA nucleotide
    *     corresponing to the row of the open-gap-matrix cell that schould
    *     be created (beeing a gap for open-mRNA-gap-matrix)
    * @param match_position matchPosition indicating whether the score
    *    of the overall-matrix cell should be weighted (Seed) or not
    *    (ThreePrime)
    *
    * @return an open-gap-matrix cell with one entry opening a new gap
    *********************************************************************/
    
    openGapMatrixCell::openGapMatrixCell(const overallMatrixCell & preceding_overall_cell, const nucleotide & miRNA_nucleotide, const nucleotide & mRNA_nucleotide, matchPosition match_position) {
       /*****************************************************************\ 
      | Calculate the best score of an alignment up to the cell opening a |
      | new gap and create a new entry for it before updating the score:  |
       \*****************************************************************/
      alignmentColumn gap_open_column(mRNA_nucleotide,miRNA_nucleotide,match_position,Open);
      alignmentScore gap_open_score(gap_open_column.get_match().get_score() + preceding_overall_cell.get_score());
      entries.push_back(alignmentMatrixCellEntry(gap_open_column,&preceding_overall_cell));
      score = gap_open_score;
}

    /*****************************************************************//**
    * @brief standard constructor
    *
    * This is used to construct an object of the class openGapMatrixCell.
    * It is not intended to bo used directly. It is only provided to
    * allow using(two-dimensional) arrays (beeing the open gap score
    * alignment matrices) containing objects of this class.
    *
    * @return an open-gap-matrix cell with no entries and a score of 0.
    *********************************************************************/
    openGapMatrixCell::openGapMatrixCell() {
}

    /*****************************************************************//**
    * @brief constructor
    *
    * This is used to construct the cells of the overall score matrix
    * (except the first row and column).
    * It will have all entrys entries out of those of the given open gap
    * matrix cells and the one matching the given nucleotides that lead to
    * the best score.
    *
    * @param upper_left_overall_cell const reference to the overall-matrix
    *     cell one row up and one column left from the overall score
    *     matrix cell to create.
    * @param miRNA_open_gap_cell const reference to the open-miRNA-gap-
    *     matrix cell with the same coordinates as the overall-matrix cell
    *     to create
    * @param miRNA_open_gap_cell const reference to the open-miRNA-gap-
    *     matrix cell with the same coordinates as the overall-matrix cell
    *     to create
    * @param miRNA_nucleotide const reference to the nucleotide of the
    *    miRNA that corresponds to the row of the overall-matrix cell to
    *    create
    * @param mRNA_nucleotide const reference to the nucleotide of the mRNA
    *    that corresponds to the column of the overall-matrix cell to
    *    create
    * @param match_position matchPosition indicating whether the score
    *    of the overall-matrix cell should be weighted (Seed) or not
    *    (ThreePrime)
    *
    * @return an overall-matrix cell with the optimal score and all
    *     entries leading to it
    *********************************************************************/
    
    overallMatrixCell::overallMatrixCell(const overallMatrixCell & upper_left_overall_cell, const openGapMatrixCell & miRNA_open_gap_cell, const openGapMatrixCell & mRNA_open_gap_cell, const nucleotide & miRNA_nucleotide, const nucleotide & mRNA_nucleotide, matchPosition match_position) {
       /****************************************************************\ 
      | Calculate the best score of an alignment ending on a gap and the |
      | best score of an alignment ending with a (mis-)match:            |
       \****************************************************************/
      alignmentScore best_score(std::max(miRNA_open_gap_cell.get_score(),mRNA_open_gap_cell.get_score()));
      alignmentColumn match_column(mRNA_nucleotide,miRNA_nucleotide,match_position);
      alignmentScore match_score(match_column.get_match().get_score() + upper_left_overall_cell.get_score());
       /*****************************************************************\ 
      | Check whether the best alignment ending with a (mis-)match scores |
      | at least as good as the best alignment ending on a gap and if so  |
      | update the best score of any alignment an create a new entry for  |
      | the nucleotide match:                                             |
       \*****************************************************************/
      if(match_score >= best_score)
      {
        best_score = match_score;
        entries.push_back(alignmentMatrixCellEntry(match_column,&upper_left_overall_cell));
      }
       /*****************************************************************\ 
      | Check whether there are alignments ending on a gap that score as  |
      | good as the best alignment overall and if so update include those |
      | entries:                                                          |
       \*****************************************************************/
      if(miRNA_open_gap_cell.get_score() == best_score)
      {
        entries.insert(entries.begin(),miRNA_open_gap_cell.begin(),miRNA_open_gap_cell.end());
      }
      if(mRNA_open_gap_cell.get_score() == best_score)
      {
        entries.insert(entries.begin(),mRNA_open_gap_cell.begin(),mRNA_open_gap_cell.end());
      }
       /****************************************\ 
      | Set the score to the overall best score: |
       \****************************************/
      score = best_score;
}

    /*****************************************************************//**
    * @brief first row and column constructor
    *
    * This is used to construct the cells in the first row and column of
    * the overall score matrix (except the upper left corner).
    * It will have the same entries and score as the given open gap matrix
    * cell (which should have only one entry by the way).
    *
    * @param open_gap_cell const reference to the open-gap-matrix cell
    *     with the same coordinates as the overall-matrix cell to create
    *     (open-mRNA-gap-matrix for the first row and open-miRNA-gap-
    *     matrix for the first column)
    *
    * @return an overall-matrix cell with the same entries and score as
    *     the given open gap matrix cell
    *********************************************************************/
    
    overallMatrixCell::overallMatrixCell(const openGapMatrixCell & open_gap_cell)
    : alignmentMatrixCell(std::vector<alignmentMatrixCellEntry>(open_gap_cell.begin(),open_gap_cell.end()),open_gap_cell.get_score()) {
      
}

    /*****************************************************************//**
    * @brief upper left corner constructor
    *
    * This is used to construct the upper left corner of the overall score
    * matrix.
    * It will have one entry aligning the given nucleotides without taking
    * the match score into account (Since the miRNA's 3' end is needed for
    * RISC (RNA-induced silencing complex) binding.
    *
    * @param miRNA_five_prime const reference to the nucleotide at the 5'
    *     end of the miRNA
    * @param mRNA_three_prime const reference to the nucleotide at the 3'
    *     end of the mRNA
    *
    * @return an overall-matrix cell with one entry aligning the given
    *     nucleotides and a score of 0.
    *********************************************************************/
    overallMatrixCell::overallMatrixCell(const nucleotide & miRNA_five_prime, const nucleotide & mRNA_three_prime) {
      entries.push_back(alignmentColumn(mRNA_three_prime,miRNA_five_prime));
}

    /*****************************************************************//**
    * @brief standard constructor
    *
    * This is used to construct an object of the class overallMatrixCell.
    * It is not intended to bo used directly. It is only provided to
    * allow using a (two-dimensional) array (beeing the overall score
    * alignment matrices) containing objects of this class.
    *
    * @return an overall-matrix cell with no entries and a score of 0.
    *********************************************************************/
    overallMatrixCell::overallMatrixCell() {
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
        openGapMatrixCell mRNA_gap_matrix[mRNA_length][miRNA_length];
        openGapMatrixCell miRNA_gap_matrix[mRNA_length][miRNA_length];
        overallMatrixCell score_matrix[mRNA_length][miRNA_length];
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
    alignmentScore optimalAlignmentList::fill_matrices(openGapMatrixCell * matrix_mRNA_gap, openGapMatrixCell * matrix_miRNA_gap, overallMatrixCell * matrix_overall, const mRNA & the_mRNA, const miRNA & the_miRNA)
    {
       /****************************************************************\ 
      | Define function wide constants, get sequence lengths, initialize |
      | return value and check if the sequences actually contain         |
      | nucleotides (i.e. if a non-empty alignment is queried):          |
       \****************************************************************/
      const sequencePosition seed_start = 2;  // position 1 won't be aligned
      const sequencePosition seed_end = 8;    // weighting until (inclusive) that position
      sequenceLength mRNA_length(the_mRNA.get_length());
      sequenceLength miRNA_length(the_miRNA.get_length());
      alignmentScore best_overall_score(0);   // score of empty alignment
      if(mRNA_length != 0 && miRNA_length != 0)
      {
         /******************************************************************\ 
        | Iterate over the miRNA (5' to 3') to fill the matrices linewise    |
        | calculating the loop wide constant stating whether the current row |
        | should be weighted or not:                                         |
         \******************************************************************/
        unsigned short int column = 0;
        for(sequence::const_iterator miRNA_it(the_miRNA.begin());miRNA_it!=the_miRNA.end();++miRNA_it,++column)
        {
          const matchPosition match_pos (seed_start-1 <= column && column < seed_end ? Seed : ThreePrime);
           /******************************************************************\ 
          | Iterate over the mRNA (3' to 5') to fill the current row column by |
          | column calculating loop wide constants representing the indices of |
          | the relevant cells in the given arrays (which represent arrays in  |
          | two dimensions but because the matrix width - beeing the length of |
          | the miRNA - is unkwon at compilation time it is not possible to    |
          | use the two-dimensional subscript operator directly which is why   |
          | we need to calculate the offset from the [0][0] pointer - note: If |
          | the given matrices have smaller dimensions than the sequence's     |
          | length the bevavior of this implementation is undefined!) and the  |
          | gap-nucleotides that might be inserted in this iteration loop:     |
           \******************************************************************/
          unsigned short int row = 0;
          for(sequence::const_iterator mRNA_it(the_mRNA.end()-1);mRNA_it>=the_mRNA.begin();--mRNA_it,++row)
          {
            const unsigned short int index = row * miRNA_length + column;
            const unsigned short int index_left = index - 1;
            const unsigned short int index_up = index - miRNA_length;
            const unsigned short int index_upleft = index_up - 1;
            const nucleotide mRNA_gap(Gap,mRNA_it->get_sequence_position(),mRNA_it->get_chromosome_position());
            const nucleotide miRNA_gap(Gap,miRNA_it->get_sequence_position(),miRNA_it->get_chromosome_position());
             /******************************************************************\ 
            | While the central iteration of the linear programming algorithm is |
            | quite clear there are many different cases since the recursive     |
            | formulars rely on precalculated values (as it is the nature of     |
            | linear programming) that might not exist for the cells in the      |
            | first two rows and columns.                                        |
            | The first special case is the first matrix element at all: the     |
            | upper-left corner representing the alignment position given by the |
            | miRNA target prediction tool. Note that there are no starting gap  |
            | row and column since the reported mRNA position has to be aligned  |
            | with the miRNA's 5' end in any case to match the prediction:       |
             \******************************************************************/
            if(row==0) // row == 0
            {
              if(column==0) // row == 0 && column == 0
              {
                 /******************************************************************\ 
                | Since there is no preceeding alignment part, there is no recursion |
                | to apply (i.e. the base case of the overall score recursive        |
                | formular) and there cannot be any open gap (i.e. undefined for the |
                | open gap score recursive formulars). Note that the match score of  |
                | the first pair of nucleotides is not considered, since we already  |
                | know that the miRNA's 5' end will not bind to the mRNA in any case |
                | because it is needed for binding the RISC (RNA-induced silencing   |
                | complex):                                                          |
                 \******************************************************************/
                matrix_overall[index] = overallMatrixCell(*miRNA_it,*mRNA_it);
              } // row == 0 && col == 0
               /**************************************************************\ 
              | The second special case is the second column of the first row: |
               \**************************************************************/
              else if(column==1) // row == 0 && column == 1
              {
                 /*****************************************************************\ 
                | Since we are still in the first row, there cannot be any gap in   |
                | the miRNA (i.e. undefined in the open miRNA gap score formular)   |
                | and since we are in the second column there cannot be any         |
                | preceeding gap in the mRNA but there must be a preceeding match   |
                | (i.e. the base case of the open mRNA gap score formular).         |
                | Note that there is no need to consider the preceeding score since |
                | it has to be 0.                                                   |
                | Since we are still in the first row there is no other option than |
                | inserting a mRNA gap and thus the overall score equals the open   |
                | mRNA gap score in any case:                                       |
                 \*****************************************************************/
                matrix_mRNA_gap[index] = openGapMatrixCell(matrix_overall[index_left],*miRNA_it,mRNA_gap,match_pos);
                matrix_overall[index] = overallMatrixCell(matrix_mRNA_gap[index]);
              } // row == 0 && column == 1
               /*********************************************************\ 
              | The third special case is the remainder of the first row: |
               \*********************************************************/
              else // row == 0 && column > 1
              {
                 /*****************************************************************\ 
                | Since we are still in the first row, there cannot be any gap in   |
                | the miRNA (i.e. undefined in the open miRNA gap score formular)   |
                | and since we are behind the second column there must be at least  |
                | one gap in the mRNA (i.e. the gap open case of the open mRNA gap  |
                | score formular can be omitted).                                   |
                | Since we are still in the first row there is no other option than |
                | inserting a mRNA gap and thus the overall score equals the open   |
                | mRNA gap score in any case:                                       |
                 \*****************************************************************/
                matrix_mRNA_gap[index] = openGapMatrixCell(matrix_mRNA_gap[index_left],*miRNA_it,mRNA_gap,match_pos);
                matrix_overall[index] = overallMatrixCell(matrix_mRNA_gap[index]);
              } // row == 0 && column > 1
               /******************************************************************\ 
              | Since we now have completed the first row of the matrices we have  |
              | calculated the first overall alignment score which is yet the best |
              | score. Note that an optimal alignment can end anywhere in the last |
              | row because we want to align the whole miRNA while we do not care  |
              | about how much of the mRNA is covered by the alignment:            |
               \******************************************************************/
              if(column==miRNA_length-1) // first overall score
              {
                best_overall_score = matrix_overall[index].get_score();
              } // first overall score
            } // row==0
             /**************************************************************\ 
            | The fourth special case is the first column of the second row: |
             \**************************************************************/
            else if(row==1) // row == 1
            {
              if(column==0) // row == 1 && column == 0
              {
                 /*****************************************************************\ 
                | Since we are in the first column, there cannot be any gap in the  |
                | mRNA (i.e. undefined in the open mRNA gap score formular) and     |
                | since we are still in the second row there cannot be any          |
                | preceeding gap in the miRNA but there must be a preceeding match  |
                | (i.e. the base case of the open miRNA gap score formular).        |
                | Note that there is no need to consider the preceeding score since |
                | it has to be 0.                                                   |
                | Since we are in the first row there is no other option than       |
                | inserting a miRNA gap and thus the overall score equals the open  |
                | miRNA gap score in any case:                                      |
                 \*****************************************************************/
                matrix_miRNA_gap[index] = openGapMatrixCell(matrix_overall[index_up],miRNA_gap,*mRNA_it,match_pos);
                matrix_overall[index] = overallMatrixCell(matrix_miRNA_gap[index]);
              } // row == 1 && column == 0
               /**************************************************************\ 
              | The fifth special case is the second column of the second row: |
               \**************************************************************/
              else if(column==1) // row == 1 && column == 1
              {
                 /******************************************************************\ 
                | Since we are in the second column, there cannot be a be any        |
                | preceeding gap in the mRNA but there has must be a preceeding gap  |
                | in the miRNA (i.e. the base case of the open mRNA gap score        |
                | formular) and since we are still in the second row there cannot be |
                | any preceeding gap in the miRNA but there must be a preceeding gap |
                | in the mRNA (i.e. the base case of the open miRNA gap score        |
                | formular).                                                         |
                | Since we are neither in the first row nor in the first column we   |
                | can try to match the two nucleotides without adding a gap in one   |
                | of the sequences. After calculating all possible scores we add all |
                | those predecessors leading to the maximum to the overall score     |
                | matrix cell (i.e. the full recursive step of the overall score     |
                | formular):                                                         |
                 \******************************************************************/
                matrix_mRNA_gap[index] = openGapMatrixCell(matrix_overall[index_left],*miRNA_it,mRNA_gap,match_pos);
                matrix_miRNA_gap[index] = openGapMatrixCell(matrix_overall[index_up],miRNA_gap,*mRNA_it,match_pos);
                matrix_overall[index] = overallMatrixCell(matrix_overall[index_upleft],matrix_miRNA_gap[index],
                                                          matrix_mRNA_gap[index],*miRNA_it,*mRNA_it,match_pos);
              } // row == 1 && column == 1
               /**********************************************************\ 
              | The sixth special case is the remainder of the second row: |
               \**********************************************************/
              else // row == 1 && column > 1
              {
                 /******************************************************************\ 
                | Since we are behind the second column, there could be a preceeding |
                | gap in the mRNA thus we need to check whether it is better to open |
                | a new one or to extend the existing one (i.e. the full recursive   |
                | step of the open mRNA gap score formular) but since we are still   |
                | in the second row there cannot be any preceeding gap in the miRNA  |
                | while there must be a preceeding gap in the mRNA or a preceeding   |
                | match (i.e. the base case of the open miRNA gap score formular).   |
                | Since we are neither in the first row nor in the first column we   |
                | can try to match the two nucleotides without adding a gap in one   |
                | of the sequences. After calculating all possible scores we add all |
                | those predecessors leading to the maximum to the overall score     |
                | matrix cell (i.e. the full recursive step of the overall score     |
                | formular):                                                         |
                 \******************************************************************/
                matrix_mRNA_gap[index] = openGapMatrixCell(matrix_mRNA_gap[index_left],matrix_overall[index_left],
                                                           *miRNA_it,mRNA_gap,match_pos);
                matrix_miRNA_gap[index] = openGapMatrixCell(matrix_overall[index_up],miRNA_gap,*mRNA_it,match_pos);
                matrix_overall[index] = overallMatrixCell(matrix_overall[index_upleft],matrix_miRNA_gap[index],
                                                          matrix_mRNA_gap[index],*miRNA_it,*mRNA_it,match_pos);
              } // row ==1 && column > 1
               /******************************************************************\ 
              | Since we now have completed the second row of the matrices we have |
              | calculated the second overall alignment score which might be yet   |
              | the best score:                                                    |
               \******************************************************************/
              if(column==miRNA_length-1 && matrix_overall[index].get_score()>best_overall_score) // new best overall score
              {
                best_overall_score = matrix_overall[index].get_score();
              } // new best overall score
             /**************************************************************\ 
            | The seventh special case is the remainder of the first column: |
             \**************************************************************/
            } // row == 1
            else // row > 1
            {
              if(column==0) // row > 1 && column == 0
              {
                 /*****************************************************************\ 
                | Since we are in the first column, there cannot be any gap in the  |
                | mRNA (i.e. undefined in the open mRNA gap score formular) and     |
                | since we are behind the second row there must be at least one gap |
                | in the miRNA (i.e. the gap open case of the open miRNA gap score  |
                | formular can be omitted).                                         |
                | Since we are in the first row there is no other option than       |
                | inserting a miRNA gap and thus the overall score equals the open  |
                | miRNA gap score in any case:                                      |
                 \*****************************************************************/
                matrix_miRNA_gap[index] = openGapMatrixCell(matrix_miRNA_gap[index_up],*miRNA_it,mRNA_gap,match_pos);
                matrix_overall[index] = overallMatrixCell(matrix_miRNA_gap[index]);
              } // row > 1 && column == 0
               /**************************************************************\ 
              | The eight and last special case is the remainder of the second |
              | column:                                                        |
               \**************************************************************/
              else if(column==1) // row > 1 && column == 1
              {
                 /******************************************************************\ 
                | Since we are in the second column, there cannot be a be any        |
                | preceeding gap in the mRNA but there has must be a preceeding gap  |
                | in the miRNA or a preceeding match (i.e. the base case of the open |
                | mRNA gap score formular) but since we are behind the second row,   |
                | there could be a preceeding gap in the miRNA thus we need to check |
                | whether it is better to open  a new one or to extend the existing  |
                | one (i.e. the full recursive step of the open miRNA gap score      |
                | formular).                                                         |
                | Since we are neither in the first row nor in the first column we   |
                | can try to match the two nucleotides without adding a gap in one   |
                | of the sequences. After calculating all possible scores we add all |
                | those predecessors leading to the maximum to the overall score     |
                | matrix cell (i.e. the full recursive step of the overall score     |
                | formular):                                                         |
                 \******************************************************************/
                matrix_mRNA_gap[index] = openGapMatrixCell(matrix_overall[index_left],*miRNA_it,mRNA_gap,match_pos);
                matrix_miRNA_gap[index] = openGapMatrixCell(matrix_miRNA_gap[index_up],matrix_overall[index_up],
                                                           miRNA_gap,*mRNA_it,match_pos);
                matrix_overall[index] = overallMatrixCell(matrix_overall[index_upleft],matrix_miRNA_gap[index],
                                                          matrix_mRNA_gap[index],*miRNA_it,*mRNA_it,match_pos);
              } // row > 1 && column == 1
               /************************************************************\ 
              | The last case is the default one and is applied to the whole |
              | matrices except the first to columns and rows:               |
               \************************************************************/
              else // row > 1 && column > 1
              {
                 /******************************************************************\ 
                | Since we are behind the second column, there could be a preceeding |
                | gap in the mRNA thus we need to check whether it is better to open |
                | a new one or to extend the existing one (i.e. the full recursive   |
                | step of the open mRNA gap score formular) and since we are behind  |
                | the second row, the same argumentation holds for the miRNA.        |
                | Since we are neither in the first row nor in the first column we   |
                | can try to match the two nucleotides without adding a gap in one   |
                | of the sequences. After calculating all possible scores we add all |
                | those predecessors leading to the maximum to the overall score     |
                | matrix cell (i.e. the full recursive step of the overall score     |
                | formular).                                                         |
                 \******************************************************************/
                matrix_mRNA_gap[index] = openGapMatrixCell(matrix_mRNA_gap[index_left],matrix_overall[index_left],
                                                           *miRNA_it,mRNA_gap,match_pos);
                matrix_miRNA_gap[index] = openGapMatrixCell(matrix_miRNA_gap[index_up],matrix_overall[index_up],
                                                           miRNA_gap,*mRNA_it,match_pos);
                matrix_overall[index] = overallMatrixCell(matrix_overall[index_upleft],matrix_miRNA_gap[index],
                                                          matrix_mRNA_gap[index],*miRNA_it,*mRNA_it,match_pos);
              } // row > 1 && column > 1
               /*****************************************************************\ 
              | Each time we completed another row of the matrices we have        |
              | calculated another overall alignment score which might be yet the |
              | the best score:                                                   |
               \*****************************************************************/
              if(column==miRNA_length-1 && matrix_overall[index].get_score()>best_overall_score) // new best overall score
              {
                best_overall_score = matrix_overall[index].get_score();
              } // new best overall score
            } // row > 1
          }  // mRNA_it
        } // miRNA_it
      } // mRNA_length != 0 && miRNA_length != 0
       /*******************************************************\ 
      | In the end we return the score of an optimal alignment: |
       \*******************************************************/
      return best_overall_score;
      
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
    
    void optimalAlignmentList::backtrace_alignments(const alignmentMatrixCell * cell, std::vector<alignment> & alignment_vector, bool first_call, alignmentScore the_score, std::vector<alignmentColumn> * postfix)
    {
       /******************************************************************\ 
      | Checkwhether there is no more cell to trace move through and if so |
      | create the alignment from the postfix vector. Note that we collect |
      | the postfix values in the inverse order because vectors have good  |
      | performance if the insert and remove operations take place at the  |
      | end and there will be a lot more insert and remove operations than |
      | alignment creations (in fact that more often as the average        |
      | alignment has columns) justifying this approach:                   |
       \******************************************************************/
      if(cell == NULL) // upper-left corner
      {
        std::vector<alignmentColumn> alignment_columns(postfix->size());
        std::reverse_copy(postfix->begin(),postfix->end(),alignment_columns.begin());
        alignment_vector.push_back(alignment(alignment_columns,the_score));
      } // upper-left corner
       /******************************************************************\ 
      | Check whether there is no entry in the cell, if so state an error: |
       \******************************************************************/
      else if(cell->begin()==cell->end()) // undefined
      {
        std::cerr << "microSNPscore::optimalAlignmentList::backtrace_alignments\n";
        std::cerr << " ==> unitialized alignment matrix cell (call fill_matrices before)\n";
        std::cerr << "  --> no further alignment will be added\n";
      } // undefined
      else // somewhere in the matrix
      {
         /******************************************************************\ 
        | Check whether this is the initial call and if so adjust the score: |
         \******************************************************************/
        if(first_call) // right-most column
        {
          the_score = cell->get_score();
        } // right-most column
         /******************************************************************\ 
        | Iterate entries adding their alignment column to the posfix (mind  |
        | the inverse order) recursive calling the traceback from its        |
        | predecessor. Note that we need to remove the added coloumn         |
        | afterwards because the instances of this recursive function share  |
        | one prefix vector to avoid the need to copy it each time:          |
         \******************************************************************/
        for(alignmentMatrixCell::const_iterator entry_it(cell->begin());entry_it!=cell->end();++entry_it)
        {
          postfix->push_back(entry_it->get_column());
          backtrace_alignments(entry_it->get_predecessor(),alignment_vector,false,the_score,postfix);
          postfix->pop_back();
        }
      } // somewhere in the matrix
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

/*****************************************************************//**
* @brief output stream alignment insertion operator
*
* This operator is used to insert an alignment to an output stream
* (e.g. to print it on screen).
* The alignment will be represented by the chromosome positions of the
* first and last mRNA nucleotide and its seed type and score followed
* by the miRNA (5' to 3') over the match symbols and the mRNA (3' to
* 5').
*
* @param the_stream output stream the seed type should be inserted in
* @param the_alignment alignment to be inserted in the output stream
*
* @return output stream with the inserted alignment
*********************************************************************/
std::ostream & operator<<(std::ostream & the_stream, const alignment & the_alignment)
{
   /*****************************************************************\ 
  | Insert range of first and last mRNA nucleotide as well as the     |
  | alignments seed type and overall score for non-empty alignments   |
  | and iterate over the alignment for every line inserting the miRNA |
  | nucleotides, the match types (except for the first column because |
  | the miRNA 3' end will not bind to the mRNA since it is needed for |
  | RISC (RNA-induced silencing complex) binding) and the mRNA        |
  | nucleotides, respectively:                                        |
   \*****************************************************************/
  if(the_alignment.begin() != the_alignment.end())
  {
    the_stream << "\nmRNA range: ";
    the_stream << the_alignment.begin()->get_mRNA_nucleotide().get_chromosome_position();
    the_stream << " - ";
    the_stream << (the_alignment.end()-1)->get_mRNA_nucleotide().get_chromosome_position();
    the_stream << "\nseed type: ";
    the_stream << the_alignment.get_seed_type();
    the_stream << "\nscore: ";
    the_stream << the_alignment.get_score();
    the_stream << "\n\nmiRNA\t5'    ";
    for(alignment::const_iterator column_it(the_alignment.begin());column_it!=the_alignment.end();++column_it)
    {
      the_stream << column_it->get_miRNA_nucleotide();
    }
    the_stream << "    3'\n\t       ";
    for(alignment::const_iterator column_it(the_alignment.begin()+1);column_it!=the_alignment.end();++column_it)
    {
      the_stream << column_it->get_match();
    }
    the_stream << "\nmRNA\t3' ...";
    for(alignment::const_iterator column_it(the_alignment.begin());column_it!=the_alignment.end();++column_it)
    {
      the_stream << column_it->get_mRNA_nucleotide();
    }
    the_stream << "... 5'\n";
  }
  else
  {
    the_stream << "\n[empty alignment]\n";
  }
  return the_stream;
}

/*****************************************************************//**
* @brief output stream optimalAlignmentList insertion operator
*
* This operator is used to insert an optimal alignment list to an
* output stream (e.g. to print it on screen).
* The alignment list will be represented by its alignments separated
* by an empty line.
*
* @param the_stream output stream the seed type should be inserted in
* @param alignment_list optimalAlignmentList to be inserted in the
*     output stream
*
* @return output stream with the inserted optimal alignment list
*********************************************************************/
std::ostream & operator<<(std::ostream & the_stream, const optimalAlignmentList & alignment_list)
{
   /***************************************************************\ 
  | Insert the first alignment for non-empty lists and iterate over |
  | the list inserting the alignments:                              |
   \***************************************************************/
  if(alignment_list.begin() != alignment_list.end())
  {
    the_stream << *alignment_list.begin();
    for(optimalAlignmentList::const_iterator alignment_it(alignment_list.begin()+1);alignment_it!=alignment_list.end();++alignment_it)
    {
      the_stream << "\n" << std::endl;
      the_stream << *alignment_it;
    }
  }
  else
  {
    the_stream << "\n[empty alignment list]\n";
  }
  return the_stream;
}

} // namespace microSNPscore
