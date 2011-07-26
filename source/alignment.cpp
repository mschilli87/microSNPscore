
#include <algorithm>
// for std::max (alignment score calculation)
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
          const matchPosition match_pos (seed_start <= column && column <= seed_end ? Seed : ThreePrime);
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
          for(sequence::const_iterator mRNA_it(the_mRNA.end());mRNA_it>=the_mRNA.begin();--mRNA_it,++row)
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
            | But in all cases we need a vector of const alignment matrix cell   |
            | pointers as predecessor parameter for the alignment matrix cell's  |
            | constructor since we will create alignment matrix cells in any     |
            | case:                                                              |
             \******************************************************************/
            std::vector<const alignmentMatrixCell *> predecessors;
             /******************************************************************\ 
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
                matrix_mRNA_gap[index] = alignmentMatrixCell(&*mRNA_it,&*miRNA_it,0,predecessors); // dummy
                matrix_miRNA_gap[index] = alignmentMatrixCell(&*mRNA_it,&*miRNA_it,0,predecessors); // dummy
                matrix_overall[index] = alignmentMatrixCell(&*mRNA_it,&*miRNA_it,0,predecessors);
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
                matrix_miRNA_gap[index] = alignmentMatrixCell(&*mRNA_it,&*miRNA_it,0,predecessors); // dummy
                predecessors.push_back(&matrix_overall[index_left]);
                matrix_mRNA_gap[index] = alignmentMatrixCell(&*mRNA_it,&*miRNA_it,
                                                             miRNA_it->get_match(mRNA_gap,match_pos,Open).get_score(),
                                                             predecessors);
                matrix_overall[index] = matrix_mRNA_gap[index];
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
                matrix_miRNA_gap[index] = alignmentMatrixCell(&*mRNA_it,&*miRNA_it,0,predecessors); // dummy
                predecessors.push_back(&matrix_mRNA_gap[index_left]);
                matrix_mRNA_gap[index] = alignmentMatrixCell(&*mRNA_it,&*miRNA_it,
                                                             matrix_mRNA_gap[index_left].get_score()+
                                                             miRNA_it->get_match(mRNA_gap,match_pos,Extend).get_score(),
                                                             predecessors);
                matrix_overall[index] = matrix_mRNA_gap[index];
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
                matrix_mRNA_gap[index] = alignmentMatrixCell(&*mRNA_it,&*miRNA_it,0,predecessors); // dummy
                predecessors.push_back(&matrix_overall[index_up]);
                matrix_miRNA_gap[index] = alignmentMatrixCell(&*miRNA_it,&*miRNA_it,
                                                              mRNA_it->get_match(miRNA_gap,match_pos,Open).get_score(),
                                                              predecessors);
                matrix_overall[index] = matrix_miRNA_gap[index];
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
                | formular).                                                         |
                | Note that we need to remove the predecessors from our vector       |
                | before we can re-use it for another alignment matrix cell:         |
                 \******************************************************************/
                predecessors.push_back(&matrix_overall[index_left]);
                matrix_mRNA_gap[index] = alignmentMatrixCell(&*mRNA_it,&*miRNA_it,
                                                             matrix_overall[index_left].get_score()+
                                                             miRNA_it->get_match(mRNA_gap,match_pos,Open).get_score(),
                                                             predecessors);
                predecessors.pop_back();
                predecessors.push_back(&matrix_overall[index_up]);
                matrix_miRNA_gap[index] = alignmentMatrixCell(&*miRNA_it,&*miRNA_it,
                                                              matrix_overall[index_up].get_score()+
                                                              mRNA_it->get_match(miRNA_gap,match_pos,Open).get_score(),
                                                              predecessors);
                predecessors.pop_back();
                alignmentScore best_score(std::max(matrix_mRNA_gap[index].get_score(),matrix_miRNA_gap[index].get_score()));
                alignmentScore match_score(matrix_overall[index_upleft].get_score()+
                                           mRNA_it->get_match(*miRNA_it,match_pos).get_score());
                if(match_score >= best_score)
                {
                  best_score = match_score;
                  predecessors.push_back(&matrix_overall[index_upleft]);
                }
                if(matrix_mRNA_gap[index].get_score() == best_score)
                {
                  predecessors.push_back(&matrix_overall[index_left]);
                }
                if(matrix_miRNA_gap[index].get_score() == best_score)
                {
                  predecessors.push_back(&matrix_overall[index_up]);
                }
                matrix_overall[index] = alignmentMatrixCell(&*mRNA_it,&*miRNA_it,best_score,predecessors);
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
                | formular).                                                         |
                | Note that we need to remove the predecessors from our vector       |
                | before we can re-use it for another alignment matrix cell:         |
                 \******************************************************************/
                alignmentScore open_score(matrix_overall[index_left].get_score()+
                                          miRNA_it->get_match(mRNA_gap,match_pos,Open).get_score());
                alignmentScore extend_score(matrix_mRNA_gap[index_left].get_score()+
                                            miRNA_it->get_match(mRNA_gap,match_pos,Extend).get_score());
                alignmentScore best_score(std::max(open_score,extend_score));
                if(open_score==best_score)
                {
                  predecessors.push_back(&matrix_overall[index_left]);
                }
                if(extend_score==best_score)
                {
                  predecessors.push_back(&matrix_mRNA_gap[index_left]);
                }
                matrix_mRNA_gap[index] = alignmentMatrixCell(&*mRNA_it,&*miRNA_it,best_score,predecessors);
                predecessors.erase(predecessors.begin(),predecessors.end());
                predecessors.push_back(&matrix_overall[index_up]);
                matrix_miRNA_gap[index] = alignmentMatrixCell(&*miRNA_it,&*miRNA_it,
                                                              matrix_overall[index_up].get_score()+
                                                              mRNA_it->get_match(miRNA_gap,match_pos,Open).get_score(),
                                                              predecessors);
                predecessors.pop_back();
                best_score = std::max(best_score,matrix_miRNA_gap[index].get_score());
                alignmentScore match_score(matrix_overall[index_upleft].get_score()+
                                           mRNA_it->get_match(*miRNA_it,match_pos).get_score());
                if(match_score >= best_score)
                {
                  best_score = match_score;
                  predecessors.push_back(&matrix_overall[index_upleft]);
                }
                if(open_score == best_score)
                {
                  predecessors.push_back(&matrix_overall[index_left]);
                }
                if(extend_score == best_score)
                {
                  predecessors.push_back(&matrix_mRNA_gap[index_left]);
                }
                if(matrix_miRNA_gap[index].get_score() == best_score)
                {
                  predecessors.push_back(&matrix_overall[index_up]);
                }
                matrix_overall[index] = alignmentMatrixCell(&*mRNA_it,&*miRNA_it,best_score,predecessors);
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
                | in the mRNA (i.e. the gap open case of the open mRNA gap score    |
                | formular can be omitted).                                         |
                | Since we are in the first row there is no other option than       |
                | inserting a miRNA gap and thus the overall score equals the open  |
                | miRNA gap score in any case:                                      |
                 \*****************************************************************/
                matrix_mRNA_gap[index] = alignmentMatrixCell(&*mRNA_it,&*miRNA_it,0,predecessors); // dummy
                predecessors.push_back(&matrix_miRNA_gap[index_up]);
                matrix_miRNA_gap[index] = alignmentMatrixCell(&*miRNA_it,&*miRNA_it,
                                                              matrix_miRNA_gap[index_up].get_score()+
                                                              mRNA_it->get_match(miRNA_gap,match_pos,Extend).get_score(),
                                                              predecessors);
                matrix_overall[index] = matrix_miRNA_gap[index];
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
                | one (i.e. the full recursive step of the open mRNA gap score       |
                | formular).                                                         |
                | Since we are neither in the first row nor in the first column we   |
                | can try to match the two nucleotides without adding a gap in one   |
                | of the sequences. After calculating all possible scores we add all |
                | those predecessors leading to the maximum to the overall score     |
                | matrix cell (i.e. the full recursive step of the overall score     |
                | formular).                                                         |
                | Note that we need to remove the predecessors from our vector       |
                | before we can re-use it for another alignment matrix cell:         |
                 \******************************************************************/
                predecessors.push_back(&matrix_overall[index_left]);
                matrix_mRNA_gap[index] = alignmentMatrixCell(&*miRNA_it,&*miRNA_it,
                                                             matrix_overall[index_left].get_score()+
                                                             miRNA_it->get_match(mRNA_gap,match_pos,Open).get_score(),
                                                             predecessors);
                predecessors.pop_back();
                alignmentScore open_score(matrix_overall[index_up].get_score()+
                                          mRNA_it->get_match(miRNA_gap,match_pos,Open).get_score());
                alignmentScore extend_score(matrix_miRNA_gap[index_up].get_score()+
                                            mRNA_it->get_match(miRNA_gap,match_pos,Extend).get_score());
                alignmentScore best_score(std::max(open_score,extend_score));
                if(open_score==best_score)
                {
                  predecessors.push_back(&matrix_overall[index_up]);
                }
                if(extend_score==best_score)
                {
                  predecessors.push_back(&matrix_miRNA_gap[index_up]);
                }
                matrix_miRNA_gap[index] = alignmentMatrixCell(&*mRNA_it,&*miRNA_it,best_score,predecessors);
                predecessors.erase(predecessors.begin(),predecessors.end());
                best_score = std::max(best_score,matrix_mRNA_gap[index].get_score());
                alignmentScore match_score(matrix_overall[index_upleft].get_score()+
                                           mRNA_it->get_match(*miRNA_it,match_pos).get_score());
                if(match_score >= best_score)
                {
                  best_score = match_score;
                  predecessors.push_back(&matrix_overall[index_upleft]);
                }
                if(open_score == best_score)
                {
                  predecessors.push_back(&matrix_overall[index_up]);
                }
                if(extend_score == best_score)
                {
                  predecessors.push_back(&matrix_miRNA_gap[index_up]);
                }
                if(matrix_mRNA_gap[index].get_score() == best_score)
                {
                  predecessors.push_back(&matrix_overall[index_left]);
                }
                matrix_overall[index] = alignmentMatrixCell(&*mRNA_it,&*miRNA_it,best_score,predecessors);
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
                | Note that we need to remove the predecessors from our vector       |
                | before we can re-use it for another alignment matrix cell:         |
                 \******************************************************************/
                alignmentScore mRNA_open_score(matrix_overall[index_left].get_score()+
                                               miRNA_it->get_match(mRNA_gap,match_pos,Open).get_score());
                alignmentScore mRNA_extend_score(matrix_mRNA_gap[index_left].get_score()+
                                                 miRNA_it->get_match(mRNA_gap,match_pos,Extend).get_score());
                alignmentScore best_score(std::max(mRNA_open_score,mRNA_extend_score));
                if(mRNA_open_score==best_score)
                {
                  predecessors.push_back(&matrix_overall[index_left]);
                }
                if(mRNA_extend_score==best_score)
                {
                  predecessors.push_back(&matrix_mRNA_gap[index_left]);
                }
                matrix_mRNA_gap[index] = alignmentMatrixCell(&*mRNA_it,&*miRNA_it,best_score,predecessors);
                predecessors.erase(predecessors.begin(),predecessors.end());
                alignmentScore miRNA_open_score(matrix_overall[index_up].get_score()+
                                                mRNA_it->get_match(miRNA_gap,match_pos,Open).get_score());
                alignmentScore miRNA_extend_score(matrix_miRNA_gap[index_up].get_score()+
                                                  mRNA_it->get_match(miRNA_gap,match_pos,Extend).get_score());
                best_score = std::max(miRNA_open_score,miRNA_extend_score);
                if(miRNA_open_score==best_score)
                {
                  predecessors.push_back(&matrix_overall[index_up]);
                }
                if(miRNA_extend_score==best_score)
                {
                  predecessors.push_back(&matrix_miRNA_gap[index_up]);
                }
                matrix_miRNA_gap[index] = alignmentMatrixCell(&*mRNA_it,&*miRNA_it,best_score,predecessors);
                predecessors.erase(predecessors.begin(),predecessors.end());
                best_score = std::max(best_score,matrix_mRNA_gap[index].get_score());
                alignmentScore match_score(matrix_overall[index_upleft].get_score()+
                                           mRNA_it->get_match(*miRNA_it,match_pos).get_score());
                if(match_score >= best_score)
                {
                  best_score = match_score;
                  predecessors.push_back(&matrix_overall[index_upleft]);
                }
                if(mRNA_open_score == best_score)
                {
                  predecessors.push_back(&matrix_overall[index_left]);
                }
                if(mRNA_extend_score == best_score)
                {
                  predecessors.push_back(&matrix_mRNA_gap[index_left]);
                }
                if(miRNA_open_score == best_score)
                {
                  predecessors.push_back(&matrix_overall[index_up]);
                }
                if(miRNA_extend_score == best_score)
                {
                  predecessors.push_back(&matrix_miRNA_gap[index_up]);
                }
                matrix_overall[index] = alignmentMatrixCell(&*mRNA_it,&*miRNA_it,best_score,predecessors);
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
    
    void optimalAlignmentList::backtrace_alignments(alignmentMatrixCell * cell, std::vector<alignment> & alignment_vector, bool first_call, alignmentScore the_score, std::vector<alignmentColumn> * postfix)
    {
      if(first_call) // right-most column
      {
        the_score = cell->get_score();
      } // right-most column
      if(cell->predecessors_begin()==cell->predecessors_end()) // upper-left corner
      {
      //  postfix->push_back(alignmentColumn(*(cell->get_mRNA_nucleotide()),*(cell->get_miRNA_nucleotide())));
        // alignment_vector.push_back(alignment(std::vector(postfix->rbegin(),postfix->rend()),the_score));
      //  postfix->pop_back();
      } // upper-left corner
      else // somewhere in the middle
      {
        for(alignmentMatrixCell::const_iterator predecessor_it(cell->predecessors_begin());predecessor_it!=cell->predecessors_end();++predecessor_it)
        {
          const bool mRNA_move = cell->get_mRNA_nucleotide()->get_sequence_position() != (*predecessor_it)->get_mRNA_nucleotide()->get_sequence_position();
          const bool miRNA_move = cell->get_miRNA_nucleotide()->get_sequence_position() != (*predecessor_it)->get_miRNA_nucleotide()->get_sequence_position();
          if(mRNA_move && miRNA_move)
          {
      //      postfix->push_back(alignmentColumn(*(cell->get_mRNA_nucleotide()),*(cell->get_miRNA_nucleotide())));
          }
          else if(mRNA_move)
          {
          }
          else if(miRNA_move)
          {
          }
          else
          {
          }
        }
      } // somewhere in the middle
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
