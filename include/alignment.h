#ifndef MICROSNPSCORE_ALIGNMENT_H
#define MICROSNPSCORE_ALIGNMENT_H


#include "nucleotide.h"
#include <vector>

namespace microSNPscore { class mRNA; } 
namespace microSNPscore { class miRNA; } 

namespace microSNPscore {

class alignmentMatrixCell; // to allow reverse dependency with alignmentMatrixCellEntry
class overallMatrixCell; // to allow reverse dependency with openGapMatrixCell
/*****************************************************************//**
* @brief Alignment column class
*
* This represent a column of an alignment of mRNA and miRNA.
*********************************************************************/
class alignmentColumn {
  public:
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
    
    alignmentColumn(const nucleotide & the_mRNA_nucleotide = nucleotide(), const nucleotide & the_miRNA_nucleotide = nucleotide(), matchPosition position = ThreePrime, IndelType indel_type = Open);

    /*****************************************************************//**
    * @brief get method for messenger RNA nucleotide attribute
    *
    * This method is used to access the nucleotide of the messenger RNA
    * that is aligned in that alignment column.
    *
    * @return the mRNA nucleotide aligned in that alignment column
    *********************************************************************/
    inline const nucleotide get_mRNA_nucleotide() const;

    /*****************************************************************//**
    * @brief get method for microRNA nucleotide attribute
    *
    * This method is used to access the nucleotide of the microRNA
    * that is aligned in that alignment column.
    *
    * @return the miRNA nucleotide aligned in that alignment column
    *********************************************************************/
    inline const nucleotide get_miRNA_nucleotide() const;

    /*****************************************************************//**
    * @brief get method for microRNA nucleotide attribute
    *
    * This method is used to access the nucleotide of the microRNA
    * that is aligned in that alignment column.
    *
    * @return the miRNA nucleotide aligned in that alignment column
    *********************************************************************/
    inline const matchType get_match() const;


  private:
    /*****************************************************************//**
    * @brief messenger RNA nucleotide
    *
    * This is the nucleotide of the messenger RNA that is aligned in that
    * alignment column.
    * It should be const but because alignmentColumns shall be used in a
    * vector and std::vector tries to assign its elements to an internal
    * array it needs a working assignment operator which has to change the
    * object's members and therefore they cannot be declared const.
    * Nevertheless this attribute is not intended to be changed in any
    * other context than assigning an initialized object to an unitialized
    * one produced by the standard constructor which is in fact not
    * designed to be used directly but only provided to allow array
    * allocation which is needed to create containers, too.
    *
    * @see alignmentColumn()
    *********************************************************************/
    nucleotide mRNA_nucleotide;

    /*****************************************************************//**
    * @brief microRNA nucleotide
    *
    * This is the nucleotide of the microRNA that is aligned in that
    * alignment column.
    * It should be const but because alignmentColumns shall be used in a
    * vector and std::vector tries to assign its elements to an internal
    * array it needs a working assignment operator which has to change the
    * object's members and therefore they cannot be declared const.
    * Nevertheless this attribute is not intended to be changed in any
    * other context than assigning an initialized object to an unitialized
    * one produced by the standard constructor which is in fact not
    * designed to be used directly but only provided to allow array
    * allocation which is needed to create containers, too.
    *
    * @see alignmentColumn()
    *********************************************************************/
    nucleotide miRNA_nucleotide;

    /*****************************************************************//**
    * @brief match state
    *
    * This is the match state (IndelOpen, IndelExtend, Mismatch, Masked,
    * Wobble, Match) of the alignment column.
    * It should be const but because alignmentColumns shall be used in a
    * vector and std::vector tries to assign its elements to an internal
    * array it needs a working assignment operator which has to change the
    * object's members and therefore they cannot be declared const.
    * Nevertheless this attribute is not intended to be changed in any
    * other context than assigning an initialized object to an unitialized
    * one produced by the standard constructor which is in fact not
    * designed to be used directly but only provided to allow array
    * allocation which is needed to create containers, too.
    *
    * @see alignmentColumn()
    *********************************************************************/
    
    matchType match;

};
    /*****************************************************************//**
    * @brief get method for messenger RNA nucleotide attribute
    *
    * This method is used to access the nucleotide of the messenger RNA
    * that is aligned in that alignment column.
    *
    * @return the mRNA nucleotide aligned in that alignment column
    *********************************************************************/
    inline const nucleotide alignmentColumn::get_mRNA_nucleotide() const {
      return mRNA_nucleotide;
    }

    /*****************************************************************//**
    * @brief get method for microRNA nucleotide attribute
    *
    * This method is used to access the nucleotide of the microRNA
    * that is aligned in that alignment column.
    *
    * @return the miRNA nucleotide aligned in that alignment column
    *********************************************************************/
    inline const nucleotide alignmentColumn::get_miRNA_nucleotide() const {
      return miRNA_nucleotide;
    }

    /*****************************************************************//**
    * @brief get method for microRNA nucleotide attribute
    *
    * This method is used to access the nucleotide of the microRNA
    * that is aligned in that alignment column.
    *
    * @return the miRNA nucleotide aligned in that alignment column
    *********************************************************************/
    inline const matchType alignmentColumn::get_match() const {
      return match;
    }

/*****************************************************************//**
* @brief alignment score type
*
* This represents the overall score of an alignment.
*********************************************************************/
typedef short alignmentScore;
/*****************************************************************//**
* @brief seed type type
*
* This represents the seed type (sixMer, sevenMerAOne, sevenMerMEight,
* eightMer) of a mRNA:miRNA alignement.
*********************************************************************/

enum seedType {
  sixMer,
  sevenMerAOne,
  sevenMerMEight,
  eightMer

};
/*****************************************************************//**
* @brief mRNA:miRNA alignment class
*
* This represents an alignment of mRNA and miRNA.
*********************************************************************/
class alignment {
  public:
    /*****************************************************************//**
    * @brief const iterartor type
    *
    * This type is used to access the alignment's columns.
    *********************************************************************/
    typedef std::vector<alignmentColumn>::const_iterator const_iterator;

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
    alignment(const std::vector<alignmentColumn> & the_columns, alignmentScore the_score);

    /*****************************************************************//**
    * @brief alignment begin
    *
    * This is used to get the first column (beeing the 5' end of the miRNA
    * and thus the 3' end of the mRNA) of the alignment.
    *
    * @return const_iterator pointing to the first alignment column
    *********************************************************************/
    inline const_iterator begin() const;

    /*****************************************************************//**
    * @brief alignment begin
    *
    * This is used to get the end of the alignment's column vector (beeing
    * the 3' end of the miRNA and thus the 5' end of the mRNA).
    *
    * @return const_iterator pointing behind the last alignment column
    *********************************************************************/
    inline const_iterator end() const;

    /*****************************************************************//**
    * @brief get method for score attribute
    *
    * This method is used to access the alignment's score
    *
    * @return the score of the alignment
    *********************************************************************/
    inline const alignmentScore get_score() const;

    /*****************************************************************//**
    * @brief get method for seed type attribute
    *
    * This method is used to access the alignment's seed type (sixMer,
    * sevenMerAOne, sevenMerMEight, eightMer).
    *
    * @return the seed type of the alignment
    *********************************************************************/
    inline const seedType get_seed_type() const;


  private:
    /*****************************************************************//**
    * @brief alignment columns
    *
    * A vector containing the alignment's columns from miRNA 5' to 3'
    *********************************************************************/
    std::vector<alignmentColumn> columns;

    alignmentScore score;

    /*****************************************************************//**
    * @brief seed type
    *
    * This is the seed type (sixMer, sevenMerAOne, sevenMerMEight,
    * EightMer) of the alignment.
    *********************************************************************/
    seedType seed_type;

};
    /*****************************************************************//**
    * @brief alignment begin
    *
    * This is used to get the first column (beeing the 5' end of the miRNA
    * and thus the 3' end of the mRNA) of the alignment.
    *
    * @return const_iterator pointing to the first alignment column
    *********************************************************************/
    inline alignment::const_iterator alignment::begin() const {
      return columns.end();
}

    /*****************************************************************//**
    * @brief alignment begin
    *
    * This is used to get the end of the alignment's column vector (beeing
    * the 3' end of the miRNA and thus the 5' end of the mRNA).
    *
    * @return const_iterator pointing behind the last alignment column
    *********************************************************************/
    inline alignment::const_iterator alignment::end() const {
      return columns.begin();
}

    /*****************************************************************//**
    * @brief get method for score attribute
    *
    * This method is used to access the alignment's score
    *
    * @return the score of the alignment
    *********************************************************************/
    inline const alignmentScore alignment::get_score() const {
      return score;
    }

    /*****************************************************************//**
    * @brief get method for seed type attribute
    *
    * This method is used to access the alignment's seed type (sixMer,
    * sevenMerAOne, sevenMerMEight, eightMer).
    *
    * @return the seed type of the alignment
    *********************************************************************/
    inline const seedType alignment::get_seed_type() const {
      return seed_type;
    }

/*****************************************************************//**
* @brief alignment matrix cell class
*
* This represents a cell in an alignment matrix.
*********************************************************************/
class alignmentMatrixCellEntry {
  public:
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
    
    alignmentMatrixCellEntry(const alignmentColumn & the_column = alignmentColumn(), const alignmentMatrixCell * the_predecessor = NULL);

    /*****************************************************************//**
    * @brief get method for alignment column attribute
    *
    * This method is used to access the alignment column corresponding to
    * the entry.
    *
    * @return the entry's alignment cloumn
    *********************************************************************/
    inline const alignmentColumn get_column() const;

    /*****************************************************************//**
    * @brief get method for predecessor attribute
    *
    * This method is used to access the cell that predecesses the entry in
    * an optimal alignment up to its cell.
    *
    * @return a pointer to the entry's predecessor
    *********************************************************************/
    inline const alignmentMatrixCell * get_predecessor() const;


  private:
    /*****************************************************************//**
    * @brief alignment column
    *
    * This is the alignment column corresponding to this entry.
    *********************************************************************/
    alignmentColumn column;

    /*****************************************************************//**
    * @brief predecessor vector
    *
    * This is a pointers to the cell that predecesses this entry in an
    * optimal alignment up to this cell.
    *********************************************************************/
    const alignmentMatrixCell * predecessor;

};
    /*****************************************************************//**
    * @brief get method for alignment column attribute
    *
    * This method is used to access the alignment column corresponding to
    * the entry.
    *
    * @return the entry's alignment cloumn
    *********************************************************************/
    inline const alignmentColumn alignmentMatrixCellEntry::get_column() const {
      return column;
    }

    /*****************************************************************//**
    * @brief get method for predecessor attribute
    *
    * This method is used to access the cell that predecesses the entry in
    * an optimal alignment up to its cell.
    *
    * @return a pointer to the entry's predecessor
    *********************************************************************/
    inline const alignmentMatrixCell * alignmentMatrixCellEntry::get_predecessor() const {
      return predecessor;
    }

/*****************************************************************//**
* @brief alignment matrix cell class
*
* This represents a cell in an alignment matrix it is thought as base
* class for the cells of the two different kinds of alignment matrices
* (the one containing overall scores and those storing alignments that
* end with a gap in their associated sequence).
*********************************************************************/
class alignmentMatrixCell {
  public:
    /*****************************************************************//**
    * @brief const iterator type
    *
    * This type is used to access the entries of the alignment matrix
    * cell.
    *********************************************************************/
    
    typedef std::vector<alignmentMatrixCellEntry>::const_iterator const_iterator;

    /*****************************************************************//**
    * @brief entry begin
    *
    * This is used to get the first entry of the cell.
    *
    * @return const_iterator pointing to the first entry of the cell
    *********************************************************************/
    inline const_iterator begin() const;

    /*****************************************************************//**
    * @brief entry end
    *
    * This is used to get the end of the cell's entry vector.
    *
    * @return const_iterator pointing behind last entry of the cell
    *********************************************************************/
    inline const_iterator end() const;

    /*****************************************************************//**
    * @brief get method for alignment score attribute
    *
    * This method is used to access the score of an optimal aligment up to
    * the cell.
    *
    * @return the alignment score of the cell
    *********************************************************************/
    inline const alignmentScore get_score() const;


  protected:
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
    
    alignmentMatrixCell(const std::vector<alignmentMatrixCellEntry> & the_entries = std::vector<alignmentMatrixCellEntry>(), const alignmentScore the_score = 0);

    /*****************************************************************//**
    * @brief cell entry vector
    *
    * This is a vector containing the cell's entries.
    *********************************************************************/
    std::vector<alignmentMatrixCellEntry> entries;

    /*****************************************************************//**
    * @brief alignment score
    *
    * This is the score of an optimal aligment up to this cell.
    *********************************************************************/
    alignmentScore score;

};
    /*****************************************************************//**
    * @brief entry begin
    *
    * This is used to get the first entry of the cell.
    *
    * @return const_iterator pointing to the first entry of the cell
    *********************************************************************/
    inline alignmentMatrixCell::const_iterator alignmentMatrixCell::begin() const {
      return entries.begin();
}

    /*****************************************************************//**
    * @brief entry end
    *
    * This is used to get the end of the cell's entry vector.
    *
    * @return const_iterator pointing behind last entry of the cell
    *********************************************************************/
    inline alignmentMatrixCell::const_iterator alignmentMatrixCell::end() const {
      return entries.end();
}

    /*****************************************************************//**
    * @brief get method for alignment score attribute
    *
    * This method is used to access the score of an optimal aligment up to
    * the cell.
    *
    * @return the alignment score of the cell
    *********************************************************************/
    inline const alignmentScore alignmentMatrixCell::get_score() const {
      return score;
    }

/*****************************************************************//**
* @brief overall-matrix cell class
*
* This represents a cell in one of the open-gap-matrices (those
* containing the score of the best alignment up to each cell that end
* with a gap in their associated sequence).
*********************************************************************/

class openGapMatrixCell : public alignmentMatrixCell {
  public:
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
    
    openGapMatrixCell(const openGapMatrixCell & preceding_open_gap_cell, const overallMatrixCell & preceding_overall_cell, const nucleotide & miRNA_nucleotide, const nucleotide & mRNA_nucleotide, matchPosition match_position);

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
    openGapMatrixCell(const openGapMatrixCell & preceding_open_gap_cell, const nucleotide & miRNA_nucleotide, const nucleotide & mRNA_nucleotide, matchPosition match_position);

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
    
    openGapMatrixCell(const overallMatrixCell & preceding_overall_cell, const nucleotide & miRNA_nucleotide, const nucleotide & mRNA_nucleotide, matchPosition match_position = ThreePrime);

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
    openGapMatrixCell();

};
/*****************************************************************//**
* @brief overall-matrix cell class
*
* This represents a cell in the overall-matrix (the one containing
* the score of the overall best alignment up to each cell).
*********************************************************************/

class overallMatrixCell : public alignmentMatrixCell {
  public:
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
    
    overallMatrixCell(const overallMatrixCell & upper_left_overall_cell, const openGapMatrixCell & miRNA_open_gap_cell, const openGapMatrixCell & mRNA_open_gap_cell, const nucleotide & miRNA_nucleotide, const nucleotide & mRNA_nucleotide, matchPosition match_position);

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
    
    overallMatrixCell(const openGapMatrixCell & open_gap_cell);

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
    overallMatrixCell(const nucleotide & miRNA_five_prime, const nucleotide & mRNA_three_prime);

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
    overallMatrixCell();

};
/*****************************************************************//**
* @brief optimal alignment list class
*
* This represents a list of all optimal alignments between a mRNA and
* a miRNA.
* It does NOT implement the std::list<> interface.
*********************************************************************/

class optimalAlignmentList {
  public:
    /*****************************************************************//**
    * @brief const iterator type
    *
    * This type is used to access the optimal alignments in the list.
    *********************************************************************/
    typedef std::vector<alignment>::const_iterator const_iterator;

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
    optimalAlignmentList(const mRNA & the_mRNA, const miRNA & the_miRNA);

    /*****************************************************************//**
    * @brief alignment list begin
    *
    * This is used to get the first optimal alignment in the list.
    *
    * @return const_iterator pointing to the first optimal alignment
    *********************************************************************/
    inline const_iterator begin() const;

    /*****************************************************************//**
    * @brief alignment list end
    *
    * This is used to get the end of the optimal alignment list.
    *
    * @return const_iterator pointing behind the last optimal alignment
    *********************************************************************/
    inline const_iterator end() const;


  private:
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
    static alignmentScore fill_matrices(openGapMatrixCell * matrix_mRNA_gap, openGapMatrixCell * matrix_miRNA_gap, overallMatrixCell * matrix_overall, const mRNA & the_mRNA, const miRNA & the_miRNA);

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
    
    static void backtrace_alignments(const alignmentMatrixCell * cell, std::vector<alignment> & alignment_vector, bool first_call = true, alignmentScore the_score = 0, std::vector<alignmentColumn> * postfix = new std::vector<alignmentColumn>);

    /*****************************************************************//**
    * @brief alignment vector
    *
    * This is a vector containing the optimal alignments.
    *********************************************************************/
    std::vector<alignment> alignments;

};
    /*****************************************************************//**
    * @brief alignment list begin
    *
    * This is used to get the first optimal alignment in the list.
    *
    * @return const_iterator pointing to the first optimal alignment
    *********************************************************************/
    inline optimalAlignmentList::const_iterator optimalAlignmentList::begin() const {
      return alignments.begin();
}

    /*****************************************************************//**
    * @brief alignment list end
    *
    * This is used to get the end of the optimal alignment list.
    *
    * @return const_iterator pointing behind the last optimal alignment
    *********************************************************************/
    inline optimalAlignmentList::const_iterator optimalAlignmentList::end() const {
      return alignments.end();
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
std::ostream & operator<<(std::ostream & the_stream, const seedType & seed_type);

} // namespace microSNPscore
#endif
