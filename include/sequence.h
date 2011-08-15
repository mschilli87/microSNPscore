#ifndef MICROSNPSCORE_SEQUENCE_H
#define MICROSNPSCORE_SEQUENCE_H


#include <iostream>
//for std::ostream (operator<<)
#include <string>
#include "nucleotide.h"
#include <vector>

namespace microSNPscore { class nucleotide; } 
namespace microSNPscore { class conservationList; } 
namespace microSNPscore { class SNP; } 


namespace microSNPscore {

/*****************************************************************//**
* @brief sequence ID type
*
* This represents the ID of a sequence.
*********************************************************************/
typedef std::string sequenceID;
/*****************************************************************//**
* @brief chromosome type
*
* This represents a chromosome. It is defined as string to handle
* different notations (like "chr1" or "1" and "MIT" or "24") (but this
* is done without any consistency checking) and special 'chromosomes'
* (like "HSCHR12_3_CTG2_1" and "GL000195.1").
*********************************************************************/
typedef std::string chromosomeType;
/*****************************************************************//**
* @brief strand type
*
* This represents the strand of a sequence (Plus or Minus).
*********************************************************************/
enum strandType {
  Plus,
  Minus

};
typedef unsigned short sequenceLength;
/*****************************************************************//**
* @brief exon class
*
* This represents an exon as a pair of positions (start/end).
* All exons are defined by their position on the + strand in 5' --> 3'
* direction, no matter which strand they are placed on.
*********************************************************************/
class exon {
  public:
    /*****************************************************************//**
    * @brief constructor - Do not call without parameter values!
    *
    * This is used to create an instance of the exon class.
    * If @p end_position < @p start_position then a zero-sized exon
    * reaching from @p start_position to @p start_position - 1 is
    * created.
    * The default values are not intended to be used directly.
    * They are only provided to allow array allocation but you will need
    * to assign a valid object created by giving those parameters a value
    * to actually use it. This is done by containers like std::vector and
    * the reason for providing those default values is to allow using
    * containers containing objects of this class.
    *
    * @param start_position (pseudo-optional) chromosomePosition that
    *     represents the start position of the exon on the chromosome -
    *     Defaults to 0
    * @param end_position (pseudo-optional) chromosomePosition that
    *     represents the end position of the exon on the chromosome -
    *     Defaults to 0
    *
    * @return an exon located at the given positions on chromosome
    *********************************************************************/
    
    exon(chromosomePosition start_position = 0, chromosomePosition end_position = 0);

    /*****************************************************************//**
    * @brief get method for start position attribute
    *
    * This method is used to access the start position of the exon on the
    * chromosome (i.e. the end with the smaller distance to the chromosome
    * start beeing the 5' end of the + strand and accordingly the 3' end
    * of the - strand), the 5' end of the + strand (i.e. the 3' end of
    * the - strand) beeing position 1.
    *
    * @return the start position of the exon on the chromosome
    *********************************************************************/
    inline const chromosomePosition get_start() const;

    /*****************************************************************//**
    * @brief get method for end position attribute
    *
    * This method is used to access the end position of the exon on the
    * chromosome (i.e. the end with the smaller distance to the chromosome
    * end beeing the 3' end of the + strand and accordingly the 5' end
    * of the - strand), the 5' end of the + strand (i.e. the 3' end of
    * the - strand) beeing position 1.
    *
    * @return the end position of the exon on the chromosome
    *********************************************************************/
    inline const chromosomePosition get_end() const;

    /*****************************************************************//**
    * @brief length calculation
    *
    * This method is used to calculate the length of the exon.
    *
    * @return  the exon's length
    *********************************************************************/
    sequenceLength get_length() const;


  private:
    /*****************************************************************//**
    * @brief start position
    *
    * This is the start position of the exon on the chromosome (i.e. the
    * end with the smaller distance to the chromosome start beeing the 5'
    * end of the + strand and accordingly the 3' end of the - strand), the
    * 5' end of the + strand (i.e. the 3' end of the - strand) beeing
    * position 1.
    * It should be const but because exons shall be used in a vector and
    * std::vector tries to assign its elements to an internal array it
    * needs a working assignment operator which has to change the object's
    * members and therefore they cannot be declared const.
    * Nevertheless this attribute is not intended to be changed in any
    * other context than assigning an initialized object to an unitialized
    * one produced by the standard constructor which is in fact not
    * designed to be used directly but only provided to allow array
    * allocation which is needed to create containers, too.
    *
    * @see exon()
    *********************************************************************/
    chromosomePosition start;

    /*****************************************************************//**
    * @brief start position
    *
    * This is the end position of the exon on the chromosome (i.e. the
    * end with the smaller distance to the chromosome end beeing the 3'
    * end of the + strand and accordingly the 5' end of the - strand), the
    * 5' end of the + strand (i.e. the 3' end of the - strand) beeing
    * position 1.
    * It should be const but because exons shall be used in a vector
    * and std::vector tries to assign its elements to an internal array it
    * needs a working assignment operator which has to change the object's
    * members and therefore they cannot be declared const.
    * Nevertheless this attribute is not intended to be changed in any
    * other context than assigning an initialized object to an unitialized
    * one produced by the standard constructor which is in fact not
    * designed to be used directly but only provided to allow array
    * allocation which is needed to create containers, too.
    *
    * @see exon()
    *********************************************************************/
    chromosomePosition end;

};
    /*****************************************************************//**
    * @brief get method for start position attribute
    *
    * This method is used to access the start position of the exon on the
    * chromosome (i.e. the end with the smaller distance to the chromosome
    * start beeing the 5' end of the + strand and accordingly the 3' end
    * of the - strand), the 5' end of the + strand (i.e. the 3' end of
    * the - strand) beeing position 1.
    *
    * @return the start position of the exon on the chromosome
    *********************************************************************/
    inline const chromosomePosition exon::get_start() const {
      return start;
    }

    /*****************************************************************//**
    * @brief get method for end position attribute
    *
    * This method is used to access the end position of the exon on the
    * chromosome (i.e. the end with the smaller distance to the chromosome
    * end beeing the 3' end of the + strand and accordingly the 5' end
    * of the - strand), the 5' end of the + strand (i.e. the 3' end of
    * the - strand) beeing position 1.
    *
    * @return the end position of the exon on the chromosome
    *********************************************************************/
    inline const chromosomePosition exon::get_end() const {
      return end;
    }

/*****************************************************************//**
* @brief sequence class
*
* This represents a RNA sequence and is used as base class for the
* mRNA and miRNA classes.
*
* @see mRNA
* @see miRNA
*********************************************************************/

class sequence {
  public:
    /*****************************************************************//**
    * @brief const iterartor type
    *
    * This type is used to access the sequence's nucleotide.
    *********************************************************************/
    typedef std::vector<nucleotide>::const_iterator const_iterator;

    /*****************************************************************//**
    * @brief const iterartor type
    *
    * This type is used to access the sequence's exons.
    *********************************************************************/
    typedef std::vector<exon>::const_iterator const_exon_iterator;

    /*****************************************************************//**
    * @brief constructor
    *
    * This is used to create an instance of the class sequence.
    * Lowercase letters are treated as uppercase ones.
    * T is understood as Thymine and is treated as Uracil (simulating
    * transscription) raising an error message.
    * Dashes (-) are understood as Gaps and are omitted.
    * Other characters than A,a,C,c,G,g,U,u,T,t,X,x or - raise an error
    * and are treated as Mask.
    * The ordering of exon starts and ends does not matter.
    * Overlapping exons are be merged (reporting an error).
    * If the count of exon starts does not match the count of exon ends
    * an error is raised and the additional starts or ends are omitted.
    * If the calculated sequence length (as defined by the exons) does
    * not match the count of nucleotides an error message is raised and
    * the additional nucleotides are omitted or the missing nucleotides
    * are treated as masked, respectively.
    *
    * @param the_ID sequenceID representing the ID of the sequence
    * @param sequence_string String representing the nucleotide sequence
    *     (Adenine: A, Cytosine: C, Guanine: G, Uracil: U, Mask: X)
    * @param the_chromosome chromosomeType representing the chromosome the
    *     sequence is located on
    * @param the_strand strandType representing the strand (Plus/Minus) on
    *     which the sequence is located
    * @param exon_starts: String representing the start positions (i.e.
    *     the end with the smaller distance to the chromosome start beeing
    *     the 5' end of the + strand and accordingly the 3' end of
    *     the - strand) of the exons containing the sequence as
    *     comma-separated list.
    * @param exon_ends: String representing the end positions (i.e.
    *     the end with the smaller distance to the chromosome end beeing
    *     the 3' end of the + strand and accordingly the 5' end of
    *     the - strand) of the exons containing the sequence as
    *     comma-separated list.
    * @param conservations conservationList containing the conservaton
    *     ranges for the sequence
    *
    * @return a sequence containing the given nucleotides located on the
    *     given chromosome, strand and positions with the given
    *     conservation scores.
    *********************************************************************/
    sequence(sequenceID the_ID, std::string sequence_string, chromosomeType the_chromosome, strandType the_strand, std::string exon_starts, std::string exon_ends, const conservationList & conservations);

    /*****************************************************************//**
    * @brief standard constructor - do not use directly
    *
    * This is only provided to allow array allocation but you will need
    * to assign a valid object created by the parameterized constructor
    * to actually use it. This is done by containers like std::vector and
    * the reason for providing those default values is to allow using
    * containers containing objects of this class.
    *
    * @return: an unitialized sequence object
    *********************************************************************/
    sequence();

    /*****************************************************************//**
    * @brief get method for ID attribute
    *
    * This method is used to access the ID of the sequence.
    *
    * @return the ID of the sequence
    *********************************************************************/
    inline const sequenceID get_ID() const;

    /*****************************************************************//**
    * @brief get method for chromosome attribute
    *
    * This method is used to access the chromosome the sequence is
    * located on.
    *
    * @return the chromosome of the sequence
    *********************************************************************/
    inline const chromosomeType get_chromosome() const;

    /*****************************************************************//**
    * @brief get method for strand attribute
    *
    * This method is used to access the strand of the chromosome the
    * sequence is on.
    *
    * @return the strand of the sequence (Plus or Minus)
    *********************************************************************/
    inline const strandType get_strand() const;

    /*****************************************************************//**
    * @brief get method for length attribute
    *
    * This method is used to access the length of the sequence.
    *
    * @return the length of the sequence.
    *********************************************************************/
    inline const sequenceLength get_length() const;

    /*****************************************************************//**
    * @brief nucleotide vector begin
    *
    * This is used to get the first nucleotide of the sequence's vector.
    *
    * @return const_iterator pointing to the first nucleotide
    *********************************************************************/
    inline const_iterator begin() const;

    /*****************************************************************//**
    * @brief nucleotide vector end
    *
    * This is used to get the end of the sequence's nucleotide vector.
    *
    * @return const_iterator pointing behind the last nucleotide
    *********************************************************************/
    inline const_iterator end() const;

    /*****************************************************************//**
    * @brief index operator
    *
    * This method is used to access a nucleotide by sequence position, the
    * 5' end beeing position 1.
    * If the position is not part of the sequence, the sequence's end is
    * returned.
    *
    * @param position position of the queried nucleotide in the sequence
    *
    * @return const_iterator pointing to the queried sequence if it exists
    *     or behind the last nucleotide otherwise
    *********************************************************************/
    inline const_iterator operator[](const sequencePosition & position) const;

    /*****************************************************************//**
    * @brief get nucleotide by sequence position
    *
    * This method is used to access a nucleotide by sequence position, the
    * 5' end beeing position 1.
    * If the position is not part of the sequence, the sequence's end is
    * returned.
    *
    * @param position position of the queried nucleotide in the sequence
    *
    * @return const_iterator pointing to the queried nucleotide if it
    *     exists or behind the last nucleotide otherwise
    *********************************************************************/
    
    inline const_iterator get_nucleotide(const sequencePosition & position) const;

    /*****************************************************************//**
    * @brief get nucleotide by chromosome position
    *
    * This method is used to access a nucleotide by chromosome position,
    * the 5' end  of the + strand (i.e. the 3' end of the - strand) beeing
    * position 1.
    * If the position is not part of the sequence, the sequence's end is
    * returned.
    *
    * @param position position of the queried nucleotide on chromosome
    *
    * @return const_iterator pointing to the queried nucleotide if its
    *     exists or behind the last nucleotide otherwise
    *********************************************************************/
    
    inline const_iterator get_nucleotide_chr(const chromosomePosition & position) const;

    /*****************************************************************//**
    * @brief exon vector begin
    *
    * This is used to get the first exon of the sequence.
    *
    * @return const_iterator pointing to the first exon
    *********************************************************************/
    inline const_exon_iterator exons_begin() const;

    /*****************************************************************//**
    * @brief exon vector end
    *
    * This is used to get the end of the sequence's exon vector.
    *
    * @return const_iterator pointing behind the last exon
    *********************************************************************/
    inline const_exon_iterator exons_end() const;

    /*****************************************************************//**
    * @brief get subsequence from sequence position
    *
    * This method can be used to extract a subsequence of a given length
    * starting (i.e. 5' end) at a given position from the sequence.
    * If the length is too high so that the queried subsequence would
    * reach over the (3') end of the sequence, a shorter subsequence
    * starting at the disired start position and ending at the (3') end of
    * the sequence will be returned.
    * If the given position is not part of the sequence, an empty sequence
    * will be returned.
    *
    * @param from the start position in the sequence of the subsequence
    * @param len the (maximal) length of the subsequence
    * @return the subsequence starting at the given position and ending
    *         after the given length (5' to 3') or at the end of the
    *         sequence
    *********************************************************************/
    inline sequence get_subsequence_from(sequencePosition from, const sequenceLength & len) const;

    /*****************************************************************//**
    * @brief get subsequence to sequence position
    *
    * This method can be used to extract a subsequence of a given length
    * ending (i.e. 3' end) at a given position from the sequence.
    * If the length is too high so that the queried subsequence would
    * reach over the (5') end of the sequence, a shorter subsequence
    * ending at the disired end position and starting at the (5') end of
    * the sequence will be returned.
    * If the given position is not part of the sequence, an empty sequence
    * will be returned.
    *
    * @param to the end position in the sequence of the subsequence
    * @param len the (maximal) length of the subsequence
    * @return the subsequence ending at the given position and starting
    *         after the given length (3' to 5') or at the start of the
    *         sequence
    *********************************************************************/
    inline sequence get_subsequence_to(sequencePosition to, const sequenceLength & len) const;

    /*****************************************************************//**
    * @brief get subsequence between sequence positions
    *
    * This method can be used to extract a subsequence starting (i.e. 5'
    * end) end ending (i.e. 3' end) at given positions from the sequence.
    * If one of the given positions is not part of the sequence the,
    * corresponding end will be used instead to delimit the subsequence.
    * If from >= to an empty sequence will be returned.
    *
    * @param from the start position in the sequence of the subsequence
    * @param to the end position in the sequence of the subsequence
    * @return the subsequence starting and ending at the given positions
    *********************************************************************/
    sequence get_subsequence_from_to(sequencePosition from, sequencePosition to) const;

    /*****************************************************************//**
    * @brief get subsequence from chromosome position
    *
    * This method can be used to extract a subsequence of a given length
    * starting (i.e. 5' end) at a given chromosome position from the
    * sequence.
    * If the length is too high so that the queried subsequence would
    * reach over the (3') end of the sequence, a shorter subsequence
    * starting at the disired start position and ending at the (3') end of
    * the sequence will be returned.
    * If the given chromosome position is not part of the sequence, an
    * empty sequence will be returned.
    *
    * @param from the start position on the chromosome of the subsequence
    * @param len the (maximal) length of the subsequence
    * @return the subsequence starting at the given position and ending
    *         after the given length (5' to 3') or at the end of the
    *         sequence
    *********************************************************************/
    virtual sequence get_subsequence_chr_from(chromosomePosition from, const sequenceLength & len) const;

    /*****************************************************************//**
    * @brief get subsequence to chromosome position
    *
    * This method can be used to extract a subsequence of a given length
    * ending (i.e. 3' end) at a given chromosome position from the
    * sequence.
    * If the length is too high so that the queried subsequence would
    * reach over the (5') end of the sequence, a shorter subsequence
    * ending at the disired end position and starting at the (5') end of
    * the sequence will be returned.
    * If the given chromosome position is not part of the sequence, an
    * empty sequence will be returned.
    *
    * @param to the end position on the chromosome of the subsequence
    * @param len the (maximal) length of the subsequence
    * @return the subsequence ending at the given position and starting
    *         after the given length (3' to 5') or at the start of the
    *         sequence
    *********************************************************************/
    inline sequence get_subsequence_chr_to(chromosomePosition to, const sequenceLength & len) const;

    /*****************************************************************//**
    * @brief get subsequence between chromosome positions
    *
    * This method can be used to extract a subsequence starting (i.e. 5'
    * end) end ending (i.e. 3' end) at given chromosome positions from the
    * sequence.
    * If one of the given positions is not part of the sequence the,
    * corresponding end will be used instead to delimit the subsequence.
    * If from >= to an empty sequence will be returned.
    *
    * @param from the start position on the chromosome of the subsequence
    * @param to the end position on the chromosome of the subsequence
    * @return the subsequence starting and ending at the given positions
    *********************************************************************/
    sequence get_subsequence_chr_from_to(chromosomePosition from, chromosomePosition to) const;

    /*****************************************************************//**
    * @brief apply SNP on sequence
    *
    * This method is used to change a mRNA as described by a given SNP.
    * If the SNP information do not match those in the sequence an
    * unchanged copy of the seqeunce is returned (after stating an error).
    *
    * @param the_SNP SNP containing the information how the mRNA should be
    *     changed
    *
    * @return a copy of the sequence with the changes defined by the SNP
    *********************************************************************/
    sequence mutate(const SNP & the_SNP) const;


  private:
    /*****************************************************************//**
    * @brief internal constructor
    *
    * This is used to create an instance of the class sequence from
    * another one (e.g. to construct subsequences).
    *
    * @param the_ID sequenceID representing the ID of the sequence
    * @param the_chromosome chromosomeType representing the chromosome the
    *     sequence is located on
    * @param the_exons std::vector<exon> representing the exons on
    *     which the sequence is located
    * @param the_length sequenceLength representing the length of the
    *    sequence
    * @param the_nucleotides: std::vector<nucleotide> representing the
    *     sequence's nucleotides
    *
    * @return a sequence with the given attributes
    *********************************************************************/
    sequence(sequenceID the_ID, chromosomeType the_chromosome, strandType the_strand, std::vector<exon> the_exons, sequenceLength the_length, const std::vector<nucleotide> & the_nucleotides);

    /*****************************************************************//**
    * @brief exon initialisation
    *
    * This method is used to create the exon vector of a sequence.
    * The ordering of exon starts and ends does not matter.
    * Overlapping exons are be merged (reporting an error).
    * If the count of exon starts does not match the count of exon ends
    * an error is raised and the additional starts or ends are omitted.
    *
    * @param starts String representing the start positions (i.e.
    *     the end with the smaller distance to the chromosome start beeing
    *     the 5' end of the + strand and accordingly the 3' end of
    *     the - strand) of the exons containing the sequence as
    *     comma-separated list.
    * @param ends String representing the end positions (i.e.
    *     the end with the smaller distance to the chromosome end beeing
    *     the 3' end of the + strand and accordingly the 5' end of
    *     the - strand) of the exons containing the sequence as
    *     comma-separated list.
    *
    * @return a vector containing exons with the given coordinates
    *********************************************************************/
    static std::vector<exon> initialize_exons(const std::vector<chromosomePosition> & starts, const std::vector<chromosomePosition> & ends);

    /*****************************************************************//**
    * @brief length calculation
    *
    * This method is used to calculate the length of a sequence.
    *
    * @param begin_of_exons const_exon_iterator pointing the sequence's
    *     first exon
    * @param end_of_exons const_exon_iterator pointing behind the
    *     sequence's exon vector
    *
    * @return the sequence's length
    *********************************************************************/
    
    static sequenceLength initialize_length(const const_exon_iterator & begin_of_exons, const const_exon_iterator & end_of_exons);

    /*****************************************************************//**
    * @brief nucleotide initialization
    *
    * This method is used to calculate the nucleotide vector of a sequence.
    * Lowercase letters are treated as uppercase ones.
    * T is understood as Thymine and is treated as Uracil (simulating
    * transscription) raising an error message.
    * Dashes (-) are understood as Gaps and are omitted.
    * Other characters than A,a,C,c,G,g,U,u,T,t,X,x or - raise an error
    * and are treated as Mask.
    * If the given sequence length does
    * not match the count of nucleotides an error message is raised and
    * the additional nucleotides are omitted or the missing nucleotides
    * are treated as masked, respectively.
    *
    * @param the_sequence String representing the nucleotide sequence
    *     (Adenine: A, Cytosine: C, Guanine: G, Uracil: U, Mask: X)
    * @param the_chromosome chromosomeType representing the chromosome the
    *     sequence is located on
    * @param the_strand strandType representing the strand (Plus/Minus) on
    *     which the sequence is located
    * @param begin_of_exons const_exon_iterator pointing the sequence's
    *     first exon
    * @param end_of_exons const_exon_iterator pointing behind the
    *     sequence's exon vector
    * @param the_length the requested length of the sequence
    * @param conservations conservationList containing the conservaton
    *     ranges for the sequence
    *
    * @return a vector containing the sequence's nucleotides
    *********************************************************************/
    static std::vector<nucleotide> initialize_nucleotides(const std::string & the_sequence, chromosomeType the_chromosome, strandType the_strand
    , const const_exon_iterator & begin_of_exons, const const_exon_iterator & end_of_exons, sequenceLength the_length, const conservationList & conservations);

    /*****************************************************************//**
    * @brief string to position vector conversion
    *
    * This method is used to convert a string containing chromosome
    * positions as comma-separated list to a sorted vector of
    * chromosomePositions.
    * Illegal values raise an error message and are omitted.
    *
    * @param string_list: String representing positions on a chromosome as
    *     comma-separated list.
    *
    * @return a sorted vector containing the converted positions
    *********************************************************************/
    
    static std::vector<chromosomePosition> position_string_to_vector(std::string string_list);


  protected:
    /*****************************************************************//**
    * @brief chromosome position to sequence position conversion
    *
    * This method is used to convert a position on chromosome to the
    * corresponding position in the sequence.
    * If the given position is not part of the sequence, 0 is returned.
    *
    * @param chromosome_position the position on chromosome to convert
    *
    * @return the position in the sequence that corresponds to the given
    *     position on chromosome
    *********************************************************************/
    sequencePosition chromosome_position_to_sequence_position(chromosomePosition chromosome_position) const;


  private:
    /*****************************************************************//**
    * @brief sequence ID
    *
    * This is the ID of the sequence.
    *********************************************************************/
    const sequenceID ID;

    /*****************************************************************//**
    * @brief chromosome
    *
    * This is the chromosome the sequence is located on.
    *********************************************************************/
    const chromosomeType chromosome;

    /*****************************************************************//**
    * @brief strand on chromosome
    *
    * This is the chromosome strand (Plus or Minus) the sequence is on.
    *********************************************************************/
    const strandType strand;

    /*****************************************************************//**
    * @brief exon segmentation
    *
    * A vector containing the sequence's exons sorted by start position
    * (i.e. the end with the smaller distance to the chromosome start
    * beeing the 5' end of the + strand and accordingly the 3' end of
    * the - strand) on the chromosome
    *********************************************************************/
    const std::vector<exon> exons;

    /*****************************************************************//**
    * @brief sequence length
    *
    * This is the length of the sequence.
    *********************************************************************/
    const sequenceLength length;

    /*****************************************************************//**
    * @brief nucleotide sequence
    *
    * A vector containing the sequence's nucleotides from 5' to 3'
    *********************************************************************/
    const std::vector<nucleotide> nucleotides;

};
    /*****************************************************************//**
    * @brief get method for ID attribute
    *
    * This method is used to access the ID of the sequence.
    *
    * @return the ID of the sequence
    *********************************************************************/
    inline const sequenceID sequence::get_ID() const {
      return ID;
    }

    /*****************************************************************//**
    * @brief get method for chromosome attribute
    *
    * This method is used to access the chromosome the sequence is
    * located on.
    *
    * @return the chromosome of the sequence
    *********************************************************************/
    inline const chromosomeType sequence::get_chromosome() const {
      return chromosome;
    }

    /*****************************************************************//**
    * @brief get method for strand attribute
    *
    * This method is used to access the strand of the chromosome the
    * sequence is on.
    *
    * @return the strand of the sequence (Plus or Minus)
    *********************************************************************/
    inline const strandType sequence::get_strand() const {
      return strand;
    }

    /*****************************************************************//**
    * @brief get method for length attribute
    *
    * This method is used to access the length of the sequence.
    *
    * @return the length of the sequence.
    *********************************************************************/
    inline const sequenceLength sequence::get_length() const {
      return length;
    }

    /*****************************************************************//**
    * @brief nucleotide vector begin
    *
    * This is used to get the first nucleotide of the sequence's vector.
    *
    * @return const_iterator pointing to the first nucleotide
    *********************************************************************/
    inline sequence::const_iterator sequence::begin() const {
      return nucleotides.begin();
}

    /*****************************************************************//**
    * @brief nucleotide vector end
    *
    * This is used to get the end of the sequence's nucleotide vector.
    *
    * @return const_iterator pointing behind the last nucleotide
    *********************************************************************/
    inline sequence::const_iterator sequence::end() const {
      return nucleotides.end();
}

    /*****************************************************************//**
    * @brief index operator
    *
    * This method is used to access a nucleotide by sequence position, the
    * 5' end beeing position 1.
    * If the position is not part of the sequence, the sequence's end is
    * returned.
    *
    * @param position position of the queried nucleotide in the sequence
    *
    * @return const_iterator pointing to the queried sequence if it exists
    *     or behind the last nucleotide otherwise
    *********************************************************************/
    inline sequence::const_iterator sequence::operator[](const sequencePosition & position) const {
      return get_nucleotide(position);
}

    /*****************************************************************//**
    * @brief get nucleotide by sequence position
    *
    * This method is used to access a nucleotide by sequence position, the
    * 5' end beeing position 1.
    * If the position is not part of the sequence, the sequence's end is
    * returned.
    *
    * @param position position of the queried nucleotide in the sequence
    *
    * @return const_iterator pointing to the queried nucleotide if it
    *     exists or behind the last nucleotide otherwise
    *********************************************************************/
    
    inline sequence::const_iterator sequence::get_nucleotide(const sequencePosition & position) const {
      return position <= get_length() ? begin()+(position-1) : end();
}

    /*****************************************************************//**
    * @brief get nucleotide by chromosome position
    *
    * This method is used to access a nucleotide by chromosome position,
    * the 5' end  of the + strand (i.e. the 3' end of the - strand) beeing
    * position 1.
    * If the position is not part of the sequence, the sequence's end is
    * returned.
    *
    * @param position position of the queried nucleotide on chromosome
    *
    * @return const_iterator pointing to the queried nucleotide if its
    *     exists or behind the last nucleotide otherwise
    *********************************************************************/
    
    inline sequence::const_iterator sequence::get_nucleotide_chr(const chromosomePosition & position) const {
      return get_nucleotide(chromosome_position_to_sequence_position(position));
}

    /*****************************************************************//**
    * @brief exon vector begin
    *
    * This is used to get the first exon of the sequence.
    *
    * @return const_iterator pointing to the first exon
    *********************************************************************/
    inline sequence::const_exon_iterator sequence::exons_begin() const {
      return exons.begin();
}

    /*****************************************************************//**
    * @brief exon vector end
    *
    * This is used to get the end of the sequence's exon vector.
    *
    * @return const_iterator pointing behind the last exon
    *********************************************************************/
    inline sequence::const_exon_iterator sequence::exons_end() const {
      return exons.end();
}

/*****************************************************************//**
* @brief get subsequence from sequence position
*
* This method can be used to extract a subsequence of a given length
* starting (i.e. 5' end) at a given position from the sequence.
* If the length is too high so that the queried subsequence would
* reach over the (3') end of the sequence, a shorter subsequence
* starting at the disired start position and ending at the (3') end of
* the sequence will be returned.
* If the given position is not part of the sequence, an empty sequence
* will be returned.
*
* @param from the start position in the sequence of the subsequence
* @param len the (maximal) length of the subsequence
* @return the subsequence starting at the given position and ending
*         after the given length (5' to 3') or at the end of the
*         sequence
*********************************************************************/
inline sequence sequence::get_subsequence_from(sequencePosition from, const sequenceLength & len) const {
return get_subsequence_from_to(from,from+len-1);
}

/*****************************************************************//**
* @brief get subsequence to sequence position
*
* This method can be used to extract a subsequence of a given length
* ending (i.e. 3' end) at a given position from the sequence.
* If the length is too high so that the queried subsequence would
* reach over the (5') end of the sequence, a shorter subsequence
* ending at the disired end position and starting at the (5') end of
* the sequence will be returned.
* If the given position is not part of the sequence, an empty sequence
* will be returned.
*
* @param to the end position in the sequence of the subsequence
* @param len the (maximal) length of the subsequence
* @return the subsequence ending at the given position and starting
*         after the given length (3' to 5') or at the start of the
*         sequence
*********************************************************************/
inline sequence sequence::get_subsequence_to(sequencePosition to, const sequenceLength & len) const {
return get_subsequence_from_to(to-len+1,to);
}

/*****************************************************************//**
* @brief get subsequence to chromosome position
*
* This method can be used to extract a subsequence of a given length
* ending (i.e. 3' end) at a given chromosome position from the
* sequence.
* If the length is too high so that the queried subsequence would
* reach over the (5') end of the sequence, a shorter subsequence
* ending at the disired end position and starting at the (5') end of
* the sequence will be returned.
* If the given chromosome position is not part of the sequence, an
* empty sequence will be returned.
*
* @param to the end position on the chromosome of the subsequence
* @param len the (maximal) length of the subsequence
* @return the subsequence ending at the given position and starting
*         after the given length (3' to 5') or at the start of the
*         sequence
*********************************************************************/
inline sequence sequence::get_subsequence_chr_to(chromosomePosition to, const sequenceLength & len) const {
return get_subsequence_to(chromosome_position_to_sequence_position(to),len);
}

/*****************************************************************//**
* @brief output stream sequence insertion operator
*
* This operator is used to insert a sequence to an output stream (e.g.
* to print it on screen).
* The sequence will be represented by its nucleotides (Adenine: A,
* Cytosine: C, Guanine: G, Uracil: U, Masked: X) from 5' to 3'.
*
* @param the_stream output stream the sequence should be inserted in
* @param the_sequence sequence to be inserted in the output stream
*
* @return output stream with the inserted sequence
*********************************************************************/
std::ostream & operator<<(std::ostream & the_stream, const sequence & the_sequence);

/*****************************************************************//**
* @brief output stream strand insertion operator
*
* This operator is used to insert a strand to an output stream (e.g.
* to print it on screen).
* The strand will be represented by its sign (+ or -).
*
* @param the_stream output stream the strand should be inserted in
* @param the_strand strandType to be inserted in the output stream
*
* @return output stream with the inserted strand
*********************************************************************/
std::ostream & operator<<(std::ostream & the_stream, const strandType & the_strand);

} // namespace microSNPscore
#endif
