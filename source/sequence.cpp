
#include <iostream>
//for std::cerr and std::endl (error stating)
#include <sstream>
//for std::istringstream (type conversion)
#include <algorithm>
//for std::sort (exon sorting)
#include "sequence.h"

namespace microSNPscore {

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
    
    exon::exon(sequencePosition start_position, sequencePosition end_position)
    :start(start_position),end(end_position) {
       /*********************************************************\ 
      | Check for negative length and if so set to zero length by |
      | changing end position (after reporting the error):        |
       \*********************************************************/
      if (end_position<start_position)
      {
        std::cerr << "microSNPscore::exon::exon\n";
        std::cerr << " ==> negative length exon range: ";
        std::cerr << start_position << "-" << end_position << std::endl;
        std::cerr << "  --> setting to zero-length: ";
        std::cerr << start_position << "-" << start_position-1 << std::endl;
        end=start_position-1;
      }
}

    /*****************************************************************//**
    * @brief length calculation
    *
    * This method is used to calculate the length of the exon.
    *
    * @return  the exon's length
    *********************************************************************/
    sequenceLength exon::get_length() const {
       /******************************************************\ 
      | Calculate the length (start end end position included) |
       \******************************************************/
      return (get_end()-get_start()+1);
}

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
    * @return a sequence containing the given nucleotides located on the
    *     given chromosome, strand and positions.
    *********************************************************************/
    sequence::sequence(std::string sequence_string, const chromosomeType & the_chromosome, strandType the_strand, std::string exon_starts, std::string exon_ends)
    :chromosome(the_chromosome),strand(the_strand),exons(initialize_exons(position_string_to_vector(exon_starts),position_string_to_vector(exon_ends)))
    ,length(initialize_length(exons)),nucleotides(initialize_nucleotides(sequence_string,the_chromosome,the_strand,exons,length)) {
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
sequence sequence::get_subsequence_from(sequencePosition from, const sequenceLength & len) const {
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
sequence sequence::get_subsequence_to(sequencePosition to, const sequenceLength & len) const {
}

/*****************************************************************//**
* @brief get subsequence between sequence positions
*
* This method can be used to extract a subsequence starting (i.e. 5'
* end) end ending (i.e. 3' end) at given positions from the sequence.
* If at least one of the given positions is not part of the sequence,
* an empty sequence will be returned.
*
* @param from the start position in the sequence of the subsequence
* @param to the end position in the sequence of the subsequence
* @return the subsequence starting and ending at the given positions
*********************************************************************/
sequence sequence::get_subsequence_from_to(sequencePosition from, sequencePosition to) const {
}

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
sequence sequence::get_subsequence_chr_from(chromosomePosition from, const sequenceLength & len) const {
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
sequence sequence::get_subsequence_chr_to(chromosomePosition to, const sequenceLength & len) const {
}

/*****************************************************************//**
* @brief get subsequence between chromosome positions
*
* This method can be used to extract a subsequence starting (i.e. 5'
* end) end ending (i.e. 3' end) at given chromosome positions from the
* sequence.
* If at least one of the given chromosome positions is not part of the
* sequence, an empty sequence will be returned.
*
* @param from the start position on the chromosome of the subsequence
* @param to the end position on the chromosome of the subsequence
* @return the subsequence starting and ending at the given positions
*********************************************************************/
sequence sequence::get_subsequence_chr_from_to(chromosomePosition from, chromosomePosition to) const {
}

    /*****************************************************************//**
    * @brief standard constructor
    *
    * This is used to create an instance of the class sequence.
    * It is just created as dummy for mRNA and miRNA.
    *
    * @return an empty sequence
    *
    * @todo delete when mRNA and miRNA constructors are done.
    *********************************************************************/
    sequence::sequence()
    :strand(Plus),chromosome(""),length(0) {
}

    /*****************************************************************//**
    * @brief exon initialisation
    *
    * This method is used to create the exon vector of a sequence.
    * The ordering of exon starts and ends does not matter.
    * Overlapping exons are be merged (reporting an error).
    * If the count of exon starts does not match the count of exon ends
    * an error is raised and the additional starts or ends are omitted.
    *
    * @param starts: String representing the start positions (i.e.
    *     the end with the smaller distance to the chromosome start beeing
    *     the 5' end of the + strand and accordingly the 3' end of
    *     the - strand) of the exons containing the sequence as
    *     comma-separated list.
    * @param ends: String representing the end positions (i.e.
    *     the end with the smaller distance to the chromosome end beeing
    *     the 3' end of the + strand and accordingly the 5' end of
    *     the - strand) of the exons containing the sequence as
    *     comma-separated list.
    *
    * @return a vector containing exons with the given coordinates
    *********************************************************************/
    std::vector<exon> sequence::initialize_exons(const std::vector<chromosomePosition> & starts, const std::vector<chromosomePosition> & ends)
    {
       /***************************************************************\ 
      | Pair each start with the next unpaired end (skipping those less |
      | than the start) to build the first exon vector and ther iterate |
      | over that vector merging overlapping exons to build the final   |
      | exon vector:                                                    |
       \***************************************************************/
      std::vector<chromosomePosition>::const_iterator start_it(starts.begin());
      std::vector<chromosomePosition>::const_iterator end_it(ends.begin());
      std::vector<exon> unmerged_exon_vector;
      while(start_it != starts.end() && end_it != ends.end())
      {
        while(*end_it < *start_it && end_it != ends.end()) ++end_it;
        unmerged_exon_vector.push_back(exon(*start_it,*end_it));
        ++start_it;
        ++end_it;
      }
      std::vector<exon> exon_vector;
      for(std::vector<exon>::const_iterator exon_it(unmerged_exon_vector.begin());exon_it != unmerged_exon_vector.end();)
      {
        chromosomePosition start(exon_it->get_start());
        chromosomePosition end(exon_it->get_end());
        ++exon_it;
        while (exon_it->get_start() <= end && exon_it != unmerged_exon_vector.end())
        {
          end=exon_it->get_end();
          ++exon_it;
        }
        exon_vector.push_back(exon(start,end));
      }
      return exon_vector;
}

    /*****************************************************************//**
    * @brief length calculation
    *
    * This method is used to calculate the length of a sequence.
    *
    * @param exon_vector const std::vector<exon> reference to a vector
    *     containing the sequence's exons
    *
    * @return the sequence's length
    *********************************************************************/
    sequenceLength sequence::initialize_length(const std::vector<exon> & exon_vector)
    {
       /**********************************************\ 
      | Iterate over the exons summing up the lengths: |
       \**********************************************/ 
      sequenceLength the_length;
      for(std::vector<exon>::const_iterator exon_it(exon_vector.begin());exon_it != exon_vector.end();++exon_it)
      {
        the_length += exon_it->get_length();
      }
      return the_length;
}

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
    * @param the_exons a vector containing the sequence's exons
    * @param the_length the requested length of the sequence
    *
    * @return a vector containing the sequence's nucleotides
    *********************************************************************/
    std::vector<nucleotide> sequence::initialize_nucleotides(const std::string & the_sequence, chromosomeType the_chromosome, strandType the_strand, const std::vector<exon> & the_exons, sequenceLength the_length)
    {
       /**************************************************************\ 
      | Initialize empty nucleotide vector, counters and iterators and |
      | loop up to requested length:                                   |
       \**************************************************************/
      std::vector<nucleotide> nucleotide_vector;
      std::string::const_iterator sequence_it(the_strand == Plus ?
                                              the_sequence.begin() :
                                              the_sequence.end());
      std::vector<exon>::const_iterator exon_it(the_exons.begin());
      chromosomePosition position_on_chromosome(exon_it != the_exons.end() ? 
                                                exon_it->get_start() :
                                                0);
      sequenceLength length_of_sequence(0);
      sequenceLength length_of_exon(0);
      char the_base_char(the_sequence.begin() != the_sequence.end() ?
                         *sequence_it :
                         '\0');
      while(length_of_sequence != the_length &&
           (length_of_exon != exon_it->get_length() || exon_it != the_exons.end()))
      {
         /**************************************************************\ 
        | For every nucleobase that is no gap: add nucleotide, increment |
        | counters and if needed move on to next exon:                   |
         \**************************************************************/
        if(the_base_char != '-')
        {
          nucleoBase nucleo_base(Adenine);
          switch(the_base_char)
          {
            case 'a':
            case 'A':
              break;
            case 't':
            case 'T':
              std::cerr << "microSNPscore::sequence::initialize_nucleotides\n";
              std::cerr << " ==> illegal nucleo base character: \n";
              std::cerr << the_base_char << std::endl;
              std::cerr << "  --> assuming Uracil\n";
            case 'u':
            case 'U':
              nucleo_base=Uracil;
              break;
            case 'c':
            case 'C':
              nucleo_base=Cytosine;
              break;
            case 'g':
            case 'G':
              nucleo_base=Guanine;
              break;
            case '\0':
              std::cerr << "microSNPscore::sequence::initialize_nucleotides\n";
              std::cerr << " ==> missing nucleo base character\n";
              std::cerr << "  --> assuming Mask\n";
            case 'x':
            case 'X':
              nucleo_base=Mask;
              break;
            default:
              std::cerr << "microSNPscore::sequence::initialize_nucleotides\n";
              std::cerr << " ==> illegal nucleo base character: \n";
              std::cerr << the_base_char << std::endl;
              std::cerr << "  --> assuming Mask\n";
              nucleo_base=Mask;
          } // switch(the_base_char)
          nucleotide_vector.push_back(nucleotide(nucleo_base,length_of_sequence++,position_on_chromosome));
          if(position_on_chromosome != exon_it->get_end())
          {
            ++position_on_chromosome;
            ++length_of_exon;
          }
          else
          {
            ++exon_it;
            position_on_chromosome = exon_it->get_start();
            length_of_exon = 0;
          }
        }  // if(the_base_char != '-')
        else
        {
          std::cerr << "microSNPscore::sequence::initialize_nucleotides\n";
          std::cerr << " ==> illegal nucleo base character: \n";
          std::cerr << the_base_char << std::endl;
          std::cerr << "  --> assuming Gap --> omitting\n";
        }
         /**************************************************************\ 
        | Move on in given sequence and a append zero-chars if the given |
        | sequece is too short:                                          |
         \**************************************************************/
        if(the_base_char != '\0')
        {
          ++sequence_it;
        }
        the_base_char = (((sequence_it != the_sequence.begin() || the_strand == Plus) && (sequence_it != the_sequence.end() || the_strand == Minus)) ?
                         *sequence_it :
                         '\0');
      } // while-loop
       /*****************************************************************\ 
      | Check for additional characters in given sequence that have to be |
      | omitted and return constructed nucleotide vector:                 |
       \*****************************************************************/
      if((the_strand == Plus && sequence_it != the_sequence.end()) ||
         (the_strand == Minus && sequence_it != the_sequence.begin()))
      {
            std::cerr << "microSNPscore::sequence::initialize_nucleotides\n";
            std::cerr << " ==> additional nucleo base characters: \n";
            std::cerr << (the_strand == Plus ?
                          the_sequence.substr(sequence_it - the_sequence.begin(),the_sequence.end() - sequence_it +1) :
                          the_sequence.substr(0,sequence_it - the_sequence.begin() + 1)) << std::endl;
            std::cerr << "  --> omitting\n";
      }
      return nucleotide_vector;
}

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
    
    std::vector<chromosomePosition> sequence::position_string_to_vector(std::string string_list)
    {
       /****************************************************************\ 
      | Put the string in a stream, split at comma, put each part in its |
      | own stream, try to read a chromosomePosition from that stream,   |
      | append to vector and sort the result:                            |
       \****************************************************************/
      std::vector<chromosomePosition> position_vector;
      char string_to_split_at=',';
      std::istringstream string_to_be_splitted(string_list);
      std::string splitted_string;
      while(std::getline(string_to_be_splitted,splitted_string,string_to_split_at))
      {
        std::istringstream string_to_be_converted(splitted_string);
        chromosomePosition converted_string;
        string_to_be_converted >> converted_string;
        if(! converted_string)
        {
          std::cerr << "microSNPscore::sequence::position_string_to_position_vector\n";
          std::cerr << " ==> illegal chromosome position: ";
          std::cerr << splitted_string << std::endl;
          std::cerr << "  --> omitting position\n";
        }
        else
        {
           position_vector.push_back(converted_string);
        }
      }
      std::sort(position_vector.begin(),position_vector.end());
      return position_vector;
}


} // namespace microSNPscore
