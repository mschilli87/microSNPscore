
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
    
    exon::exon(const chromosomePosition & start_position, const chromosomePosition & end_position)
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
    * @return a sequence containing the given nucleotides located on the
    *     given chromosome, strand and positions.
    *********************************************************************/
    sequence::sequence(sequenceID the_ID, std::string sequence_string, chromosomeType the_chromosome, strandType the_strand, std::string exon_starts, std::string exon_ends)
    :ID(the_ID),chromosome(the_chromosome),strand(the_strand),exons(initialize_exons(position_string_to_vector(exon_starts),position_string_to_vector(exon_ends)))
    ,length(initialize_length(exons.begin(),exons.end())),nucleotides(initialize_nucleotides(sequence_string,the_chromosome,the_strand,exons.begin(),exons.end(),length)) {
}

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
sequence sequence::get_subsequence_from_to(sequencePosition from, sequencePosition to) const {
 /*************************************************************\ 
| Only extract non-empty subsequences from non-empty sequences: |
 \*************************************************************/
std::vector<exon> exon_vector;
sequenceLength sequence_length(0);
std::vector<nucleotide> nucleotide_vector;
if(get_length()!=0 && from<to)
{
   /****************************************************************\ 
  | Check borders and move them to the respective ends if not valid: |
   \****************************************************************/
  if(from>get_length() || from<1)
  {
    from=1;
  }
  if(to>get_length() || to<1)
  {
    to=get_length();
  }
   /******************************************************************\ 
  | Iterate forward over + strand sequences or backward over - strand  |
  | sequences looking for non-successive chromosome positions to build |
  | up the exon vector:                                                |
   \******************************************************************/
  chromosomePosition exon_start(get_nucleotide(get_strand()==Plus ? from : to)->get_chromosome_position());
  chromosomePosition exon_end(exon_start);
  for(const_iterator nucleotide_it(get_strand()==Plus ?
                                   get_nucleotide(from+1) :
                                   get_nucleotide(to-1));get_strand()==Plus ?
                                                         nucleotide_it<=get_nucleotide(to) :
                                                         nucleotide_it>=get_nucleotide(from);nucleotide_it += get_strand()==Plus ?
                                                                                                              1 :
                                                                                                              -1)
  {
    if(nucleotide_it->get_chromosome_position()!=++exon_end)
    {
      exon_vector.push_back(exon(exon_start,exon_end-1));
      exon_start=nucleotide_it->get_chromosome_position();
      exon_end=exon_start;
    }
  }
   /***************************************************************\ 
  | Iterate forward over sequence to build up the nucleotide vector |
  | while counting the subsequence's length:                        |
   \***************************************************************/
  exon_vector.push_back(exon(exon_start,exon_end));
  for(const_iterator nucleotide_it(get_nucleotide(from));nucleotide_it<=get_nucleotide(to);++nucleotide_it)
  {
    nucleotide_vector.push_back(nucleotide(nucleotide_it->get_base(),++sequence_length,nucleotide_it->get_chromosome_position()));
  }
   /*************************************************************\ 
  | Return a sequence on the same chromosome and strand, with the |
  | calculated exon and nucleotide vectors and length:            |
   \*************************************************************/
}  // if(get_length()!=0 && from<to)
return sequence(get_ID(),get_chromosome(),get_strand(),exon_vector,sequence_length,nucleotide_vector);
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
return get_subsequence_from(chromosome_position_to_sequence_position(from),len);
}

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
sequence sequence::get_subsequence_chr_from_to(chromosomePosition from, chromosomePosition to) const {
 /******************************************************************\ 
| Invert borders for - strand sequences and map chromosome positions |
| to sequence positions:                                             |
 \******************************************************************/
if(get_strand()==Minus)
{
  chromosomePosition tmp(from);
  from=to;
  to=tmp;
}
return get_subsequence_from_to(chromosome_position_to_sequence_position(from),
                               chromosome_position_to_sequence_position(to));
}

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
    * @return a sequence with the given attributes
    *********************************************************************/
    sequence::sequence(sequenceID the_ID, chromosomeType the_chromosome, strandType the_strand, std::vector<exon> the_exons, sequenceLength the_length, const std::vector<nucleotide> & the_nucleotides)
    :ID(the_ID),chromosome(the_chromosome),strand(the_strand),exons(the_exons),length(the_length),nucleotides(the_nucleotides) {
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
    * @param begin_of_exons const_exon_iterator pointing the sequence's
    *     first exon
    * @param end_of_exons const_exon_iterator pointing behind the
    *     sequence's exon vector
    *
    * @return the sequence's length
    *********************************************************************/
    
    sequenceLength sequence::initialize_length(const sequence::const_exon_iterator & begin_of_exons, const sequence::const_exon_iterator & end_of_exons)
    {
       /**********************************************\ 
      | Iterate over the exons summing up the lengths: |
       \**********************************************/ 
      sequenceLength the_length(0);
      for(const_exon_iterator exon_it(begin_of_exons);exon_it != end_of_exons;++exon_it)
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
    * @param begin_of_exons const_exon_iterator pointing the sequence's
    *     first exon
    * @param end_of_exons const_exon_iterator pointing behind the
    *     sequence's exon vector
    * @param the_length the requested length of the sequence
    *
    * @return a vector containing the sequence's nucleotides
    *********************************************************************/
    std::vector<nucleotide> sequence::initialize_nucleotides(const std::string & the_sequence, chromosomeType the_chromosome, strandType the_strand
    , const sequence::const_exon_iterator & begin_of_exons, const sequence::const_exon_iterator & end_of_exons, sequenceLength the_length)
    {
       /**************************************************************\ 
      | Initialize empty nucleotide vector, counters and iterators and |
      | loop up to requested length where the order the exons are      |
      | iterated in depends on the strand (+: forward / -: backward):  |
       \**************************************************************/
      std::vector<nucleotide> nucleotide_vector;
      std::string::const_iterator sequence_it(the_sequence.begin());
      const_exon_iterator exon_it(the_strand == Plus ?
                                  begin_of_exons :
                                  end_of_exons - (begin_of_exons != end_of_exons ?
                                                  1 :
                                                  0));
      chromosomePosition position_on_chromosome(exon_it != end_of_exons ? 
                                                (the_strand == Plus ?
                                                 exon_it->get_start() :
                                                 exon_it->get_end()) :
                                                0);
      sequenceLength length_of_sequence(0);
      char the_base_char(sequence_it != the_sequence.end() ?
                         *sequence_it :
                         '\0');
      while(length_of_sequence != the_length &&
           ((the_strand == Plus && (position_on_chromosome <= exon_it->get_end() || exon_it != end_of_exons)) ||
            (the_strand == Minus && (position_on_chromosome >= exon_it->get_start() || exon_it != begin_of_exons))))
      {
         /***************************************************\ 
        | Add nucleotide for every nucleobase that is no gap: |
         \***************************************************/
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
          nucleotide_vector.push_back(nucleotide(nucleo_base,++length_of_sequence,position_on_chromosome));
           /***************************************************************\ 
          | Increment counters and if needed move on to next exon where the |
          | direction depends on the strand (+: forward / -:backward):      |
           \***************************************************************/
          if(the_strand == Plus)
          {
            if(position_on_chromosome != exon_it->get_end())
            {
              ++position_on_chromosome;
            }
            else
            {
              ++exon_it;
              position_on_chromosome = exon_it->get_start();
            }
          }
          else
          {
            if(position_on_chromosome != exon_it->get_start())
            {
              --position_on_chromosome;
            }
            else
            {
              --exon_it;
              position_on_chromosome = exon_it->get_end();
            }
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
        if(sequence_it != the_sequence.end()) // bases remaining?
        {
          ++sequence_it;
        }
        the_base_char = (sequence_it != the_sequence.end() ? *sequence_it : '\0');
      } // while-loop
       /*****************************************************************\ 
      | Check for additional characters in given sequence that have to be |
      | omitted and return constructed nucleotide vector:                 |
       \*****************************************************************/
      if(sequence_it != the_sequence.end())
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
    sequencePosition sequence::chromosome_position_to_sequence_position(chromosomePosition chromosome_position) const {
       /**************************************************************\ 
      | Iterating over the exons, finding the one containing the given |
      | position while summing up the lengths and calculate the        |
      | corresponding position in the found exon.                      |
      | On + stranded sequences that's the predessecor's position, on  |
      | - stranded sequences that's the distance to the last position: |
       \**************************************************************/
      sequenceLength prefix_length(0);
      for(const_exon_iterator exon_it(exons_begin());exon_it!=exons_end();++exon_it)
      {
        if(exon_it->get_start()<=chromosome_position)
        {
          if(exon_it->get_end()>=chromosome_position)
          {
            prefix_length += chromosome_position - exon_it->get_start();
          }
          else
          {
            prefix_length += exon_it->get_length();
          }
        }
        else
        {
          break;
        }
      }
      return get_strand() == Plus ?
             prefix_length + 1 :
             get_length() - prefix_length;
}

/*****************************************************************//**
* @brief sequence output stream insertion operator
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
std::ostream & operator<<(std::ostream & the_stream, const sequence & the_sequence)
{
   /*********************************************************\ 
  | Copy the stream, iterate over the sequence inserting each |
  | nucleotide and return the result:                         |
   \*********************************************************/
  for(sequence::const_iterator sequence_it(the_sequence.begin());sequence_it!=the_sequence.end();++sequence_it)
  {
    the_stream << *sequence_it;
  }
  return the_stream;
}

} // namespace microSNPscore
