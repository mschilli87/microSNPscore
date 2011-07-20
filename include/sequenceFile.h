#ifndef MICROSNPSCORE_SEQUENCEFILE_H
#define MICROSNPSCORE_SEQUENCEFILE_H


#include <string>
#include "sequence.h"
#include <vector>

namespace microSNPscore {

/*****************************************************************//**
* @brief file path type
*
* This represents a file path.
*********************************************************************/
typedef std::string filePath;
class sequenceFileEntry {
  public:
    /*****************************************************************//**
    * @brief constructor
    *
    * This is used to create a sequenceFileEntry from a FASTA entry.
    * If the given string is not a valid FASTA entry the empty default
    * entry is created (see default value).
    * A valid FASTA entry starts with the > sign followed by the sequence
    * ID, a pipe (|) character, a comma-separated list of exon starts,
    * another pipe character, a comma-separated list of exon ends, another
    * pipe character, 1 or -1 representing the strand (+ or -,
    * respectively), a last pipe character, the chromosome name and
    * linebreak followed by the sequence which may contain newlines but
    * must end with a newline.
    * If the header line containes more than four pipe characters, the
    * additional one are taken as part of the sequence ID.
    * The default value is not intended to be used directly.
    * It only provided to allow array allocation but you will need
    * to assign a valid object created by giving those parameters a value
    * to actually use it. This is done by containers like std::vector and
    * the reason for providing those default values is to allow using
    * containers containing objects of this class.
    *
    * @param FASTA_entry (pseudo-optional) std::string beeing a valid
    *    FASTA entry - Defaults to ">|||1|\n\n"
    * @return sequenceFileEntry corresponding to the given FASTA entry
    *********************************************************************/
    
    sequenceFileEntry(std::string FASTA_entry);

    /*****************************************************************//**
    * @brief get method for ID attribute
    *
    * This method is used to access the ID of the sequence the entry is
    * for.
    *
    * @return the ID of the sequence the entry is for
    *
    * @see sequence::get_id()
    *********************************************************************/
    
    inline const sequenceID get_ID() const;

    /*****************************************************************//**
    * @brief get method for chromosome attribute
    *
    * This method is used to access the chromosome of the sequence
    * the entry is for.
    *
    * @return the chromosome of the sequence the entry is for
    *
    * @see sequence::get_chromosome()
    *********************************************************************/
    inline const chromosomeType get_chromosome() const;

    /*****************************************************************//**
    * @brief get method for strand attribute
    *
    * This method is used to access the strand of the sequence this entry
    * is for.
    *
    * @return the strand of the sequence the entry is for
    *
    * @see sequence::get_strand()
    *********************************************************************/
    inline const strandType get_strand() const;

    /*****************************************************************//**
    * @brief get method for exon starts attribute
    *
    * This method is used to access the starts of the exons containing the
    * sequence this entry is for.
    *
    * @return std::string representing a comma-separated list
    * containing the starts of the exons containing the sequence this
    * entry is for
    *********************************************************************/
    inline const std::string get_exon_starts() const;

    /*****************************************************************//**
    * @brief get method for exon ends attribute
    *
    * This method is used to access the ends of the exons containing the
    * sequence this entry is for.
    *
    * @return std::string representing a comma-separated list
    * containing the ends of the exons containing the sequence this
    * entry is for
    *********************************************************************/
    inline const std::string get_exon_ends() const;

    inline const std::string get_nucleotide_sequence() const;


  private:
    /*****************************************************************//**
    * @brief sequence ID
    *
    * This is the ID of the sequence this entry is for.
    *
    * @see sequence::ID
    *********************************************************************/
    sequenceID ID;

    /*****************************************************************//**
    * @brief chromosome
    *
    * This is the chromosome of the sequence this entry is for.
    *
    * @see sequence::chromosome
    *********************************************************************/
    chromosomeType chromosome;

    /*****************************************************************//**
    * @brief strand
    *
    * This is the strand of the sequence this entry is for.
    *
    * @see sequence::strand
    *********************************************************************/
    strandType strand;

    /*****************************************************************//**
    * @brief exon start positions
    *
    * This is a comma-separated list containing the starts of the exons
    * containing the sequence this entry is for.
    *********************************************************************/
    std::string exon_starts;

    /*****************************************************************//**
    * @brief exon end positions
    *
    * This is a comma-separated list containing the ends of the exons
    * containing the sequence this entry is for.
    *********************************************************************/
    std::string exon_ends;

    /*****************************************************************//**
    * @brief letter code nucleotide string
    *
    * This is a string containing the nucleotides of the sequence this
    * entry is for.
    *********************************************************************/
    
    std::string nucleotide_sequence;

};
    /*****************************************************************//**
    * @brief get method for ID attribute
    *
    * This method is used to access the ID of the sequence the entry is
    * for.
    *
    * @return the ID of the sequence the entry is for
    *
    * @see sequence::get_id()
    *********************************************************************/
    
    inline const sequenceID sequenceFileEntry::get_ID() const {
      return ID;
    }

    /*****************************************************************//**
    * @brief get method for chromosome attribute
    *
    * This method is used to access the chromosome of the sequence
    * the entry is for.
    *
    * @return the chromosome of the sequence the entry is for
    *
    * @see sequence::get_chromosome()
    *********************************************************************/
    inline const chromosomeType sequenceFileEntry::get_chromosome() const {
      return chromosome;
    }

    /*****************************************************************//**
    * @brief get method for strand attribute
    *
    * This method is used to access the strand of the sequence this entry
    * is for.
    *
    * @return the strand of the sequence the entry is for
    *
    * @see sequence::get_strand()
    *********************************************************************/
    inline const strandType sequenceFileEntry::get_strand() const {
      return strand;
    }

    /*****************************************************************//**
    * @brief get method for exon starts attribute
    *
    * This method is used to access the starts of the exons containing the
    * sequence this entry is for.
    *
    * @return std::string representing a comma-separated list
    * containing the starts of the exons containing the sequence this
    * entry is for
    *********************************************************************/
    inline const std::string sequenceFileEntry::get_exon_starts() const {
      return exon_starts;
    }

    /*****************************************************************//**
    * @brief get method for exon ends attribute
    *
    * This method is used to access the ends of the exons containing the
    * sequence this entry is for.
    *
    * @return std::string representing a comma-separated list
    * containing the ends of the exons containing the sequence this
    * entry is for
    *********************************************************************/
    inline const std::string sequenceFileEntry::get_exon_ends() const {
      return exon_ends;
    }

    inline const std::string sequenceFileEntry::get_nucleotide_sequence() const {
      return nucleotide_sequence;
    }

class sequenceFile {
  private:
    /*****************************************************************//**
    * @brief file path
    *
    * This is the file path of the sequence file.
    *********************************************************************/
    filePath path;

    /*****************************************************************//**
    * @brief sequence entries
    *
    * This is a vector containing the sequence entries corresponding to
    * the FASTA entries in the sequence file.
    *********************************************************************/
    std::vector<sequenceFileEntry> entries;

};

} // namespace microSNPscore
#endif
