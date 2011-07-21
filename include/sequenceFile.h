#ifndef MICROSNPSCORE_SEQUENCEFILE_H
#define MICROSNPSCORE_SEQUENCEFILE_H


#include <iostream>
//for std::ostream (operator<<)
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
/*****************************************************************//**
* @brief sequence file entry class
*
* This represents a sequence file entry and is convertable to a
* sequence object and to a FASTA entry making it possible to read from
* and write to the file.
*
* @see sequenceFile
* @see sequence
*********************************************************************/

class sequenceFileEntry {
  public:
    /*****************************************************************//**
    * @brief constructor from FASTA entry
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
    * @brief constructor from sequence object
    *
    * This is used to create a sequenceFileEntry from a sequence object.
    *
    * @param the_sequence const sequence reference to the sequence the
    *     entry should be created for
    * @return sequenceFileEntry corresponding to the sequence
    *********************************************************************/
    sequenceFileEntry(const sequence & the_sequence);

    /*****************************************************************//**
    * @brief sequence object creation
    *
    * This method is used to create a sequence object corresponding to the
    * sequence file entry.
    *
    * @return sequence object corresponding to the sequence file entry
    *********************************************************************/
    inline sequence get_sequence() const;

    /*****************************************************************//**
    * @brief FASTA entry creation
    *
    * This method is used to create a FASTA entry corresponding to the
    * sequence file entry.
    *
    * @return string containing a FASTA entry corresponding to the
    *     sequence file entry
    *********************************************************************/
    std::string get_FASTA() const;


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
    * @brief sequence object creation
    *
    * This method is used to create a sequence object corresponding to the
    * sequence file entry.
    *
    * @return sequence object corresponding to the sequence file entry
    *********************************************************************/
    inline sequence sequenceFileEntry::get_sequence() const {
      return sequence(ID,nucleotide_sequence,chromosome,strand,exon_starts,exon_ends);
}

/*****************************************************************//**
* @brief sequence file class
*
* This represents a sequence file containing FASTA entries. It is
* capable of reading from and writing to the file.
*
* @see sequenceFileEntry
*********************************************************************/

class sequenceFile {
  public:
    /*****************************************************************//**
    * @brief const entry iterator
    *
    * This is used to iterate over sequence file's entries
    *********************************************************************/
    typedef std::vector<sequenceFileEntry>::const_iterator const_iterator;

    /*****************************************************************//**
    * @brief constructor
    *
    * This is used to create a sequenceFile with a given path.
    * The path cannot be changed after the creation.
    * Even if the given file path exists and contains valid FASTA
    * entries the sequence file won't contain any entries until you call
    * the @p read method.
    * If the file path is not valid the sequence file will be created but
    * useless since @p read won't be able to read from the file and
    * @p write won't be able to write to the file.
    *
    * @param the_path filePath to the sequence file
    *
    * @return sequenceFile corresponding to the given file path
    *
    * @see read()
    * @see write()
    *********************************************************************/
    
    sequenceFile();

    /*****************************************************************//**
    * @brief read entries from existing file
    *
    * This method is used to read valid FASTA entries from an existing
    * file to the sequence file's entries.
    * If the file does not exist or does not contain any valid FASTA
    * entires an error is raised and no sequence file entry is created.
    * A valid FASTA entry starts with the > sign followed by the sequence
    * ID, a pipe (|) character, a comma-separated list of exon starts,
    * another pipe character, a comma-separated list of exon ends, another
    * pipe character, 1 or -1 representing the strand (+ or -,
    * respectively), a last pipe character, the chromosome name and
    * linebreak followed by the sequence which may contain newlines but
    * must end with a newline.
    * If the header line containes more than four pipe characters, the
    * additional one are taken as part of the sequence ID.
    *********************************************************************/
    void read();

    /*****************************************************************//**
    * @brief entry vector begin
    *
    * This is used to get the first entry of the sequence file.
    *
    * @return const_iterator pointing to the sequence file's first entry
    *********************************************************************/
    inline const_iterator begin() const;

    /*****************************************************************//**
    * @brief entry vector end
    *
    * This is used to get the end sequence file's entry vector.
    *
    * @return const_iterator pointing behind the last entry of the file
    *********************************************************************/
    inline const_iterator end() const;

    /*****************************************************************//**
    * @brief sequence entry creation
    *
    * This method is used to create an entry for a given sequence in the
    * sequence file.
    * The file won't actually change until you call the @p write method.
    *
    * @param the_sequence const sequence reference to the sequence an
    *     entry should be made for
    *
    * @see write()
    *********************************************************************/
    inline void add_sequence(const sequence & the_sequence);

    /*****************************************************************//**
    * @brief write entries to file
    *
    * This method is used to write FASTA entries from the sequence file's
    * entries to the file.
    * If the file path is not valid (or not accessable), an error is
    * raised.
    * If the path is valid (and accessable) and exists its current content
    * will be overwritten with the current sequence file content.
    * If the path is valid (and accessable) and does not exist it is
    * created and the current sequence file content is written to it.
    * To write a given sequence to the sequence file you first have to add
    * it to the sequnece file by calling the @p add_sequence method.
    *
    * @see add_sequence()
    *********************************************************************/
    
    void write();


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
    /*****************************************************************//**
    * @brief entry vector begin
    *
    * This is used to get the first entry of the sequence file.
    *
    * @return const_iterator pointing to the sequence file's first entry
    *********************************************************************/
    inline sequenceFile::const_iterator sequenceFile::begin() const {
      return entries.begin();
}

    /*****************************************************************//**
    * @brief entry vector end
    *
    * This is used to get the end sequence file's entry vector.
    *
    * @return const_iterator pointing behind the last entry of the file
    *********************************************************************/
    inline sequenceFile::const_iterator sequenceFile::end() const {
      return entries.end();
}

    /*****************************************************************//**
    * @brief sequence entry creation
    *
    * This method is used to create an entry for a given sequence in the
    * sequence file.
    * The file won't actually change until you call the @p write method.
    *
    * @param the_sequence const sequence reference to the sequence an
    *     entry should be made for
    *
    * @see write()
    *********************************************************************/
    inline void sequenceFile::add_sequence(const sequence & the_sequence) {
      entries.push_back(sequenceFileEntry(the_sequence));
}

/*****************************************************************//**
* @brief output stream sequence file entry insertion operator
*
* This operator is used to insert a sequence file entry to an output
* stream (e.g. to print it on screen).
* The sequenceFileEntry will be represented as FASTA entry.
*
* @param the_stream output stream the sequence file entry should be
*     inserted in
* @param the_entry sequenceFileEntry to be inserted in the output
*     stream
*
* @return output stream with the inserted sequence file entry
*********************************************************************/
inline std::ostream & operator<<(std::ostream & the_stream, const sequenceFileEntry & the_entry) {
  return the_stream << the_entry.get_FASTA();
}

} // namespace microSNPscore
#endif
