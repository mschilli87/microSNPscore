
#include "sequenceFile.h"
#include "sequence.h"

namespace microSNPscore {

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
    
    sequenceFileEntry::sequenceFileEntry(std::string FASTA_entry) {
}

    /*****************************************************************//**
    * @brief constructor from sequence object
    *
    * This is used to create a sequenceFileEntry from a sequence object.
    *
    * @param the_sequence const sequence reference to the sequence the
    *     entry should be created for
    * @return sequenceFileEntry corresponding to the sequence
    *********************************************************************/
    sequenceFileEntry::sequenceFileEntry(const sequence & the_sequence) {
}

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
    
    sequenceFile::sequenceFile() {
}

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
    void sequenceFile::read() {
}

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
    
    void sequenceFile::write() {
}

/*****************************************************************//**
* @brief output stream insertion operator
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
std::ostream & operator<<(std::ostream & the_stream, const sequenceFileEntry & the_entry)
{
}

} // namespace microSNPscore
