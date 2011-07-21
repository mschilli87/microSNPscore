
#include<regex.h>
// for regex_t, regmatch_t, regcomp, regexec (regular expressions)
#include <sstream>
//for std::istringstream (FASTA string composition)
#include "sequenceFile.h"

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
    
    sequenceFileEntry::sequenceFileEntry(std::string FASTA_entry):
    ID(""),chromosome(""),strand(Plus),exon_starts(""),exon_ends(""),nucleotide_sequence("") {
       /******************************************************************\ 
      | Try to initialize extended regular expression matching valid FASTA |
      | entry stating error in case of failure:                            |
       \******************************************************************/
      regex_t regex_FASTA;
      char pattern_FASTA[] = "^>([^|\n]*)[|]([^\n]*)[|]([^|\n]*)[|](-?1)[|]([^|\n]*)\n(.*)$";
      if(regcomp(&regex_FASTA,pattern_FASTA,REG_EXTENDED) != 0)
      {
        std::cerr << "microSNPscore::sequenceFileEntry::sequenceFileEntry\n";
        std::cerr << " ==> compiling FASTA regular expression failed\n";
        std::cerr << "  --> creating empty default sequence\n";
      }
      else
      {
         /**************************************************************\ 
        | Try to match the initialized regular expression on given FASTA |
        | entry stating error in case of failure:                        |
         \**************************************************************/
        size_t nmatch_FASTA(0);
        regmatch_t pmatch_FASTA[7];
        const char * string_FASTA = FASTA_entry.c_str();
        if(regexec(&regex_FASTA,string_FASTA,nmatch_FASTA,pmatch_FASTA,0) != 0)
        {
          std::cerr << "microSNPscore::sequenceFileEntry::sequenceFileEntry\n";
          std::cerr << " ==> no valid FASTA entry:\n";
          std::cerr << FASTA_entry << std::endl;
          std::cerr << "  --> creating empty default sequence\n";
        }
        else
        {
           /*****************************************************************\ 
          | Extract the subsequences matching the regular expression's groups |
          | from the given FASTA entry assigning them to the corresponding    |
          | attributes:                                                       |
           \*****************************************************************/
          ID = FASTA_entry.substr(pmatch_FASTA[1].rm_so,pmatch_FASTA[1].rm_eo - pmatch_FASTA[1].rm_so);
          exon_starts = FASTA_entry.substr(pmatch_FASTA[2].rm_so,pmatch_FASTA[2].rm_eo - pmatch_FASTA[2].rm_so);
          exon_ends = FASTA_entry.substr(pmatch_FASTA[3].rm_so,pmatch_FASTA[3].rm_eo - pmatch_FASTA[3].rm_so);
          strand = FASTA_entry.substr(pmatch_FASTA[4].rm_so,pmatch_FASTA[4].rm_eo - pmatch_FASTA[4].rm_so) == "1" ? Plus : Minus;
          chromosome = FASTA_entry.substr(pmatch_FASTA[5].rm_so,pmatch_FASTA[5].rm_eo - pmatch_FASTA[5].rm_so);
          nucleotide_sequence = FASTA_entry.substr(pmatch_FASTA[6].rm_so,pmatch_FASTA[6].rm_eo - pmatch_FASTA[6].rm_so);
        } // regexec(&regex_FASTA,FASTA_entry,nmatch_FASTA,pmatch_FASTA_FASTA) == 0
      } // regcomp(&regex_FASTA,patter_FASTA,REG_EXTEND) == 0
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
    * @brief FASTA entry creation
    *
    * This method is used to create a FASTA entry corresponding to the
    * sequence file entry.
    *
    * @return string containing a FASTA entry corresponding to the
    *     sequence file entry
    *********************************************************************/
    std::string sequenceFileEntry::get_FASTA() const {
      std::istringstream sequence_stream(nucleotide_sequence);
      std::ostringstream FASTA_stream;
      FASTA_stream << ">" << ID << "|" << exon_starts << "|" << exon_ends << "|";
      FASTA_stream << (strand == Plus ? "1" : "-1") << "|" << chromosome << std::endl;
      char sequence_line[61];
      while(sequence_stream.read(sequence_line,60).good())
      {
        FASTA_stream << sequence_line;
      }
      return FASTA_stream.str();
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


} // namespace microSNPscore
