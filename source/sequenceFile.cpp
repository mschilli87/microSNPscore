
#include <regex.h>
// for regex_t, regmatch_t, regcomp and regexec (regular expressions)
#include <sstream>
//for std::istringstream (newline omitting & FASTA string composition)
#include <iostream>
//for std::cerr and std::endl (error stating)
#include <fstream>
//for std::ofstream and std::ifstream (file access)
#include "sequenceFile.h"
#include "conservationList.h"

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
      char pattern_FASTA[] = "^>([^\n]*)[|]([^|\n]*)[|]([^|\n]*)[|](-?1)[|]([^|\n]*)\n(.*)$";
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
        size_t nmatch_FASTA(7);
        regmatch_t pmatch_FASTA[nmatch_FASTA];
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
          | attributes (omitting newlines in the nucleotide sequence):        |
           \*****************************************************************/
          ID = FASTA_entry.substr(pmatch_FASTA[1].rm_so,pmatch_FASTA[1].rm_eo - pmatch_FASTA[1].rm_so);
          exon_starts = FASTA_entry.substr(pmatch_FASTA[2].rm_so,pmatch_FASTA[2].rm_eo - pmatch_FASTA[2].rm_so);
          exon_ends = FASTA_entry.substr(pmatch_FASTA[3].rm_so,pmatch_FASTA[3].rm_eo - pmatch_FASTA[3].rm_so);
          strand = FASTA_entry.substr(pmatch_FASTA[4].rm_so,pmatch_FASTA[4].rm_eo - pmatch_FASTA[4].rm_so) == "1" ? Plus : Minus;
          chromosome = FASTA_entry.substr(pmatch_FASTA[5].rm_so,pmatch_FASTA[5].rm_eo - pmatch_FASTA[5].rm_so);
          std::istringstream nucleotide_stream(FASTA_entry.substr(pmatch_FASTA[6].rm_so,pmatch_FASTA[6].rm_eo - pmatch_FASTA[6].rm_so));
          std::string nucleotide_line;
          while(std::getline(nucleotide_stream,nucleotide_line))
          {
            nucleotide_sequence += nucleotide_line;
          }
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
    sequenceFileEntry::sequenceFileEntry(const sequence & the_sequence)
    :ID(the_sequence.get_ID()),chromosome(the_sequence.get_chromosome()),strand(the_sequence.get_strand()),exon_starts(""),exon_ends(""),nucleotide_sequence("") {
       /****************************************************************\ 
      | If required create output string streams for the exon starts and |
      | end lists, iterate over the sequnece's exons inserting their     |
      | start and end to the stream, get the nucleotide sequence by      |
      | inserting the sequence to another output string stream and get   |
      | the final string values from the streams:                        |
       \****************************************************************/
      if(the_sequence.get_length() != 0)
      {
        std::ostringstream exon_starts_stream;
        std::ostringstream exon_ends_stream;
        exon_starts_stream << the_sequence.exons_begin()->get_start();
        exon_ends_stream << the_sequence.exons_begin()->get_end();
        for(sequence::const_exon_iterator exon_it(the_sequence.exons_begin()+1);exon_it!=the_sequence.exons_end();++exon_it)
        {
          exon_starts_stream << "," << exon_it->get_start();
          exon_ends_stream << "," << exon_it->get_end();
        }
        std::ostringstream nucleotide_stream;
        nucleotide_stream << the_sequence;
        exon_starts = exon_starts_stream.str();
        exon_ends = exon_ends_stream.str();
        nucleotide_sequence = nucleotide_stream.str();
      }
}

    /*****************************************************************//**
    * @brief FASTA entry creation
    *
    * This method is used to create a FASTA entry corresponding to the
    * sequence file entry.
    *
    * @param nucleotides_per_line (optional) sequenceLength of one line
    * in the FASTA output (a newline will be insterted after every that
    * number of nucleotides) - Defaults to 60
    * @return string containing a FASTA entry corresponding to the
    *     sequence file entry
    *********************************************************************/
    
    std::string sequenceFileEntry::get_FASTA(sequenceLength nucleotides_per_line) const {
       /******************************************************************\ 
      | Create an output string stream for the FASTA entry composition and |
      | insert the header informations into it. Then create an input       |
      | string stream from the nucleotide sequence to read lines of the    |
      | desired length from it inserting them (followed by a newline) into |
      | the output stream before converting it to a string for return:     |
       \******************************************************************/
      std::ostringstream FASTA_stream;
      FASTA_stream << ">" << ID << "|" << exon_starts << "|" << exon_ends << "|";
      FASTA_stream << (strand == Plus ? "1" : "-1") << "|" << chromosome << std::endl;
      std::istringstream sequence_stream(nucleotide_sequence);
      char sequence_line[nucleotides_per_line+2];
      sequence_line[nucleotides_per_line] = '\n';
      sequence_line[nucleotides_per_line+1] = '\0';
      while(sequence_stream.read(sequence_line,nucleotides_per_line).good())
      {
        FASTA_stream << sequence_line;
      }
      sequenceLength last_line_length(sequence_stream.gcount());
      if(last_line_length > 0)
      {
        sequence_line[last_line_length] = '\n';
        sequence_line[++last_line_length] = '\0';
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
    
    sequenceFile::sequenceFile(filePath the_path)
    :path(the_path),entries(std::vector<sequenceFileEntry>()) {
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
       /**********************************************************\ 
      | Try to open an input file stream associated to the file    |
      | corresponding to the sequence file's path stating an error |
      | in the case of failure:                                    |
       \**********************************************************/
      std::ifstream FASTA_file(path.c_str());
      if(FASTA_file.fail())
      {
        std::cerr << "microSNPscore::sequenceFile::sequenceFile\n";
        std::cerr << " ==> Cannot open file to read from: ";
        std::cerr << path << std::endl;
        std::cerr << "  --> no sequences will be read from the file\n";
      }
      else
      {
         /******************************************************************\ 
        | Try to initialize extended regular expression matching a valid     |
        | FASTA file stating error in case of failure:                       |
         \******************************************************************/
        regex_t regex_FASTA_file;
        char pattern_FASTA_file[] = "^(>[^\n]*[|][^|\n]*[|][^|\n]*[|]-?1[|][^|\n]*\n[^>]*)+$";
        if(regcomp(&regex_FASTA_file,pattern_FASTA_file,REG_EXTENDED|REG_NOSUB) != 0)
        {
          std::cerr << "microSNPscore::sequenceFile::read\n";
          std::cerr << " ==> compiling FASTA file regular expression failed\n";
          std::cerr << "  --> no sequences will be read from the file\n";
        }
        else
        {
           /****************************************************************\ 
          | Read the content of the file into string stream, convert it to a |
          | string and a c-type string and try to match the initialized      |
          | regular expression stating error in case of failure:             |
           \****************************************************************/
          std::stringstream FASTA_file_stream;
          FASTA_file_stream << FASTA_file.rdbuf();
          FASTA_file.close();
          std::string FASTA_file_string(FASTA_file_stream.str());
          const char * FASTA_file_cstring = FASTA_file_string.c_str();
          if(regexec(&regex_FASTA_file,FASTA_file_cstring,(size_t) 0,NULL,0) != 0)
          {
            std::cerr << "microSNPscore::sequenceFile::read\n";
            std::cerr << " ==> no valid FASTA file:\n";
            std::cerr << FASTA_file_string << std::endl;
            std::cerr << "  --> creating empty default sequence\n";
          }
          else
          {
             /******************************************************************\ 
            | Try to initialize extended regular expression matching a valid     |
            | FASTA entry stating error in case of failure:                      |
             \******************************************************************/
            regex_t regex_FASTA_entry;
            char pattern_FASTA_entry[] = "^>[^\n]*[|][^|\n]*[|][^|\n]*[|]-?1[|][^|\n]*\n[^>]*";
            if(regcomp(&regex_FASTA_entry,pattern_FASTA_entry,REG_EXTENDED) != 0)
            {
              std::cerr << "microSNPscore::sequenceFile::read\n";
              std::cerr << " ==> compiling FASTA entry regular expression failed\n";
              std::cerr << "  --> no sequences will be read from the file\n";
            }
            else
            {
               /*************************************************************\ 
              | Match the initialized regular expression as often as possible |
              | (at least one time because it matches valid FAST file regex)  |
              | inserting a sequence file entry for every match found:        |
               \*************************************************************/
              size_t nmatch_FASTA_entry(1);
              regmatch_t pmatch_FASTA_entry[nmatch_FASTA_entry];
              regoff_t start(0);
              regoff_t end(0);
              regexec(&regex_FASTA_entry,FASTA_file_cstring,nmatch_FASTA_entry,pmatch_FASTA_entry,0);
              do
              {
                start = end + pmatch_FASTA_entry[0].rm_so;
                end += pmatch_FASTA_entry[0].rm_eo;
                entries.push_back(sequenceFileEntry(FASTA_file_string.substr(start,end-start)));
              }
              while(regexec(&regex_FASTA_entry,FASTA_file_cstring+end,nmatch_FASTA_entry,pmatch_FASTA_entry,0) == 0);
            } // regcomp(&regex_FASTA_entry,pattern_FASTA_entry,REG_EXTENDED) == 0
          } // regexec(&regex_FASTA_file,FASTA_file_cstring,(size_t) 0,NULL,0) == 0
        } // regcomp(&regex_FASTA_file,pattern_FASTA_file,REG_EXTENDED|REG_NOSUB) == 0
      } // ! FASTA_file.fail()
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
       /*****************************************************************\ 
      | Try to open an outut file stream associated to the file           |
      | corresponding to the sequence file's path stating an error in the |
      | case of failure and iterate over the entries inserting each one   |
      | into the output stream before closing the output file stream:     |
       \*****************************************************************/
      std::ofstream the_file(path.c_str());
      if(the_file.fail())
      {
        std::cerr << "microSNPscore::sequenceFile::sequenceFile\n";
        std::cerr << " ==> Cannot open file to write to: ";
        std::cerr << path << std::endl;
        std::cerr << "  --> sequences won't be written to the file\n";
      }
      else
      {
        for(const_iterator entry_it(begin());entry_it!=end();++entry_it)
        {
          the_file << *entry_it;
        }
        the_file.close();
      }
}


} // namespace microSNPscore
