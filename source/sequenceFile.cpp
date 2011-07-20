
#include "sequenceFile.h"

namespace microSNPscore {

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
    
    sequenceFileEntry::sequenceFileEntry(std::string FASTA_entry) {
}


} // namespace microSNPscore
