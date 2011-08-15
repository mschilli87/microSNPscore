
#include "mRNA.h"
#include "conservationList.h"
#include "SNP.h"

namespace microSNPscore {

    /*****************************************************************//**
    * @brief constructor
    *
    * This is used to create an instance of the class mRNA.
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
    * @param the_ID sequenceID representing the ID of the mRNA
    * @param sequence_string String representing the nucleotide sequence
    *     (Adenine: A, Cytosine: C, Guanine: G, Uracil: U, Mask: X)
    * @param the_chromosome chromosomeType representing the chromosome the
    *     mRNA is located on
    * @param the_strand strandType representing the strand (Plus/Minus) on
    *     which the mRNA is located
    * @param exon_starts: String representing the start positions (i.e.
    *     the end with the smaller distance to the chromosome start beeing
    *     the 5' end of the + strand and accordingly the 3' end of
    *     the - strand) of the exons containing the mRNA as
    *     comma-separated list.
    * @param exon_ends: String representing the end positions (i.e.
    *     the end with the smaller distance to the chromosome end beeing
    *     the 3' end of the + strand and accordingly the 5' end of
    *     the - strand) of the exons containing the mRNA as
    *     comma-separated list.
    * @param conservations conservationList containing the conservaton
    *     ranges for the sequence
    *
    * @return a mRNA containing the given nucleotides located on the
    *     given chromosome, strand and positions.
    *********************************************************************/
    mRNA::mRNA(const sequenceID the_ID, std::string sequence_string, chromosomeType the_chromosome, strandType the_strand, std::string exons_starts, std::string exon_ends, const conservationList & conservations)
    :sequence(the_ID,sequence_string,the_chromosome,the_strand,exons_starts,exon_ends,conservations) {
}

    /*****************************************************************//**
    * @brief standard constructor - do not use directly
    *
    * This is only provided to allow array allocation but you will need
    * to assign a valid object created by the parameterized constructor
    * to actually use it. This is done by containers like std::vector and
    * the reason for providing those default values is to allow using
    * containers containing objects of this class.
    *
    * @return: an unitialized mRNA object
    *********************************************************************/
    mRNA::mRNA()
    :sequence() {
}

    /*****************************************************************//**
    * @brief extract subsequence relevant for accessability calculations
    *
    * This method is used to query the subsequence from the whole mRNA
    * sequence ('the whole mRNA' actually means 'only' the 3'UTR which is
    * everything considered by microSNPscore at all) that is used to
    * calculate accessability scores (i.e. +/- 80 nucleotides around seed
    * match region end).
    *
    * @param predicted_miRNA_three_prime_position the position on the
    *     chromosome the target prediction algorithm has predicted to be
    *     the position the 3' end of the miRNA would be aligned to (if it
    *     would actually bind) (i.e. one base downstream from the seed 
    *     matching region)
    *
    * @return mRNA containing the subsequence relevant for accessability
    *    score calculation
    *********************************************************************/
    mRNA mRNA::get_subsequence_for_accessability(chromosomePosition predicted_miRNA_three_prime_position) const {
       /******************************************************************\ 
      | Calculating sequence position before subsequence querying to avoid |
      | caring about strand:                                               |
       \******************************************************************/
      sequenceLength window_size(80);
      sequencePosition sequence_position(chromosome_position_to_sequence_position(predicted_miRNA_three_prime_position));
      return get_subsequence_from_to(sequence_position - window_size,sequence_position + window_size);
}

    /*****************************************************************//**
    * @brief extract subsequence relevant for downstream AU content 
    *    calculation
    *
    * This method is used to query the subsequence from the whole mRNA
    * sequence ('the whole mRNA' actually means 'only' the 3'UTR which is
    * everything considered by microSNPscore at all) that is used to
    * calculate the downstream (i.e. +30 nucleotides behind seed
    * 8mer match region).
    *
    * @param predicted_miRNA_three_prime_position the position on the
    *     chromosome the target prediction algorithm has predicted to be
    *     the position the 3' end of the miRNA would be aligned to (if it
    *     would actually bind) (i.e. one base downstream from the seed 
    *     matching region)
    *
    * @return mRNA containing the subsequence relevant for downstream AU
    *    content calculation
    *********************************************************************/
    mRNA mRNA::get_subsequence_for_downstream_AU_content(chromosomePosition predicted_miRNA_three_prime_position) const {
       /******************************************************************\ 
      | Calculating sequence position before subsequence querying to avoid |
      | caring about strand:                                               |
       \******************************************************************/
      sequenceLength window_size(30);
      sequencePosition sequence_position(chromosome_position_to_sequence_position(predicted_miRNA_three_prime_position) + 1);
      return get_subsequence_from_to(sequence_position,sequence_position + window_size);
}

    /*****************************************************************//**
    * @brief extract subsequence relevant for upstream AU content 
    *    calculation
    *
    * This method is used to query the subsequence from the whole mRNA
    * sequence ('the whole mRNA' actually means 'only' the 3'UTR which is
    * everything considered by microSNPscore at all) that is used to
    * calculate the upstream (i.e. -30 nucleotides before seed
    * 8mer match region).
    *
    * @param predicted_miRNA_three_prime_position the position on the
    *     chromosome the target prediction algorithm has predicted to be
    *     the position the 3' end of the miRNA would be aligned to (if it
    *     would actually bind) (i.e. one base downstream from the seed 
    *     matching region)
    *
    * @return mRNA containing the subsequence relevant for downstream AU
    *    content calculation
    *********************************************************************/
    mRNA mRNA::get_subsequence_for_upstream_AU_content(chromosomePosition predicted_miRNA_three_prime_position) const {
       /******************************************************************\ 
      | Calculating sequence position before subsequence querying to avoid |
      | caring about strand:                                               |
       \******************************************************************/
      sequenceLength window_size(30);
      sequencePosition sequence_position(chromosome_position_to_sequence_position(predicted_miRNA_three_prime_position) - 8);
      return get_subsequence_from_to(sequence_position - window_size,sequence_position);
}

    /*****************************************************************//**
    * @brief internal constructor
    *
    * This method is used to convert a sequence to a mRNA.
    *
    * @param the_sequence const sequence reference to the sequence that
    *     should become an mRNA
    *
    * @return mRNA with the same attributes as the given sequence
    *********************************************************************/
    mRNA::mRNA(const sequence & the_sequence)
    :sequence(the_sequence) {
}


} // namespace microSNPscore
