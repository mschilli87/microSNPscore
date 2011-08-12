#ifndef MICROSNPSCORE_MRNA_H
#define MICROSNPSCORE_MRNA_H


#include "sequence.h"
#include <string>
#include "nucleotide.h"

namespace microSNPscore { class conservationList; } 
namespace microSNPscore { class SNP; } 

namespace microSNPscore {

/*****************************************************************//**
* @brief messenger RNA class
*
* This represents a messenger RNA sequence.
* It is a specialisation of the sequence class.
*********************************************************************/
class mRNA : public sequence {
  public:
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
    mRNA(sequenceID the_ID, std::string sequence_string, const chromosomeType & the_chromosome, strandType the_strand, std::string exons_starts, std::string exon_ends, const conservationList & conservations);

    /*****************************************************************//**
    * @brief extract subsequence relevant for mRNA:miRNA alignment
    *
    * This method is used to query the subsequence from the whole mRNA
    * sequence ('the whole mRNA' actually means 'only' the 3'UTR which is
    * everything considered by microSNPscore at all) that is expected to
    * have any relevance in aligning the miRNA to the mRNA.
    *
    * @param predicted_miRNA_three_prime_position the position on the
    *     chromosome the target prediction algorithm has predicted to be
    *     the position the 3' end of the miRNA would be aligned to (if it
    *     would actually bind) (i.e. one base downstream from the seed 
    *     matching region)
    * @param len the length of the part that is expected to be relevant
    *    for the alignment - Defaults to 30 (since miRNAs are about 21-22
    *    nucleotides long and we assume that there would be no big loops
    *    in the optimal alignment)
    *
    * @return mRNA containing the subsequence relevant for the alignment
    *********************************************************************/
    
    inline mRNA get_subsequence_for_alignment(chromosomePosition predicted_miRNA_three_prime_position, sequenceLength len = 30) const;

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
    mRNA get_subsequence_for_accessability(chromosomePosition predicted_miRNA_three_prime_position) const;

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
    mRNA get_subsequence_for_downstream_AU_content(chromosomePosition predicted_miRNA_three_prime_position) const;

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
    mRNA get_subsequence_for_upstream_AU_content(chromosomePosition predicted_miRNA_three_prime_position) const;

    /*****************************************************************//**
    * @brief apply SNP on mRNA
    *
    * This method is used to change a mRNA as described by a given SNP.
    * If the SNP information do not match those in the mRNA an unchanged
    * copy of the mRNA is returned (after stating an error).
    *
    * @param the_SNP SNP containing the information how the mRNA should be
    *     changed
    *
    * @return a copy of the mRNA with the changes defined by the SNP
    *********************************************************************/
    inline mRNA mutate(const SNP & the_SNP) const;


  private:
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
    mRNA(const sequence & the_sequence);

};
    /*****************************************************************//**
    * @brief extract subsequence relevant for mRNA:miRNA alignment
    *
    * This method is used to query the subsequence from the whole mRNA
    * sequence ('the whole mRNA' actually means 'only' the 3'UTR which is
    * everything considered by microSNPscore at all) that is expected to
    * have any relevance in aligning the miRNA to the mRNA.
    *
    * @param predicted_miRNA_three_prime_position the position on the
    *     chromosome the target prediction algorithm has predicted to be
    *     the position the 3' end of the miRNA would be aligned to (if it
    *     would actually bind) (i.e. one base downstream from the seed 
    *     matching region)
    * @param len the length of the part that is expected to be relevant
    *    for the alignment - Defaults to 30 (since miRNAs are about 21-22
    *    nucleotides long and we assume that there would be no big loops
    *    in the optimal alignment)
    *
    * @return mRNA containing the subsequence relevant for the alignment
    *********************************************************************/
    
    inline mRNA mRNA::get_subsequence_for_alignment(chromosomePosition predicted_miRNA_three_prime_position, sequenceLength len) const {
      return get_subsequence_chr_to(predicted_miRNA_three_prime_position,len);
}

    /*****************************************************************//**
    * @brief apply SNP on mRNA
    *
    * This method is used to change a mRNA as described by a given SNP.
    * If the SNP information do not match those in the mRNA an unchanged
    * copy of the mRNA is returned (after stating an error).
    *
    * @param the_SNP SNP containing the information how the mRNA should be
    *     changed
    *
    * @return a copy of the mRNA with the changes defined by the SNP
    *********************************************************************/
    inline mRNA mRNA::mutate(const SNP & the_SNP) const {
      return mRNA(sequence::mutate(the_SNP));
}


} // namespace microSNPscore
#endif
