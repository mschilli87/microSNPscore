#ifndef MICROSNPSCORE_MIRNA_H
#define MICROSNPSCORE_MIRNA_H


#include "sequence.h"
#include <string>
#include "nucleotide.h"

namespace microSNPscore { class mRNA; } 
namespace microSNPscore { class alignment; } 
namespace microSNPscore { class SNP; } 

namespace microSNPscore {

/*****************************************************************//**
* @brief downregulation score
*
* This represents the score measuring how much the translation of a
* mRNA is downregulated induced by a miRNA.
*********************************************************************/

typedef double downregulation_score;
/*****************************************************************//**
* @brief microRNA class
*
* This represents a microRNA sequence.
* It is a specialisation of the sequence class.
*********************************************************************/
class miRNA : public sequence {
  public:
    /*****************************************************************//**
    * @brief constructor
    *
    * This is used to create an instance of the class miRNA.
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
    * @param the_ID sequenceID representing the ID of the miRNA
    * @param sequence_string String representing the nucleotide sequence
    *     (Adenine: A, Cytosine: C, Guanine: G, Uracil: U, Mask: X)
    * @param the_chromosome chromosomeType representing the chromosome the
    *     miRNA is located on
    * @param the_strand strandType representing the strand (Plus/Minus) on
    *     which the miRNA is located
    * @param exon_starts: String representing the start positions (i.e.
    *     the end with the smaller distance to the chromosome start beeing
    *     the 5' end of the + strand and accordingly the 3' end of
    *     the - strand) of the exons containing the miRNA as
    *     comma-separated list.
    * @param exon_ends: String representing the end positions (i.e.
    *     the end with the smaller distance to the chromosome end beeing
    *     the 3' end of the + strand and accordingly the 5' end of
    *     the - strand) of the exons containing the miRNA as
    *     comma-separated list.
    * @return a miRNA containing the given nucleotides located on the
    *     given chromosome, strand and positions.
    *********************************************************************/
    miRNA(sequenceID the_id, std::string sequence_string, chromosomeType the_chromosome, strandType the_strand, std::string exon_starts, std::string exon_ends);

    /*****************************************************************//**
    * @brief calculate downregulation score
    *
    * This method calculates the sore measuring how much the translation of
    * a given mRNA will be downregulated by this miRNA through a target site
    * starting (miRNA 5' or mRNA 3') at a given position (reported from a
    * target prediction tool).
    * It implements the mirSVR algorithm.
    *
    * @param the_mRNA  mRNA that is predicted to be downregulated by the
    *     miRNA
    * @param predicted_three_prime_position position on chromosome (the
    *     5' end of the + strand (i.e. the 3' end of the - strand) beeing
    *     position 1) that is predicted to be the mRNA nucleotide that
    *     would bind the miRNA 5' end (if it would bind) (i.e. one base
    *     downstream (3') from the seed match region)
    *
    * @return the downregulation score for the target site of the miRNA
    *     starting at the given position in the given mRNA
    *********************************************************************/
    downregulation_score get_downregulation_score(const mRNA & the_mRNA, const chromosomePosition & predicted_three_prime_position) const;


  private:
    /*****************************************************************//**
    * @brief calculate deregulation score for one possible alignment
    *
    * This method calculates the sore measuring how much the translation of
    * a given mRNA will be downregulated by this miRNA through a target site
    * starting (miRNA 5' or mRNA 3') at a given position (reported from a
    * target prediction tool) considering a given alignment.
    * It implements the mirSVR algorithm.
    *
    * @param the_mRNA  mRNA that is predicted to be downregulated by the
    *     miRNA
    * @param predicted_three_prime_position position on chromosome (the
    *     5' end of the + strand (i.e. the 3' end of the - strand) beeing
    *     position 1) that is predicted to be the mRNA nucleotide that
    *     would bind the miRNA 5' end (if it would bind) (i.e. one base
    *     downstream (3') from the seed match region)
    * @param the_alignment an alignment that is considered to be the best
    *     one for the miRNA-induced downregulation
    *
    * @return the downregulation score for the target site of the miRNA
    *     starting at the given position in the given mRNA considering the
    *     given alignment
    *********************************************************************/
    downregulation_score downregulation_score_candidate(const mRNA & the_mRNA, const chromosomePosition & predicted_three_prime_position, const alignment & the_alignment) const;


  public:
    inline miRNA mutate(const SNP & the_SNP) const;


  private:
    /*****************************************************************//**
    * @brief internal constructor
    *
    * This method is used to convert a sequence to a miRNA.
    *
    * @param the_sequence const sequence reference to the sequence that
    *     should become an miRNA
    *
    * @return miRNA with the same attributes as the given sequence
    *********************************************************************/
    miRNA(const sequence & the_sequence);

};
    inline miRNA miRNA::mutate(const SNP & the_SNP) const {
      return miRNA(sequence::mutate(the_SNP));
}


} // namespace microSNPscore
#endif
