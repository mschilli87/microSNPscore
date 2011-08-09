#ifndef MICROSNPSCORE_MIRNA_H
#define MICROSNPSCORE_MIRNA_H


#include "sequence.h"
#include <string>
#include "nucleotide.h"
#include "alignment.h"

namespace microSNPscore { class SNP; } 
namespace microSNPscore { class mRNA; } 
namespace microSNPscore { class conservationList; } 
namespace microSNPscore { class alignment; } 

namespace microSNPscore {

/*****************************************************************//**
* @brief downregulation score
*
* This represents the score measuring how much the translation of a
* mRNA is downregulated induced by a miRNA.
*********************************************************************/

typedef double downregulationScore;
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

    inline miRNA mutate(const SNP & the_SNP) const;

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
    * @param conservations conservationList containing the conservation
    *     ranges to use for the conservation scoring
    *
    * @return the downregulation score for the target site of the miRNA
    *     starting at the given position in the given mRNA
    *
    * @see SNP::get_deregulation_score()
    *********************************************************************/
    downregulationScore get_downregulation_score(const mRNA & the_mRNA, chromosomePosition predicted_three_prime_position, const conservationList & conservations) const;


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
    * @param conservations conservationList containing the conservation
    *     ranges to use for the conservation scoring
    *
    * @return the downregulation score for the target site of the miRNA
    *     starting at the given position in the given mRNA considering the
    *     given alignment
    *
    * @see get_downregulation_score()
    *********************************************************************/
    static downregulationScore downregulation_score_candidate(const mRNA & the_mRNA, chromosomePosition predicted_three_prime_position, const alignment & the_alignment, const conservationList & conservations);

    /*****************************************************************//**
    * @brief secondary structure features calculation
    *
    * This method is used to calculate the accessability features for the
    * downregulation score calculation.
    * It uses rnafold.
    *
    * @param features array the twenty raw feature scores should be
    *     stored in
    * @param mRNA_subsequence the part of the mRNA that should be
    *     considerd for the secondary structure prediction
    * @param predicted_three_prime_position position on chromosome (the
    *     5' end of the + strand (i.e. the 3' end of the - strand) beeing
    *     position 1) that is predicted to be the mRNA nucleotide that
    *     would bind the miRNA 5' end (if it would bind) (i.e. one base
    *     downstream (3') from the seed match region)
    *
    * @see downregulation_score_candidate()
    *********************************************************************/
    static void calculate_accessability_features(downregulationScore features[], const mRNA & mRNA_subsequence, chromosomePosition predicted_three_prime_position);

    /*****************************************************************//**
    * @brief target site conservation feature calculation
    *
    * This method is used to calculate the conservation feature for the
    * downregulation score calculation.
    *
    * @param the_alignment an alignment that is considered to be the best
    *     one for the miRNA-induced downregulation
    * @param mRNA_strand the strand (Plus or Minus) the mRNA is transcribed
    *     from
    * @param mRNA_chromosome the chromosome the mRNA is located on
    * @param conservations conservationList containing the conservation
    *     ranges to use for the scoring
    *
    * @return conservation score for the target site represented by the
    *     given alignment
    *
    * @see downregulation_score_candidate()
    *********************************************************************/
    static downregulationScore calculate_conservation_feature(const alignment & the_alignment, strandType mRNA_strand, const chromosomeType & mRNA_chromosome, const conservationList & conservations);

    /*****************************************************************//**
    * @brief local AU content feature calculation
    *
    * This method is used to calculate the AU content feature for the
    * downregulation score calculation.
    *
    * @param downstream_mRNA_subsequence the part of the mRNA that should
    *     be considerd for the downstream AU content calculation
    * @param upstream_mRNA_subsequence the part of the mRNA that should be
    *     considerd for the upstream AU content calculation
    * @param seed_type the seed type (sixMer, sevenMerAOne,
    *     sevenMerMEight, EightMer) of the alignment that is considered to
    *     be the best one for the miRNA-induced downregulation
    *
    * @return AU content score for the target site of the given seed type
    *     with the given flanking sequences
    *
    * @see downregulation_score_candidate()
    *********************************************************************/
    static downregulationScore calculate_AU_content_feature(const mRNA & downstream_mRNA_subsequence, const mRNA & upstream_mRNA_subsequence, seedType seed_type);

    /*****************************************************************//**
    * @brief 3' match feature calculation
    *
    * This method is used to calculate the 3'-score feature for the
    * downregulation score calculation.
    *
    * @param the_alignment an alignment that is considered to be the best
    *     one for the miRNA-induced downregulation
    *
    * @return 3' feature score for the target site represented by the
    *     given alignment
    *
    * @see downregulation_score_candidate()
    *********************************************************************/
    static downregulationScore calculate_three_prime_feature(const alignment & the_alignment);

    /*****************************************************************//**
    * @brief UTR distance feature calculation
    *
    * This method is used to calculate the UTR distance feature for the
    * downregulation score calculation.
    *
    * @param the_mRNA mRNA sequence the target site is on
    * @param predicted_three_prime_position position on chromosome (the
    *     5' end of the + strand (i.e. the 3' end of the - strand) beeing
    *     position 1) that is predicted to be the mRNA nucleotide that
    *     would bind the miRNA 5' end (if it would bind) (i.e. one base
    *     downstream (3') from the seed match region)
    * @param seed_type the seed type (sixMer, sevenMerAOne,
    *     sevenMerMEight, EightMer) of the alignment that is considered to
    *     be the best one for the miRNA-induced downregulation
    *
    * @return UTR distance score for a target site of the given site
    *     located at the given position in the given mRNA
    *
    * @see downregulation_score_candidate()
    *********************************************************************/
    static downregulationScore calculate_UTR_dist_feature(const mRNA & the_mRNA, chromosomePosition predicted_three_prime_position, seedType seed_type);

    /*****************************************************************//**
    * @brief seed match features calculation
    *
    * This method is used to calculate the seed match features for the
    * downregulation score calculation.
    *
    * @param features array the ninefeature scores should be stored in
    * @param the_alignment an alignment that is considered to be the best
    *     one for the miRNA-induced downregulation
    *
    * @see downregulation_score_candidate()
    *********************************************************************/
    static void calculate_seed_match_features(downregulationScore features[], const alignment & the_alignment);

};
    inline miRNA miRNA::mutate(const SNP & the_SNP) const {
      return miRNA(sequence::mutate(the_SNP));
}


} // namespace microSNPscore
#endif
