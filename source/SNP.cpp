
#include "SNP.h"
#include "miRNA.h"
#include "mRNA.h"

namespace microSNPscore {

    /*****************************************************************//**
    * @brief constructor
    *
    * This is used to create an instance of the class SNP.
    * Lowercase letters are treated as uppercase ones.
    * T is understood as Thymine and is treated as Uracil (simulating
    * transscription).
    * Dashes (-) are understood as Gaps and are omitted.
    * Other characters than A,a,C,c,G,g,U,u,T,t,X,x or - raise an error
    * and are treated as Mask.
    *
    * @param the_ID SNPID representing the ID of the SNP
    * @param reference_string String representing the nucleotide sequence
    *     of the SNP's reference sequence (Adenine: A, Cytosine: C,
    *     Guanine: G, Uracil: U, Mask: X)
    * @param alternative_string String representing the nucleotide sequence
    *     of the SNP's alternative sequence (Adenine: A, Cytosine: C,
    *     Guanine: G, Uracil: U, Mask: X)
    * @param the_chromosome chromosomeType representing the chromosome the
    *     SNP is located on
    * @param the_strand strandType representing the strand (Plus/Minus) on
    *     which the SNP is defined
    * @param the_position: chromosomePosition representing the position of
    *     the 5' end of the reference sequence
    * @return a SNP containing the given nucleotides located on the given
    *     chromosome, strand and position.
    *********************************************************************/
    
    SNP::SNP(SNPID the_ID, std::string reference_string, std::string alternative_string, chromosomeType the_chromosome, strandType the_strand, const chromosomePosition the_position) {
}

    /*****************************************************************//**
    * @brief get method for position attribute
    *
    * This method is used to access the  position on the chromosome (the
    * 5' end of the + strand (i.e. the 3' end of the - strand) beeing
    * position 1) of the 5' end of the SNP's reference sequence on the
    * given strand.
    *
    * @param the_strand strand (Plus or Minus) the SNP should be
    *     evaluated on
    *
    * @return the reference sequence's 5' end position on chromosome
    *********************************************************************/
    chromosomePosition SNP::get_position(strandType the_strand) const {
}

    /*****************************************************************//**
    * @brief compare sequence information
    *
    * This method is used to check whether the SNP would influence a given
    * sequence and if so whether the information about the reference
    * sequence stored in the SNP match those stored in the sequence (if
    * not an error is stated and @p false is returned).
    *
    * @param the_sequence sequence the SNP should be mapped on
    *
    * @return @p true if the SNP is on the sequence and the reference
    *     matches, @p false otherwise
    *********************************************************************/
    bool SNP::matches(const sequence & the_sequence) const {
}

    /*****************************************************************//**
    * @brief calculate deregulation score
    *
    * This method calculates the sore measuring measuring how much the
    * downregulation of the translation of a given mRNA induced by a given
    * miRNA binding at a target site starting at a given position is
    * changed by the SNP.
    *
    * @param the_miRNA  miRNA that is predicted to downregulate by the
    *     mRNA
    * @param the_mRNA  mRNA that is predicted to be downregulated by the
    *     miRNA
    * @param predicted_three_prime_position position on chromosome (the
    *     5' end of the + strand (i.e. the 3' end of the - strand) beeing
    *     position 1) that is predicted to be the mRNA nucleotide that
    *     would bind the miRNA 5' end (if it would bind) (i.e. one base
    *     downstream (3') from the seed match region)
    *
    * @return the deregulation score of the SNP for the target site of
    *     the miRNA starting at the given position in the given mRNA
    *********************************************************************/
    
    deregulation_score SNP::get_deregulation_score(const miRNA & the_miRNA, const mRNA & the_mRNA, chromosomePosition predicted_three_prime_position) const {
}


} // namespace microSNPscore
