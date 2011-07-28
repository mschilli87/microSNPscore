
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
    
    SNP::SNP(SNPID the_ID, std::string reference_string, std::string alternative_string, chromosomeType the_chromosome, strandType the_strand, const chromosomePosition the_position)
    :ID(the_ID),chromosome(the_chromosome),position(the_position),shift(alternative_string.length()-reference_string.length()) {
       /******************************************************************\ 
      | Convert the strings to nucleo base vectors and invert them for the |
      | sequence not matching the strand and move the position if needed:  |
       \******************************************************************/
      std::vector<nucleoBase> reference;
      for(std::string::const_iterator base_it(reference_string.begin());base_it!=reference_string.end();++base_it)
      {
        reference.push_back(make_base(*base_it));
      }
      std::vector<nucleoBase> alternative;
      for(std::string::const_iterator base_it(alternative_string.begin());base_it!=alternative_string.end();++base_it)
      {
        alternative.push_back(make_base(*base_it));
      }
      if(the_strand == Plus)
      {
        reference_plus=reference;
        alternative_plus=alternative;
        reference_minus=invert(reference);
        alternative_minus=invert(alternative);
      }
      else
      {
        position+=reference_string.length()-1;
        reference_minus=reference;
        alternative_minus=alternative;
        reference_plus=invert(reference);
        alternative_plus=invert(alternative);
      }
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

    /*****************************************************************//**
    * @brief convert char to nucleo base
    *
    * This method 'translates' chars to nucleo bases.
    * Lowercase letters are treated as uppercase ones.
    * T is understood as Thymine and is treated as Uracil (simulating
    * transscription).
    * Dashes (-) are understood as Gaps and are omitted.
    * Other characters than A,a,C,c,G,g,U,u,T,t,X,x or - raise an error
    * and are treated as Mask.
    *
    * @param the_char character telling which base to return
    *
    * @return nucleo base corresponding to the given character
    *********************************************************************/
    nucleoBase SNP::make_base(char the_char)
    {
}

    /*****************************************************************//**
    * @brief invert nucleo base sequence
    *
    * This method creates the reverse complement of a given sequence to
    * switch strands.
    *
    * @param vector containing the nucleo bases of a sequence on one
    * strand from 5' to 3'
    *
    * @return vector containing the nucleo bases of the given sequence on
    *     the other strand from 5' to 3'
    *********************************************************************/
    std::vector<nucleoBase> SNP::invert(const std::vector<nucleoBase> & the_bases)
    {
}


} // namespace microSNPscore
