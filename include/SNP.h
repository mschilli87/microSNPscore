#ifndef MICROSNPSCORE_SNP_H
#define MICROSNPSCORE_SNP_H


#include <string>
#include <vector>
#include "nucleotide.h"
#include "sequence.h"

namespace microSNPscore { class sequence; } 
namespace microSNPscore { class miRNA; } 
namespace microSNPscore { class mRNA; } 

namespace microSNPscore {

/*****************************************************************//**
* @brief deregulation score
*
* This represents the score measuring how much the downregulation of
* the translation of a mRNA induced by a miRNA is changed by a SNP.
*********************************************************************/

typedef double deregulationScore;
/*****************************************************************//**
* @brief SNP ID type
*
* This represents the ID of a SNP.
*********************************************************************/
typedef std::string SNPID;
/*****************************************************************//**
* @brief SNP class
*
* This represents a SNP (meaning any variant - even longer indels).
*********************************************************************/
class SNP {
  public:
    /*****************************************************************//**
    * @brief const iterartor type
    *
    * This type is used to access the reference and alternative
    * sequences of the SNP.
    *********************************************************************/
    typedef std::vector<nucleoBase>::const_iterator const_iterator;

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
    
    SNP(SNPID the_ID, std::string reference_string, std::string alternative_string, chromosomeType the_chromosome, strandType the_strand, const chromosomePosition the_position);

    /*****************************************************************//**
    * @brief get method for ID attribute
    *
    * This method is used to access the ID of the SNP.
    *
    * @return the ID of the SNP
    *********************************************************************/
    inline const SNPID get_ID() const;

    /*****************************************************************//**
    * @brief get method for chromosome attribute
    *
    * This method is used to access the chromosome the SNP is located on.
    *
    * @return the chromosome of the SNP
    *********************************************************************/
    inline const chromosomeType get_chromosome() const;

    /*****************************************************************//**
    * @brief get method for position attribute
    *
    * This method is used to access the number of nucleotides the
    * alternative sequence is longer than the reference sequence.
    *
    * @return the position shift of the SNP
    *********************************************************************/
    inline const short get_shift() const;

    /*****************************************************************//**
    * @brief get method for position attribute
    *
    * This method is used to access the position on the chromosome (the
    * 5' end of the + strand (i.e. the 3' end of the - strand) beeing
    * position 1) of the 5' end of the SNP's reference sequence on the
    * given strand.
    *
    * @param the_strand strand (Plus or Minus) the SNP should be
    *     evaluated on
    *
    * @return the reference sequence's 5' end position on chromosome
    *********************************************************************/
    inline chromosomePosition get_position(strandType the_strand) const;

    /*****************************************************************//**
    * @brief reference vector begin
    *
    * This method is used to access the first nucleo base of the reference
    * sequence of the SNP.
    *
    * @param the_strand strand (Plus or Minus) the SNP should be
    *     evaluated on
    *
    * @return const iterator pointing to the first nucleo base of the
    *     SNP's reference sequence
    *********************************************************************/
    inline const_iterator reference_begin(strandType the_strand) const;

    /*****************************************************************//**
    * @brief reference vector end
    *
    * This method is used to access the end of the SNP's reference
    * sequence nucleo base vector.
    *
    * @param the_strand strand (Plus or Minus) the SNP should be
    *     evaluated on
    *
    * @return const iterator pointing behind the last nucleo base of the
    *     SNP's reference sequence
    *********************************************************************/
    inline const_iterator reference_end(strandType the_strand) const;

    /*****************************************************************//**
    * @brief alternative vector begin
    *
    * This method is used to access the first nucleo base of the
    * alternative sequence of the SNP.
    *
    * @param the_strand strand (Plus or Minus) the SNP should be
    *     evaluated on
    *
    * @return const iterator pointing to the first nucleo base of the
    *     SNP's alternative sequence
    *********************************************************************/
    inline const_iterator alternative_begin(strandType the_strand) const;

    /*****************************************************************//**
    * @brief alternative vector end
    *
    * This method is used to access the end of the SNP's
    * alternative sequence nucleo base vector.
    *
    * @param the_strand strand (Plus or Minus) the SNP should be
    *     evaluated on
    *
    * @return const iterator pointing behind the last nucleo base of the
    *     SNP's alternative sequence
    *********************************************************************/
    inline const_iterator alternative_end(strandType the_strand) const;

    /*****************************************************************//**
    * @brief compare sequence information
    *
    * This method is used to check whether the SNP would influence a given
    * sequence and if so whether the information about the reference
    * sequence stored in the SNP match those stored in the sequence (if
    * not an error is stated and @p false is returned).
    * A SNP is said to match a sequence if its whole reference sequence is
    * located on a single exon of the sequence and the nucleo bases in
    * that region are the same.
    *
    * @param the_sequence sequence the SNP should be mapped on
    *
    * @return @p true if the SNP is on the sequence and the reference
    *     matches, @p false otherwise
    *********************************************************************/
    
    bool matches(const sequence & the_sequence) const;

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
    
    deregulationScore get_deregulation_score(const miRNA & the_miRNA, const mRNA & the_mRNA, chromosomePosition predicted_three_prime_position) const;


  private:
    /*****************************************************************//**
    * @brief convert char to nucleo base
    *
    * This method 'translates' chars to nucleo bases.
    * Lowercase letters are treated as uppercase ones.
    * T is understood as Thymine and is treated as Uracil (simulating
    * transscription).
    * Other characters than A,a,C,c,G,g,U,u,T,t,X or x raise an error and
    * are treated as Mask.
    *
    * @param the_char character telling which base to return
    *
    * @return nucleo base corresponding to the given character
    *********************************************************************/
    
    static nucleoBase make_base(char the_char);

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
    static std::vector<nucleoBase> invert(const std::vector<nucleoBase> & the_bases);

    /*****************************************************************//**
    * @brief SNP ID
    *
    * This is the ID of this SNP.
    *********************************************************************/
    SNPID ID;

    /*****************************************************************//**
    * @brief chromosome
    *
    * This is the chromosome the SNP is located on.
    *********************************************************************/
    chromosomeType chromosome;

    /*****************************************************************//**
    * @brief position on chromosome
    *
    * This is the position on the chromosome (the 5' end of the + strand
    * (i.e. the 3' end of the - strand) beeing position 1) of the 5' end
    * of the SNP's reference sequence on the + strand.
    *********************************************************************/
    chromosomePosition position;

    /*****************************************************************//**
    * @brief reference sequence on + strand
    *
    * This is a vector containing the nucleo bases of the SNP's reference
    * sequence (5' to 3' on the + strand).
    *********************************************************************/
    std::vector<nucleoBase> reference_plus;

    /*****************************************************************//**
    * @brief alternative sequence on + strand
    *
    * This is a vector containing the nucleo bases of the SNP's alternative
    * sequence (5' to 3' on the - strand).
    *********************************************************************/
    std::vector<nucleoBase> alternative_plus;

    /*****************************************************************//**
    * @brief reference sequence on - strand
    *
    * This is a vector containing the nucleo bases of the SNP's reference
    * sequence (5' to 3' on the - strand).
    *********************************************************************/
    std::vector<nucleoBase> reference_minus;

    /*****************************************************************//**
    * @brief alternative sequence on - strand
    *
    * This is a vector containing the nucleo bases of the SNP's alternative
    * sequence.
    *********************************************************************/
    std::vector<nucleoBase> alternative_minus;

    /*****************************************************************//**
    * @brief length difference between reference and alternative
    *
    * This is the number of nucleotides the alternative sequence is longer
    * than the reference sequence.
    *********************************************************************/
    short shift;

};
    /*****************************************************************//**
    * @brief get method for ID attribute
    *
    * This method is used to access the ID of the SNP.
    *
    * @return the ID of the SNP
    *********************************************************************/
    inline const SNPID SNP::get_ID() const {
      return ID;
    }

    /*****************************************************************//**
    * @brief get method for chromosome attribute
    *
    * This method is used to access the chromosome the SNP is located on.
    *
    * @return the chromosome of the SNP
    *********************************************************************/
    inline const chromosomeType SNP::get_chromosome() const {
      return chromosome;
    }

    /*****************************************************************//**
    * @brief get method for position attribute
    *
    * This method is used to access the number of nucleotides the
    * alternative sequence is longer than the reference sequence.
    *
    * @return the position shift of the SNP
    *********************************************************************/
    inline const short SNP::get_shift() const {
      return shift;
    }

    /*****************************************************************//**
    * @brief get method for position attribute
    *
    * This method is used to access the position on the chromosome (the
    * 5' end of the + strand (i.e. the 3' end of the - strand) beeing
    * position 1) of the 5' end of the SNP's reference sequence on the
    * given strand.
    *
    * @param the_strand strand (Plus or Minus) the SNP should be
    *     evaluated on
    *
    * @return the reference sequence's 5' end position on chromosome
    *********************************************************************/
    inline chromosomePosition SNP::get_position(strandType the_strand) const {
      return position + (the_strand == Plus ? 0 : (reference_end(Plus)-reference_begin(Plus)-1));
}

    /*****************************************************************//**
    * @brief reference vector begin
    *
    * This method is used to access the first nucleo base of the reference
    * sequence of the SNP.
    *
    * @param the_strand strand (Plus or Minus) the SNP should be
    *     evaluated on
    *
    * @return const iterator pointing to the first nucleo base of the
    *     SNP's reference sequence
    *********************************************************************/
    inline SNP::const_iterator SNP::reference_begin(strandType the_strand) const {
      return (the_strand == Plus ? reference_plus : reference_minus).begin();
}

    /*****************************************************************//**
    * @brief reference vector end
    *
    * This method is used to access the end of the SNP's reference
    * sequence nucleo base vector.
    *
    * @param the_strand strand (Plus or Minus) the SNP should be
    *     evaluated on
    *
    * @return const iterator pointing behind the last nucleo base of the
    *     SNP's reference sequence
    *********************************************************************/
    inline SNP::const_iterator SNP::reference_end(strandType the_strand) const {
      return (the_strand == Plus ? reference_plus : reference_minus).end();
}

    /*****************************************************************//**
    * @brief alternative vector begin
    *
    * This method is used to access the first nucleo base of the
    * alternative sequence of the SNP.
    *
    * @param the_strand strand (Plus or Minus) the SNP should be
    *     evaluated on
    *
    * @return const iterator pointing to the first nucleo base of the
    *     SNP's alternative sequence
    *********************************************************************/
    inline SNP::const_iterator SNP::alternative_begin(strandType the_strand) const {
      return (the_strand == Plus ? alternative_plus : alternative_minus).begin();
}

    /*****************************************************************//**
    * @brief alternative vector end
    *
    * This method is used to access the end of the SNP's
    * alternative sequence nucleo base vector.
    *
    * @param the_strand strand (Plus or Minus) the SNP should be
    *     evaluated on
    *
    * @return const iterator pointing behind the last nucleo base of the
    *     SNP's alternative sequence
    *********************************************************************/
    inline SNP::const_iterator SNP::alternative_end(strandType the_strand) const {
      return (the_strand == Plus ? alternative_plus : alternative_minus).end();
}


} // namespace microSNPscore
#endif
