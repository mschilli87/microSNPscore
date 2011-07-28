
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
       /*****************************************************************\ 
      | Convert the strings to nucleo base vectors (ignoring gaps) and    |
      | invert them for the sequence not matching the strand and move the | 
      | position if needed:                                               |
       \*****************************************************************/
      std::vector<nucleoBase> reference;
      for(std::string::const_iterator base_it(reference_string.begin());base_it!=reference_string.end();++base_it)
      {
        if(*base_it != '-')
        {
          reference.push_back(make_base(*base_it));
        }
      }
      std::vector<nucleoBase> alternative;
      for(std::string::const_iterator base_it(alternative_string.begin());base_it!=alternative_string.end();++base_it)
      {
        if(*base_it != '-')
        {
          alternative.push_back(make_base(*base_it));
        }
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
    * A SNP is said to match a sequence if its whole reference sequence is
    * located on a single exon of the sequence and the nucleo bases in
    * that region are the same.
    *
    * @param the_sequence sequence the SNP should be mapped on
    *
    * @return @p true if the SNP is on the sequence and the reference
    *     matches, @p false otherwise
    *********************************************************************/
    
    bool SNP::matches(const sequence & the_sequence) const {
       /*****************************************************************\ 
      | Compare the chromosomes and if they are the same search for an    |
      | exon starting before and ending after the reference and if one is |
      | found, compare the sequence's bases with those of the reference   |
      | sequence (on the corresponding strand):                           |
       \*****************************************************************/
      if(get_chromosome() != the_sequence.get_chromosome())
      {
        return false;
      }
      else // same chromosome
      {
        const sequenceLength length = reference_end(Plus)-reference_begin(Plus);
        sequence::const_exon_iterator exon_it(the_sequence.exons_begin());
        for(;exon_it != the_sequence.exons_end() && exon_it->get_start() > get_position(Plus);++exon_it ) {/* nothing */}
        if(exon_it == the_sequence.exons_end() || exon_it->get_end() < (get_position(Plus) + length))
        {
          return false;
        }
        else // whole reference on one exon of the sequence
        {
          const strandType the_strand(the_sequence.get_strand());
          const const_iterator begin(reference_begin(the_strand));
          const const_iterator end(reference_end(the_strand));
          const_iterator bases_it(begin);
          for(sequence::const_iterator sequence_it(the_sequence.get_nucleotide_chr(get_position(the_strand)));
              bases_it!=end && sequence_it->get_base() == *bases_it;++bases_it,++sequence_it) {/* nothing */}
          return bases_it == end;
        }
      }
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
    
    deregulationScore SNP::get_deregulation_score(const miRNA & the_miRNA, const mRNA & the_mRNA, chromosomePosition predicted_three_prime_position) const {
       /*************************************************************\ 
      | Verify that the SNP may have influence on the downregulation  |
      | score and if so return the score difference between reference |
      | and alternative shifting the predicted 3' position if needed: |
       \*************************************************************/
      const bool SNP_on_miRNA = matches(the_miRNA);
      const bool SNP_on_mRNA = matches(the_mRNA);
      if(!SNP_on_miRNA && !SNP_on_mRNA)
      {
        return 0;
      }
      else
      {
        return the_miRNA.get_downregulation_score(the_mRNA,predicted_three_prime_position) -
               (SNP_on_miRNA ?
                the_miRNA.mutate(*this).get_downregulation_score(the_mRNA,predicted_three_prime_position) :
                the_miRNA.get_downregulation_score(the_mRNA.mutate(*this),predicted_three_prime_position +
                                                                          (predicted_three_prime_position < (get_position(Plus) +
                                                                                                             reference_end(Plus) -
                                                                                                             reference_begin(Plus)) ?
                                                                           0 : get_shift())));
      }
}

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
    
    nucleoBase SNP::make_base(char the_char)
    {
       /******************************************************************\ 
      | Return the nucleo base corresponding to the char or state an error |
      | and return Mask. Because return exists the function there is no    |
      | need for break statements:                                         |
       \******************************************************************/
      switch(the_char)
      {
        case 'a':
        case 'A': return Adenine;
        case 't':
        case 'T':
        case 'u':
        case 'U': return Uracil;
        case 'c':
        case 'C': return Cytosine;
        case 'g':
        case 'G': return Guanine;
        case 'x':
        case 'X': return Mask;
        default:
          std::cerr << "microSNPscore::SNP::make_base\n";
          std::cerr << " ==> illegal nucleo base character: \n";
          std::cerr << the_char << std::endl;
          std::cerr << "  --> assuming Mask\n";
          return Mask;
      }
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
       /***********************************************************\ 
      | Create an empty nucleo base vector and iterate over the the |
      | givenbases (if there are some) in reverse order appending   |
      | their complement or Mask (stating an error) to the vector   |
      | before returning it:                                        |
       \***********************************************************/
      std::vector<nucleoBase> complement;
      if(the_bases.begin() != the_bases.end())
      {
        for(const_iterator base_it(the_bases.end()-1);base_it>=the_bases.begin();--base_it)
        {
          nucleoBase the_base(Adenine);
          switch(the_base)
          {
            case Uracil:
              break;
            case Adenine:
              the_base = Uracil;
              break;
            case Guanine:
              the_base = Cytosine;
              break;
            case Cytosine:
              the_base = Guanine;
              break;
            case Gap:
              the_base = Gap;
              break;
            case Mask:
              the_base = Mask;
              break;
            default:
              std::cerr << "microSNPscore::SNP::invert\n";
              std::cerr << " ==> Undefined nucleo base: ";
              std::cerr << the_base << std:: endl;
              std::cerr << "  --> assuming Masked\n";
              the_base = Mask;
          }
          complement.push_back(the_base);
        }
      }
      return complement;
}


} // namespace microSNPscore
