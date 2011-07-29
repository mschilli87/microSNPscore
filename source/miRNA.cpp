
#include <algorithm>
// for std::max (downregulation score calculation)
#include <cmath>
// for exp (downregulation score sigmoid function)
#include "miRNA.h"
#include "mRNA.h"
#include "alignment.h"
#include "SNP.h"

namespace microSNPscore {

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
    miRNA::miRNA(sequenceID the_id, std::string sequence_string, chromosomeType the_chromosome, strandType the_strand, std::string exon_starts, std::string exon_ends)
    :sequence(the_id,sequence_string,the_chromosome,the_strand,exon_starts,exon_ends) {
}

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
    *
    * @see SNP::get_deregulation_score()
    *********************************************************************/
    downregulationScore miRNA::get_downregulation_score(const mRNA & the_mRNA, chromosomePosition predicted_three_prime_position) const {
       /*********************************************************\ 
      | Calculate downregulation score candidates for all optimal |
      | alignments and return the maximum:                        |
       \*********************************************************/
      optimalAlignmentList alignments(the_mRNA.get_subsequence_for_alignment(predicted_three_prime_position),*this);
      if(alignments.begin() == alignments.end())
      {
        return 0;
      }
      else
      {
        downregulationScore downregulation_score(downregulation_score_candidate(the_mRNA,predicted_three_prime_position,*alignments.begin()));
        for(optimalAlignmentList::const_iterator alignment_it(alignments.begin()+1);alignment_it!=alignments.end();++alignment_it)
        {
          downregulation_score = std::max(downregulation_score,downregulation_score_candidate(the_mRNA,predicted_three_prime_position,*alignment_it));
        }
        return downregulation_score;
      }
}

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
    *
    * @see get_downregulation_score()
    *********************************************************************/
    downregulationScore miRNA::downregulation_score_candidate(const mRNA & the_mRNA, chromosomePosition predicted_three_prime_position, const alignment & the_alignment)
    {
       /****************************************\ 
      | Define number and names of the features: |
       \****************************************/
      const unsigned short int feature_count = 34;
      
      const unsigned short int    UTRLength =  0;
      const unsigned short int         SS01 =  1;/*
      const unsigned short int         SS02 =  2;
      const unsigned short int         SS03 =  3;
      const unsigned short int         SS04 =  4;
      const unsigned short int         SS05 =  5;
      const unsigned short int         SS06 =  6;
      const unsigned short int         SS07 =  7;
      const unsigned short int         SS08 =  8;
      const unsigned short int         SS09 =  9;
      const unsigned short int         SS10 = 10;
      const unsigned short int         SS11 = 11;
      const unsigned short int         SS12 = 12;
      const unsigned short int         SS13 = 13;
      const unsigned short int         SS14 = 14;
      const unsigned short int         SS15 = 15;
      const unsigned short int         SS16 = 16;
      const unsigned short int         SS17 = 17;
      const unsigned short int         SS18 = 18;
      const unsigned short int         SS19 = 19;
      const unsigned short int         SS20 = 20;*/
      const unsigned short int conservation = 21;
      const unsigned short int   AU_content = 22;
      const unsigned short int  three_prime = 23;
      const unsigned short int     UTR_dist = 24;
      const unsigned short int           A1 = 25;/*
      const unsigned short int           m2 = 26;
      const unsigned short int           m3 = 27;
      const unsigned short int           m4 = 28;
      const unsigned short int           m5 = 29;
      const unsigned short int           m6 = 30;
      const unsigned short int           m7 = 31;
      const unsigned short int           m8 = 32;
      const unsigned short int           m9 = 33;*/
      
       /************************************************************\ 
      | Define the means of the features as in mirSVR (seed features |
      | where not z-transformed in mirSVR - simulated by mean 0):    |
       \************************************************************/
      const downregulationScore feature_means[feature_count] = {
                                             /*    UTRLength */  1007.05587
                                             /*         SS01 */,    6.42691
                                             /*         SS02 */,    6.36598
                                             /*         SS03 */,    6.27593
                                             /*         SS04 */,    6.18784
                                             /*         SS05 */,    5.92570
                                             /*         SS06 */,    5.91493
                                             /*         SS07 */,    6.07210
                                             /*         SS08 */,    6.17630
                                             /*         SS09 */,    6.20022
                                             /*         SS10 */,    6.18562
                                             /*         SS11 */,    6.20394
                                             /*         SS12 */,    6.36488
                                             /*         SS13 */,    6.58276
                                             /*         SS14 */,    7.05142
                                             /*         SS15 */,    7.05196
                                             /*         SS16 */,    6.99292
                                             /*         SS17 */,    7.00693
                                             /*         SS18 */,    7.14608
                                             /*         SS19 */,    7.06098
                                             /*         SS20 */,    6.94628
                                             /* conservation */,    0.57633
                                             /*   AU_content */,    0.58134
                                             /*  three_prime */,    2.32868
                                             /*     UTR_dist */,  233.56983
                                             /*           1A */,    0
                                             /*           m2 */,    0
                                             /*           m3 */,    0
                                             /*           m4 */,    0
                                             /*           m5 */,    0
                                             /*           m6 */,    0
                                             /*           m7 */,    0
                                             /*           m8 */,    0
                                             /*           m9 */,    0
                                                               };
      
       /*************************************************************\ 
      | Define the sigmas of the features as in mirSVR (seed features |
      | where not z-transformed in mirSVR - simulated by sigma 1):    |
       \*************************************************************/
      const downregulationScore feature_sigmas[feature_count] = {
                                              /*    UTRLength */  694.87216
                                              /*         SS01 */,   3.32566
                                              /*         SS02 */,   3.21532
                                              /*         SS03 */,   3.17074
                                              /*         SS04 */,   2.96199
                                              /*         SS05 */,   2.71384
                                              /*         SS06 */,   2.47704
                                              /*         SS07 */,   2.59893
                                              /*         SS08 */,   2.58327
                                              /*         SS09 */,   2.67000
                                              /*         SS10 */,   2.63056
                                              /*         SS11 */,   2.64559
                                              /*         SS12 */,   2.68505
                                              /*         SS13 */,   2.95636
                                              /*         SS14 */,   3.37445
                                              /*         SS15 */,   3.53977
                                              /*         SS16 */,   3.59986
                                              /*         SS17 */,   3.65355
                                              /*         SS18 */,   3.76228
                                              /*         SS19 */,   3.76480
                                              /*         SS20 */,   3.71696
                                              /* conservation */,   0.07974
                                              /*   AU_content */,   0.14308
                                              /*  three_prime */,   0.73776
                                              /*     UTR_dist */, 227.99686
                                              /*           1A */,   1
                                              /*           m2 */,   1
                                              /*           m3 */,   1
                                              /*           m4 */,   1
                                              /*           m5 */,   1
                                              /*           m6 */,   1
                                              /*           m7 */,   1
                                              /*           m8 */,   1
                                              /*           m9 */,   1
                                                                };
      
       /************************************************\ 
      | Define the weights of the features as in mirSVR: |
       \************************************************/
      const downregulationScore feature_weights[feature_count] = {
                                               /*    UTRLength */   0.042589787580216
                                               /*         SS01 */,  0.003374756562456
                                               /*         SS02 */, -0.003802322920479
                                               /*         SS03 */, -0.008904121278094
                                               /*         SS04 */,  0.019566604120901
                                               /*         SS05 */, -0.020852554122244
                                               /*         SS06 */,  0.032737102268958
                                               /*         SS07 */, -0.029488605561073
                                               /*         SS08 */,  0.011107428379054
                                               /*         SS09 */,  0.013359252892328
                                               /*         SS10 */, -0.028984578128118
                                               /*         SS11 */,  0.062695067984425
                                               /*         SS12 */, -0.021906784041198
                                               /*         SS13 */, -0.009700573458590
                                               /*         SS14 */,  0.005674996762089
                                               /*         SS15 */, -0.011251607785794
                                               /*         SS16 */, -0.010795836406302
                                               /*         SS17 */,  0.006186147032099
                                               /*         SS18 */,  0.059234976420213
                                               /*         SS19 */, -0.031649418324568
                                               /*         SS20 */, -0.005148762246524
                                               /* conservation */, -0.046310553500448
                                               /*   AU_content */, -0.112958003179723
                                               /*  three_prime */,  0.001502205540968
                                               /*     UTR_dist */,  0.001660551262213
                                               /*           1A */, -0.069506963524673
                                               /*           m2 */, -0.370993899195771
                                               /*           m3 */, -0.464640990134142
                                               /*           m4 */, -0.548318765529733
                                               /*           m5 */, -0.492871203746197
                                               /*           m6 */, -0.481534081295122
                                               /*           m7 */, -0.512792852776873
                                               /*           m8 */, -0.350587322571024
                                               /*           m9 */,  0.007880767476343
                                                                 };
      
       /****************************************************\ 
      | Define the sigmoid function parameters as in mirSVR: |
       \****************************************************/
      const downregulationScore score_sigmoid_alpha = 10.6915;
      const downregulationScore score_sigmoid_beta = 4.2222;
      const downregulationScore score_sigmoid_C = 1.3681;
      
       /****************************************************************\ 
      | Define the the score bias as in mirSVR (they defined it negative |
      | and substracted it from the score, we define it positive and add |
      | the score to it):                                                |
       \****************************************************************/
      const downregulationScore score_bias = 2.8827265367288781;
      
       /*********************************\ 
      | Calculate the raw feature scores: |
       \*********************************/
      downregulationScore features[feature_count];
      features[UTRLength]=the_mRNA.get_length();
      calculate_accessability_features(&features[SS01],the_mRNA.get_subsequence_for_accessability(predicted_three_prime_position),predicted_three_prime_position);
      features[conservation]=calculate_conservation_feature(the_alignment,the_mRNA.get_strand());
      features[AU_content]=calculate_AU_content_feature(the_mRNA.get_subsequence_for_downstream_AU_content(predicted_three_prime_position),
                                                        the_mRNA.get_subsequence_for_upstream_AU_content(predicted_three_prime_position),
                                                        the_alignment.get_seed_type());
      features[three_prime]=calculate_three_prime_feature(the_alignment);
      features[UTR_dist]=calculate_UTR_dist_feature(the_mRNA,predicted_three_prime_position,the_alignment.get_seed_type());
      calculate_seed_match_features(&features[A1],the_alignment);
      
       /***************************************************************\ 
      | Initialize raw score with bias and add sum of the z-transformed |
      | feature scores:                                                 |
       \***************************************************************/
      downregulationScore the_score(score_bias);
      for(unsigned short int feature_number=0;feature_number<feature_count;++feature_number)
      {
        the_score += (features[feature_number] - feature_means[feature_number]) / feature_sigmas[feature_number] * feature_weights[feature_number];
      }
      
       /********************************************************\
| Apply sigmoid function before returning the final score: |
       \********************************************************/
      return score_sigmoid_C / (1 + exp(score_sigmoid_alpha * the_score + score_sigmoid_beta));
}

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
    miRNA::miRNA(const sequence & the_sequence)
    :sequence(the_sequence) {
}

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
    void miRNA::calculate_accessability_features(downregulationScore features[], const mRNA & mRNA_subsequence, chromosomePosition predicted_three_prime_position)
    {
      features[0]=6.42691;
      features[1]=6.36598;
      features[2]=6.27593;
      features[3]=6.18784;
      features[4]=5.92570;
      features[5]=5.91493;
      features[6]=6.07210;
      features[7]=6.17630;
      features[8]=6.20022;
      features[9]=6.18562;
      features[10]=6.20394;
      features[11]=6.36488;
      features[12]=6.58276;
      features[13]=7.05142;
      features[14]=7.05196;
      features[15]=6.99292;
      features[16]=7.00693;
      features[17]=7.14608;
      features[18]=7.06098;
      features[19]=6.94628;
}

    /*****************************************************************//**
    * @brief target site conservation feature calculation
    *
    * This method is used to calculate the conservation feature for the
    * downregulation score calculation.
    *
    * @param the_alignment an alignment that is considered to be the best
    *     one for the miRNA-induced downregulation
    * @param the_strand the strand (Plus or Minus) the mRNA is transcribed
    *     from
    *
    * @return conservation score for the target site represented by the
    *     given alignment
    *
    * @see downregulation_score_candidate()
    *********************************************************************/
    downregulationScore miRNA::calculate_conservation_feature(const alignment & the_alignment, strandType the_strand)
    {
      return 0.57633;
}

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
    downregulationScore miRNA::calculate_AU_content_feature(const mRNA & downstream_mRNA_subsequence, const mRNA & upstream_mRNA_subsequence, seedType seed_type)
    {
      return 0.58134;
}

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
    downregulationScore miRNA::calculate_three_prime_feature(const alignment & the_alignment)
    {
      return 2.32868;
}

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
    downregulationScore miRNA::calculate_UTR_dist_feature(const mRNA & the_mRNA, chromosomePosition predicted_three_prime_position, seedType seed_type)
    {
       /**********************************************\
| Define the regression parameters as in mirSVR: |
       \**********************************************/
      const downregulationScore intercept = -0.07;
      const downregulationScore     slope =  0.000172;
       /*****************************************\
| Define the distance cutoff to not penalty |
      | long UTRs too much as in  mirSVR:         |
       \*****************************************/
      const sequenceLength distance_cutoff = 1500;
       /*****************************************************\ 
      | Get sequence position of predicted 3' target site end |
      | and calculate distances to both sequence ends:        |
       \*****************************************************/
      const sequencePosition three_prime_position = the_mRNA.get_nucleotide_chr(predicted_three_prime_position)->get_sequence_position();
      const sequenceLength    five_prime_distance = the_mRNA.get_length()-three_prime_position;
      const sequenceLength   three_prime_distance = three_prime_position-((seed_type==eightMer || seed_type==sevenMerMEight) ? 9 : 8);
       /***************************************************************\ 
      | Calculate the minimal distance (or cut of) and apply regression |
      | function before returning the final feature score:              |
       \***************************************************************/
      const sequenceLength min_distance(std::min(three_prime_distance,five_prime_distance));
      return std::min(min_distance,distance_cutoff) * slope + intercept;
}

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
    void miRNA::calculate_seed_match_features(downregulationScore features[], const alignment & the_alignment)
    {
       /******************************************************************\ 
      | The first feature (1A) is 1 if the mRNA 3' end of the site is      |
      | Adenine, the other 8 features (m2-m9) are 1 if their corresponding |
      | alignment column (counting from miRNA 5' or mRNA 3', respectively) |
      | contain a match. Otherwise the features are 0:                     |
       \******************************************************************/
      features[0]=the_alignment.begin()->get_mRNA_nucleotide().get_base()==Adenine;
      for(unsigned int alignment_position=1;alignment_position<10;++alignment_position)
      {
        features[alignment_position]=(the_alignment.begin()+alignment_position)->get_match().get_identifier()==Match;
      }
}


} // namespace microSNPscore
