#include <string>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include "mRNA.h"
#include "miRNA.h"
#include "sequenceFile.h"
#include "alignment.h"
#include "SNP.h"
#include "conservationList.h"
#include "filepath.h"

using namespace microSNPscore;

int main(){
  miRNA miR195("hsa-miR195","aGCTTCCCUGGCUCUAGCAGCACAGAAAUAUUGGCACAGGGAAGCGAGUCUGCCAAUAUUGGCUGUGCUGCUCCAGGCAGGGUGGUG","chr17",Plus,"6920934","6921020");
  mRNA BCL2("BCL2","AGUCAACAUGCCUGCCCCAAACAAAUAUGCAAAAGGUUCACUAAAGCAGUAGAAAUAAUAUGCAUUGUCAGUGAUGUACCAUGAAACAAAGCUGCAGGCUGUUUAAGAAAAAAUAACACACAUAUAAACAU\
CACACACACAGACAGACACACACACACACAACAAUUAACAGUCUUCAGGCAAAACGUCGAAUCAGCUAUUUACUGCCAAAGGGAAAUAUCAUUUAUUUUUUACAUUAUUAAGAAAAAAAGAUUUAUUUAUUUAAGACAGUCCCAUCAAAAC\
UCCUGUCUUUGGAAAUCCGACCACUAAUUGCCAAGCACCGCUUCGUGUGGCUCCACCUGGAUGUUCUGUGCCUGUAAACAUAGAUUCGCUUUCCAUGUUGUUGGCCGGAUCACCAUCUGAAGAGCAGACGGAUGGAAAAAGGACCUGAUCA\
UUGGGGAAGCUGGCUUUCUGGCUGCUGGAGGCUGGGGAGAAGGUGUUCAUUCACUUGCAUUUCUUUGCCCUGGGGGCUGUGAUAUUAACAGAGGGAGGGUUCCUGUGGGGGGAAGUCCAUGCCUCCCUGGCCUGAAGAAGAGACUCUUUGC\
AUAUGACUCACAUGAUGCAUACCUGGUGGGAGGAAAAGAGUUGGGAACUUCAGAUGGACCUAGUACCCACUGAGAUUUCCACGCCGAAGGACAGCGAUGGGAAAAAUGCCCUUAAAUCAUAGGAAAGUAUUUUUUUAAGCUACCAAUUGUG\
CCGAGAAAAGCAUUUUAGCAAUUUAUACAAUAUCAUCCAGUACCUUAAGCCCUGAUUGUGUAUAUUCAUAUAUUUUGGAUACGCACCCCCCAACUCCCAAUACUGGCUCUGUCUGAGUAAGAAACAGAAUCCUCUGGAACUUGAGGAAGUG\
AACAUUUCGGUGACUUCCGCAUCAGGAAGGCUAGAGUUACCCAGAGCAUCAGGCCGCCACAAGUGCCUGCUUUUAGGAGACCGAAGUCCGCAGAACCUGCCUGUGUCCCAGCUUGGAGGCCUGGUCCUGGAACUGAGCCGGGGCCCUCACU\
GGCCUCCUCCAGGGAUGAUCAACAGGGCAGUGUGGUCUCCGAAUGUCUGGAAGCUGAUGGAGCUCAGAAUUCCACUGUCAAGAAAGAGCAGUAGAGGGGUGUGGCUGGGCCUGUCACCCUGGGGCCCUCCAGGUAGGCCCGUUUUCACGUG\
GAGCAUGGGAGCCACGACCCUUCUUAAGACAUGUAUCACUGUAGAGGGAAGGAACAGAGGCCCUGGGCCCUUCCUAUCAGAAGGACAUGGUGAAGGCUGGGAACGUGAGGAGAGGCAAUGGCCACGGCCCAUUUUGGCUGUAGCACAUGGC\
ACGUUGGCUGUGUGGCCUUGGCCCACCUGUGAGUUUAAAGCAAGGCUUUAAAUGACUUUGGAGAGGGUCACAAAUCCUAAAAGAAGCAUUGAAGUGAGGUGUCAUGGAUUAAUUGACCCCUGUCUAUGGAAUUACAUGUAAAACAUUAUCU\
UGUCACUGUAGUUUGGUUUUAUUUGAAAACCUGACAAAAAAAAAGUUCCAGGUGUGGAAUAUGGGGGUUAUCUGUACAUCCUGGGGCAUUAAAAAAAAAAUCAAUGGUGGGGAACUAUAAAGAAGUAACAAAAGAAGUGACAUCUUCAGCA\
AAUAAACUAGGAAAUUUUUUUUUCUUCCAGUUUAGAAUCAGCCUUGAAACAUUGAUGGAAUAACUCUGUGGCAUUAUUGCAUUAUAUACCAUUUAUCUGUAUUAACUUUGGAAUGUACUCUGUUCAAUGUUUAAUGCUGUGGUUGAUAUUU\
CGAAAGCUGCUUUAAAAAAAUACAUGCAUCUCAGCGUUUUUUUGUUUUUAAUUGUAUUUAGUUAUGGCCUAUACACUAUUUGUGAGCAAAGGUGAUCGUUUUCUGUUUGAGAUUUUUAUCUCUUGAUUCUUCAAAAGCAUUCUGAGAAGGU\
GAGAUAAGCCCUGAGUCUCAGCUACCUAAGAAAAACCUGGAUGUCACUGGCCACUGAGGAGCUUUGUUUCAACCAAGUCAUGUGCAUUUCCACGUCAACAGAAUUGUUUAUUGUGACAGUUAUAUCUGUUGUCCCUUUGACCUUGUUUCUU\
GAAGGUUUCCUCGUCCCUGGGCAAUUCCGCAUUUAAUUCAUGGUAUUCAGGAUUACAUGCAUGUUUGGUUAAACCCAUGAGAUUCAUUCAGUUAAAAAUCCAGAUGGCAAAUGACCAGCAGAUUCAAAUCUAUGGUGGUUUGACCUUUAGA\
GAGUUGCUUUACGUGGCCUGUUUCAACACAGACCCACCCAGAGCCCUCCUGCCCUCCUUCCGCGGGGGCUUUCUCAUGGCUGUCCUUCAGGGUCUUCCUGAAAUGCAGUGGUGCUUACGCUCCACCAAGAAAGCAGGAAACCUGUGGUAUG\
AAGCCAGACCUCCCCGGCGGGCCUCAGGGAACAGAAUGAUCAGACCUUUGAAUGAUUCUAAUUUUUAAGCAAAAUAUUAUUUUAUGAAAGGUUUACAUUGUCAAAGUGAUGAAUAUGGAAUAUCCAAUCCUGUGCUGCUAUCCUGCCAAAA\
UCAUUUUAAUGGAGUCAGUUUGCAGUAUGCUCCACGUGGUAAGAUCCUCCAAGCUGCUUUAGAAGUAACAAUGAAGAACGUGGACGUUUUUAAUAUAAAGCCUGUUUUGUCUUUUGUUGUUGUUCAAACGGGAUUCACAGAGUAUUUGAAA\
AAUGUAUAUAUAUUAAGAGGUCACGGGGGCUAAUUGCUGGCUGGCUGCCUUUUGCUGUGGGGUUUUGUUACCUGGUUUUAAUAACAGUAAAUGUGCCCAGCCUCUUGGCCCCAGAACUGUACAGUAUUGUGGCUGCACUUGCUCUAAGAGU\
AGUUGAUGUUGCAUUUUCCUUAUUGUUAAAAACAUGUUAGAAGCAAUGAAUGUAUAUAAAAGCCUCAACUAGUCAUUUUUUUCUCCUCUUCUUUUUUUUCAUUAUAUCUAAUUAUUUUGCAGUUGGGCAACAGAGAACCAUCCCUAUUUUG\
UAUUGAAGAGGGAUUCACAUCUGCAUCUUAACUGCUCUUUAUGAAUGAAAAAACAGUCCUCUGUAUGUACUCCUCUUUACACUGGCCAGGGUCAGAGUUAAAUAGAGUAUAUGCACUUUCCAAAUUGGGGACAAGGGCUCUAAAAAAAGCC\
CCAAAAGGAGAAGAACAUCUGAGAACCUCCUCGGCCCUCCCAGUCCCUCGCUGCACAAAUACUCCGCAAGAGAGGCCAGAAUGACAGCUGACAGGGUCUAUGGCCAUCGGGUCGUCUCCGAAGAUUUGGCAGGGGCAGAAAACUCUGGCAG\
GCUUAAGAUUUGGAAUAAAGUCACAGAAUUAAGGAAGCACCUCAAUUUAGUUCAAACAAGACGCCAACAUUCUCUCCACAGCUCACUUACCUCUCUGUGUUCAGAUGUGGCCUUCCAUUUAUAUGUGAUCUUUGUUUUAUUAGUAAAUGCU\
UAUCAUCUAAAGAUGUAGCUCUGGCCCAGUGGGAAAAAUUAGGAAGUGAUUAUAAAUCGAGAGGAGUUAUAAUAAUCAAGAUUAAAUGUAAAUAAUCAGGGCAAUCCCAACACAUGUCUAGCUUUCACCUCCAGGAUCUAUUGAGUGAACA\
GAAUUGCAAAUAGUCUCUAUUUGUAAUUGAACUUAUCCUAAAACAAAUAGUUUAUAAAUGUGAACUUAAACUCUAAUUAAUUCCAACUGUACUUUUAAGGCAGUGGCUGUUUUUAGACUUUCUUAUCACUUAUAGUUAGUAAUGUACACCU\
ACUCUAUCAGAGAAAAACAGGAAAGGCUCGAAAUACAAGCCAUUCUAAGGAAAUUAGGGAGUCAGUUGAAAUUCUAUUCUGAUCUUAUUCUGUGGUGUCUUUUGCAGCCCAGACAAAUGUGGUUACACACUUUUUAAGAAAUACAAUUCUA\
CAUUGUCAAGCUUAUGAAGGUUCCAAUCAGAUCUUUAUUGUUAUUCAAUUUGGAUCUUUCAGGGAUUUUUUUUUUAAAUUAUUAUGGGACAAAGGACAUUUGUUGGAGGGGUGGGAGGGAGGAAGAAUUUUUAAAUGUAAAACAUUCCCAA\
GUUUGGAUCAGGGAGUUGGAAGUUUUCAGAAUAACCAGAACUAAGGGUAUGAAGGACCUGUAUUGGGGUCGAUGUGAUGCCUCUGCGAAGAACCUUGUGUGACAAAUGAGAAACAUUUUGAAGUUUGUGGUACGACCUUUAGAUUCCAGAG\
ACAUCAGCAUGGCUCAAAGUGCAGCUCCGUUUGGCAGUGCAAUGGUAUAAAUUUCAAGCUGGAUAUGUCUAAUGGGUAUUUAAACAAUAAAUGUGCAGUUUUAACUAACAGGAUAUUUAAUGACAACCUUCUGGUUGGUAGGGACAUCUGU\
UUCUAAAUGUUUAUUAUGUACAAUACAGAAAAAAAUUUUAUAAAAUUAAGCAAUGUGAAACUGAAUUGGAGAGUGAUAAUACAAGUCCUUUAGUCUUACCCAGUGAAUCAUUCUGUUCCAUGUCUUUGGACAACCAUGACCUUGGACAAUC\
AUGAAAUAUGCAUCUCACUGGAUGCAAAGAAAAUCAGAUGGAGCAUGAAUGGUACUGUACCGGUUCAUCUGGACUGCCCCAGAAAAAUAACUUCAAGCAAACAUCCUAUCAACAACAAGGUUGUUCUGCAUACCAAGCUGAGCACAGAAGA\
UGGGAACACUGGUGGAGGAUGGAAAGGCUCGCUCAAUCAAGAAAAUUCUGAGACUAUUAAUAAAUAAGACUGUAGUGUAGAUACUGAGUAAAUCCAUGCACCUAAACCUUUUGGAAAAUCUGCCGUGGGCCCUCCAGAUAGCUCAUUUCAU\
UAAGUUUUUCCCUCCAAGGUAGAAUUUGCAAGAGUGACAGUGGAUUGCAUUUCUUUUGGGGAAGCUUUCUUUUGGUGGUUUUGUUUAUUAUACCUUCUUAAGUUUUCAACCAAGGUUUGCUUUUGUUUUGAGUUACUGGGGUUAUUUUUGU\
UUUAAAUAAAAAUAAGUGUACAAUAAGUGUUUUUGUAUUGAAAGCUUUUGUUAUCAAGAUUUUCAUACUUUUACCUUCCAUGGCUCUUUUUAAGAUUGAUACUUUUAAGAGGUGGCUGAUAUUCUGCAACACUGUACACAUAAAAAAUACG\
GUAAGGAUACUUUACAUGGUUAAGGUAAAGUAAGUCUCCAGUUGGCCACCAUUAGCUAUAAUGGCACUUUGUUUGUGUUGUUGGAAAAAGUCACAUUGCCAUUAAACUUUCCUUGUCUGUCUAGUUAAUAUUGUGAAGAAAAAUAAAGUAC\
AGUGUGAGAUACUG","chr18",Minus,"60790579","60795857");
  std::cout << "Length of miR195: " << miR195.get_length() << std::endl;
  std::cout << "BCL2 is transcribed from the " << BCL2.get_strand() << " strand.\n";
  std::cout << "First nucleotide of BCL2 is " << BCL2.begin()->get_base() << " (" << *(BCL2.begin());
  std::cout << ") and is located at position " << BCL2.begin()->get_chromosome_position();
  std::cout << " on the chromosome " << BCL2.get_chromosome() << " which corresponds to sequence position ";
  std::cout << BCL2.get_nucleotide_chr(BCL2[1]->get_chromosome_position())->get_sequence_position() << ".\n";
  std::cout << "Under no special conditions aligning it with the third nucleotide of miR195, which is a ";
  std::cout << miR195[3]->get_base() << ", would result in a score of ";
  std::cout << BCL2.begin()->get_match(*(miR195[3])).get_score() << " because it is a ";
  std::cout << BCL2.begin()->get_match(*(miR195[3])).get_identifier() << ".\n";
  std::cout << "The 1st nucleotide of miR195 should be on chromosome position 6920934 and is on position ";
  std::cout << miR195.begin()->get_chromosome_position() << ".\n";
  std::cout << "And of course chromosome position 6920934 should be position 1 in miR195 and is ";
  std::cout << miR195.get_nucleotide_chr(6920934)->get_sequence_position() << ".\n";
  std::cout << "At least chromosome position 6920934 should be position 6920934 on chromosome and is ";
  std::cout << miR195.get_nucleotide_chr(6920934)->get_chromosome_position() << ".\n";
  std::cout << "The first exon of miR195 starts at chromosome position ";
  std::cout << miR195.exons_begin()->get_start() << ".\n";
  std::cout << "The 27th nucleotide of BCL2 is at chromosome position ";
  std::cout << BCL2[27]->get_chromosome_position() << " and on sequence position ";
  std::cout << BCL2[27]->get_sequence_position() << ".\n";
  std::cout << "The first nucleotide of the 10 nucleotide BCL2 subsequence starting at that chromosome position is on chromosome position ";
  std::cout << BCL2.get_subsequence_chr_from(BCL2[27]->get_chromosome_position(),10).begin()->get_chromosome_position() << " and on position ";
  std::cout << BCL2.get_subsequence_chr_from(BCL2[27]->get_chromosome_position(),10).begin()->get_sequence_position() << " in sequence.\n";
  std::cout << "This sequence has a length of ";
  std::cout << BCL2.get_subsequence_chr_from(BCL2[27]->get_chromosome_position(),10).get_length() << " and its last exon ends at position ";
  std::cout << (BCL2.get_subsequence_chr_from(BCL2[27]->get_chromosome_position(),10).exons_end()-1)->get_end() << " on chromosome which is the ";
  std::cout << BCL2.get_subsequence_chr_from(BCL2[27]->get_chromosome_position(),10).get_nucleotide_chr((BCL2.get_subsequence_chr_from(BCL2[27]->get_chromosome_position(),10).exons_end()-1)->get_end())->get_sequence_position();
  std::cout << "th nucleotide of the subsequence while its first exon starts at chromosome position ";
  std::cout << BCL2.get_subsequence_chr_from(BCL2[27]->get_chromosome_position(),10).exons_begin()->get_start() << ".\n";
  std::cout << "If we take position 60793322 on chromosome to be reported from the target prediction tool we get ";
  std::cout << BCL2.get_subsequence_for_alignment(60793322) << " to be the subsequence of interest for the alignment, ";
  std::cout << BCL2.get_subsequence_for_accessability(60793322) << " to be the subsequence of interest for the accessability score calculation, ";
  std::cout << BCL2.get_subsequence_for_downstream_AU_content(60793322) << " to be the subsequence for the downstream AU content calculation and ";
  std::cout << BCL2.get_subsequence_for_upstream_AU_content(60793322) << " to be the subsequence for the upstream AU content calculation.\n";
  std::cout << "The ID of the alignment relevant BCL2 subsequence is " << BCL2.get_subsequence_for_alignment(60793322).get_ID() << ".\n";
  std::string FASTA_test(">test|no1|103,160,150|165,135,270|-1|my_chrom\n\
AAAAAAACAUACdACGaaaUAAAAAAAAAAAAACCCCCCCCCCCCCCCUUUUUGGXGGGC\n\
UAGUCGUACCAGUCAUGACGUGUCUGCAGUUUCACCGUCGUACUACGUACGUCUGCUGCU\n\
caugaucGUCGAUACCAGUAXGGGGGGGGGGG\n");
  std::cout << "If we put the following FASTA entry:\n" << FASTA_test << " into a sequence file entry and print it we get:\n";
  std::cout << sequenceFileEntry(FASTA_test);
  mRNA test(sequenceFileEntry(FASTA_test).get_mRNA());
  std::cout << "If we convert that sequence file entry to a mRNA object and then convert that object back to a sequence file entry we get:\n";
  std::cout << sequenceFileEntry(test);
  std::cout << "The next test is to create a sequence file \"test.seq\"";
  sequenceFile test_file("test.seq");
  std::cout << ", append BCL2";
  test_file.add_sequence(BCL2);
  std::cout << " and miR195";
  test_file.add_sequence(miR195);
  std::cout << " and write it to the hard disk.\n";
  test_file.write();
  std::cout << "Now we create another sequence file with the same file name";
  sequenceFile scnd_file("test.seq");
  std::cout << " and add the test mRNA from above";
  scnd_file.add_sequence(test);
  std::cout << " before we read the entries already on disk";
  scnd_file.read();
  std::cout << " and remove the file (no test - just cleanup).\n";
  remove("test.seq");
  std::cout << "If we would write the second sequence file's content to the file it would look like this:\n";
  for(sequenceFile::const_iterator entry_it(scnd_file.begin());entry_it!=scnd_file.end();++entry_it)
  {
    std::cout << *entry_it;
  }
  std::cout << "If we align the relevant BCL2 substring with the mature miR195 we get:\n";
  miRNA mature195(sequenceFileEntry(">mature195|6920949|6920969|1|17\nUAGCAGCACAGAAAUAUUGGC\n").get_miRNA());
  std::cout << optimalAlignmentList(BCL2.get_subsequence_for_alignment(60793322),mature195);
  std::cout << "\nNow it is time to introduce the SNP rs4987856.\n";
  SNP rs4987856("rs4987856","G","A","chr18",Minus,60793494);
  std::cout << "It is located on chromosome " << rs4987856.get_chromosome() << " at position " << rs4987856.get_position(Minus);
  std::cout << " and has a shift length of " << rs4987856.get_shift() << ".\n";
  std::cout << "The BCL nucleotide at that position is " << *BCL2.get_nucleotide_chr(rs4987856.get_position(BCL2.get_strand())) << ".\n";
  std::cout << "rs4987856 will turn " << *rs4987856.reference_begin(BCL2.get_strand());
  std::cout <<  " into " << *rs4987856.alternative_begin(BCL2.get_strand()) << ".\n";
  std::cout << "Thus rs4987856 " << (rs4987856.matches(BCL2) ? "matches" : "does not match") << " BCL2.\n";
  std::cout << "And if we mutate BCL2 with rs4987856 we get:\n";
  std::cout << sequenceFileEntry(BCL2.mutate(rs4987856));
  std::cout << "As you might see (or not ^^) the position of interest now contains ";
  std::cout << *BCL2.mutate(rs4987856).get_nucleotide_chr(rs4987856.get_position(BCL2.mutate(rs4987856).get_strand())) << ".\n";
//  std::cout << "The deregulation score of rs4987856 for the miR195 target site in BCL2 starting at position 60793322 on the chromosome is ";
//  std::cout << rs4987856.get_deregulation_score(mature195,BCL2,60793322,/* conservations */) << ".\n";
  std::cout << "\n\n====================================================================================================\n\n";
  std::cout << "For a more comprehensive test of the deregulation score calculation we need more examples testing more cases:\n\n";
  std::cout << "We start with a made up miRNA we call miR-0815:\n\n";
  std::string FASTA0815(">miR-0815|11|31|1|42\naACCGUGAagucaucaccagc\n");
  std::cout << FASTA0815;
  miRNA miR0815(sequenceFileEntry(FASTA0815).get_miRNA());
  std::cout << "\nThe seed is marked in upper case\n";
  std::cout << "\nIn addition to that miRNA we make up a mRNA we call gene-4711:\n\n";
  std::string FASTA4711(">gene-4711|21,242|129,360|-1|007\n\
acguagcgguacgucacacacacugguacuacuacagugcacacaguguacugaugcuac\n\
cggggaucguacuacuacugacuuacgacuacgacguguacggcuagcauccccgugauu\n\
uucacgguAcaucagucuagcgcgagagagaucuucucagcuagcugacuagcugaucgu\n\
agcuagcugacuagcguagcuacguagcuagucagucgaugcuagcga\n");
  std::cout << FASTA4711;
  mRNA gene4711(sequenceFileEntry(FASTA4711).get_mRNA());
  std::cout << "\nThe base marked in upper case is chromosome position 120 and we take that position to be reported as 3' end of the target site.\n";
  std::cout << "\nThis leads to the following optimal alignment(s) of miR-0815 and gene-4711:\n";
  std::cout << optimalAlignmentList(gene4711.get_subsequence_for_alignment(120),miR0815);
  std::vector<SNP> SNPs;
  std::cout << "\nThe first SNP we use for our test is a big deletion from position 25 to 100 we call del-25-100.\n";
  SNP del25to100("del-25-100","agaucuucucagcuagcugacuagcugaucguagcuagcugacuagcguagcuacguagcuagucagucgaugcua","","007",Minus,100);
  SNPs.push_back(del25to100);
  std::cout << "This SNP should simulate the impact of the distance to the UTR 5' end on the downregulation.\n";
  std::cout << "\nEven though we expect no changes we align miR-0815 and gene-4711:del-25-100 to test the 3' position shift:\n";
  std::cout << optimalAlignmentList(gene4711.mutate(del25to100).get_subsequence_for_alignment(120+del25to100.get_shift()),miR0815);
  std::cout << "\nThe next SNP we use is a big insertion of 200 Guanine (+ strand) after position 300 we call ins-200G-300.\n";
  SNP ins200Gat300("ins-200G-300","T","TGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG","007",Plus,300);
  std::cout << "This SNP should simulate the impact of the distance to the UTR 3' end on the downregulation.\n";
  SNPs.push_back(ins200Gat300);
  std::cout << "\nAgain we expect no changes for the alignment(s) for miR-0815 and gene-4711:ins-200G-300 but this time there should be no shift:\n";
  std::cout << optimalAlignmentList(gene4711.mutate(ins200Gat300).get_subsequence_for_alignment(120),miR0815);
  std::cout << "\nThe next SNP we use is a point mutation from Adenine to Uracil (- strand) at position 120 we call mut-A->U-120.\n";
  SNP mutAtoUat120("mut-A->U-120","A","U","007",Minus,120);
  SNPs.push_back(mutAtoUat120);
  std::cout << "This SNP should simulate the impact of the UTR 3' end of the target site beeing Adenin to the UTR 3' end ";
  std::cout << "on the downregulation as well as the lack of impact of a match at that position.\n";
  std::cout << "\nAgain the optimal alignment(s) for miR-0815 and gene-4711:mut-A->U-120:\n";
  std::cout << optimalAlignmentList(gene4711.mutate(mutAtoUat120).get_subsequence_for_alignment(120),miR0815);
  std::cout << "\nThe next SNP we use is a point mutation from Adenine to Uracil (- strand) at position 125 we call mut-A->U-125.\n";
  SNP mutAtoUat125("mut-A->U-125","A","U","007",Minus,125);
  SNPs.push_back(mutAtoUat125);
  std::cout << "This SNP should simulate the impact of a seed mismatch on the downregulation.\n";
  std::cout << "\nAgain the optimal alignment(s) for miR-0815 and gene-4711:mut-A->U-125:\n";
  std::cout << optimalAlignmentList(gene4711.mutate(mutAtoUat125).get_subsequence_for_alignment(120),miR0815);
  std::cout << "\nThe next SNP we use is a point mutation from Thymine to Cytosine (+ strand) at position 125 we call mut-T->C-125.\n";
  SNP mutTtoCat125("mut-T->C-125","T","C","007",Plus,125);
  SNPs.push_back(mutTtoCat125);
  std::cout << "This SNP should simulate the impact of a seed wobble pair on the downregulation.\n";
  std::cout << "\nAgain the optimal alignment(s) for miR-0815 and gene-4711:mut-T->C-125:\n";
  std::cout << optimalAlignmentList(gene4711.mutate(mutTtoCat125).get_subsequence_for_alignment(120),miR0815);
  std::cout << "\nThe next SNP we use is a copy number variation from one to three Adenine (+ strand) at position 128 we call cnv-A-1->3-128.\n";
  SNP cnvA1to3at128("cnv-A-1->3-128","A","AAA","007",Plus,128);
  SNPs.push_back(cnvA1to3at128);
  std::cout << "This SNP should simulate the impact of removing the gap behind the seed match on the downregulation.\n";
  std::cout << "\nAgain the optimal alignment(s) for miR-0815 and gene-4711:cnv-A-1->3-128:\n";
  std::cout << optimalAlignmentList(gene4711.mutate(cnvA1to3at128).get_subsequence_for_alignment(120),miR0815);
  std::cout << "\nThe next SNP we use is a point deletion at position 129 we call del-129.\n";
  SNP del129("del-129","T","","007",Minus,129);
  SNPs.push_back(del129);
  std::cout << "This SNP should simulate the impact of moving the 3' match.\n";
  std::cout << "\nAgain the optimal alignment(s) for miR-0815 and gene-4711:del-129:\n";
  std::cout << optimalAlignmentList(gene4711.mutate(del129).get_subsequence_for_alignment(120),miR0815);
  std::cout << "\nThe next SNP we use is a mutation from GAU to G (just to test a weird way to describe a deletion) (- strand) ";
  std::cout << "at position 244 we call mut-GAU->G-244.\n";
  SNP mutGAUtoGat244("mut-GAU->G-244","GAU","G","007",Minus,244);
  SNPs.push_back(mutGAUtoGat244);
  std::cout << "This SNP should simulate the impact of shortening the 3' match.\n";
  std::cout << "\nAgain the optimal alignment(s) for miR-0815 and gene-4711:mut-GAU->G-244:\n";
  std::cout << optimalAlignmentList(gene4711.mutate(mutGAUtoGat244).get_subsequence_for_alignment(120),miR0815);
  std::cout << "\nThe last SNP we use is a mutation from GG to CA (+ strand) at position 247 we call mut-GG->CA-247.\n";
  SNP mutGGtoCAat247("mut-GG->CA-247","GG","CA","007",Plus,247);
  SNPs.push_back(mutGGtoCAat247);
  std::cout << "This SNP should simulate the impact of extending the 3' match.\n";
  std::cout << "\nAgain the optimal alignment(s) for miR-0815 and gene-4711:mut-GG->CA-247:\n";
  std::cout << optimalAlignmentList(gene4711.mutate(mutGGtoCAat247).get_subsequence_for_alignment(120),miR0815);
  std::cout << "\nIn addition to that data we need conservation information for our scoring algorithm.\n";
  std::cout << "Therefor we create a conservation ranges file with the following content:\n";
  std::string conservation_file_content("\
001	1	0.12\n\
001	3	0.4\n\
001	1213	0.82\n\
001	1278	0.32\n\
007	1	0.32\n\
007	56	0.45\n\
007	120	0.83\n\
007	130	0.53\n\
007	200	0.23\n\
007	240	0.001\n\
007	242	0.93\n\
007	251	0.12\n\
chr1	1	0.2\n\
chr1	56	0.645\n");
  filePath conservation_file_path("test.conservations");
  std::ofstream conservation_file_stream(conservation_file_path.c_str());
  conservation_file_stream << conservation_file_content;
  conservation_file_stream.close();
  conservationList conservation_list(conservation_file_path);
  std::cout << std::endl << conservation_list << std::endl;
  std::cout << "\nNow we calculate the deregulation score for those SNPs:\n\n\tSNP\t\t|\tscore\n------------------------+------------------------\n";
  for(std::vector<SNP>::const_iterator SNP_it(SNPs.begin());SNP_it!=SNPs.end();++SNP_it)
  {
    std::cout << "\t" << SNP_it->get_ID() << (SNP_it->get_ID().length() < 8 ? "\t" : "") << "\t|\t";
    std::cout << SNP_it->get_deregulation_score(miR0815,gene4711,120,conservation_list) << std::endl;
  }
  remove(conservation_file_path.c_str());
}
