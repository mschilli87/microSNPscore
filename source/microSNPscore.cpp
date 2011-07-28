#include <string>
#include <iostream>
#include <stdio.h>
#include "mRNA.h"
#include "miRNA.h"
#include "sequenceFile.h"
#include "alignment.h"
#include "SNP.h"

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
  std::cout << "The deregulation score of rs4987856 for the miR195 target site in BCL2 starting at position 60793322 on the chromosome is ";
  std::cout << rs4987856.get_deregulation_score(mature195,BCL2,60793322) << ".\n";
}
