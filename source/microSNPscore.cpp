#include <string>
#include <iostream>
#include "mRNA.h"
#include "miRNA.h"

using namespace microSNPscore;

std::string base2str(nucleoBase b)
{
  switch(b)
  {
    case Adenine: return "Adenine";
    case Uracil: return "Uracil";
    case Cytosine: return "Cytosine";
    case Guanine: return "Guanine";
    case Gap: return "Gap";
    case Mask: return "Mask";
    default:
      std::cerr << "base2str\n";
      std::cerr << " ==> unkown base: ";
      std::cerr << b << std::endl;
      std::cerr << "  --> assuming Mask\n";
      return "Mask";
  }
}

char base2chr(nucleoBase b)
{
  switch(b)
  {
    case Adenine: return 'A';
    case Uracil: return 'U';
    case Cytosine: return 'C';
    case Guanine: return 'G';
    case Gap: return '-';
    case Mask: return 'X';
    default:
      std::cerr << "base2chr\n";
      std::cerr << " ==> unkown base: ";
      std::cerr << b << std::endl;
      std::cerr << "  --> assuming Mask\n";
      return 'X';
  }
}

std::ostream & operator<<(std::ostream & the_stream,const nucleoBase & the_base)
{
  return the_stream << base2str(the_base);
}

int main(){
  miRNA miR195("aGCTTCCCUGGCUCUAGCAGCACAGAAAUAUUGGCACAGGGAAGCGAGUCUGCCAAUAUUGGCUGUGCUGCUCCAGGCAGGGUGGUG","chr17",Plus,"6920934","6921020");
  mRNA BCL2("AGUCAACAUGCCUGCCCCAAACAAAUAUGCAAAAGGUUCACUAAAGCAGUAGAAAUAAUAUGCAUUGUCAGUGAUGUACCAUGAAACAAAGCUGCAGGCUGUUUAAGAAAAAAUAACACACA\
UAUAAACAUCACACACACAGACAGACACACACACACACAACAAUUAACAGUCUUCAGGCAAAACGUCGAAUCAGCUAUUUACUGCCAAAGGGAAAUAUCAUUUAUUUUUUACAUUAUUAAGAAAAAAAGAUUUAUUUAU\
UUAAGACAGUCCCAUCAAAACUCCUGUCUUUGGAAAUCCGACCACUAAUUGCCAAGCACCGCUUCGUGUGGCUCCACCUGGAUGUUCUGUGCCUGUAAACAUAGAUUCGCUUUCCAUGUUGUUGGCCGGAUCACCAUCU\
GAAGAGCAGACGGAUGGAAAAAGGACCUGAUCAUUGGGGAAGCUGGCUUUCUGGCUGCUGGAGGCUGGGGAGAAGGUGUUCAUUCACUUGCAUUUCUUUGCCCUGGGGGCUGUGAUAUUAACAGAGGGAGGGUUCCUGU\
GGGGGGAAGUCCAUGCCUCCCUGGCCUGAAGAAGAGACUCUUUGCAUAUGACUCACAUGAUGCAUACCUGGUGGGAGGAAAAGAGUUGGGAACUUCAGAUGGACCUAGUACCCACUGAGAUUUCCACGCCGAAGGACAG\
CGAUGGGAAAAAUGCCCUUAAAUCAUAGGAAAGUAUUUUUUUAAGCUACCAAUUGUGCCGAGAAAAGCAUUUUAGCAAUUUAUACAAUAUCAUCCAGUACCUUAAGCCCUGAUUGUGUAUAUUCAUAUAUUUUGGAUAC\
GCACCCCCCAACUCCCAAUACUGGCUCUGUCUGAGUAAGAAACAGAAUCCUCUGGAACUUGAGGAAGUGAACAUUUCGGUGACUUCCGCAUCAGGAAGGCUAGAGUUACCCAGAGCAUCAGGCCGCCACAAGUGCCUGC\
UUUUAGGAGACCGAAGUCCGCAGAACCUGCCUGUGUCCCAGCUUGGAGGCCUGGUCCUGGAACUGAGCCGGGGCCCUCACUGGCCUCCUCCAGGGAUGAUCAACAGGGCAGUGUGGUCUCCGAAUGUCUGGAAGCUGAU\
GGAGCUCAGAAUUCCACUGUCAAGAAAGAGCAGUAGAGGGGUGUGGCUGGGCCUGUCACCCUGGGGCCCUCCAGGUAGGCCCGUUUUCACGUGGAGCAUGGGAGCCACGACCCUUCUUAAGACAUGUAUCACUGUAGAG\
GGAAGGAACAGAGGCCCUGGGCCCUUCCUAUCAGAAGGACAUGGUGAAGGCUGGGAACGUGAGGAGAGGCAAUGGCCACGGCCCAUUUUGGCUGUAGCACAUGGCACGUUGGCUGUGUGGCCUUGGCCCACCUGUGAGU\
UUAAAGCAAGGCUUUAAAUGACUUUGGAGAGGGUCACAAAUCCUAAAAGAAGCAUUGAAGUGAGGUGUCAUGGAUUAAUUGACCCCUGUCUAUGGAAUUACAUGUAAAACAUUAUCUUGUCACUGUAGUUUGGUUUUAU\
UUGAAAACCUGACAAAAAAAAAGUUCCAGGUGUGGAAUAUGGGGGUUAUCUGUACAUCCUGGGGCAUUAAAAAAAAAAUCAAUGGUGGGGAACUAUAAAGAAGUAACAAAAGAAGUGACAUCUUCAGCAAAUAAACUAG\
GAAAUUUUUUUUUCUUCCAGUUUAGAAUCAGCCUUGAAACAUUGAUGGAAUAACUCUGUGGCAUUAUUGCAUUAUAUACCAUUUAUCUGUAUUAACUUUGGAAUGUACUCUGUUCAAUGUUUAAUGCUGUGGUUGAUAU\
UUCGAAAGCUGCUUUAAAAAAAUACAUGCAUCUCAGCGUUUUUUUGUUUUUAAUUGUAUUUAGUUAUGGCCUAUACACUAUUUGUGAGCAAAGGUGAUCGUUUUCUGUUUGAGAUUUUUAUCUCUUGAUUCUUCAAAAG\
CAUUCUGAGAAGGUGAGAUAAGCCCUGAGUCUCAGCUACCUAAGAAAAACCUGGAUGUCACUGGCCACUGAGGAGCUUUGUUUCAACCAAGUCAUGUGCAUUUCCACGUCAACAGAAUUGUUUAUUGUGACAGUUAUAU\
CUGUUGUCCCUUUGACCUUGUUUCUUGAAGGUUUCCUCGUCCCUGGGCAAUUCCGCAUUUAAUUCAUGGUAUUCAGGAUUACAUGCAUGUUUGGUUAAACCCAUGAGAUUCAUUCAGUUAAAAAUCCAGAUGGCAAAUG\
ACCAGCAGAUUCAAAUCUAUGGUGGUUUGACCUUUAGAGAGUUGCUUUACGUGGCCUGUUUCAACACAGACCCACCCAGAGCCCUCCUGCCCUCCUUCCGCGGGGGCUUUCUCAUGGCUGUCCUUCAGGGUCUUCCUGA\
AAUGCAGUGGUGCUUACGCUCCACCAAGAAAGCAGGAAACCUGUGGUAUGAAGCCAGACCUCCCCGGCGGGCCUCAGGGAACAGAAUGAUCAGACCUUUGAAUGAUUCUAAUUUUUAAGCAAAAUAUUAUUUUAUGAAA\
GGUUUACAUUGUCAAAGUGAUGAAUAUGGAAUAUCCAAUCCUGUGCUGCUAUCCUGCCAAAAUCAUUUUAAUGGAGUCAGUUUGCAGUAUGCUCCACGUGGUAAGAUCCUCCAAGCUGCUUUAGAAGUAACAAUGAAGA\
ACGUGGACGUUUUUAAUAUAAAGCCUGUUUUGUCUUUUGUUGUUGUUCAAACGGGAUUCACAGAGUAUUUGAAAAAUGUAUAUAUAUUAAGAGGUCACGGGGGCUAAUUGCUGGCUGGCUGCCUUUUGCUGUGGGGUUU\
UGUUACCUGGUUUUAAUAACAGUAAAUGUGCCCAGCCUCUUGGCCCCAGAACUGUACAGUAUUGUGGCUGCACUUGCUCUAAGAGUAGUUGAUGUUGCAUUUUCCUUAUUGUUAAAAACAUGUUAGAAGCAAUGAAUGU\
AUAUAAAAGCCUCAACUAGUCAUUUUUUUCUCCUCUUCUUUUUUUUCAUUAUAUCUAAUUAUUUUGCAGUUGGGCAACAGAGAACCAUCCCUAUUUUGUAUUGAAGAGGGAUUCACAUCUGCAUCUUAACUGCUCUUUA\
UGAAUGAAAAAACAGUCCUCUGUAUGUACUCCUCUUUACACUGGCCAGGGUCAGAGUUAAAUAGAGUAUAUGCACUUUCCAAAUUGGGGACAAGGGCUCUAAAAAAAGCCCCAAAAGGAGAAGAACAUCUGAGAACCUC\
CUCGGCCCUCCCAGUCCCUCGCUGCACAAAUACUCCGCAAGAGAGGCCAGAAUGACAGCUGACAGGGUCUAUGGCCAUCGGGUCGUCUCCGAAGAUUUGGCAGGGGCAGAAAACUCUGGCAGGCUUAAGAUUUGGAAUA\
AAGUCACAGAAUUAAGGAAGCACCUCAAUUUAGUUCAAACAAGACGCCAACAUUCUCUCCACAGCUCACUUACCUCUCUGUGUUCAGAUGUGGCCUUCCAUUUAUAUGUGAUCUUUGUUUUAUUAGUAAAUGCUUAUCA\
UCUAAAGAUGUAGCUCUGGCCCAGUGGGAAAAAUUAGGAAGUGAUUAUAAAUCGAGAGGAGUUAUAAUAAUCAAGAUUAAAUGUAAAUAAUCAGGGCAAUCCCAACACAUGUCUAGCUUUCACCUCCAGGAUCUAUUGA\
GUGAACAGAAUUGCAAAUAGUCUCUAUUUGUAAUUGAACUUAUCCUAAAACAAAUAGUUUAUAAAUGUGAACUUAAACUCUAAUUAAUUCCAACUGUACUUUUAAGGCAGUGGCUGUUUUUAGACUUUCUUAUCACUUA\
UAGUUAGUAAUGUACACCUACUCUAUCAGAGAAAAACAGGAAAGGCUCGAAAUACAAGCCAUUCUAAGGAAAUUAGGGAGUCAGUUGAAAUUCUAUUCUGAUCUUAUUCUGUGGUGUCUUUUGCAGCCCAGACAAAUGU\
GGUUACACACUUUUUAAGAAAUACAAUUCUACAUUGUCAAGCUUAUGAAGGUUCCAAUCAGAUCUUUAUUGUUAUUCAAUUUGGAUCUUUCAGGGAUUUUUUUUUUAAAUUAUUAUGGGACAAAGGACAUUUGUUGGAG\
GGGUGGGAGGGAGGAAGAAUUUUUAAAUGUAAAACAUUCCCAAGUUUGGAUCAGGGAGUUGGAAGUUUUCAGAAUAACCAGAACUAAGGGUAUGAAGGACCUGUAUUGGGGUCGAUGUGAUGCCUCUGCGAAGAACCUU\
GUGUGACAAAUGAGAAACAUUUUGAAGUUUGUGGUACGACCUUUAGAUUCCAGAGACAUCAGCAUGGCUCAAAGUGCAGCUCCGUUUGGCAGUGCAAUGGUAUAAAUUUCAAGCUGGAUAUGUCUAAUGGGUAUUUAAA\
CAAUAAAUGUGCAGUUUUAACUAACAGGAUAUUUAAUGACAACCUUCUGGUUGGUAGGGACAUCUGUUUCUAAAUGUUUAUUAUGUACAAUACAGAAAAAAAUUUUAUAAAAUUAAGCAAUGUGAAACUGAAUUGGAGA\
GUGAUAAUACAAGUCCUUUAGUCUUACCCAGUGAAUCAUUCUGUUCCAUGUCUUUGGACAACCAUGACCUUGGACAAUCAUGAAAUAUGCAUCUCACUGGAUGCAAAGAAAAUCAGAUGGAGCAUGAAUGGUACUGUAC\
CGGUUCAUCUGGACUGCCCCAGAAAAAUAACUUCAAGCAAACAUCCUAUCAACAACAAGGUUGUUCUGCAUACCAAGCUGAGCACAGAAGAUGGGAACACUGGUGGAGGAUGGAAAGGCUCGCUCAAUCAAGAAAAUUC\
UGAGACUAUUAAUAAAUAAGACUGUAGUGUAGAUACUGAGUAAAUCCAUGCACCUAAACCUUUUGGAAAAUCUGCCGUGGGCCCUCCAGAUAGCUCAUUUCAUUAAGUUUUUCCCUCCAAGGUAGAAUUUGCAAGAGUG\
ACAGUGGAUUGCAUUUCUUUUGGGGAAGCUUUCUUUUGGUGGUUUUGUUUAUUAUACCUUCUUAAGUUUUCAACCAAGGUUUGCUUUUGUUUUGAGUUACUGGGGUUAUUUUUGUUUUAAAUAAAAAUAAGUGUACAAU\
AAGUGUUUUUGUAUUGAAAGCUUUUGUUAUCAAGAUUUUCAUACUUUUACCUUCCAUGGCUCUUUUUAAGAUUGAUACUUUUAAGAGGUGGCUGAUAUUCUGCAACACUGUACACAUAAAAAAUACGGUAAGGAUACUU\
UACAUGGUUAAGGUAAAGUAAGUCUCCAGUUGGCCACCAUUAGCUAUAAUGGCACUUUGUUUGUGUUGUUGGAAAAAGUCACAUUGCCAUUAAACUUUCCUUGUCUGUCUAGUUAAUAUUGUGAAGAAAAAUAAAGUAC\
AGUGUGAGAUACUG","chr18",Minus,"60790579","60795857");
  std::cout << "Length of miR195: " << miR195.get_length() << std::endl;
  std::cout << "BCL2 is transcribed from the " << (BCL2.get_strand() == Plus ? "+" : "-") << " strand.\n";
  std::cout << "First nucleotide of BCL2 is " << BCL2.begin()->get_base();
  std::cout << " and is located at position " << BCL2.begin()->get_chromosome_position();
  std::cout << " on the chromosome " << BCL2.get_chromosome() << " which corresponds to sequence position ";
  std::cout << BCL2.get_nucleotide_chr(BCL2[1]->get_chromosome_position())->get_sequence_position() << ".\n";
  std::cout << "Under no special conditions aligning it with the third nucleotide of miR195, which is a ";
  std::cout << miR195[3]->get_base() << ", would result in a score of ";
  std::cout << BCL2.begin()->get_match(*(miR195[3])).get_score() << ".\n";
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
  mRNA BCL2_alignment_part(BCL2.get_subsequence_for_alignment(60793322));
  for(sequence::const_iterator sequence_it(BCL2_alignment_part.begin());sequence_it!=BCL2_alignment_part.end();++sequence_it)
  {
    std::cout << base2chr(sequence_it->get_base());
  }
  std::cout << " to be the subsequence of interest for the alignment and ";
  mRNA BCL2_accessability_part(BCL2.get_subsequence_for_accessability(60793322));
  for(sequence::const_iterator sequence_it(BCL2_accessability_part.begin());sequence_it!=BCL2_accessability_part.end();++sequence_it)
  {
    std::cout << base2chr(sequence_it->get_base());
  }
  std::cout << " to be the subsequence of interest for the accessability score calculation.\n";
}
