#include<string>
#include<iostream>
#include "sequence.h"

using namespace microSNPscore;

int main(){
  sequence miR195("aGCTTCCCUGGCUCUAGCAGCACAGAAAUAUUGGCACAGGGAAGCGAGUCUGCCAAUAUUGGCUGUGCUGCUCCAGGCAGGGUGGUG","chr17",Plus,"6920934","6921020");
  std::string BCL2_seq("AGUCAACAUGCCUGCCCCAAACAAAUAUGCAAAAGGUUCACUAAAGCAGUAGAAAUAAUAUGCAUUGUCAGUGAUGUACCAUGAAACAAAGCUGCAGGCUGUUUAAGAAAAAAUA\
ACACACAUAUAAACAUCACACACACAGACAGACACACACACACACAACAAUUAACAGUCUUCAGGCAAAACGUCGAAUCAGCUAUUUACUGCCAAAGGGAAAUAUCAUUUAUUUUUUACAUUAUUAAGAAAAAAAGAUU\
UAUUUAUUUAAGACAGUCCCAUCAAAACUCCUGUCUUUGGAAAUCCGACCACUAAUUGCCAAGCACCGCUUCGUGUGGCUCCACCUGGAUGUUCUGUGCCUGUAAACAUAGAUUCGCUUUCCAUGUUGUUGGCCGGAUC\
ACCAUCUGAAGAGCAGACGGAUGGAAAAAGGACCUGAUCAUUGGGGAAGCUGGCUUUCUGGCUGCUGGAGGCUGGGGAGAAGGUGUUCAUUCACUUGCAUUUCUUUGCCCUGGGGGCUGUGAUAUUAACAGAGGGAGGG\
UUCCUGUGGGGGGAAGUCCAUGCCUCCCUGGCCUGAAGAAGAGACUCUUUGCAUAUGACUCACAUGAUGCAUACCUGGUGGGAGGAAAAGAGUUGGGAACUUCAGAUGGACCUAGUACCCACUGAGAUUUCCACGCCGA\
AGGACAGCGAUGGGAAAAAUGCCCUUAAAUCAUAGGAAAGUAUUUUUUUAAGCUACCAAUUGUGCCGAGAAAAGCAUUUUAGCAAUUUAUACAAUAUCAUCCAGUACCUUAAGCCCUGAUUGUGUAUAUUCAUAUAUUU\
UGGAUACGCACCCCCCAACUCCCAAUACUGGCUCUGUCUGAGUAAGAAACAGAAUCCUCUGGAACUUGAGGAAGUGAACAUUUCGGUGACUUCCGCAUCAGGAAGGCUAGAGUUACCCAGAGCAUCAGGCCGCCACAAG\
UGCCUGCUUUUAGGAGACCGAAGUCCGCAGAACCUGCCUGUGUCCCAGCUUGGAGGCCUGGUCCUGGAACUGAGCCGGGGCCCUCACUGGCCUCCUCCAGGGAUGAUCAACAGGGCAGUGUGGUCUCCGAAUGUCUGGA\
AGCUGAUGGAGCUCAGAAUUCCACUGUCAAGAAAGAGCAGUAGAGGGGUGUGGCUGGGCCUGUCACCCUGGGGCCCUCCAGGUAGGCCCGUUUUCACGUGGAGCAUGGGAGCCACGACCCUUCUUAAGACAUGUAUCAC\
UGUAGAGGGAAGGAACAGAGGCCCUGGGCCCUUCCUAUCAGAAGGACAUGGUGAAGGCUGGGAACGUGAGGAGAGGCAAUGGCCACGGCCCAUUUUGGCUGUAGCACAUGGCACGUUGGCUGUGUGGCCUUGGCCCACC\
UGUGAGUUUAAAGCAAGGCUUUAAAUGACUUUGGAGAGGGUCACAAAUCCUAAAAGAAGCAUUGAAGUGAGGUGUCAUGGAUUAAUUGACCCCUGUCUAUGGAAUUACAUGUAAAACAUUAUCUUGUCACUGUAGUUUG\
GUUUUAUUUGAAAACCUGACAAAAAAAAAGUUCCAGGUGUGGAAUAUGGGGGUUAUCUGUACAUCCUGGGGCAUUAAAAAAAAAAUCAAUGGUGGGGAACUAUAAAGAAGUAACAAAAGAAGUGACAUCUUCAGCAAAU\
AAACUAGGAAAUUUUUUUUUCUUCCAGUUUAGAAUCAGCCUUGAAACAUUGAUGGAAUAACUCUGUGGCAUUAUUGCAUUAUAUACCAUUUAUCUGUAUUAACUUUGGAAUGUACUCUGUUCAAUGUUUAAUGCUGUGG\
UUGAUAUUUCGAAAGCUGCUUUAAAAAAAUACAUGCAUCUCAGCGUUUUUUUGUUUUUAAUUGUAUUUAGUUAUGGCCUAUACACUAUUUGUGAGCAAAGGUGAUCGUUUUCUGUUUGAGAUUUUUAUCUCUUGAUUCU\
UCAAAAGCAUUCUGAGAAGGUGAGAUAAGCCCUGAGUCUCAGCUACCUAAGAAAAACCUGGAUGUCACUGGCCACUGAGGAGCUUUGUUUCAACCAAGUCAUGUGCAUUUCCACGUCAACAGAAUUGUUUAUUGUGACA\
GUUAUAUCUGUUGUCCCUUUGACCUUGUUUCUUGAAGGUUUCCUCGUCCCUGGGCAAUUCCGCAUUUAAUUCAUGGUAUUCAGGAUUACAUGCAUGUUUGGUUAAACCCAUGAGAUUCAUUCAGUUAAAAAUCCAGAUG\
GCAAAUGACCAGCAGAUUCAAAUCUAUGGUGGUUUGACCUUUAGAGAGUUGCUUUACGUGGCCUGUUUCAACACAGACCCACCCAGAGCCCUCCUGCCCUCCUUCCGCGGGGGCUUUCUCAUGGCUGUCCUUCAGGGUC\
UUCCUGAAAUGCAGUGGUGCUUACGCUCCACCAAGAAAGCAGGAAACCUGUGGUAUGAAGCCAGACCUCCCCGGCGGGCCUCAGGGAACAGAAUGAUCAGACCUUUGAAUGAUUCUAAUUUUUAAGCAAAAUAUUAUUU\
UAUGAAAGGUUUACAUUGUCAAAGUGAUGAAUAUGGAAUAUCCAAUCCUGUGCUGCUAUCCUGCCAAAAUCAUUUUAAUGGAGUCAGUUUGCAGUAUGCUCCACGUGGUAAGAUCCUCCAAGCUGCUUUAGAAGUAACA\
AUGAAGAACGUGGACGUUUUUAAUAUAAAGCCUGUUUUGUCUUUUGUUGUUGUUCAAACGGGAUUCACAGAGUAUUUGAAAAAUGUAUAUAUAUUAAGAGGUCACGGGGGCUAAUUGCUGGCUGGCUGCCUUUUGCUGU\
GGGGUUUUGUUACCUGGUUUUAAUAACAGUAAAUGUGCCCAGCCUCUUGGCCCCAGAACUGUACAGUAUUGUGGCUGCACUUGCUCUAAGAGUAGUUGAUGUUGCAUUUUCCUUAUUGUUAAAAACAUGUUAGAAGCAA\
UGAAUGUAUAUAAAAGCCUCAACUAGUCAUUUUUUUCUCCUCUUCUUUUUUUUCAUUAUAUCUAAUUAUUUUGCAGUUGGGCAACAGAGAACCAUCCCUAUUUUGUAUUGAAGAGGGAUUCACAUCUGCAUCUUAACUG\
CUCUUUAUGAAUGAAAAAACAGUCCUCUGUAUGUACUCCUCUUUACACUGGCCAGGGUCAGAGUUAAAUAGAGUAUAUGCACUUUCCAAAUUGGGGACAAGGGCUCUAAAAAAAGCCCCAAAAGGAGAAGAACAUCUGA\
GAACCUCCUCGGCCCUCCCAGUCCCUCGCUGCACAAAUACUCCGCAAGAGAGGCCAGAAUGACAGCUGACAGGGUCUAUGGCCAUCGGGUCGUCUCCGAAGAUUUGGCAGGGGCAGAAAACUCUGGCAGGCUUAAGAUU\
UGGAAUAAAGUCACAGAAUUAAGGAAGCACCUCAAUUUAGUUCAAACAAGACGCCAACAUUCUCUCCACAGCUCACUUACCUCUCUGUGUUCAGAUGUGGCCUUCCAUUUAUAUGUGAUCUUUGUUUUAUUAGUAAAUG\
CUUAUCAUCUAAAGAUGUAGCUCUGGCCCAGUGGGAAAAAUUAGGAAGUGAUUAUAAAUCGAGAGGAGUUAUAAUAAUCAAGAUUAAAUGUAAAUAAUCAGGGCAAUCCCAACACAUGUCUAGCUUUCACCUCCAGGAU\
CUAUUGAGUGAACAGAAUUGCAAAUAGUCUCUAUUUGUAAUUGAACUUAUCCUAAAACAAAUAGUUUAUAAAUGUGAACUUAAACUCUAAUUAAUUCCAACUGUACUUUUAAGGCAGUGGCUGUUUUUAGACUUUCUUA\
UCACUUAUAGUUAGUAAUGUACACCUACUCUAUCAGAGAAAAACAGGAAAGGCUCGAAAUACAAGCCAUUCUAAGGAAAUUAGGGAGUCAGUUGAAAUUCUAUUCUGAUCUUAUUCUGUGGUGUCUUUUGCAGCCCAGA\
CAAAUGUGGUUACACACUUUUUAAGAAAUACAAUUCUACAUUGUCAAGCUUAUGAAGGUUCCAAUCAGAUCUUUAUUGUUAUUCAAUUUGGAUCUUUCAGGGAUUUUUUUUUUAAAUUAUUAUGGGACAAAGGACAUUU\
GUUGGAGGGGUGGGAGGGAGGAAGAAUUUUUAAAUGUAAAACAUUCCCAAGUUUGGAUCAGGGAGUUGGAAGUUUUCAGAAUAACCAGAACUAAGGGUAUGAAGGACCUGUAUUGGGGUCGAUGUGAUGCCUCUGCGAA\
GAACCUUGUGUGACAAAUGAGAAACAUUUUGAAGUUUGUGGUACGACCUUUAGAUUCCAGAGACAUCAGCAUGGCUCAAAGUGCAGCUCCGUUUGGCAGUGCAAUGGUAUAAAUUUCAAGCUGGAUAUGUCUAAUGGGU\
AUUUAAACAAUAAAUGUGCAGUUUUAACUAACAGGAUAUUUAAUGACAACCUUCUGGUUGGUAGGGACAUCUGUUUCUAAAUGUUUAUUAUGUACAAUACAGAAAAAAAUUUUAUAAAAUUAAGCAAUGUGAAACUGAA\
UUGGAGAGUGAUAAUACAAGUCCUUUAGUCUUACCCAGUGAAUCAUUCUGUUCCAUGUCUUUGGACAACCAUGACCUUGGACAAUCAUGAAAUAUGCAUCUCACUGGAUGCAAAGAAAAUCAGAUGGAGCAUGAAUGGU\
ACUGUACCGGUUCAUCUGGACUGCCCCAGAAAAAUAACUUCAAGCAAACAUCCUAUCAACAACAAGGUUGUUCUGCAUACCAAGCUGAGCACAGAAGAUGGGAACACUGGUGGAGGAUGGAAAGGCUCGCUCAAUCAAG\
AAAAUUCUGAGACUAUUAAUAAAUAAGACUGUAGUGUAGAUACUGAGUAAAUCCAUGCACCUAAACCUUUUGGAAAAUCUGCCGUGGGCCCUCCAGAUAGCUCAUUUCAUUAAGUUUUUCCCUCCAAGGUAGAAUUUGC\
AAGAGUGACAGUGGAUUGCAUUUCUUUUGGGGAAGCUUUCUUUUGGUGGUUUUGUUUAUUAUACCUUCUUAAGUUUUCAACCAAGGUUUGCUUUUGUUUUGAGUUACUGGGGUUAUUUUUGUUUUAAAUAAAAAUAAGU\
GUACAAUAAGUGUUUUUGUAUUGAAAGCUUUUGUUAUCAAGAUUUUCAUACUUUUACCUUCCAUGGCUCUUUUUAAGAUUGAUACUUUUAAGAGGUGGCUGAUAUUCUGCAACACUGUACACAUAAAAAAUACGGUAAG\
GAUACUUUACAUGGUUAAGGUAAAGUAAGUCUCCAGUUGGCCACCAUUAGCUAUAAUGGCACUUUGUUUGUGUUGUUGGAAAAAGUCACAUUGCCAUUAAACUUUCCUUGUCUGUCUAGUUAAUAUUGUGAAGAAAAAU\
AAAGUACAGUGUGAGAUACUG");
std::cout << BCL2_seq.length() << std::endl;
sequence BCL2(BCL2_seq,"chr18",Minus,"60790579","60795857");
}
