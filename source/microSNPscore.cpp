#include <iostream>
#include "sequence.h"
#include "nucleotide.h"

using namespace microSNPscore;

int main(){
   sequence testRNA;
   nucleotide testNuc(Uracil,5,1);
   std::cout << "Test: " << testNuc.get_chromosome_position() << "/" << testNuc.get_sequence_position() << "-" << testNuc.get_base() << std::endl;
}