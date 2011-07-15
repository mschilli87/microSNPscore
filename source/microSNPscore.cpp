#include <iostream>
#include <vector>
#include "nucleotide.h"
#include "alignment.h"

using namespace microSNPscore;

int main(){
  nucleotide testNuc(Uracil,5,1);
  std::vector<nucleotide> test_match;
  test_match.push_back(nucleotide(Adenine,15,3));
  test_match.push_back(nucleotide(Cytosine,7,3));
  test_match.push_back(nucleotide(Guanine,17,1));
  test_match.push_back(nucleotide(Uracil,12,23232));
  test_match.push_back(nucleotide(Gap,16,2));
  test_match.push_back(nucleotide(Mask,17,1));
  for(std::vector<nucleotide>::const_iterator match_it=test_match.begin();match_it!=test_match.end();++match_it)
  {
    nucleoBase test_base(testNuc.get_base());
    nucleoBase match_base(match_it->get_base());
    matchIdentifier the_match(alignmentColumn(testNuc,*match_it,Seed,Open).get_match().get_identifier());
    matchScore the_score(alignmentColumn(testNuc,*match_it,Seed,Open).get_match().get_score());
    std::cout << "Aligning " << test_base << " and " << match_base  << " in seed with no open indel: " << the_match << ":" << the_score << std::endl;
    the_match=alignmentColumn(testNuc,*match_it,Seed,Extend).get_match().get_identifier();
    the_score=alignmentColumn(testNuc,*match_it,Seed,Extend).get_match().get_score();
    std::cout << "Aligning " << test_base << " and " << match_base  << " in seed with open indel: " << the_match << ":" << the_score << std::endl;
    the_match=alignmentColumn(testNuc,*match_it,ThreePrime,Open).get_match().get_identifier();
    the_score=alignmentColumn(testNuc,*match_it,ThreePrime,Open).get_match().get_score();
    std::cout << "Aligning " << test_base << " and " << match_base  << " in miRNA 3' with no open indel: " << the_match << ":" << the_score << std::endl;
    the_match=alignmentColumn(testNuc,*match_it,ThreePrime,Extend).get_match().get_identifier();
    the_score=alignmentColumn(testNuc,*match_it,ThreePrime,Extend).get_match().get_score();
    std::cout << "Aligning " << test_base << " and " << match_base  << " in miRNA 3' with open indel: " << the_match << ":" << the_score << std::endl;
  }
}
