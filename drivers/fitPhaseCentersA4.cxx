#include "AcclaimPhaseCenterFitter.h"
#include <iostream>

int main(int argc, char* argv[]){

  Acclaim::PhaseCenterFitter fitter("AcclaimCorrelationSummary_*.root");
  fitter.SetFitParameterSpace(Acclaim::PhaseCenterFitter::ParameterSpace::PitchRoll);
  
  auto& results = fitter.fit();

  std::cout << "Results = ";
  for(auto param : results){
    std::cout << param << ", ";
  }
  std::cout << std::endl;
  return 0;
}
