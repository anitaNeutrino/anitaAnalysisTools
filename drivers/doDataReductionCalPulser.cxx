#include "SummaryDraw.h"
#include "SumTreeReductionSelector.h"
#include "SummarySet.h"
#include "OutputConvention.h"
#include <iostream>
#include "TSystem.h"

using namespace Acclaim;

int main(int argc, char* argv[]){

  if(argc != 2){
    std::cerr << argv[0] << " 'glob' " << std::endl;
    return 1;
  }
  OutputConvention oc(1, argv);  
  const char* glob = argv[1];

  // Cuts::setMode(Cuts::kAcclaimAnalysis); // not that this should matter

  SummarySet ss(glob);
  ss.SetUseProof(true);

  SumTreeReductionSelector reduce(oc.getOutputFileName(), "sumTree");
  
  reduce.addEventSelectionCut(&Cuts::isTaggedAsWaisPulser);
  reduce.addEventSelectionCut(&Cuts::closeToWais);

  ss.Process(&reduce);
  
  gSystem->Exit(0); /// Must be called after initializing PROOF!

  return 0;
}
