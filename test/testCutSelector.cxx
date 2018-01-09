#include "SummaryDraw.h"
#include "CutTreeSelector.h"
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

  SummarySet ss(glob);
  ss.SetUseProof(true);

  CutTreeSelector cts(oc.getOutputFileName(), "cutTree");
  std::vector<const char*> fs;
  fs.push_back("run");
  fs.push_back("eventNumber");  
  fs.push_back("sum.trainingPeak().value");
  cts.setFormulaStrings(fs);

  cts.addCut(&SumTree::isAboveHorizontal);
  
  ss.Process(&cts);

  gSystem->Exit(0); /// Must be called after initializing PROOF!

  return 0;
}
