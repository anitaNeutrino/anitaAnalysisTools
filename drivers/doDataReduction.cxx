#include "AnalysisCuts.h"
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

  // AnalysisCuts::setMode(AnalysisCuts::kAcclaimAnalysis);

  SummarySet ss(glob);
  ss.SetUseProof(true);

  SumTreeReductionSelector reduce(oc.getOutputFileName(), "sumTree");
  
  reduce.addEventSelectionCut(&AnalysisCuts::isNotTaggedAsPulser);
  reduce.addEventSelectionCut(&AnalysisCuts::isGood);
  reduce.addEventSelectionCut(&AnalysisCuts::smallDeltaRough);
  reduce.addEventSelectionCut(&AnalysisCuts::goodGPS);
  reduce.addEventSelectionCut(&AnalysisCuts::realSNR);
  reduce.addEventSelectionCut(&AnalysisCuts::isRfTrigger);
  reduce.addEventSelectionCut(&AnalysisCuts::higherHilbertPeakAfterDedispersion);
  reduce.addEventSelectionCut(&AnalysisCuts::higherImpulsivityMeasureAfterDedispersion);
  reduce.addEventSelectionCut(&AnalysisCuts::lowerFracPowerWindowGradientAfterDedispersion);
  // reduce.addEventSelectionCut(&AnalysisCuts::fisherScoreAboveThreshold); ///@todo fix this

  // reduce.addEventSelectionCut(&AnalysisCuts::dedispersedFracPowerWindowGradientBelowThreshold);

  ss.Process(&reduce);
  
  gSystem->Exit(0); /// Must be called after initializing PROOF!

  return 0;
}
