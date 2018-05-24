#include "DrawStrings.h"
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

  // SumTree::setMode(SumTree::kAcclaimAnalysis);

  SummarySet ss(glob);
  ss.SetUseProof(true);

  SumTreeReductionSelector reduce(oc.getOutputFileName(), "sumTree");

  const int n = 15; //16;
  const TCut* selection[n] = {&SumTree::highestPeak,
			      &SumTree::isRfTrigger,
			      &SumTree::isNotTaggedAsPulser,
			      &SumTree::newPayloadBlastCutPart0A,
			      &SumTree::newPayloadBlastCutPart0B,
			      &SumTree::newPayloadBlastCutPart1,
			      &SumTree::newPayloadBlastCutPart2,
			      &SumTree::newPayloadBlastCutPart3,
			      &SumTree::smallDeltaRough,
			      &SumTree::goodGPS,
			      &SumTree::isBelowHorizontal,
			      &SumTree::reasonableHilbertPeakTimeShiftAfterDedispersion,
			      &SumTree::higherHilbertPeakAfterDedispersion,
			      &SumTree::higherImpulsivityMeasureAfterDedispersion,
			      &SumTree::lowerFracPowerWindowGradientAfterDedispersion
  };
			      // &SumTree::fisherDiscriminantAboveThreshold};
  for(int i=0; i < n; i++){
    reduce.addCut(selection[i]);
  }

  // reduce.fDoAnalysisCutTree = false;

  // reduce.addCut(&SumTree::isNotTaggedAsPulser);
  // reduce.addCut(&SumTree::isGood);
  // reduce.addCut(&SumTree::smallDeltaRough);
  // reduce.addCut(&SumTree::goodGPS);
  // // reduce.addCut(&SumTree::realSNR);
  // reduce.addCut(&SumTree::isRfTrigger);
  // reduce.addCut(&SumTree::higherHilbertPeakAfterDedispersion);
  // reduce.addCut(&SumTree::higherImpulsivityMeasureAfterDedispersion);
  // reduce.addCut(&SumTree::lowerFracPowerWindowGradientAfterDedispersion);
  // reduce.addCut(&SumTree::fisherScoreAboveThreshold); ///@todo fix this

  // reduce.addCut(&SumTree::dedispersedFracPowerWindowGradientBelowThreshold);

  ss.Process(&reduce);
  
  gSystem->Exit(0); /// Must be called after initializing PROOF!

  return 0;
}
