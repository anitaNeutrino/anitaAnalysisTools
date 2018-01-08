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

  // Cuts::setMode(Cuts::kAcclaimAnalysis);

  SummarySet ss(glob);
  ss.SetUseProof(true);

  SumTreeReductionSelector reduce(oc.getOutputFileName(), "sumTree");

  const int n = 15; //16;
  const TCut* selection[n] = {&Cuts::highestPeak,
			      &Cuts::isRfTrigger,
			      &Cuts::isNotTaggedAsPulser,
			      &Cuts::npbc0A,
			      &Cuts::npbc0B,
			      &Cuts::npbc1,
			      &Cuts::npbc2,
			      &Cuts::npbc3,
			      &Cuts::smallDeltaRough,
			      &Cuts::goodGPS,
			      &Cuts::isBelowHorizontal,
			      &Cuts::reasonableHilbertPeakTimeShiftAfterDedispersion,
			      &Cuts::higherHilbertPeakAfterDedispersion,
			      &Cuts::higherImpulsivityMeasureAfterDedispersion,
			      &Cuts::lowerFracPowerWindowGradientAfterDedispersion
  };
			      // &Cuts::fisherDiscriminantAboveThreshold};
  for(int i=0; i < n; i++){
    reduce.addCut(selection[i]);
  }

  // reduce.fDoAnalysisCutTree = false;

  // reduce.addCut(&Cuts::isNotTaggedAsPulser);
  // reduce.addCut(&Cuts::isGood);
  // reduce.addCut(&Cuts::smallDeltaRough);
  // reduce.addCut(&Cuts::goodGPS);
  // // reduce.addCut(&Cuts::realSNR);
  // reduce.addCut(&Cuts::isRfTrigger);
  // reduce.addCut(&Cuts::higherHilbertPeakAfterDedispersion);
  // reduce.addCut(&Cuts::higherImpulsivityMeasureAfterDedispersion);
  // reduce.addCut(&Cuts::lowerFracPowerWindowGradientAfterDedispersion);
  // reduce.addCut(&Cuts::fisherScoreAboveThreshold); ///@todo fix this

  // reduce.addCut(&Cuts::dedispersedFracPowerWindowGradientBelowThreshold);

  ss.Process(&reduce);
  
  gSystem->Exit(0); /// Must be called after initializing PROOF!

  return 0;
}
