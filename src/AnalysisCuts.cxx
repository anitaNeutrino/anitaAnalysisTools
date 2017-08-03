#include "AnalysisCuts.h"
#include "AnitaEventSummary.h"


int Acclaim::AnalysisCuts::isAboveHorizontal(const AnitaEventSummary* sum){
  if(!sum){
    return 2;
  }
  else{
    return sum->higherPeak().theta > 0 ? 1 : 0;
  }
}


int Acclaim::AnalysisCuts::isTaggedAsWaisPulser(const AnitaEventSummary* sum){
  if(!sum){
    return 2;
  }
  else{
    return sum->flags.pulser == AnitaEventSummary::EventFlags::WAIS;
  }
}


int Acclaim::AnalysisCuts::higherPol(const AnitaEventSummary* sum){
  if(!sum){
    return AnitaPol::kNotAPol;
  }
  else{
    return sum->higherPeakPol();
  }
}
