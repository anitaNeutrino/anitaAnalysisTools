#include "AnalysisCuts.h"
#include "AnitaEventSummary.h"
#include "RampdemReader.h"

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


int Acclaim::AnalysisCuts::hasSourceLocation(const AnitaEventSummary* sum){
  if(!sum){
    return 2;
  }
  else{
    bool didReconstruct = (sum->higherPeak().latitude < -900 || sum->higherPeak().theta_adjustment_needed > 0) ? false : true;
    return didReconstruct;
  }
}


int Acclaim::AnalysisCuts::isOnContinent(const AnitaEventSummary* sum){
  if(!sum){
    return 2;
  }
  else{
    return RampdemReader::isOnContinent(sum->higherPeak().longitude, sum->higherPeak().latitude);
  }
}
