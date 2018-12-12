#include "AcclaimPhaseCenterFitter.h"
#include <iostream>
#include "TString.h"

int main(int argc, char* argv[]){

  Acclaim::PhaseCenterFitter fitter("AcclaimCorrelationSummary_*.root");

  // fitter.setFitParameterSpace(Acclaim::PhaseCenterFitter::ParameterSpace::PitchRollHeading);
  // fitter.fit(AnitaPol::kNotAPol,  "pitchRollHeadingFit.root");
  // fitter.printResults();

  fitter.setFitParameterSpace(Acclaim::PhaseCenterFitter::ParameterSpace::RingPhiRZ);
  // fitter.setPrintOnEval(true);
  for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
    TString fName = TString::Format("RingPhiRZ_%d.root", pol);
    fitter.fit(pol, fName.Data());
    fitter.printResults();
  }


  // fitter.setFitParameterSpace(Acclaim::PhaseCenterFitter::ParameterSpace::RingR);
  // // fitter.setPrintOnEval(true);
  // for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
  //   TString fName = TString::Format("RingR_%d.root", pol);
  //   fitter.fit(pol, fName.Data());
  //   fitter.printResults();
  // }

  // fitter.setFitParameterSpace(Acclaim::PhaseCenterFitter::ParameterSpace::RingZ);
  // // fitter.setPrintOnEval(true);
  // for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
  //   TString fName = TString::Format("RingR_%d.root", pol);
  //   fitter.fit(pol, fName.Data());
  //   fitter.printResults();
  // }
  

}
