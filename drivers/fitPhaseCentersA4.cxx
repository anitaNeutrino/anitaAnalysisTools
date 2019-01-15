#include "AcclaimPhaseCenterFitter.h"
#include <iostream>
#include "TString.h"

int main(int argc, char* argv[]){

  Acclaim::PhaseCenterFit::Minimizer fitter("data/AcclaimCorrelationSummary_*.root");

  fitter.setParameterSpace(Acclaim::PhaseCenterFit::ParameterSpace::None);
  fitter.fit(AnitaPol::kNotAPol,  "None.root");
  fitter.printResults();

  // fitter.setParameterSpace(Acclaim::PhaseCenterFit::ParameterSpace::PitchRollHeading);
  // fitter.fit(AnitaPol::kNotAPol,  "pitchRollHeading.root");
  // fitter.printResults();

  // fitter.setParameterSpace(Acclaim::PhaseCenterFit::ParameterSpace::PitchRollHeading);
  // fitter.fit(AnitaPol::kNotAPol,  "pitchRollHeadingFit.root");
  // fitter.printResults();

  // fitter.setParameterSpace(Acclaim::PhaseCenterFit::ParameterSpace::RingPhiRZ);
  // // fitter.setPrintOnEval(true);
  // for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
  //   TString fName = TString::Format("RingPhiRZ_%d.root", pol);
  //   fitter.fit(pol, fName.Data());
  //   fitter.printResults();
  // }

  // fitter.setParameterSpace(Acclaim::PhaseCenterFit::ParameterSpace::RingR);
  // // fitter.setPrintOnEval(true);
  // for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
  //   TString fName = TString::Format("RingR_%d.root", pol);
  //   fitter.fit(pol, fName.Data());
  //   fitter.printResults();
  // }

  // fitter.setParameterSpace(Acclaim::PhaseCenterFit::ParameterSpace::RingZ);
  // // fitter.setPrintOnEval(true);
  // for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
  //   TString fName = TString::Format("RingR_%d.root", pol);
  //   fitter.fit(pol, fName.Data());
  //   fitter.printResults();
  // }

  // fitter.setParameterSpace(Acclaim::PhaseCenterFit::ParameterSpace::RingEllipse);
  // fitter.fit(AnitaPol::kNotAPol, "RingEllipse.root");
  // fitter.printResults();

  // fitter.setPrintOnEval(true);
  // for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
  //   TString fName = TString::Format("RingEllipse_%d.root", pol);
  //   fitter.fit(pol, fName.Data());
  //   fitter.printResults();
  // }


  // fitter.setParameterSpace(Acclaim::PhaseCenterFit::ParameterSpace::ExtraDeltaT);
  // fitter.fit(AnitaPol::kNotAPol, "ExtraDeltaT.root");
  // fitter.printResults();

}
