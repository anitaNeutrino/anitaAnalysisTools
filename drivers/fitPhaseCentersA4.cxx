#include "AcclaimPhaseCenterMinimizer.h"
#include <iostream>
#include "TString.h"

int main(int argc, char* argv[]){

  Acclaim::PhaseCenter::Minimizer fitter("data/AcclaimCorrelationSummary_*.root");

  fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::None);
  fitter.setAllowedPairs(Acclaim::PhaseCenter::Minimizer::AllowedPairs::SamePhiSector);
  fitter.fit(AnitaPol::kNotAPol);
  fitter.printResults();

  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::None);
  // fitter.setAllowedPairs(Acclaim::PhaseCenter::Minimizer::AllowedPairs::SamePhiSector);
  // fitter.fit(AnitaPol::kNotAPol, "None_SamePhiSector.root");
  // fitter.printResults();

  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::PitchRoll);
  // fitter.setAllowedPairs(Acclaim::PhaseCenter::Minimizer::AllowedPairs::SamePhiSector);
  // fitter.fit(AnitaPol::kNotAPol);
  // fitter.printResults();
  
  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::PitchRollHeading);
  // fitter.fit(AnitaPol::kNotAPol);
  // fitter.printResults();

  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::PitchRollHeading);
  // fitter.fit(AnitaPol::kNotAPol);
  // fitter.printResults();

  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::RingPhiRZ);
  // // fitter.setPrintOnEval(true);
  // for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
  //   fitter.fit(pol);
  //   fitter.printResults();
  // }

  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::RingR);
  // // fitter.setPrintOnEval(true);
  // for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
  //   fitter.fit(pol);
  //   fitter.printResults();
  // }

  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::RingZ);
  // // fitter.setPrintOnEval(true);
  // for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
  //   fitter.fit(pol);
  //   fitter.printResults();
  // }

  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::RingEllipse);
  // fitter.fit(AnitaPol::kNotAPol);
  // fitter.printResults();

  // fitter.setPrintOnEval(true);
  // for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
  //   TString fName = TString::Format("RingEllipse_%d.root", pol);
  //   fitter.fit(pol, fName.Data());
  //   fitter.printResults();
  // }


  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::ExtraDeltaT);
  // fitter.fit(AnitaPol::kNotAPol);
  // fitter.printResults();

}
