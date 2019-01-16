#include "AcclaimPhaseCenterMinimizer.h"
#include <iostream>
#include "TString.h"

int main(int argc, char* argv[]){

  Acclaim::PhaseCenter::Minimizer fitter("data/AcclaimCorrelationSummary_*.root");

  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::None);
  // fitter.setAllowedPairs(Acclaim::PhaseCenter::Minimizer::AllowedPairs::SamePhiSectorOrHorizontalNeigbour);
  // fitter.fit(AnitaPol::kNotAPol, "None.root");
  // fitter.printResults();

  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::None);
  // fitter.setAllowedPairs(Acclaim::PhaseCenter::Minimizer::AllowedPairs::SamePhiSector);
  // fitter.fit(AnitaPol::kNotAPol, "None_SamePhiSector.root");
  // fitter.printResults();

  fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::PitchRoll);
  fitter.setAllowedPairs(Acclaim::PhaseCenter::Minimizer::AllowedPairs::SamePhiSector);
  // fitter.setPrintOnEval(true);
  // fitter.fit(AnitaPol::kNotAPol,  "pitchRoll_SamePhiSectorOrHorizontalNeigbour.root");
  fitter.fit(AnitaPol::kNotAPol,  "pitchRoll_SamePhiSector.root");  
  fitter.printResults();
  
  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::PitchRollHeading);
  // fitter.fit(AnitaPol::kNotAPol,  "pitchRollHeading.root");
  // fitter.printResults();

  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::PitchRollHeading);
  // fitter.fit(AnitaPol::kNotAPol,  "pitchRollHeadingFit.root");
  // fitter.printResults();

  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::RingPhiRZ);
  // // fitter.setPrintOnEval(true);
  // for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
  //   TString fName = TString::Format("RingPhiRZ_%d.root", pol);
  //   fitter.fit(pol, fName.Data());
  //   fitter.printResults();
  // }

  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::RingR);
  // // fitter.setPrintOnEval(true);
  // for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
  //   TString fName = TString::Format("RingR_%d.root", pol);
  //   fitter.fit(pol, fName.Data());
  //   fitter.printResults();
  // }

  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::RingZ);
  // // fitter.setPrintOnEval(true);
  // for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
  //   TString fName = TString::Format("RingR_%d.root", pol);
  //   fitter.fit(pol, fName.Data());
  //   fitter.printResults();
  // }

  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::RingEllipse);
  // fitter.fit(AnitaPol::kNotAPol, "RingEllipse.root");
  // fitter.printResults();

  // fitter.setPrintOnEval(true);
  // for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
  //   TString fName = TString::Format("RingEllipse_%d.root", pol);
  //   fitter.fit(pol, fName.Data());
  //   fitter.printResults();
  // }


  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::ExtraDeltaT);
  // fitter.fit(AnitaPol::kNotAPol, "ExtraDeltaT.root");
  // fitter.printResults();

}
