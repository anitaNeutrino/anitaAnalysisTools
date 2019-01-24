#include "AcclaimPhaseCenterMinimizer.h"
#include <iostream>
#include "TString.h"

int main(int argc, char* argv[]){

  Acclaim::PhaseCenter::Minimizer fitter("data/AcclaimCorrelationSummary_*.root");

  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::None);
  // fitter.setAllowedPairs(Acclaim::PhaseCenter::Minimizer::AllowedPairs::SamePhiSector);
  // fitter.fit(AnitaPol::kNotAPol);
  // fitter.printResults();

  fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::PitchRoll);
  fitter.setAllowedPairs(Acclaim::PhaseCenter::Minimizer::AllowedPairs::SamePhiSector);
  fitter.fit(AnitaPol::kNotAPol);
  fitter.printResults();

  // pitch = -0.201669
  // roll = -0.115879

  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::RingEllipse);
  // fitter.setAllowedPairs(Acclaim::PhaseCenter::Minimizer::AllowedPairs::SamePhiSector);
  // fitter.fit(AnitaPol::kNotAPol);
  // fitter.printResults();

  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::RingR);
  // fitter.setAllowedPairs(Acclaim::PhaseCenter::Minimizer::AllowedPairs::SamePhiSector);
  // fitter.fit(AnitaPol::kNotAPol);
  // fitter.printResults();


  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::RingZ);
  // fitter.setAllowedPairs(Acclaim::PhaseCenter::Minimizer::AllowedPairs::SamePhiSector);
  // fitter.fit(AnitaPol::kNotAPol);
  // fitter.printResults();


  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::RingPhiRZ);
  // fitter.setAllowedPairs(Acclaim::PhaseCenter::Minimizer::AllowedPairs::SamePhiSector);
  // fitter.fit(AnitaPol::kNotAPol);
  // fitter.printResults();

  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::AntPhi);
  // fitter.setAllowedPairs(Acclaim::PhaseCenter::Minimizer::AllowedPairs::SamePhiSector);
  // fitter.fit(AnitaPol::kNotAPol);
  // fitter.printResults();

  // fitter.setParameterSpace(Acclaim::PhaseCenter::ParameterSpace::AntPhiRZ);
  // fitter.setAllowedPairs(Acclaim::PhaseCenter::Minimizer::AllowedPairs::SamePhiSector);
  // fitter.fit(AnitaPol::kNotAPol);
  // fitter.printResults();


}
