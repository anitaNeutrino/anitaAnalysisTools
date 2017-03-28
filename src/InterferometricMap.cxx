#include "InterferometricMap.h"
#include "InterferometricMapMaker.h" // for the geometric definitions
#include "TAxis.h"
#include "TMath.h"
#include <iostream>

InterferometricMap::InterferometricMap() : TH2D() {
  
}


void InterferometricMap::initializeInternals(){

  // generic
  bool thetaAngleInSinTheta = true;

  // funk up the theta bin spacing...  
  UInt_t nBinsTheta = GetNbinsY();
  Double_t minTheta = fYaxis.GetBinLowEdge(1);
  Double_t maxTheta = fYaxis.GetBinLowEdge(nBinsTheta+1);  
  
  // calculate the bin spaces in sin(theta)
  Double_t sinThetaMin = sin(minTheta*TMath::DegToRad());
  Double_t sinThetaMax = sin(maxTheta*TMath::DegToRad());
  Double_t sinThetaRange = sinThetaMax - sinThetaMin;
  Double_t dSinTheta = sinThetaRange/nBinsTheta;

  std::vector<Double_t> binEdges(nBinsTheta+1);

  for(unsigned bt = 0; bt <= nBinsTheta; bt++){
    Double_t thisSinTheta = sinThetaMin + bt*dSinTheta;
    Double_t thisTheta = TMath::RadToDeg()*TMath::ASin(thisSinTheta);
    binEdges.at(bt) = thisTheta;
  }
  
  // fXaxis is phi
  fYaxis.Set(nBinsTheta, &binEdges[0]);
}


InterferometricMap::InterferometricMap(TString name, TString title, Int_t nBinsPhi, Double_t phiMin, Double_t phiMax, Int_t nBinsTheta, Double_t minTheta, Double_t maxTheta)
  : TH2D(name, title, nBinsPhi, phiMin, phiMax, nBinsTheta, minTheta, maxTheta)
{

  initializeInternals();
}







InterferometricMap::InterferometricMap(TString name, TString title, Double_t phiMin)
  : TH2D(name, title, NUM_PHI*NUM_BINS_PHI, phiMin, phiMin+DEGREES_IN_CIRCLE, NUM_BINS_THETA, MIN_THETA, MAX_THETA)
{
  initializeInternals();
}




