#include "InterferometricMap.h"
#include "InterferometricMapMaker.h" // for the geometric definitions
#include "TAxis.h"
#include "TMath.h"
#include <iostream>



std::vector<Double_t> coarseBinEdgesPhi; // has size NUM_BINS_PHI+1
std::vector<Double_t> fineBinEdgesPhi; // has size NUM_BINS_PHI_ZOOM_TOTAL+1

std::vector<Double_t> coarseBinEdgesTheta; // has size NUM_BINS_THETA+1
std::vector<Double_t> fineBinEdgesTheta; // has size NUM_BINS_THETA_ZOOM_TOTAL+1
Double_t bin0PhiDeg = -9999;

// static member function
const std::vector<Double_t>& InterferometricMap::getCoarseBinEdgesTheta(){

  if(coarseBinEdgesTheta.size()==0) // then not initialized so do it here...
  {

    // funk up the theta bin spacing...  
    UInt_t nBinsTheta = NUM_BINS_THETA;
    Double_t minTheta = MIN_THETA;
    Double_t maxTheta = MAX_THETA;
  
    // calculate the bin spaces in sin(theta)
    Double_t sinThetaMin = sin(minTheta*TMath::DegToRad());
    Double_t sinThetaMax = sin(maxTheta*TMath::DegToRad());
    Double_t sinThetaRange = sinThetaMax - sinThetaMin;
    Double_t dSinTheta = sinThetaRange/nBinsTheta;

    // std::vector<Double_t> binEdges(nBinsTheta+1);
    coarseBinEdgesTheta.reserve(nBinsTheta+1);
    for(unsigned bt = 0; bt <= nBinsTheta; bt++){
      Double_t thisSinTheta = sinThetaMin + bt*dSinTheta;
      Double_t thisTheta = TMath::RadToDeg()*TMath::ASin(thisSinTheta);
      // coarseBinEdgesTheta.at(bt) = thisTheta;
      coarseBinEdgesTheta.push_back(thisTheta);
    }
  }
  return coarseBinEdgesTheta;
}


const std::vector<Double_t>& InterferometricMap::getFineBinEdgesTheta(){

  if(fineBinEdgesTheta.size()==0) // then not initialized so do it here...
  {

    // funk up the theta bin spacing...  
    UInt_t nBinsTheta = NUM_BINS_THETA_ZOOM_TOTAL;
    Double_t minTheta = MIN_THETA - THETA_RANGE_ZOOM/2;
    Double_t maxTheta = MAX_THETA + THETA_RANGE_ZOOM/2;    
    // Double_t maxTheta = MAX_THETA;
  
    // calculate the bin spaces in sin(theta)
    Double_t sinThetaMin = sin(minTheta*TMath::DegToRad());
    Double_t sinThetaMax = sin(maxTheta*TMath::DegToRad());
    Double_t sinThetaRange = sinThetaMax - sinThetaMin;
    Double_t dSinTheta = sinThetaRange/nBinsTheta;

    // std::vector<Double_t> binEdges(nBinsTheta+1);
    fineBinEdgesTheta.reserve(nBinsTheta+1);
    for(unsigned bt = 0; bt <= nBinsTheta; bt++){
      Double_t thisSinTheta = sinThetaMin + bt*dSinTheta;
      Double_t thisTheta = TMath::RadToDeg()*TMath::ASin(thisSinTheta);
      // coarseBinEdgesTheta.at(bt) = thisTheta;
      fineBinEdgesTheta.push_back(thisTheta);
    }
  }
  return fineBinEdgesTheta;
}


const std::vector<Double_t>& InterferometricMap::getCoarseBinEdgesPhi(){

  if(coarseBinEdgesPhi.size()==0) // then not initialized so do it here...
  {

    // funk up the theta bin spacing...  
    UInt_t nBinsPhi = NUM_BINS_PHI*NUM_PHI;
    Double_t minPhi = getBin0PhiDeg();
    Double_t dPhi = double(DEGREES_IN_CIRCLE)/nBinsPhi;
    // std::vector<Double_t> binEdges(nBinsTheta+1);
    coarseBinEdgesTheta.reserve(nBinsPhi+1);
    for(unsigned bp = 0; bp <= nBinsPhi; bp++){
      Double_t thisPhi = minPhi + dPhi*bp;
      coarseBinEdgesPhi.push_back(thisPhi);
    }
  }
  return coarseBinEdgesPhi;
}



Double_t InterferometricMap::getBin0PhiDeg(){

  if(bin0PhiDeg == -9999){

    AnitaGeomTool* geom = AnitaGeomTool::Instance();
    Double_t aftForeOffset = geom->aftForeOffsetAngleVertical*TMath::RadToDeg();
    
    Double_t phi0 = -aftForeOffset;
    if(phi0 < -DEGREES_IN_CIRCLE/2){
      phi0+=DEGREES_IN_CIRCLE;
    }
    else if(phi0 >= DEGREES_IN_CIRCLE/2){
      phi0-=DEGREES_IN_CIRCLE;
    }
    bin0PhiDeg = phi0 - PHI_RANGE/2;
  }
  return bin0PhiDeg;
}








// class members functions

InterferometricMap::InterferometricMap() : TH2D() {
  initializeInternals();
}




InterferometricMap::InterferometricMap(TString name, TString title, Int_t nBinsPhi, Double_t phiMin, Double_t phiMax, Int_t nBinsTheta, Double_t minTheta, Double_t maxTheta)
  : TH2D(name, title, nBinsPhi, phiMin, phiMax, nBinsTheta, minTheta, maxTheta)
{

  initializeInternals();
}






InterferometricMap::InterferometricMap(TString name, TString title, Double_t phiMin)
  : TH2D(name, title, NUM_PHI*NUM_BINS_PHI, phiMin, phiMin+DEGREES_IN_CIRCLE, getCoarseBinEdgesTheta().size()-1, &getCoarseBinEdgesTheta()[0])
{
  initializeInternals();
}










void InterferometricMap::initializeInternals(){

  // generic
  bool thetaAngleInSinTheta = true;

  // // funk up the theta bin spacing...  
  // UInt_t nBinsTheta = GetNbinsY();
  // Double_t minTheta = fYaxis.GetBinLowEdge(1);
  // Double_t maxTheta = fYaxis.GetBinLowEdge(nBinsTheta+1);  
  
  // // calculate the bin spaces in sin(theta)
  // Double_t sinThetaMin = sin(minTheta*TMath::DegToRad());
  // Double_t sinThetaMax = sin(maxTheta*TMath::DegToRad());
  // Double_t sinThetaRange = sinThetaMax - sinThetaMin;
  // Double_t dSinTheta = sinThetaRange/nBinsTheta;

  // std::vector<Double_t> binEdges(nBinsTheta+1);

  // for(unsigned bt = 0; bt <= nBinsTheta; bt++){
  //   Double_t thisSinTheta = sinThetaMin + bt*dSinTheta;
  //   Double_t thisTheta = TMath::RadToDeg()*TMath::ASin(thisSinTheta);
  //   binEdges.at(bt) = thisTheta;
  // }
  
  // // fXaxis is phi
  // fYaxis.Set(nBinsTheta, &binEdges[0]);
}
