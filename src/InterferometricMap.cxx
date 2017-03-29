#include "InterferometricMap.h"
#include "InterferometricMapMaker.h" // for the geometric definitions
#include "InterferometryCache.h"

#include "TAxis.h"
#include "TMath.h"
#include <iostream>

ClassImp(InterferometricMap)

std::vector<Double_t> coarseBinEdgesPhi; // has size NUM_BINS_PHI+1
std::vector<Double_t> fineBinEdgesPhi; // has size NUM_BINS_PHI_ZOOM_TOTAL+1

std::vector<Double_t> coarseBinEdgesTheta; // has size NUM_BINS_THETA+1
std::vector<Double_t> fineBinEdgesTheta; // has size NUM_BINS_THETA_ZOOM_TOTAL+1
Double_t bin0PhiDeg = -9999;



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
    // Double_t sinThetaMin = sin(minTheta*TMath::DegToRad());
    // Double_t sinThetaMax = sin(maxTheta*TMath::DegToRad());
    // Double_t sinThetaRange = sinThetaMax - sinThetaMin;
    // Double_t dSinTheta = sinThetaRange/nBinsTheta;


    Double_t dTheta = (maxTheta - minTheta) / nBinsTheta;
    // std::vector<Double_t> binEdges(nBinsTheta+1);
    fineBinEdgesTheta.reserve(nBinsTheta+1);
    for(unsigned bt = 0; bt <= nBinsTheta; bt++){
      // Double_t thisSinTheta = sinThetaMin + bt*dSinTheta;
      // Double_t thisTheta = TMath::RadToDeg()*TMath::ASin(thisSinTheta);
      // Double_t thisSinTheta = sinThetaMin + bt*dSinTheta;
      
      double thisTheta = minTheta + bt*dTheta;
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


const std::vector<Double_t>& InterferometricMap::getFineBinEdgesPhi(){

  if(fineBinEdgesPhi.size()==0) // then not initialized so do it here...
  {

    // funk up the theta bin spacing...  
    UInt_t nBinsPhi = NUM_BINS_PHI_ZOOM_TOTAL;
    Double_t minPhi = getBin0PhiDeg() - PHI_RANGE_ZOOM/2;
    Double_t dPhi = double(DEGREES_IN_CIRCLE)/nBinsPhi;
    // std::vector<Double_t> binEdges(nBinsTheta+1);
    fineBinEdgesTheta.reserve(nBinsPhi+1);
    for(unsigned bp = 0; bp <= nBinsPhi; bp++){
      Double_t thisPhi = minPhi + dPhi*bp;
      fineBinEdgesPhi.push_back(thisPhi);
    }
  }
  return fineBinEdgesPhi;
}























// class members functions







InterferometricMap::InterferometricMap(TString name, TString title, Double_t zoomCentrePhi, Double_t phiRange, Double_t zoomCentreTheta, Double_t thetaRange)
  : TH2D(name, title, 1, 0, 1, 1, 0, 1) // minimal allocation as we will have to reallocate
{
  isZoomMap = true;

  Double_t minPhiDesired = zoomCentrePhi - phiRange/2;
  Double_t maxPhiDesired = zoomCentrePhi + phiRange/2;
  const std::vector<double> fineBinsPhi = InterferometricMap::getFineBinEdgesPhi();
  int minPhiBin, maxPhiBin;
  getIndicesOfEdgeBins(fineBinsPhi, minPhiDesired, maxPhiDesired, minPhiBin, maxPhiBin);


  Double_t minThetaDesired = zoomCentreTheta - thetaRange/2;
  Double_t maxThetaDesired = zoomCentreTheta + thetaRange/2;
  const std::vector<double> fineBinsTheta = InterferometricMap::getFineBinEdgesTheta();
  int minThetaBin, maxThetaBin;
  getIndicesOfEdgeBins(fineBinsTheta, minThetaDesired, maxThetaDesired, minThetaBin, maxThetaBin);
  
  SetBins(maxPhiBin - minPhiBin, &fineBinsPhi[minPhiBin], maxThetaBin - minThetaBin, &fineBinsTheta[minThetaBin]); // changes also errors array (if any)  
  
  initializeInternals();
}




InterferometricMap::InterferometricMap(TString name, TString title)
  : TH2D(name, title, NUM_PHI*NUM_BINS_PHI, InterferometricMap::getBin0PhiDeg(), InterferometricMap::getBin0PhiDeg()+DEGREES_IN_CIRCLE,
	 getCoarseBinEdgesTheta().size()-1, &getCoarseBinEdgesTheta()[0])
{
  isZoomMap = false;
  initializeInternals();
}



void InterferometricMap::Fill(AnitaPol::AnitaPol_t pol, CrossCorrelator* cc, InterferometryCache* dtCache){


  if(!isZoomMap){
    std::vector<Int_t>* combosToUse = NULL;
    Int_t binsPerPhiSector = GetNbinsPhi()/NUM_PHI;
    for(Int_t phiSector = 0; phiSector < NUM_PHI; phiSector++){
      combosToUse = &cc->combosToUseGlobal[phiSector];

      Int_t startPhiBin = phiSector*binsPerPhiSector;
      Int_t endPhiBin = startPhiBin + binsPerPhiSector;

      for(UInt_t comboInd=0; comboInd<combosToUse->size(); comboInd++){
	Int_t combo = combosToUse->at(comboInd);
	if(cc->kOnlyThisCombo >= 0 && combo!=cc->kOnlyThisCombo){
	  continue;
	}
	for(Int_t phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){
	  for(Int_t thetaBin = 0; thetaBin < GetNbinsTheta(); thetaBin++){

	    Int_t bin = (thetaBin+1)*(GetNbinsPhi()+2) + phiBin+1;
	    Double_t cInterp = cc->getCrossCorrelation(pol, combo, dtCache->coarseDt(pol, combo, phiBin, thetaBin));	  
	    AddBinContent(bin,cInterp);
	    fEntries++; // otherwise drawing don't work at all
	  }
	}
      }
    }

    for(Int_t phiSector = 0; phiSector < NUM_PHI; phiSector++){
      combosToUse = &cc->combosToUseGlobal[phiSector];

      Double_t normFactor = cc->kOnlyThisCombo < 0 && combosToUse->size() > 0 ? combosToUse->size() : 1;
      // absorb the removed inverse FFT normalization
      normFactor*=(cc->numSamples*cc->numSamples);

      Int_t startPhiBin = phiSector*binsPerPhiSector;
      Int_t endPhiBin = (phiSector+1)*binsPerPhiSector;
      for(Int_t phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){
	for(Int_t thetaBin = 0; thetaBin < GetNbinsTheta(); thetaBin++){

	  Double_t val = GetBinContent(phiBin+1, thetaBin+1);
	  SetBinContent(phiBin+1, thetaBin+1, val/normFactor);
	}
      }
    }
  }
}





void InterferometricMap::findPeakValues(Int_t numPeaks, std::vector<Double_t>& peakValues, std::vector<Double_t>& phiDegs, std::vector<Double_t>& thetaDegs){

  
  // In this function I want to find numPeak peaks and set an exclusion zone around each peak
  // so that the next peak isn't just a close neighbour of a true peak.

  // Set not crazy, but still debug visible values for peak values/location
  // -DBL_MAX was causing me some while loop issues in RootTools::getDeltaAngleDeg(...)
  // which has an unrestricted while loop inside.

  peakValues.clear();
  phiDegs.clear();
  thetaDegs.clear();
  
  peakValues.resize(numPeaks, -999);
  phiDegs.resize(numPeaks, -999);
  thetaDegs.resize(numPeaks, -999);

  
  // As we start, all regions are allowed.
  // Int_t allowedBins[NUM_BINS_PHI*NUM_PHI][NUM_BINS_THETA];
  // Int_t allowedBins[NUM_BINS_PHI*NUM_PHI][NUM_BINS_THETA];
  int nTheta = GetNbinsTheta();
  int nPhi = GetNbinsTheta();  
  std::vector<short> allowedBins(nPhi*nTheta, 1);

  
  for(Int_t peakInd=0; peakInd < numPeaks; peakInd++){
    Int_t gotHere = 0;
    for(Int_t phiBin=0; phiBin<nPhi; phiBin++){
      for(Int_t thetaBin = 0; thetaBin < nTheta; thetaBin++){

	// if(allowedBins[phiBin][thetaBin] > 0){
	if(allowedBins[phiBin*nTheta + thetaBin] > 0){	  
	  // then we are in an allowed region

	  // Double_t mapValue = coarseMap[pol][phiBin][thetaBin];
	  Double_t mapValue = GetBinContent(phiBin+1, thetaBin+1);
	  if(mapValue > peakValues[peakInd]){
	    // higher than the current highest map value

	    gotHere = 1;
	    // Time for something a little tricky to read...
	    // I want to stop bins that are on the edge of allowed regions registering as peaks
	    // if the bins just inside the disallowed regions are higher valued.
	    // To ensure a local maxima so I am going to ensure that all neighbouring
	    // bins are less than the current value.
	    // Checking the 3x3 grid around the current bin (wrapping in phi).

	    // Wrap phi edges
	    Int_t lastPhiBin = phiBin != 0 ? phiBin - 1 : (NUM_BINS_PHI*NUM_PHI) - 1;
	    Int_t nextPhiBin = phiBin != (NUM_BINS_PHI*NUM_PHI) - 1 ? phiBin + 1 : 0;


	    // Is current bin higher than phi neighbours in same theta row?
	    // if(coarseMap[pol][lastPhiBin][thetaBin] < mapValue &&
	    //    coarseMap[pol][nextPhiBin][thetaBin] < mapValue){
	    if(GetBinContent(lastPhiBin+1, thetaBin+1) < mapValue &&
	       GetBinContent(nextPhiBin+1, thetaBin+1) < mapValue){

	      gotHere = 2;
	      // So far so good, now check theta neigbours
	      // This doesn't wrap, so I will allow edge cases to pass as maxima
	      Int_t nextThetaBin = thetaBin != (NUM_BINS_THETA) - 1 ? thetaBin+1 : thetaBin;
	      Int_t lastThetaBin = thetaBin != 0 ? thetaBin-1 : thetaBin;

	      // Check theta bins below
	      gotHere = 3;
	      if(thetaBin == 0 || (GetBinContent(lastPhiBin+1, lastThetaBin+1) < mapValue &&
				   GetBinContent(phiBin+1, lastThetaBin+1) < mapValue &&
				   GetBinContent(nextPhiBin+1, lastThetaBin+1) < mapValue)){

		// Check theta bins above
		gotHere = 4;
		if(thetaBin == NUM_BINS_THETA-1 || (GetBinContent(lastPhiBin+1, nextThetaBin+1) < mapValue &&
						    GetBinContent(phiBin+1, nextThetaBin+1) < mapValue &&
						    GetBinContent(nextPhiBin+1, nextThetaBin+1) < mapValue)){

		  // Okay okay okay... you're a local maxima, you can be my next peak.
		  peakValues[peakInd] = GetBinContent(phiBin+1, thetaBin+1);
		  // phiDegs[peakInd] = phiWaveLookup[phiBin]*TMath::RadToDeg();
		  // thetaDegs[peakInd] = thetaWaves[thetaBin]*TMath::RadToDeg();
		  phiDegs[peakInd] = GetPhiAxis()->GetBinLowEdge(phiBin+1);
		  thetaDegs[peakInd] = GetThetaAxis()->GetBinLowEdge(thetaBin+1);		  
		  // thetaDegs[peakInd] = InterferometricMap::getCoarseBinEdgesTheta()[thetaBin];
		}
	      }
	    }
	  }
	}
      } // thetaBin
    } // phiBin

    // Still inside the peakInd loop

    // Here I set the exclusion zone around the peak bin.
    if(peakValues[peakInd] >= -999){ // Checks that a peak was found, probably unnecessary

      for(Int_t phiBin=0; phiBin<nPhi; phiBin++){
	// Double_t phiDeg = phiWaveLookup[phiBin]*TMath::RadToDeg();
	Double_t phiDeg = GetPhiAxis()->GetBinLowEdge(phiBin+1);
	Double_t absDeltaPhi = TMath::Abs(RootTools::getDeltaAngleDeg(phiDegs[peakInd], phiDeg));
	if(absDeltaPhi < PEAK_PHI_DEG_RANGE){

	  for(Int_t thetaBin = 0; thetaBin < nTheta; thetaBin++){

	    // Double_t thetaDeg = thetaWaves[thetaBin]*TMath::RadToDeg();
	    Double_t thetaDeg = GetThetaAxis()->GetBinLowEdge(thetaBin+1);
	    Double_t absDeltaTheta = TMath::Abs(RootTools::getDeltaAngleDeg(thetaDegs[peakInd], thetaDeg));
	    if(absDeltaTheta < PEAK_THETA_DEG_RANGE){

	      // Now this region in phi/theta is disallowed on the next loop through...
	      // allowedBins[phiBin][thetaBin] = 0;
	      allowedBins[phiBin*nTheta + thetaBin] = 0;	      
	    }
	  }
	}
      }
    }
    if(peakValues[peakInd] < 0){
      std::cerr << "Peak " << peakInd << " = " << peakValues[peakInd] << "\t" << gotHere << std::endl;
      for(int pi=0; pi <= peakInd; pi++){
    	std::cerr << peakValues[pi] << "\t" << phiDegs[pi] << "\t" << thetaDegs[pi] << std::endl;
      }
    }
  }
}







void InterferometricMap::initializeInternals(){
  thetaAxisInSinTheta = true;
}



void InterferometricMap::getIndicesOfEdgeBins(const std::vector<double>& binEdges, Double_t lowVal, Double_t highVal, Int_t& lowIndex, Int_t& highIndex){

  // get last below
  lowIndex = -1;
  for(unsigned i=0; i < binEdges.size(); i++){
    if(binEdges.at(i) < lowVal){
      lowIndex = i;
    }
    else{
      break;
    }
  }

  highIndex = binEdges.size() - 1;
  // get last phi above
  for(int i=binEdges.size()-1; i >= 0; i--){
    if(binEdges.at(i) > highVal){
      highIndex = i;
    }
    else{
      break;
    }
  }
}
