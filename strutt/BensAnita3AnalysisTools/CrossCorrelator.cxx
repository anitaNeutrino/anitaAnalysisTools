/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             A Cross Correlator to interact with the ROOTified ANITA-3 data and do some interferometry. 
*************************************************************************************************************** */

#include "CrossCorrelator.h"

ClassImp(CrossCorrelator)


/************************************************************************************************************
Constructor and destructor functions
************************************************************************************************************/

/*!
  \brief Constructor
*/
CrossCorrelator::CrossCorrelator(){
  initializeVariables();
}

/*!
  \brief Destructor
*/
CrossCorrelator::~CrossCorrelator(){
  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    deleteAllWaveforms((AnitaPol::AnitaPol_t)pol);
    deleteCrossCorrelations((AnitaPol::AnitaPol_t)pol);
  }
}


/*!
  \brief Workhorse function to set internal variables.
  /param upSampleFactorTemp is passed here from the constructor.
*/
void CrossCorrelator::initializeVariables(){

  // Initialize with NULL otherwise very bad things will happen with gcc 
  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      grs[pol][ant] = NULL;
      grsInterp[pol][ant] = NULL;
      interpRMS[pol][ant] = 0;
    }
    for(int combo=0; combo<NUM_COMBOS; combo++){
      crossCorrelations[pol][combo] = NULL;
    }
    lastEventNormalized[pol] = 0;
    eventNumber[pol] = 0;    
  }


  nominalSamplingDeltaT = 1./2.6;
  upsampleFactor = 40; // 0.38ns -> ~10ps
  correlationDeltaT = nominalSamplingDeltaT/upsampleFactor;
  numSamplesUpsampled = 2*NUM_SAMPLES*upsampleFactor;

  deltaTMax = 0;
  deltaTMin = 0;

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    rArray.push_back(geom->getAntR(ant));
    zArray.push_back(geom->getAntZ(ant));
    phiArrayDeg.push_back(geom->getAntPhiPositionRelToAftFore(ant)*TMath::RadToDeg());
  }

  do5PhiSectorCombinatorics();
  fillDeltaTLookup();
  printInfo();


  mapModeNames[kGlobal] = "Global";
  mapModeNames[kTriggered] = "Triggered";
  zoomModeNames[kZoomedOut] = "";
  zoomModeNames[kZoomedIn] = "Zoom";
  
}



void CrossCorrelator::printInfo(){
  std::cerr << __PRETTY_FUNCTION__ << std::endl;
  std::cerr << "\tupsample factor = " << upsampleFactor << std::endl;
  std::cerr << "\tdeltaT max = " << deltaTMax << std::endl;
  std::cerr << "\tdeltaT min = " << deltaTMin << std::endl;
  std::cerr << "\tBin size theta (deg) = " << Double_t(THETA_RANGE)/NUM_BINS_THETA << std::endl;
  std::cerr << "\tBin size phi (deg) = " << Double_t(PHI_RANGE)/NUM_BINS_PHI << std::endl;
  std::cerr << "\tdeltaTs array size = " << sizeof(dtIndex_t)*numCombos*NUM_PHI*NUM_BINS_PHI*NUM_BINS_THETA << " bytes" << std::endl;
}





/*!
  \brief Single function to get the angle of the first bin of the interferometric histogram.
  \returns position of antenna 0 in ADU5Pat coordinates, offset by half a phi-sector.
*/

Double_t CrossCorrelator::getBin0PhiDeg(){
  Double_t phi0 = phiArrayDeg.at(0);
  if(phi0 < -180){
    phi0+=360;
  }
  else if(phi0 >= 180){
    phi0-=360;
  }
  return phi0 - PHI_RANGE/2;
}










/************************************************************************************************************
Waveform manipulation functions
************************************************************************************************************/


/*!
  \brief Loops through all waveform graphs in the UsefulAnitaEvent and makes an evenly re-sampled, normalized copy of each one.
  \param usefulEvent points to the UsefulAnitaEvent of interest
*/
void CrossCorrelator::getNormalizedInterpolatedTGraphs(UsefulAnitaEvent* usefulEvent, AnitaPol::AnitaPol_t pol){
  // Potentially needed in a few places, so it gets its own function 

  // Pretty much just for profiling 
  if(usefulEvent->eventNumber!=lastEventNormalized[pol]){

    // Delete any old waveforms (at start rather than end to leave in memory to be examined if need be)
    deleteAllWaveforms(pol);

    // Find the start time of all waveforms 
    Double_t earliestStart[NUM_POL] = {100};
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      grs[pol][ant] = usefulEvent->getGraph(ant, (AnitaPol::AnitaPol_t)pol);

      if(grs[pol][ant]->GetX()[0]<earliestStart[pol]){
	earliestStart[pol] = grs[pol][ant]->GetX()[0];
      }
    }

    // Interpolate with earliest start time 
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      grsInterp[pol][ant] = interpolateWithStartTime(grs[pol][ant], earliestStart[pol]);
      Double_t mean; // don't care about this.
      RootTools::normalize(grsInterp[pol][ant], mean, interpRMS[pol][ant]);
    }
    lastEventNormalized[pol] = usefulEvent->eventNumber;
  }
}


/*!
  \brief Wrapper function for ROOT's interpolator, can zero pad the front to start from a particular time.
  \param grIn points to the TGraph containing the waveform to interpolate
  \param startTime is the start time for interpolation: zero pads if this is earlier than the TGraph start time.
*/
TGraph* CrossCorrelator::interpolateWithStartTime(TGraph* grIn, Double_t startTime){

  std::vector<Double_t> newTimes = std::vector<Double_t>(numSamplesUpsampled, 0);
  std::vector<Double_t> newVolts = std::vector<Double_t>(numSamplesUpsampled, 0);  
  Double_t thisStartTime = grIn->GetX()[0];
  Double_t lastTime = grIn->GetX()[grIn->GetN()-1];


  // Quantizes the start and end times so data poInt_ts lie at Int_teger multiples of nominal sampling 
  startTime = correlationDeltaT*TMath::Nint(startTime/correlationDeltaT + 0.5);
  lastTime = correlationDeltaT*TMath::Nint(lastTime/correlationDeltaT - 0.5);

   //ROOT Int_terpolator object constructor takes std::vector objects
  std::vector<Double_t> tVec(grIn->GetX(), grIn->GetX() + grIn->GetN());
  std::vector<Double_t> vVec(grIn->GetY(), grIn->GetY() + grIn->GetN());
  
  // This is ROOT's Int_terpolator object
  ROOT::Math::Interpolator chanInterp(tVec,vVec,ROOT::Math::Interpolation::kAKIMA);
  
  // Put new data Int_to arrays
  Double_t time = startTime;
  for(Int_t samp = 0; samp < numSamplesUpsampled; samp++){
    newTimes.at(samp) = time;
    if(time >= thisStartTime && time <= lastTime){
      newVolts.at(samp) = chanInterp.Eval(time);
    }
    else{
      newVolts.at(samp) = 0;
    }
    time += correlationDeltaT;
  }

  return new TGraph(numSamplesUpsampled, &newTimes[0], &newVolts[0]);

}





















/************************************************************************************************************
All correlation functions
************************************************************************************************************/


/*!
  \brief Loops through the set of cross correlations and returns a vector of vectors containing the maximum times for the specified polarization.
  \returns A vector of length numCombos containing the set of times of maximum correlation.
*/
std::vector<Double_t> CrossCorrelator::getMaxCorrelationTimes(AnitaPol::AnitaPol_t pol){
  
  std::vector<Double_t> peakTimes(NUM_COMBOS, 0);

  for(Int_t combo=0; combo<NUM_COMBOS; combo++){
    Int_t maxInd = RootTools::getIndexOfMaximum(numSamplesUpsampled, crossCorrelations[pol][combo]);
    Int_t offset = maxInd >= numSamplesUpsampled/2 ? maxInd - numSamplesUpsampled : maxInd;
    peakTimes.at(combo) = offset*correlationDeltaT;      
  }
  return peakTimes;
}

/*!
  \brief Loops through the set of cross correlations and returns a vector of vectors containing the maximum times.
  \returns A vector length 2 of vector length numCombos containing the set of times of maximum correlation.
*/
std::vector<std::vector<Double_t> > CrossCorrelator::getMaxCorrelationTimes(){
  
  std::vector<std::vector<Double_t> > peakTimes;
  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    peakTimes.push_back(getMaxCorrelationTimes((AnitaPol::AnitaPol_t)pol));
  }
  return peakTimes;
}


/*!
  \brief Loops through the set of cross correlations and returns a vector of vectors containing the maximum correlation values for the specified polarization.
  \returns A vector of length numCombos containing the set of maximum correlation values.
*/
std::vector<Double_t> CrossCorrelator::getMaxCorrelationValues(AnitaPol::AnitaPol_t pol){
  
  std::vector<Double_t> peakVals(NUM_COMBOS, 0);
  for(Int_t combo=0; combo<NUM_COMBOS; combo++){
    Int_t maxInd = RootTools::getIndexOfMaximum(numSamplesUpsampled, crossCorrelations[pol][combo]);
    peakVals.at(combo) = crossCorrelations[pol][combo][maxInd];
  }
  return peakVals;
}


/*!
  \brief Loops through the set of cross correlations and returns a vector of vectors containing the maximum correlation values.
  \returns A vector length 2 of vector length numCombos containing the set of maximum correlation values.
*/
std::vector<std::vector<Double_t> > CrossCorrelator::getMaxCorrelationValues(){
  
  std::vector<std::vector<Double_t> > peakVals;
  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    peakVals.push_back(getMaxCorrelationValues((AnitaPol::AnitaPol_t)pol));    
  }
  return peakVals;
}



/*!
  \brief Correlate the event

  First step in interferometry, pass pointer to UsefulAnitaEvent in here, which populates the internal arrays of normalized, interpolated waveforms and set of cross correlations. These internal values are then used as look up tables for generating interferometic images.
*/
void CrossCorrelator::correlateEvent(UsefulAnitaEvent* usefulEvent){

  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){  
    correlateEvent(usefulEvent, (AnitaPol::AnitaPol_t)pol);
  }
}


/*!
  \brief Correlate the event

  First step in interferometry, pass pointer to UsefulAnitaEvent in here, which populates the internal arrays of normalized, interpolated waveforms and set of cross correlations. These internal values are then used as look up tables for generating interferometic images.
*/
void CrossCorrelator::correlateEvent(UsefulAnitaEvent* usefulEvent, AnitaPol::AnitaPol_t pol){

  // Read TGraphs from events into memory (also deletes old TGraphs)
  getNormalizedInterpolatedTGraphs(usefulEvent, pol);

  // Cross correlate waveforms using the normalized TGraphs.
  doAllCrossCorrelations(pol);

  // Safety check to make sure we don't do any hard work twice.
  eventNumber[pol] = usefulEvent->eventNumber;
  
}



/*!
  \brief Loop over both polarizations and all combinations, and generate set of cross correlations.
*/
void CrossCorrelator::doAllCrossCorrelations(AnitaPol::AnitaPol_t pol){

  // Delete old cross correlations first 
  deleteCrossCorrelations(pol);

  // Loop over combinations generating set of cross correlations.
  for(Int_t combo=0; combo<numCombos; combo++){
    Int_t ant1 = comboToAnt1s.at(combo);
    Int_t ant2 = comboToAnt2s.at(combo);
    crossCorrelations[pol][combo] = crossCorrelateFourier(grsInterp[pol][ant2], grsInterp[pol][ant1]);
  }
}


/*!
  \brief Interface to FFT library to generate cross correlations.
*/
Double_t* CrossCorrelator::crossCorrelateFourier(TGraph* gr1, TGraph* gr2){
  // Generate cross correlations, now using FancyFFTs 
  //  return  FFTtools::getCorrelation(gr1->GetN(), gr1->GetY(), gr2->GetY()) ;
  return FancyFFTs::crossCorrelate(gr1->GetN(), gr1->GetY(), gr2->GetY()) ;
}
















/************************************************************************************************************
Calculate deltaT between two antennas (for a plane wave unless function name says otherwise)
************************************************************************************************************/

Int_t CrossCorrelator::getDeltaTExpected(Int_t ant1, Int_t ant2, Double_t phiWave, Double_t thetaWave){
  Double_t tanThetaW = tan(-thetaWave);
  Double_t part1 = zArray.at(ant1)*tanThetaW - rArray.at(ant1)*cos(phiWave-TMath::DegToRad()*phiArrayDeg.at(ant1));
  Double_t part2 = zArray.at(ant2)*tanThetaW - rArray.at(ant2)*cos(phiWave-TMath::DegToRad()*phiArrayDeg.at(ant2));
  Double_t tdiff = 1e9*((cos(-thetaWave) * (part2 - part1))/SPEED_OF_LIGHT); // Returns time in ns

  tdiff /= correlationDeltaT;
  return TMath::Nint(tdiff);
}


Int_t CrossCorrelator::getDeltaTExpectedSpherical(Int_t ant1, Int_t ant2, Double_t phiWave, Double_t thetaWave, Double_t rWave){

  Double_t phi1 = TMath::DegToRad()*phiArrayDeg.at(ant1);
  Double_t x1 = rArray.at(ant1)*TMath::Cos(phi1);
  Double_t y1 = rArray.at(ant1)*TMath::Sin(phi1);
  Double_t z1 = zArray.at(ant1);

  Double_t phi2 = TMath::DegToRad()*phiArrayDeg.at(ant2);
  Double_t x2 = rArray.at(ant2)*TMath::Cos(phi2);
  Double_t y2 = rArray.at(ant2)*TMath::Sin(phi2);
  Double_t z2 = zArray.at(ant2);

  Double_t x0 = rWave*TMath::Cos(phiWave)*TMath::Cos(thetaWave);
  Double_t y0 = rWave*TMath::Sin(phiWave)*TMath::Cos(thetaWave);
  Double_t z0 = rWave*TMath::Sin(thetaWave);

  Double_t part1 = TMath::Sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0));
  Double_t part2 = TMath::Sqrt((x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0));

  Double_t tdiff = 1e9*(part2 - part1)/SPEED_OF_LIGHT; // 1e9 for seconds to nanoseconds 

  tdiff /= correlationDeltaT;
  return TMath::Nint(tdiff);
}









/************************************************************************************************************
Precalculate DeltaTs during initialization where appropriate
************************************************************************************************************/

void CrossCorrelator::do5PhiSectorCombinatorics(){
  // For checking later... 
  for(Int_t ant1=0; ant1 < NUM_SEAVEYS; ant1++){
    for(Int_t ant2=0; ant2 < NUM_SEAVEYS; ant2++){
      comboIndices[ant1][ant2] = -1;
    }
  }

  numCombos=0;
  for(Int_t ant1=0; ant1<NUM_SEAVEYS; ant1++){
    Double_t phiSect1 = ant1%NUM_PHI;	
    for(Int_t ant2=ant1+1; ant2<NUM_SEAVEYS; ant2++){
      Double_t phiSect2 = ant2%NUM_PHI;
      if(TMath::Abs(phiSect1 - phiSect2) <= DELTA_PHI_SECT || TMath::Abs(phiSect1 - phiSect2) >= (NUM_PHI-DELTA_PHI_SECT)){
	comboIndices[ant1][ant2] = numCombos;
	comboIndices[ant2][ant1] = numCombos;
	comboToAnt1s.push_back(ant1);
	comboToAnt2s.push_back(ant2);	    
	numCombos++;
      }
    }    
  }
  if(numCombos != NUM_COMBOS){
    std::cerr << "Warning! in = " << __PRETTY_FUNCTION__ << std::endl;
    std::cerr << "\tnumCombos = " << numCombos << "." << std::endl;
    std::cerr << "\tNUM_COMBOS = " << NUM_COMBOS << "." << std::endl;
    std::cerr << "\tExpecting NUM_COMBOS == numCombos." << std::endl;
  }
  
}

void CrossCorrelator::fillDeltaTLookup(){
  
  Double_t phi0 = getBin0PhiDeg();
  for(Int_t phiSector = 0; phiSector<NUM_PHI; phiSector++){
    for(Int_t phiInd = 0; phiInd < NUM_BINS_PHI; phiInd++){
      Int_t phiBin = phiSector*NUM_BINS_PHI + phiInd;
      Double_t phiDeg = phi0 + phiBin*Double_t(PHI_RANGE)/NUM_BINS_PHI;
      Double_t phiWave = TMath::DegToRad()*phiDeg;
      for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA; thetaBin++){
	Double_t thetaDeg = THETA_RANGE*((Double_t)thetaBin/NUM_BINS_THETA - 0.5);
	Double_t thetaWave = TMath::DegToRad()*thetaDeg;
	for(Int_t combo=0; combo<numCombos; combo++){
	  Int_t ant1 = comboToAnt1s.at(combo);
	  Int_t ant2 = comboToAnt2s.at(combo);
	  Int_t offset = getDeltaTExpected(ant1, ant2, phiWave, thetaWave);
	  if(offset > deltaTMax){
	    deltaTMax = offset;
	  }
	  if(offset < deltaTMin){
	    deltaTMin = offset;
	  }	  
	  deltaTs[phiBin][thetaBin][combo] = offset;
	}
      }
    }
  }  
}


Bool_t CrossCorrelator::useCombo(Int_t ant1, Int_t ant2, Int_t phiSector){

  Bool_t ant1InRange = TMath::Abs(phiSector - (ant1%NUM_PHI))<=DELTA_PHI_SECT;
  ant1InRange = ant1InRange || TMath::Abs(phiSector - (ant1%NUM_PHI))>=(NUM_PHI-DELTA_PHI_SECT);

  Bool_t ant2InRange = TMath::Abs(phiSector - (ant2%NUM_PHI))<=DELTA_PHI_SECT;
  ant2InRange = ant2InRange || TMath::Abs(phiSector - (ant2%NUM_PHI))>=(NUM_PHI-DELTA_PHI_SECT);  

  return (ant1InRange && ant2InRange);
  
}
    
void CrossCorrelator::fillCombosToUse(mapMode_t mapMode, UInt_t l3Trigger, std::vector<Int_t>* combosToUse){
  if(mapMode==kTriggered){
    std::vector<Int_t> combosToUseTemp;
    for(Int_t phiSector = 0; phiSector<NUM_PHI; phiSector++){
      UInt_t doPhiSector = ((l3Trigger >> phiSector) & 1);
      if(doPhiSector){
	for(Int_t combo=0; combo<numCombos; combo++){
	  Int_t ant1 = comboToAnt1s.at(combo);
	  Int_t ant2 = comboToAnt2s.at(combo);
	  if(useCombo(ant1, ant2, phiSector)){
	    if(std::find(combosToUseTemp.begin(), combosToUseTemp.end(), combo)==combosToUseTemp.end()){
	      combosToUseTemp.push_back(combo);
	    }
	  }
	}
      }
    }
    for(Int_t phiSector = 0; phiSector<NUM_PHI; phiSector++){
      combosToUse[phiSector] = combosToUseTemp;
    }
  }
  else if(mapMode==kGlobal){
    for(Int_t phiSector = 0; phiSector<NUM_PHI; phiSector++){
      for(Int_t combo=0; combo<numCombos; combo++){
        Int_t ant1 = comboToAnt1s.at(combo);
        Int_t ant2 = comboToAnt2s.at(combo);
        if(useCombo(ant1, ant2, phiSector)){
	  combosToUse[phiSector].push_back(combo);
        }
      }
    }
  }
}





























/************************************************************************************************************
Image generation functions.
************************************************************************************************************/


TH2D* CrossCorrelator::makeBlankImage(TString name, TString title){

  Double_t phiMin = getBin0PhiDeg();
  Double_t phiMax = phiMin + 360;
  Double_t thetaMin = -THETA_RANGE/2;
  Double_t thetaMax = THETA_RANGE/2;

  TH2D* hImage = new TH2D(name, title,
			  NUM_BINS_PHI*NUM_PHI, phiMin, phiMax,
			  NUM_BINS_THETA, thetaMin, thetaMax);
  hImage->GetXaxis()->SetTitle("Azimuth (Degrees)");
  hImage->GetYaxis()->SetTitle("Elevation (Degrees)");
  return hImage;

}

TH2D* CrossCorrelator::makeBlankZoomedImage(TString name, TString title,
					    Double_t zoomCenterPhiDeg,
					    Double_t zoomCenterThetaDeg){

  Double_t phiMin = zoomCenterPhiDeg - PHI_RANGE_ZOOM/2;
  Double_t phiMax = zoomCenterPhiDeg + PHI_RANGE_ZOOM/2;
  Double_t thetaMin = zoomCenterThetaDeg - THETA_RANGE_ZOOM/2;
  Double_t thetaMax = zoomCenterThetaDeg + THETA_RANGE_ZOOM/2;

  TH2D* hImage = new TH2D(name, title,
			  NUM_BINS_PHI_ZOOM, phiMin, phiMax,
			  NUM_BINS_THETA_ZOOM, thetaMin, thetaMax);
  hImage->GetXaxis()->SetTitle("Azimuth (Degrees)");
  hImage->GetYaxis()->SetTitle("Elevation (Degrees)");
  return hImage;

}


TH2D* CrossCorrelator::makeGlobalSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave){
  Double_t imagePeak, peakPhiDeg, peakThetaDeg;
  return makeGlobalSphericalImage(pol, rWave, imagePeak, peakPhiDeg, peakThetaDeg);
}

TH2D* CrossCorrelator::makeGlobalSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave, Double_t& imagePeak,
						Double_t& peakPhiDeg, Double_t& peakThetaDeg){
  return makeImage(pol, rWave, imagePeak, peakPhiDeg, peakThetaDeg, ALL_PHI_TRIGS,
		   kGlobal, kZoomedOut);
}


TH2D* CrossCorrelator::makeTriggeredSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave, UInt_t l3Trigger){
  Double_t imagePeak, peakPhiDeg, peakThetaDeg;
  return makeTriggeredSphericalImage(pol, rWave, imagePeak, peakPhiDeg, peakThetaDeg, l3Trigger);
}

TH2D* CrossCorrelator::makeTriggeredSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave,
						   Double_t& imagePeak, Double_t& peakPhiDeg,
						   Double_t& peakThetaDeg, UInt_t l3Trigger){
  return makeImage(pol, rWave, imagePeak, peakPhiDeg, peakThetaDeg, l3Trigger, kTriggered, kZoomedOut);
}


TH2D* CrossCorrelator::makeGlobalImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
				       Double_t& peakThetaDeg){
  return makeImage(pol, 0, imagePeak, peakPhiDeg, peakThetaDeg, ALL_PHI_TRIGS,
		   kGlobal, kZoomedOut);
}


TH2D* CrossCorrelator::makeTriggeredImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
					  Double_t& peakThetaDeg, UInt_t l3Trigger){
  return makeImage(pol, 0, imagePeak, peakPhiDeg, peakThetaDeg, l3Trigger, kTriggered, kZoomedOut);
}


TH2D* CrossCorrelator::makeZoomedImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
				       Double_t& peakThetaDeg, UInt_t l3Trigger, Double_t zoomCenterPhiDeg,
				       Double_t zoomCenterThetaDeg){
  return makeImage(pol, 0, imagePeak, peakPhiDeg, peakThetaDeg, l3Trigger,
		   kTriggered, kZoomedIn, zoomCenterPhiDeg, zoomCenterThetaDeg);
}

TH2D* CrossCorrelator::makeZoomedImage(AnitaPol::AnitaPol_t pol, UInt_t l3Trigger, Double_t zoomCenterPhiDeg,
				       Double_t zoomCenterThetaDeg){
  Double_t imagePeak, peakPhiDeg, peakThetaDeg;    
  return makeZoomedImage(pol, imagePeak, peakPhiDeg, peakThetaDeg,
			 l3Trigger, zoomCenterPhiDeg, zoomCenterThetaDeg);
}

void CrossCorrelator::createImageNameAndTitle(TString& name, TString& title, mapMode_t mapMode,
					      zoomMode_t zoomMode, Double_t rWave, AnitaPol::AnitaPol_t pol){

  // Generate unique name and title from eventNumber, polarization, spherical/plane wave, triggered vs global.
  name = "h" + mapModeNames[mapMode];
  name += rWave == 0 ? "" : "Spherical";
  name += zoomModeNames[zoomMode];
  name += pol == AnitaPol::kVertical ? "ImageV" : "ImageH";
  name += TString::Format("%u", eventNumber[pol]);

  title = TString::Format("Event %u ", eventNumber[pol]);
  title += (pol == AnitaPol::kVertical ? "VPOL" : "HPOL");
  title += rWave == 0 ? "" : TString::Format(" Spherical r=%4.2lf", rWave);
  title += " " + mapModeNames[mapMode];
  if(zoomMode==kZoomedIn){
    title += " Zoomed In";
  }
  title += " Map";
}

TH2D* CrossCorrelator::makeImage(AnitaPol::AnitaPol_t pol, Double_t rWave, Double_t& imagePeak,
				 Double_t& peakPhiDeg, Double_t& peakThetaDeg, UInt_t l3Trigger,
				 mapMode_t mapMode, zoomMode_t zoomMode, Double_t zoomCenterPhiDeg,
				 Double_t zoomCenterThetaDeg){

  TString name, title;
  createImageNameAndTitle(name, title, mapMode, zoomMode, rWave, pol);
  
  
  TH2D* hImage = NULL;
  if(zoomMode == kZoomedOut){
    hImage = makeBlankImage(name, title);
  }
  else if(zoomMode == kZoomedIn){
    hImage = makeBlankZoomedImage(name, title, zoomCenterPhiDeg, zoomCenterThetaDeg);
  }

  std::vector<Int_t> combosToUse[NUM_PHI];
  fillCombosToUse(mapMode, l3Trigger, combosToUse);
  
  imagePeak = DBL_MIN;
  Int_t peakPhiBin = -1;
  Int_t peakThetaBin = -1;
  for(Int_t phiInd = 0; phiInd < hImage->GetNbinsX(); phiInd++){
    Int_t phiBin = phiInd;
    Int_t phiSector = zoomMode==kZoomedIn ? 0 : phiBin/NUM_BINS_PHI;
    Double_t phiWave = hImage->GetXaxis()->GetBinLowEdge(phiBin+1)*TMath::DegToRad();
    for(Int_t thetaBin = 0; thetaBin < hImage->GetNbinsY(); thetaBin++){
      Double_t thetaWave = hImage->GetYaxis()->GetBinLowEdge(thetaBin+1)*TMath::DegToRad();	
      Double_t correlations = 0;
      Int_t contributors = 0;
      for(UInt_t comboInd=0; comboInd<combosToUse[phiSector].size(); comboInd++){
	Int_t combo = combosToUse[phiSector].at(comboInd);

	Int_t offset = 0;

	// If we are in zoomed out & plane wave mode then use the lookup table for a big speedup
	if(zoomMode==kZoomedOut){
	  offset = deltaTs[phiBin][thetaBin][combo];
	}
	// If we are in zoomed in & plane wave mode then calculate
	else if(rWave==0){
	  Int_t ant1 = comboToAnt1s.at(combo);
	  Int_t ant2 = comboToAnt2s.at(combo);
	  offset = getDeltaTExpected(ant1, ant2, phiWave, thetaWave);
	}
	// If we are in zoomed out & spherical wave mode then calculate	  
	else{
	  Int_t ant1 = comboToAnt1s.at(combo);
	  Int_t ant2 = comboToAnt2s.at(combo);
	  offset = getDeltaTExpectedSpherical(ant1, ant2, phiWave, thetaWave, rWave);
	}

	offset = offset < 0 ? offset + numSamplesUpsampled : offset;
	correlations += crossCorrelations[pol][combo][offset];
	contributors++;
      }
      if(contributors>0){
	correlations /= contributors;
      }
      hImage->SetBinContent(phiBin + 1, thetaBin + 1, correlations);
      if(correlations > imagePeak){
	imagePeak = correlations;
	peakPhiBin = phiBin;
	peakThetaBin = thetaBin;
      }
    }
  }

  peakPhiDeg = hImage->GetXaxis()->GetBinLowEdge(peakPhiBin+1);
  peakThetaDeg = hImage->GetYaxis()->GetBinLowEdge(peakThetaBin+1);

  return hImage;
}

Double_t CrossCorrelator::findImagePeak(TH2D* hist, Double_t& imagePeakTheta, Double_t& imagePeakPhi){

  Int_t nx = hist->GetNbinsX();
  Int_t ny = hist->GetNbinsY();

  Double_t maxVal = -2;
  for(Int_t by = 1; by<=ny; by++){
    for(Int_t bx = 1; bx<=nx; bx++){
      Double_t val = hist->GetBinContent(bx, by);
      if(val > maxVal){
	maxVal = val;
	imagePeakPhi = hist->GetXaxis()->GetBinLowEdge(bx);
	imagePeakTheta = hist->GetYaxis()->GetBinLowEdge(by);
      }
    }
  }

  return maxVal;
  
}





/************************************************************************************************************
Functions to delete pointers to internal variables
************************************************************************************************************/
void CrossCorrelator::deleteCrossCorrelations(AnitaPol::AnitaPol_t pol){
  for(Int_t comboInd=0; comboInd<NUM_COMBOS; comboInd++){
    if(crossCorrelations[pol][comboInd] != NULL){
      delete [] crossCorrelations[pol][comboInd];
      crossCorrelations[pol][comboInd] = NULL;
    }
  }
}

void CrossCorrelator::deleteAllWaveforms(AnitaPol::AnitaPol_t pol){
  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
    if(grs[pol][ant]){
      delete grs[pol][ant];
      grs[pol][ant] = NULL;
    }
    if(grsInterp[pol][ant]){
      delete grsInterp[pol][ant];
      grsInterp[pol][ant] = NULL;
    }
  }
}


















/************************************************************************************************************
Functions for debugging or testing
************************************************************************************************************/

void CrossCorrelator::correlateEventTest(Double_t phiDegSource, Double_t thetaDegSource, Double_t rSource){
  // Generates a set of delta function like waveforms, correlates them 

  // Delete any old waveforms (at start rather than end to leave in memory to be examined if need be)

  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){  
    deleteAllWaveforms((AnitaPol::AnitaPol_t)pol);
  }


  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      std::vector<Double_t> newVs = std::vector<Double_t>(numSamplesUpsampled, 0);
      std::vector<Double_t> newTimes = std::vector<Double_t>(numSamplesUpsampled, 0);      
      Int_t dt = numSamplesUpsampled/2 + getDeltaTExpectedSpherical(0, ant, 
								    phiDegSource*TMath::DegToRad(), 
								    thetaDegSource*TMath::DegToRad(), 
								    rSource);
      for(int samp=0; samp<numSamplesUpsampled; samp++){
	newTimes[samp] = correlationDeltaT*samp;
	if(samp==dt){
	  newVs[samp] = numSamplesUpsampled;
	}
      }
      grsInterp[pol][ant] = new TGraph(numSamplesUpsampled, &newTimes[0], &newVs[0]);    
      RootTools::normalize(grsInterp[pol][ant]);
    }
    doAllCrossCorrelations((AnitaPol::AnitaPol_t)pol);
  }
}


TGraph* CrossCorrelator::getCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2){
  // Primarily for debugging, put cross correlations in a TGraph 

  Int_t comboInd = comboIndices[ant1][ant2];
  if(comboInd < 0 ){
    return NULL;
  }
  std::vector<Double_t> offsets = std::vector<Double_t>(numSamplesUpsampled, 0);
  std::vector<Double_t> corrs = std::vector<Double_t>(numSamplesUpsampled, 0);  

  for(Int_t i=0; i<numSamplesUpsampled; i++){
    Int_t offset = (i - numSamplesUpsampled/2);
    offsets.at(i) = offset*correlationDeltaT;
    Int_t j = offset < 0 ? offset + numSamplesUpsampled : offset;
    corrs.at(i) = crossCorrelations[pol][comboInd][j];

  }


  TGraph* gr = new TGraph(numSamplesUpsampled,  &offsets[0], &corrs[0]);
  gr->SetName(TString::Format("grCorr_%d_%d", ant1, ant2));
  gr->SetTitle(TString::Format("Cross Correlation ant1 = %d, ant2 = %d", ant1, ant2));

  return gr;
}


