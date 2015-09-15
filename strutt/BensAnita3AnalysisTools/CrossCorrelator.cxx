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
  \param upsampleFactorTemp is the upsample factor, by default = 1.
*/
CrossCorrelator::CrossCorrelator(Int_t upsampleFactorTemp){
  initializeVariables(upsampleFactorTemp);
}

/*!
  \brief Destructor
*/
CrossCorrelator::~CrossCorrelator(){
  deleteAllWaveforms();
  deleteCrossCorrelations();
}


/*!
  \brief Workhorse function to set internal variables.
  /param upSampleFactorTemp is passed here from the constructor.
*/
void CrossCorrelator::initializeVariables(Int_t upSampleFactorTemp){

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
  }
  lastEventNormalized = 0;
  eventNumber = 0;

  nominalSamplingDeltaT = 1./2.6;
  upsampleFactor = upSampleFactorTemp;
  correlationDeltaT = nominalSamplingDeltaT/upsampleFactor;
  numSamplesUpsampled = 2*NUM_SAMPLES*upsampleFactor;

  // Fill geom, timing arrays and combinatorics
  // const Double_t rArrayTemp[NUM_SEAVEYS] = {0.9675,0.7402,0.9675,0.7402,0.9675,0.7402,0.9675,0.7402,
  // 					    0.9675,0.7402,0.9675,0.7402,0.9675,0.7402,0.9675,0.7402,
  // 					    2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,
  // 					    2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,
  // 					    2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,
  // 					    2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447};
  
  // const Double_t phiArrayDegTemp[NUM_SEAVEYS] = {0.0,22.5,45.0,67.5,90.0,112.5,135.0,157.5,
  // 						 180.0,202.5,225.0,247.5,270.0,292.5,315.0,337.5,
  // 						 0.0,22.5,45.0,67.5,90.0,112.5,135.0,157.5,
  // 						 180.0,202.5,225.0,247.5,270.0,292.5,315.0,337.5,
  // 						 0.0,22.5,45.0,67.5,90.0,112.5,135.0,157.5,180.0,
  // 						 202.5,225.0,247.5,270.0,292.5,315.0,337.5};

  // const Double_t zArrayTemp[NUM_SEAVEYS] = {-1.4407,-2.4135,-1.4407,-2.4135,-1.4407,-2.4135,-1.4407,-2.4135,
  // 					    -1.4407,-2.4135,-1.4407,-2.4135,-1.4407,-2.4135,-1.4407,-2.4135,
  // 					    -5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,
  // 					    -5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,
  // 					    -6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,
  // 					    -6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951};

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    rArray.push_back(geom->getAntR(ant));
    zArray.push_back(geom->getAntZ(ant));
    phiArrayDeg.push_back(geom->getAntPhiPosition(ant)*TMath::RadToDeg());
  }

  do5PhiSectorCombinatorics();

  // Here we try to read the dts from binary file (for speed)
  Int_t readFileSuccess = readDeltaTsFile();
  if(readFileSuccess != 0){
    fillDeltaTLookup();
    writeDeltaTsFile();
  }
}

















/************************************************************************************************************
Waveform manipulation functions
************************************************************************************************************/


/*!
  \brief Loops through all waveform graphs in the UsefulAnitaEvent and makes an evenly re-sampled, normalized copy of each one.
  \param usefulEvent points to the UsefulAnitaEvent of interest
*/
void CrossCorrelator::getNormalizedInterpolatedTGraphs(UsefulAnitaEvent* usefulEvent){
  // Potentially needed in a few places, so it gets its own function 

  // Pretty much just for profiling 
  if(usefulEvent->eventNumber!=lastEventNormalized){

    // Delete any old waveforms (at start rather than end to leave in memory to be examined if need be)
    deleteAllWaveforms();

    // Find the start time of all waveforms 
    Double_t earliestStart[NUM_POL] = {100};
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      grs[AnitaPol::kVertical][ant] = usefulEvent->getGraph(ant, AnitaPol::kVertical);
      grs[AnitaPol::kHorizontal][ant] = usefulEvent->getGraph(ant, AnitaPol::kHorizontal);

      for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
	if(grs[pol][ant]->GetX()[0]<earliestStart[pol]){
	  earliestStart[pol] = grs[pol][ant]->GetX()[0];
	}
      }
    }

    // Interpolate with earliest start time 
    for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
      for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
	grsInterp[pol][ant] = interpolateWithStartTime(grs[pol][ant], earliestStart[pol]);
	Double_t mean; // don't care about this.
	RootTools::normalize(grsInterp[pol][ant], mean, interpRMS[pol][ant]);
      }
    }
    lastEventNormalized = usefulEvent->eventNumber;
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
  \brief Loops through the set of cross correlations and returns a vector of vectors containing the maximum times
*/
std::vector<std::vector<Double_t> > CrossCorrelator::getMaxCorrelationTimes(){
  
  std::vector<std::vector<Double_t> > peakTimes(NUM_POL, std::vector<Double_t>(NUM_COMBOS, 0));

  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    for(Int_t combo=0; combo<NUM_COMBOS; combo++){
      Int_t maxInd = RootTools::getIndexOfMaximum(numSamplesUpsampled, crossCorrelations[pol][combo]);
      Int_t offset = maxInd >= numSamplesUpsampled/2 ? maxInd - numSamplesUpsampled : maxInd;
	//	offset = offset < 0 ? offset + numSamplesUpsampled : offset;
      peakTimes.at(pol).at(combo) = offset*correlationDeltaT;      
    }
  }
  return peakTimes;
}

/*!
  \brief Loops through the set of cross correlations and returns a vector of vectors containing the maximum correlation values
*/
std::vector<std::vector<Double_t> > CrossCorrelator::getMaxCorrelationValues(){
  
  std::vector<std::vector<Double_t> > peakVals(NUM_POL, std::vector<Double_t>(NUM_COMBOS, 0));

  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    for(Int_t combo=0; combo<NUM_COMBOS; combo++){
      Int_t maxInd = RootTools::getIndexOfMaximum(numSamplesUpsampled, crossCorrelations[pol][combo]);
      peakVals.at(pol).at(combo) = crossCorrelations[pol][combo][maxInd];
    }
  }
  return peakVals;
}




Double_t CrossCorrelator::correlationWithOffset(TGraph* gr1, TGraph* gr2, Int_t offset){
  // Time domain cross correlation, no longer used 
  Double_t correlation = 0;
  for(Int_t i=0; i<numSamplesUpsampled; i++){
    // Modulo numSamplesUpsampled wraps for circular correlation 
    Int_t j = (i + offset)%numSamplesUpsampled; 
    correlation += gr1->GetY()[i]*gr2->GetY()[j];
  }
  correlation/=numSamplesUpsampled;  
  return correlation;
}

void CrossCorrelator::correlateEvent(UsefulAnitaEvent* usefulEvent){
  // Read TGraphs from events into memory (also deletes old TGraphs) 
  getNormalizedInterpolatedTGraphs(usefulEvent);

  // Cross correlate waveforms using the TGraphs we just obtained (also deletes old cross correlations) 
  doAllCrossCorrelations();

  eventNumber = usefulEvent->eventNumber;
}



void CrossCorrelator::doAllCrossCorrelations(){

  // Delete old cross correlations first 
  deleteCrossCorrelations();


  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    for(Int_t ant1=0; ant1<NUM_SEAVEYS; ant1++){
      for(UInt_t ant2Ind=0; ant2Ind<ant2s[ant1].size(); ant2Ind++){
	Int_t ant2 = ant2s[ant1].at(ant2Ind);
	Int_t comboInd = comboIndices[ant1][ant2];
	if(!doneCrossCorrelations[pol][comboInd]){
	  crossCorrelations[pol][comboInd] = crossCorrelateFourier(grsInterp[pol][ant1], grsInterp[pol][ant2]);
	  
	  doneCrossCorrelations[pol][comboInd] = 1;
	}
      }
    }    
  }
}


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



Int_t CrossCorrelator::getDeltaTExpected(Int_t ant1, Int_t ant2, Int_t phiBin, Int_t thetaBin){
  Double_t part1 = -zArray.at(ant1)*tanThetaLookup[thetaBin] - rArray.at(ant1)*(cosPhiWaveLookup[phiBin]*cosPhiArrayLookup[ant1] + sinPhiWaveLookup[phiBin]*sinPhiArrayLookup[ant1]);
  Double_t part2 = -zArray.at(ant2)*tanThetaLookup[thetaBin] - rArray.at(ant2)*(cosPhiWaveLookup[phiBin]*cosPhiArrayLookup[ant2] + sinPhiWaveLookup[phiBin]*sinPhiArrayLookup[ant2]);
  Double_t tdiff = 1e9*((cosThetaLookup[thetaBin] * (part2 - part1))/SPEED_OF_LIGHT); // Returns time in ns
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

  Int_t numCombos=0;
  for(Int_t ant1=0; ant1 < NUM_SEAVEYS; ant1++){
    Int_t phiSect1 = ant1%NUM_PHI;
    for(Int_t deltaPhiSect=-2; deltaPhiSect<=2; deltaPhiSect++){
    //    for(Int_t deltaPhiSect=-3; deltaPhiSect<=3; deltaPhiSect++){
      for(Int_t ring=0; ring<NUM_RING; ring++){
	Int_t phiSect2 = phiSect1 + deltaPhiSect;
	phiSect2 = phiSect2 < 0 ? phiSect2 + NUM_PHI : phiSect2;
	phiSect2 = phiSect2 >= NUM_PHI ? phiSect2 - NUM_PHI : phiSect2;
	Int_t ant2 = phiSect2 + ring*NUM_PHI;
	if(ant1 != ant2){
	  ant2s[ant1].push_back(ant2);	  
	  if(comboIndices[ant1][ant2] < 0){
	    comboIndices[ant1][ant2] = numCombos;
	    comboIndices[ant2][ant1] = numCombos;
	    comboToAnt1s.push_back(ant1);
	    comboToAnt2s.push_back(ant2);	    
	    numCombos++;
	  }
	}
      }
    }
  }
  
  if(numCombos != NUM_COMBOS){
    std::cerr << "numCombos = " << numCombos
	      << ", expecting NUM_COMBOS = " << NUM_COMBOS
	      << ". Check the combinatorics... " << std::endl;
  }
  assert(numCombos==NUM_COMBOS);

}

void CrossCorrelator::fillDeltaTLookup(){

  for(Int_t thetaBin=0; thetaBin < NUM_BINS_THETA; thetaBin++){
    Double_t thetaDeg = THETA_RANGE*((Double_t)thetaBin/NUM_BINS_THETA - 0.5);
    Double_t thetaWave = TMath::DegToRad()*thetaDeg;
    tanThetaLookup[thetaBin] = TMath::Tan(thetaWave);
    cosThetaLookup[thetaBin] = TMath::Cos(thetaWave);
  }

  for(Int_t phiBin=0; phiBin < NUM_BINS_PHI*NUM_PHI; phiBin++){
    // Double_t phiDeg = -0.5*PHI_RANGE + phiBin*Double_t(PHI_RANGE)/NUM_BINS_PHI;
    Double_t phiDeg = -0.5*PHI_RANGE + phiArrayDeg.at(0) + phiBin*Double_t(PHI_RANGE)/NUM_BINS_PHI;
    Double_t phiWave = TMath::DegToRad()*phiDeg;
    //    phiWaveLookup[phiBin] = phiWave;
    cosPhiWaveLookup[phiBin] = TMath::Cos(phiWave);
    sinPhiWaveLookup[phiBin] = TMath::Sin(phiWave);
  }
  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
    Double_t phiRad = TMath::DegToRad()*phiArrayDeg.at(ant);
    cosPhiArrayLookup[ant] = TMath::Cos(phiRad);
    sinPhiArrayLookup[ant] = TMath::Sin(phiRad);
  }

  for(Int_t phiSector = 0; phiSector<NUM_PHI; phiSector++){
    for(Int_t phiInd = 0; phiInd < NUM_BINS_PHI; phiInd++){
      Int_t phiBin = phiSector*NUM_BINS_PHI + phiInd;
      // Double_t phiDeg = -0.5*PHI_RANGE + phiBin*Double_t(PHI_RANGE)/NUM_BINS_PHI;
      Double_t phiDeg = -0.5*PHI_RANGE + phiArrayDeg.at(0) + phiBin*Double_t(PHI_RANGE)/NUM_BINS_PHI;      
      Double_t phiWave = TMath::DegToRad()*phiDeg;
      for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA; thetaBin++){
	Double_t thetaDeg = THETA_RANGE*((Double_t)thetaBin/NUM_BINS_THETA - 0.5);
	Double_t thetaWave = TMath::DegToRad()*thetaDeg;
	for(Int_t ant1=0; ant1<NUM_SEAVEYS; ant1++){
	  for(UInt_t ant2Ind=0; ant2Ind<ant2s[ant1].size(); ant2Ind++){
	    Int_t ant2 = ant2s[ant1].at(ant2Ind);
	    Int_t comboInd = comboIndices[ant1][ant2];
	    Int_t offset = getDeltaTExpected(ant1, ant2, phiWave, thetaWave);
	    offset = offset < 0 ? offset + numSamplesUpsampled : offset;
	    deltaTs[comboInd][phiBin][thetaBin] = offset;
	    if(thetaBin==0 && phiBin==0 && ant1==0){
	      std::cout << thetaBin << "\t" << phiBin << "\t" << ant1 << "\t"
			<< ant2 << "\t" << offset << std::endl;
	    }

	  }
	}
      }
    }
  }  
}






























/************************************************************************************************************
Image generation functions.
************************************************************************************************************/



TH2D* CrossCorrelator::makeImage(AnitaPol::AnitaPol_t pol, UInt_t l3Trigger){
  Double_t imagePeak, peakPhiDeg, peakThetaDeg;
  return makeImage(pol, imagePeak, peakPhiDeg, peakThetaDeg, l3Trigger);
}

TH2D* CrossCorrelator::makeImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg, Double_t& peakThetaDeg, UInt_t l3Trigger){

  imagePeak = DBL_MIN;

  assert(pol == AnitaPol::kVertical || pol == AnitaPol::kHorizontal);
  TString name = pol == AnitaPol::kVertical ? "hImageV" : "hImageH";
  name += TString::Format("%u", eventNumber);

  TString title = TString::Format("Event %u ", eventNumber);
  title += (pol == AnitaPol::kVertical ? "VPOL" : "HPOL");

  TH2D* hImage = new TH2D(name, title,
			  NUM_BINS_PHI*NUM_PHI, -PHI_RANGE/2, PHI_RANGE*NUM_PHI-PHI_RANGE/2,
			  NUM_BINS_THETA, -THETA_RANGE/2, THETA_RANGE/2);
  hImage->GetXaxis()->SetTitle("Azimuth (Degrees)");
  hImage->GetYaxis()->SetTitle("Elevation (Degrees)");

  Int_t peakPhiBin = -1;
  Int_t peakThetaBin = -1;
  
  for(Int_t phiSector = 0; phiSector<NUM_PHI; phiSector++){
    UInt_t doPhiSector = ((l3Trigger >> phiSector) & 1);
    if(doPhiSector){
      for(Int_t phiInd = 0; phiInd < NUM_BINS_PHI; phiInd++){
	Int_t phiBin = phiSector*NUM_BINS_PHI + phiInd;
	for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA; thetaBin++){
	  Double_t correlations = 0;
	  Int_t contributors = 0;
	  for(Int_t ant1=phiSector; ant1<NUM_SEAVEYS; ant1+=NUM_PHI){
	  // for(Int_t ant1=1; ant1<NUM_SEAVEYS; ant1+=NUM_PHI){
	    for(UInt_t ant2Ind=0; ant2Ind<ant2s[ant1].size(); ant2Ind++){
	      Int_t ant2 = ant2s[ant1].at(ant2Ind);
	      Int_t comboInd = comboIndices[ant1][ant2];
	      if(comboInd > 0){
		Int_t offset = deltaTs[comboInd][phiBin][thetaBin];
		offset*=upsampleFactor; // Correct for variable upsample factor
		correlations += crossCorrelations[pol][comboInd][offset];
		contributors++;
	      }
	    }
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
    }
  }

  peakPhiDeg = hImage->GetXaxis()->GetBinLowEdge(peakPhiBin+1);
  peakThetaDeg = hImage->GetYaxis()->GetBinLowEdge(peakThetaBin+1);

  return hImage;
}


TH2D* CrossCorrelator::makeImageSpherical(AnitaPol::AnitaPol_t pol, Double_t rWave){
  Double_t imagePeak, peakPhiDeg, peakThetaDeg;
  return makeImageSpherical(pol, rWave, imagePeak, peakPhiDeg, peakThetaDeg);
}

TH2D* CrossCorrelator::makeImageSpherical(AnitaPol::AnitaPol_t pol, Double_t rWave, Double_t& imagePeak, Double_t& peakPhiDeg, Double_t& peakThetaDeg){

  imagePeak = -100;

  assert(pol == AnitaPol::kVertical || pol == AnitaPol::kHorizontal);
  TString name = pol == AnitaPol::kVertical ? "hImageV" : "hImageH";
  name += TString::Format("%u", eventNumber);
  name += TString::Format("_r%lf", rWave);

  TString title = TString::Format("Event %u ", eventNumber);
  title += (pol == AnitaPol::kVertical ? "VPOL" : "HPOL");  
  title += TString::Format(" (spherical #deltats r = %lf m)", rWave);

  TH2D* hImage = new TH2D(name, title,
			  NUM_BINS_PHI*NUM_PHI, phiArrayDeg.at(0)-PHI_RANGE/2,
			  phiArrayDeg.at(15)+PHI_RANGE/2,
			  NUM_BINS_THETA, -THETA_RANGE/2,
			  THETA_RANGE/2);
  hImage->GetXaxis()->SetTitle("Azimuth (Degrees)");
  hImage->GetYaxis()->SetTitle("Elevation (Degrees)");

  Int_t peakPhiBin = -1;
  Int_t peakThetaBin = -1;
  
  for(Int_t phiSector = 0; phiSector<NUM_PHI; phiSector++){
    for(Int_t phiInd = 0; phiInd < NUM_BINS_PHI; phiInd++){
      Int_t phiBin = phiSector*NUM_BINS_PHI + phiInd;
      Double_t phiWave = TMath::DegToRad()*hImage->GetXaxis()->GetBinLowEdge(phiBin+1);
      for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA; thetaBin++){
	Double_t thetaWave = TMath::DegToRad()*hImage->GetYaxis()->GetBinLowEdge(thetaBin+1);
        Double_t correlations = 0;
        Int_t contributors = 0;
	for(Int_t ant1=phiSector; ant1<NUM_SEAVEYS; ant1+=NUM_PHI){
	// for(Int_t ant1=0; ant1<NUM_SEAVEYS; ant1++){
	  for(UInt_t ant2Ind=0; ant2Ind<ant2s[ant1].size(); ant2Ind++){
	    Int_t ant2 = ant2s[ant1].at(ant2Ind);
	    Int_t comboInd = comboIndices[ant1][ant2];
	    if(comboInd < 0) continue;
	    Int_t offset = getDeltaTExpectedSpherical(ant1, ant2, phiWave, thetaWave, rWave);

	    // Not obvious, but since we correlate ant1 with ant2 but not ant2 with ant1 
	    // we need to correct the sign of deltaT in the case that ant1 < ant2.

	    if(ant1 < ant2){
	      offset *= -1;
	    }
	    offset = offset < 0 ? offset + numSamplesUpsampled : offset;
	    correlations += crossCorrelations[pol][comboInd][offset];
	    contributors++;
	  }
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
void CrossCorrelator::deleteCrossCorrelations(){
  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    for(Int_t comboInd=0; comboInd<NUM_COMBOS; comboInd++){
      doneCrossCorrelations[pol][comboInd] = 0;
      if(crossCorrelations[pol][comboInd] != NULL){
	delete [] crossCorrelations[pol][comboInd];
	crossCorrelations[pol][comboInd] = NULL;
      }
    }
  }
}

void CrossCorrelator::deleteAllWaveforms(){
  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
    for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
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
}


















/************************************************************************************************************
Functions for debugging or testing
************************************************************************************************************/

void CrossCorrelator::correlateEventTest(Double_t phiDegSource, Double_t thetaDegSource, Double_t rSource){
  // Generates a set of delta function like waveforms, correlates them 

  // Delete any old waveforms (at start rather than end to leave in memory to be examined if need be)
  deleteAllWaveforms();

  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
    for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
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
  }

  doAllCrossCorrelations();
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


void CrossCorrelator::writeDeltaTsFile(){
  // Write sets of deltaTs into a file to try and speed up class initialization by reading it 

  const char* anitaUtilEnv = "ANITA_UTIL_INSTALL_DIR";
  // std::cout << anitaUtilEnv << std::endl;
  const char* dtsDir = getenv(anitaUtilEnv);
  char dtsFileName[FILENAME_MAX];
  sprintf(dtsFileName, "%s/crossCorrelator.dts", dtsDir);

  FILE* dtsFile = fopen(dtsFileName, "w");
  fwrite(deltaTs,sizeof(unsigned char)*NUM_COMBOS*NUM_PHI*NUM_BINS_PHI*NUM_BINS_THETA, 1, dtsFile);
  fclose(dtsFile);
}

Int_t CrossCorrelator::readDeltaTsFile(){
  // Returns 0 on success 

  const char* anitaUtilEnv = "ANITA_UTIL_INSTALL_DIR";
  const char* dtsDir = getenv(anitaUtilEnv);
  char dtsFileName[FILENAME_MAX];
  sprintf(dtsFileName, "%s/crossCorrelator.dts", dtsDir);

  Int_t successState = 0;
  FILE* dtsFile = fopen(dtsFileName, "r");
  if(dtsFile != NULL){
    UInt_t numBytes = sizeof(unsigned char)*NUM_COMBOS*NUM_PHI*NUM_BINS_PHI*NUM_BINS_THETA;
    fread(deltaTs,1,numBytes,dtsFile);
    fclose(dtsFile);
  }
  else{
    successState = 1;
  }
  return successState;
}
