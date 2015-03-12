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
CrossCorrelator::CrossCorrelator(){

  /* Initialize with NULL pointers otherwise very bad things will happen with gcc */
  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      grs[pol][ant] = NULL;
      grsInterp[pol][ant] = NULL;
    }
    for(int combo=0; combo<NUM_COMBOS; combo++){
      crossCorrelations[pol][combo] = NULL;
    }
  }
  lastEventNormalized = 0;
  eventNumber = 0;
  correlationDeltaT = 1./2.6;

  /* Fill geom, timing arrays and combinatorics*/
  const Double_t rArrayTemp[NUM_SEAVEYS] = {0.9675,0.7402,0.9675,0.7402,0.9675,0.7402,0.9675,0.7402,
  					    0.9675,0.7402,0.9675,0.7402,0.9675,0.7402,0.9675,0.7402,
  					    2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,
  					    2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,
  					    2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,
  					    2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447,2.0447};
  
  const Double_t phiArrayDegTemp[NUM_SEAVEYS] = {0.0,22.5,45.0,67.5,90.0,112.5,135.0,157.5,
  						 180.0,202.5,225.0,247.5,270.0,292.5,315.0,337.5,
  						 0.0,22.5,45.0,67.5,90.0,112.5,135.0,157.5,
  						 180.0,202.5,225.0,247.5,270.0,292.5,315.0,337.5,
  						 0.0,22.5,45.0,67.5,90.0,112.5,135.0,157.5,180.0,
  						 202.5,225.0,247.5,270.0,292.5,315.0,337.5};

  const Double_t zArrayTemp[NUM_SEAVEYS] = {-1.4407,-2.4135,-1.4407,-2.4135,-1.4407,-2.4135,-1.4407,-2.4135,
  					    -1.4407,-2.4135,-1.4407,-2.4135,-1.4407,-2.4135,-1.4407,-2.4135,
  					    -5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,
  					    -5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,-5.1090,
  					    -6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,
  					    -6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951,-6.1951};
  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    rArray[ant] = rArrayTemp[ant];
    zArray[ant] = zArrayTemp[ant];
    phiArrayDeg[ant] = phiArrayDegTemp[ant];
  }

  offsetIndGPU = fillDeltaTLookupGPU();
  do5PhiSectorCombinatorics();
  fillDeltaTLookup();
}

CrossCorrelator::~CrossCorrelator(){
  free(offsetIndGPU);

  deleteAllWaveforms();
  deleteCrossCorrelations();

}
















/************************************************************************************************************
Waveform manipulation functions
************************************************************************************************************/


void CrossCorrelator::getNormalizedInterpolatedTGraphs(UsefulAnitaEvent* realEvent){
  /* Potentially needed in a few places, so it gets its own function */

  /* Pretty much just for profiling */
  if(realEvent->eventNumber!=lastEventNormalized){

    /* Delete any old waveforms (at start rather than end to leave in memory to be examined if need be)*/
    deleteAllWaveforms();

    /* Find the start time of all waveforms */
    Double_t earliestStart[NUM_POL] = {100};
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      grs[AnitaPol::kVertical][ant] = realEvent->getGraph(ant, AnitaPol::kVertical);
      grs[AnitaPol::kHorizontal][ant] = realEvent->getGraph(ant, AnitaPol::kHorizontal);

      for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
	if(grs[pol][ant]->GetX()[0]<earliestStart[pol]){
	  earliestStart[pol] = grs[pol][ant]->GetX()[0];
	}
      }
    }

    /* Interpolate with earliest start time */
    for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
      for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
	grsInterp[pol][ant] = interpolateWithStartTime(grs[pol][ant], earliestStart[pol]);
	RootTools::normalize(grsInterp[pol][ant]);
      }
    }

    lastEventNormalized = realEvent->eventNumber;
  }
}

TGraph* CrossCorrelator::interpolateWithStartTime(TGraph* grIn, Double_t startTime){

  Double_t newTimes[NUM_SAMPLES] = {0};
  Double_t newVolts[NUM_SAMPLES] = {0};
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
  for(Int_t samp = 0; samp < NUM_SAMPLES; samp++){
    newTimes[samp] = time;
    if(time >= thisStartTime && time <= lastTime){
      newVolts[samp] = chanInterp.Eval(time);
    }
    else{
      newVolts[samp] = 0;
    }
    time += correlationDeltaT;

  }

  return new TGraph(NUM_SAMPLES, newTimes, newVolts);

}





















/************************************************************************************************************
All correlation functions
************************************************************************************************************/


Double_t CrossCorrelator::correlationWithOffset(TGraph* gr1, TGraph* gr2, Int_t offset){
  /* Time domain cross correlation, no longer used */
  Double_t correlation = 0;
  for(Int_t i=0; i<NUM_SAMPLES; i++){
    /* Modulo NUM_SAMPLES wraps for circular correlation */
    Int_t j = (i + offset)%NUM_SAMPLES; 
    correlation += gr1->GetY()[i]*gr2->GetY()[j];
  }
  correlation/=NUM_SAMPLES;  
  return correlation;
}

void CrossCorrelator::correlateEvent(UsefulAnitaEvent* realEvent){
  /* Read TGraphs from events into memory (also deletes old TGraphs) */
  getNormalizedInterpolatedTGraphs(realEvent);

  /* Cross correlate waveforms using the TGraphs we just obtained (also deletes old cross correlations) */
  doAllCrossCorrelations();

  eventNumber = realEvent->eventNumber;
}


void CrossCorrelator::correlateEventGPU(UsefulAnitaEvent* realEvent){

  getNormalizedInterpolatedTGraphs(realEvent);

  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    for(Int_t phiSectorInd=0; phiSectorInd<NUM_PHI; phiSectorInd++){
      for(Int_t localCombo=0; localCombo<LOCAL_COMBOS_PER_PHI_GPU; localCombo++){

	Int_t comboLocal = localCombo;
	Int_t phiSector = (phiSectorInd + 15 + comboLocal/GLOBAL_COMBOS_PER_PHI_GPU)%NUM_PHI;
	UInt_t globalCombo = comboLocal % GLOBAL_COMBOS_PER_PHI_GPU;

	Int_t a1 = ant1Gpu[phiSectorInd][localCombo];
	Int_t a2 = ant2Gpu[phiSectorInd][localCombo];

	for(Int_t offsetInd=0; offsetInd<NUM_SAMPLES; offsetInd++){
	  correlationsGPU[pol][phiSector][globalCombo][offsetInd] = correlationWithOffset(grsInterp[pol][a1], grsInterp[pol][a2], offsetInd);
	}
      }
    }
  }
}

void CrossCorrelator::doAllCrossCorrelations(){

  /* Delete old cross correlations first */
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
  /* Generate cross correlations, now using FancyFFTs */
  //  return  FFTtools::getCorrelation(gr1->GetN(), gr1->GetY(), gr2->GetY()) ;
  return  FancyFFTs::crossCorrelate(gr1->GetN(), gr1->GetY(), gr2->GetY()) ;
}
















/************************************************************************************************************
Calculate deltaT between two antennas (for a plane wave unless function name says otherwise)
************************************************************************************************************/

Int_t CrossCorrelator::getDeltaTExpected(Int_t ant1, Int_t ant2, Double_t phiWave, Double_t thetaWave){
  Double_t tanThetaW = tan(-thetaWave);
  Double_t part1 = zArray[ant1]*tanThetaW - rArray[ant1]*cos(phiWave-TMath::DegToRad()*phiArrayDeg[ant1]);
  Double_t part2 = zArray[ant2]*tanThetaW - rArray[ant2]*cos(phiWave-TMath::DegToRad()*phiArrayDeg[ant2]);
  Double_t tdiff = 1e9*((cos(-thetaWave) * (part2 - part1))/SPEED_OF_LIGHT); // Returns time in ns

  tdiff /= correlationDeltaT;
  return TMath::Nint(tdiff);
}


Int_t CrossCorrelator::getDeltaTExpectedSpherical(Int_t ant1, Int_t ant2, Double_t phiWave, Double_t thetaWave, Double_t rWave){

  Double_t phi1 = TMath::DegToRad()*phiArrayDeg[ant1];
  Double_t x1 = rArray[ant1]*TMath::Cos(phi1);
  Double_t y1 = rArray[ant1]*TMath::Sin(phi1);
  Double_t z1 = zArray[ant1];

  Double_t phi2 = TMath::DegToRad()*phiArrayDeg[ant2];
  Double_t x2 = rArray[ant2]*TMath::Cos(phi2);
  Double_t y2 = rArray[ant2]*TMath::Sin(phi2);
  Double_t z2 = zArray[ant2];

  Double_t x0 = rWave*TMath::Cos(phiWave)*TMath::Cos(thetaWave);
  Double_t y0 = rWave*TMath::Sin(phiWave)*TMath::Cos(thetaWave);
  Double_t z0 = rWave*TMath::Sin(thetaWave);

  Double_t part1 = TMath::Sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0));
  Double_t part2 = TMath::Sqrt((x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0));

  Double_t tdiff = 1e9*(part2 - part1)/SPEED_OF_LIGHT; /* 1e9 for seconds to nanoseconds */

  tdiff /= correlationDeltaT;
  return TMath::Nint(tdiff);
}



Int_t CrossCorrelator::getDeltaTExpected(Int_t ant1, Int_t ant2, Int_t phiBin, Int_t thetaBin){
  Double_t part1 = -zArray[ant1]*tanThetaLookup[thetaBin] - rArray[ant1]*(cosPhiWaveLookup[phiBin]*cosPhiArrayLookup[ant1] + sinPhiWaveLookup[phiBin]*sinPhiArrayLookup[ant1]);
  Double_t part2 = -zArray[ant2]*tanThetaLookup[thetaBin] - rArray[ant2]*(cosPhiWaveLookup[phiBin]*cosPhiArrayLookup[ant2] + sinPhiWaveLookup[phiBin]*sinPhiArrayLookup[ant2]);
  Double_t tdiff = 1e9*((cosThetaLookup[thetaBin] * (part2 - part1))/SPEED_OF_LIGHT); // Returns time in ns
  tdiff /= correlationDeltaT;
  return TMath::Nint(tdiff);
}













/************************************************************************************************************
Precalculate DeltaTs during initialization where appropriate
************************************************************************************************************/

void CrossCorrelator::do5PhiSectorCombinatorics(){
  /* For checking later... */
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
    Double_t phiDeg = -0.5*PHI_RANGE + phiBin*Double_t(PHI_RANGE)/NUM_BINS_PHI;
    Double_t phiWave = TMath::DegToRad()*phiDeg;
    //    phiWaveLookup[phiBin] = phiWave;
    cosPhiWaveLookup[phiBin] = TMath::Cos(phiWave);
    sinPhiWaveLookup[phiBin] = TMath::Sin(phiWave);
  }
  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
    Double_t phiRad = TMath::DegToRad()*phiArrayDeg[ant];
    cosPhiArrayLookup[ant] = TMath::Cos(phiRad);
    sinPhiArrayLookup[ant] = TMath::Sin(phiRad);
  }

  for(Int_t phiSector = 0; phiSector<NUM_PHI; phiSector++){
    for(Int_t phiInd = 0; phiInd < NUM_BINS_PHI; phiInd++){
      Int_t phiBin = phiSector*NUM_BINS_PHI + phiInd;
      Double_t phiDeg = -0.5*PHI_RANGE + phiBin*Double_t(PHI_RANGE)/NUM_BINS_PHI;
      Double_t phiWave = TMath::DegToRad()*phiDeg;
      for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA; thetaBin++){
	Double_t thetaDeg = THETA_RANGE*((Double_t)thetaBin/NUM_BINS_THETA - 0.5);
	Double_t thetaWave = TMath::DegToRad()*thetaDeg;
	for(Int_t ant1=0; ant1<NUM_SEAVEYS; ant1++){
	  for(UInt_t ant2Ind=0; ant2Ind<ant2s[ant1].size(); ant2Ind++){
	    Int_t ant2 = ant2s[ant1].at(ant2Ind);
	    Int_t comboInd = comboIndices[ant1][ant2];
	    Int_t offset = getDeltaTExpected(ant1, ant2, phiWave, thetaWave);
	    offset = offset < 0 ? offset + NUM_SAMPLES : offset;
	    deltaTs[comboInd][phiBin][thetaBin] = offset;
	  }
	}
      }
    }
  }  
}

short* CrossCorrelator::fillDeltaTLookupGPU(){
  /* This function is mostly copied from the ANITA-3 Prioritizerd */

  offsetIndGPU = (short*) malloc(sizeof(short)*NUM_PHI*LOCAL_COMBOS_PER_PHI_GPU*NUM_BINS_PHI*NUM_BINS_THETA);

  for(Int_t phiSectorInd=0; phiSectorInd < NUM_PHI; phiSectorInd++){
    for(Int_t localCombo=0; localCombo < LOCAL_COMBOS_PER_PHI_GPU; localCombo++){
      // phiSector - 1... self + right
      if(localCombo < 3){
	ant1Gpu[phiSectorInd][localCombo] = (phiSectorInd+15)%16 + localCombo*NUM_PHI;
	ant2Gpu[phiSectorInd][localCombo] = (ant1Gpu[phiSectorInd][localCombo] + NUM_PHI)%NUM_SEAVEYS;
      }
      else if(localCombo < GLOBAL_COMBOS_PER_PHI_GPU){
	ant1Gpu[phiSectorInd][localCombo] = (phiSectorInd+15)%NUM_PHI + NUM_PHI*((localCombo-3)/3);
	ant2Gpu[phiSectorInd][localCombo] = phiSectorInd + (localCombo%3)*NUM_PHI;
      }
      else{
	Int_t localInd = localCombo % GLOBAL_COMBOS_PER_PHI_GPU;
	ant1Gpu[phiSectorInd][localCombo] = ant1Gpu[phiSectorInd][localInd] + localCombo/GLOBAL_COMBOS_PER_PHI_GPU;
	if(ant1Gpu[phiSectorInd][localInd]/16 - ant1Gpu[phiSectorInd][localCombo]/16 != 0){
	  ant1Gpu[phiSectorInd][localCombo] -= 16;
	}
	ant2Gpu[phiSectorInd][localCombo] = ant2Gpu[phiSectorInd][localInd] + localCombo/GLOBAL_COMBOS_PER_PHI_GPU;
	if(ant2Gpu[phiSectorInd][localInd]/16 - ant2Gpu[phiSectorInd][localCombo]/16 != 0){
	  ant2Gpu[phiSectorInd][localCombo] -= 16;
	}
      }
    }
  }

  for(Int_t phiSector = 0; phiSector < NUM_PHI; phiSector++){
    for(Int_t combo=0; combo < LOCAL_COMBOS_PER_PHI_GPU; combo++){
      Int_t a1 = ant1Gpu[phiSector][combo];
      Int_t a2 = ant2Gpu[phiSector][combo];
      for(Int_t phiInd = 0; phiInd < NUM_BINS_PHI; phiInd++){
	Double_t phiDeg = phiArrayDeg[phiSector] + PHI_RANGE*((Double_t)phiInd/NUM_BINS_PHI - 0.5);
	Double_t phiWave = phiDeg*TMath::DegToRad();
	for(Int_t thetaInd = 0; thetaInd < NUM_BINS_THETA; thetaInd++){
	  Double_t thetaDeg = THETA_RANGE*((Double_t)thetaInd/NUM_BINS_THETA - 0.5);
	  Double_t thetaWave = thetaDeg*TMath::DegToRad();
	  short offset = getDeltaTExpected(a1, a2, phiWave, thetaWave);

	  if(offset < 0){
	    offset += NUM_SAMPLES;
	  }
	  UInt_t offsetIndInd = (phiSector*LOCAL_COMBOS_PER_PHI_GPU*NUM_BINS_PHI*NUM_BINS_THETA
			       + combo*NUM_BINS_PHI*NUM_BINS_THETA
			       + phiInd*NUM_BINS_THETA
			       + thetaInd);
	  offsetIndGPU[offsetIndInd] = offset;
	}
      }
    }
  }
  return offsetIndGPU;

}





























/************************************************************************************************************
Image generation functions.
************************************************************************************************************/



TH2D* CrossCorrelator::makeImage(AnitaPol::AnitaPol_t pol){
  Double_t imagePeak, peakPhiDeg, peakThetaDeg;
  return makeImage(pol, imagePeak, peakPhiDeg, peakThetaDeg);
}

TH2D* CrossCorrelator::makeImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg, Double_t& peakThetaDeg){

  imagePeak = -100;

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
    for(Int_t phiInd = 0; phiInd < NUM_BINS_PHI; phiInd++){
      Int_t phiBin = phiSector*NUM_BINS_PHI + phiInd;
      for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA; thetaBin++){
        Double_t correlations = 0;
        Int_t contributors = 0;
	for(Int_t ant1=phiSector; ant1<NUM_SEAVEYS; ant1+=NUM_PHI){
	  for(UInt_t ant2Ind=0; ant2Ind<ant2s[ant1].size(); ant2Ind++){
	    Int_t ant2 = ant2s[ant1].at(ant2Ind);
	    Int_t comboInd = comboIndices[ant1][ant2];
	    Int_t offset = deltaTs[comboInd][phiBin][thetaBin];
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
			  NUM_BINS_PHI*NUM_PHI, phiArrayDeg[0]-PHI_RANGE/2, phiArrayDeg[15]+PHI_RANGE/2,
			  NUM_BINS_THETA, -THETA_RANGE/2, THETA_RANGE/2);
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

	    /* 
	       Not obvious, but since we correlate ant1 with ant2 but not ant2 with ant1 
	       we need to correct the sign of deltaT in the case that ant1 < ant2.
	    */
	    if(ant1 < ant2){
	      offset *= -1;
	    }
	    offset = offset < 0 ? offset + NUM_SAMPLES : offset;
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


TH2D* CrossCorrelator::makeImageGPU(AnitaPol::AnitaPol_t pol){
  UInt_t l3Trigger = 0xffff; // is equal to pow(2, 16)-1 i.e. 16 bits set.
  return makeImageGPU(pol, l3Trigger);
}


TH2D* CrossCorrelator::makeImageGPU(AnitaPol::AnitaPol_t pol, UInt_t l3Trigger){

  assert(pol == AnitaPol::kVertical || pol == AnitaPol::kHorizontal);
  TString name = pol == AnitaPol::kVertical ? "hImageV" : "hImageH";

  TH2D* hImage = new TH2D(name, name,
			  NUM_BINS_PHI*NUM_PHI, -PHI_RANGE/2, PHI_RANGE*NUM_PHI-PHI_RANGE/2,
			  NUM_BINS_THETA, -THETA_RANGE/2, THETA_RANGE/2);

  for(Int_t phiSectorInd = 0; phiSectorInd < NUM_PHI; phiSectorInd++){
    UInt_t doPhiSector = ((l3Trigger >> phiSectorInd) & 1);
    //    std::cout << phiSectorInd << " " << doPhiSector << std::endl;
    if(doPhiSector){
      for(Int_t phiInd = 0; phiInd < NUM_BINS_PHI; phiInd++){
	for(Int_t thetaInd = 0; thetaInd < NUM_BINS_THETA; thetaInd++){
	  Double_t correlations = 0;
	  Int_t contributors = 0;
	  for(Int_t localCombo=0; localCombo<LOCAL_COMBOS_PER_PHI_GPU; localCombo++){

	    Int_t comboLocal = localCombo;
	    Int_t phiSector = (phiSectorInd + 15 + comboLocal/GLOBAL_COMBOS_PER_PHI_GPU)%NUM_PHI;
	    Int_t globalCombo = comboLocal % GLOBAL_COMBOS_PER_PHI_GPU;

	    short offset = offsetIndGPU[phiSectorInd*LOCAL_COMBOS_PER_PHI_GPU*NUM_BINS_PHI*NUM_BINS_THETA
					+ localCombo*NUM_BINS_PHI*NUM_BINS_THETA
					+ phiInd*NUM_BINS_THETA
					+ thetaInd];

	    Double_t correlation = correlationsGPU[pol][phiSector][globalCombo][offset];

	    // sum correlations
	    correlations += correlation;
	    contributors++;

	  }
	  if(contributors>0){
	    correlations /= contributors;
	  }
	  hImage->SetBinContent(phiSectorInd*NUM_BINS_PHI + phiInd + 1, thetaInd + 1, correlations);
	}
      }
    }
  }

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
  /* Generates a set of delta function like waveforms, correlates them */

  /* Delete any old waveforms (at start rather than end to leave in memory to be examined if need be)*/
  deleteAllWaveforms();

  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
    for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
      Double_t newVs[NUM_SAMPLES] = {0};
      Double_t newTimes[NUM_SAMPLES] = {0};
      Int_t dt = NUM_SAMPLES/2 + getDeltaTExpectedSpherical(0, ant, 
						  phiDegSource*TMath::DegToRad(), 
						  thetaDegSource*TMath::DegToRad(), 
						  rSource);
      for(int samp=0; samp<NUM_SAMPLES; samp++){
	newTimes[samp] = correlationDeltaT*samp;
	if(samp==dt){
	  newVs[samp] = NUM_SAMPLES;
	}
      }
      grsInterp[pol][ant] = new TGraph(NUM_SAMPLES, newTimes, newVs);    
      RootTools::normalize(grsInterp[pol][ant]);
    }
  }

  doAllCrossCorrelations();
}


TGraph* CrossCorrelator::getCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2){
  /* Primarily for debugging, put cross correlations in a TGraph */

  Int_t comboInd = comboIndices[ant1][ant2];
  if(comboInd < 0 ){
    return NULL;
  }
  Double_t offsets[NUM_SAMPLES];
  Double_t corrs[NUM_SAMPLES];

  for(Int_t i=0; i<NUM_SAMPLES; i++){
    Int_t offset = (i - NUM_SAMPLES/2);
    offsets[i] = offset*correlationDeltaT;
    Int_t j = offset < 0 ? offset + NUM_SAMPLES : offset;
    corrs[i] = crossCorrelations[pol][comboInd][j];

  }


  TGraph* gr = new TGraph(NUM_SAMPLES,  offsets, corrs);
  gr->SetName(TString::Format("grCorr_%d_%d", ant1, ant2));
  gr->SetTitle(TString::Format("Cross Correlation ant1 = %d, ant2 = %d", ant1, ant2));

  return gr;
}