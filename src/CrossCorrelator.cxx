#include "CrossCorrelator.h"

CrossCorrelator::CrossCorrelator(){
  initializeVariables();
}


CrossCorrelator::~CrossCorrelator(){
}



void CrossCorrelator::initializeVariables(){

  kDeltaPhiSect = 2;
  multiplyTopRingByMinusOne = 0;


  // Initialize with NULL otherwise very bad things will happen with gcc
  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      interpRMS[pol][ant] = 0;
      interpRMS2[pol][ant] = 0;      
    }
  }

  kOnlyThisCombo = -1;
  numSamples = PAD_FACTOR*NUM_SAMPLES; // Factor of two for padding. Turns circular xcor into linear xcor.
  numSamplesUpsampled = numSamples*UPSAMPLE_FACTOR; // For upsampling

  nominalSamplingDeltaT = NOMINAL_SAMPLING_DELTAT;
  correlationDeltaT = nominalSamplingDeltaT/UPSAMPLE_FACTOR;

  do5PhiSectorCombinatorics();
}


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
    std::cerr << "\tSeriously bad things are probably about to happen." << std::endl;
    std::cerr << "\tCheck the combinatorics!" << std::endl;
  }

  fillCombosToUse();
}

// void CrossCorrelator::printInfo(){
//   std::cerr << __PRETTY_FUNCTION__ << std::endl;
//   std::cerr << "\tupsample factor = " << UPSAMPLE_FACTOR << std::endl;
//   // std::cerr << "\tBin size theta (deg) = " << Double_t(THETA_RANGE)/NUM_BINS_THETA << std::endl;
//   std::cerr << "\tBin size theta (deg) = " << Double_t(MAX_THETA - MIN_THETA)/NUM_BINS_THETA << std::endl;
//   std::cerr << "\tBin size phi (deg) = " << Double_t(PHI_RANGE)/NUM_BINS_PHI << std::endl;
//   std::cerr << "\tdeltaTs array size = "
// 	    << sizeof(Double_t)*NUM_POL*numCombos*NUM_PHI*NUM_BINS_PHI*NUM_BINS_THETA
// 	    << " bytes" << std::endl;
// }





TGraph* CrossCorrelator::interpolateWithStartTimeAndZeroMean(TGraph* grIn, Double_t startTime,
							     Double_t dt, Int_t nSamp){

  // I want to interpolate the waveform such that the mean of the post interpolation waveform is zero...

  std::vector<Double_t> newTimes = std::vector<Double_t>(nSamp, 0); // for the interp volts
  std::vector<Double_t> newVolts = std::vector<Double_t>(nSamp, 0); // for the interp times
  std::vector<Int_t> isPadding = std::vector<Int_t>(nSamp, 1); // whether or not the value is padding
  Double_t thisStartTime = grIn->GetX()[0];
  Double_t lastTime = grIn->GetX()[grIn->GetN()-1];

  // // Quantizes the start and end times so data points lie at integer multiples of nominal sampling
  // // perhaps not really necessary if all waveforms have the same start time...
  // startTime = dt*TMath::Nint(startTime/dt + 0.5);
  // lastTime = dt*TMath::Nint(lastTime/dt - 0.5);

   //ROOT interpolator object constructor takes std::vector objects
  std::vector<Double_t> tVec(grIn->GetX(), grIn->GetX() + grIn->GetN());
  std::vector<Double_t> vVec(grIn->GetY(), grIn->GetY() + grIn->GetN());

  // This is ROOT's interpolator object

  const int minGlsInterpPoints = 5;
  if(tVec.size() < minGlsInterpPoints || vVec.size() < minGlsInterpPoints){
    std::cerr << "Warning! Very short waveform! Make sure you have a data quality cut on number of samples" << std::endl;
    // not really random, but a hack to make sure GSL doesn't crash
    // this is why you should do data quality checks BEFORE reconstruction...
    const int randomPadLength = 20;
    for(int i=0; i < randomPadLength; i++){
      vVec.push_back(0);
      tVec.push_back(nominalSamplingDeltaT + (tVec.at(tVec.size()-1)));
    }
  }

  ROOT::Math::Interpolator chanInterp(tVec,vVec,ROOT::Math::Interpolation::kAKIMA);

  // Put new data into arrays
  Double_t sum = 0;
  Double_t time = startTime;
  Int_t meanCount = 0;
  for(Int_t samp = 0; samp < nSamp; samp++){
    newTimes.at(samp) = time;
    if(time >= thisStartTime && time <= lastTime){
      Double_t V = chanInterp.Eval(time);
      newVolts.at(samp) = V;
      sum += V;
      isPadding.at(samp) = 0;
      meanCount++;
    }
    else{
      newVolts.at(samp) = 0;
    }
    time += dt;
  }

  // so now we have the new arrays, time to normalize these bad boys
  Double_t mean = sum/meanCount; // the mean of the non-padded values
  for(Int_t samp=0; samp < nSamp; samp++){
    if(isPadding.at(samp)==0){
      newVolts.at(samp) -= mean;
    }
  }

  return new TGraph(nSamp, &newTimes[0], &newVolts[0]);

}






Double_t CrossCorrelator::getCrossCorrelation(AnitaPol::AnitaPol_t pol, Int_t combo, Double_t deltaT){

  Int_t ant1 = comboToAnt1s[combo];
  Int_t ant2 = comboToAnt2s[combo];
  // std::cerr << ant1 << "\t" << ant2 << "\t" << deltaT << "\t";

  deltaT += startTimes[pol][ant1];
  deltaT -= startTimes[pol][ant2];  
  // deltaT += startTimes[pol][ant2];
  // deltaT -= startTimes[pol][ant1];  

  // std::cerr << deltaT << std::endl;
  Int_t offsetLow = floor(deltaT/nominalSamplingDeltaT);
  Double_t dt1 = offsetLow*nominalSamplingDeltaT;
  Double_t interpPrefactor = (deltaT - dt1)/nominalSamplingDeltaT;
  offsetLow += numSamples/2; // account for time ordering correlations in internal memory (i.e. dt=0 is half way thrugh the array)

  Double_t c1 = crossCorrelations[pol][combo][offsetLow];
  Double_t c2 = crossCorrelations[pol][combo][offsetLow+1];
  Double_t cInterp = interpPrefactor*(c2 - c1) + c1;
  
  return cInterp;
  // return c1;  
  
}



void CrossCorrelator::getNormalizedInterpolatedTGraphs(FilteredAnitaEvent* fEv,
						       AnitaPol::AnitaPol_t pol){

  for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){
    const AnalysisWaveform* wf = fEv->getFilteredGraph(ant, pol);    
    const TGraphAligned* gr = wf->even();
    startTimes[pol][ant] = gr->GetX()[0];
    int n = gr->GetN();

    Double_t invertFactor = multiplyTopRingByMinusOne > 0 && ant < NUM_PHI ? -1 : 1;

    Double_t sumOfV = 0;
    for(int samp=0; samp < n; samp++){
      Double_t V = invertFactor*gr->GetY()[samp];
      fVolts[pol][ant][samp] = V;      
      sumOfV += V;
    }

    Double_t meanV = sumOfV/n;    

    Double_t sumOfVSquared = 0;    
    for(int samp=0; samp < n; samp++){
      fVolts[pol][ant][samp] -= meanV;
      sumOfVSquared += fVolts[pol][ant][samp]*fVolts[pol][ant][samp];
    }

    for(int samp=n; samp < numSamples; samp++){
      fVolts[pol][ant][samp] = 0;
    }
    
    interpRMS[pol][ant] = TMath::Sqrt(sumOfVSquared/numSamples);
    interpRMS2[pol][ant] = TMath::Sqrt(sumOfVSquared/numSamplesUpsampled);
  }
}




void CrossCorrelator::doFFTs(AnitaPol::AnitaPol_t pol){

  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
    FancyFFTs::doFFT(numSamples, fVolts[pol][ant], ffts[pol][ant]);
    // renormalizeFourierDomain(pol, ant); 
  }
}









void CrossCorrelator::correlateEvent(FilteredAnitaEvent* fEv){

  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    eventNumber[pol] = fEv->eventNumber;    
    correlateEvent(fEv, (AnitaPol::AnitaPol_t)pol);
  }
}


void CrossCorrelator::correlateEvent(FilteredAnitaEvent* fEv, AnitaPol::AnitaPol_t pol){

  // Read TGraphs from events into memory (also deletes old TGraphs)
  // getFftsAndStartTimes(fEv, pol);
  getNormalizedInterpolatedTGraphs(fEv, pol);
  doFFTs(pol);  
  
  // std::cout << "here" << std::endl;
  doAllCrossCorrelations(pol);

}

void CrossCorrelator::doAllCrossCorrelations(AnitaPol::AnitaPol_t pol){

  // Set variable for use in threads

  Double_t* ccInternalArray = FancyFFTs::getRealArray(std::pair<int, int>(numSamples, 0));

  for(int combo=0; combo<NUM_COMBOS; combo++){
    Int_t ant1 = comboToAnt1s.at(combo);
    Int_t ant2 = comboToAnt2s.at(combo);
    FancyFFTs::crossCorrelatePadded(numSamples,
				    1,
				    ffts[pol][ant2],
				    ffts[pol][ant1],
				    crossCorrelations[pol][combo],
				    false,
				    0,
				    false);
    const int offset = numSamples/2;

    // copies first half of original array (times >= 0) into second half of internal storage
    for(Int_t samp=0; samp < numSamples/2; samp++){
      crossCorrelations[pol][combo][samp+offset] = ccInternalArray[samp];
      // ptr->crossCorrelations[pol][combo][samp+offset] = stash[samp];
    }
    // copies second half of original array (times < 0) into first half of internal storage
    for(Int_t samp=numSamples/2; samp < numSamples; samp++){
      crossCorrelations[pol][combo][samp-offset] = ccInternalArray[samp];
      // ptr->crossCorrelations[pol][combo][samp-offset] = stash[samp];
    }
  }
  
}


void CrossCorrelator::doUpsampledCrossCorrelations(AnitaPol::AnitaPol_t pol, Int_t phiSector){


  const std::vector<Int_t>* combosToUse = &combosToUseGlobal[phiSector];
  Int_t numCombosToUpsample = combosToUse->size();

  Double_t* ccInternalArray = FancyFFTs::getRealArray(std::pair<int, int>(numSamplesUpsampled, 0));

  for(int comboInd=0; comboInd<numCombosToUpsample; comboInd++){
    Int_t combo = combosToUse->at(comboInd);

    Int_t ant1 = comboToAnt1s.at(combo);
    Int_t ant2 = comboToAnt2s.at(combo);

    FancyFFTs::crossCorrelatePadded(numSamples,
				    UPSAMPLE_FACTOR,
				    ffts[pol][ant2],
				    ffts[pol][ant1],
				    crossCorrelationsUpsampled[pol][combo],
				    false,
				    0, false);

    const int offset = numSamplesUpsampled/2;

    // copies first half of original array (times >= 0) into second half of internal storage
    for(Int_t samp=0; samp < numSamplesUpsampled/2; samp++){
      crossCorrelationsUpsampled[pol][combo][samp+offset] = ccInternalArray[samp];
    }
    // copies second half of original array (times < 0) into first half of internal storage
    for(Int_t samp=numSamplesUpsampled/2; samp < numSamplesUpsampled; samp++){
      crossCorrelationsUpsampled[pol][combo][samp-offset] = ccInternalArray[samp];
    }
  }
}


void CrossCorrelator::getMaxCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t combo,
						 Double_t& time, Double_t& value){

  Int_t maxInd = RootTools::getIndexOfMaximum(numSamples, crossCorrelations[pol][combo]);
  Int_t offset = maxInd - numSamples/2;
  value = crossCorrelations[pol][combo][maxInd];
  time = offset*nominalSamplingDeltaT;
}


void CrossCorrelator::getMaxCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
						 Double_t& time, Double_t& value){
  Int_t combo = comboIndices[ant1][ant2];
  getMaxCorrelationTimeValue(pol, combo, time, value);
}


void CrossCorrelator::getMaxUpsampledCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t combo,
							  Double_t& time, Double_t& value){

  Int_t maxInd = RootTools::getIndexOfMaximum(numSamplesUpsampled, crossCorrelationsUpsampled[pol][combo]);
  Int_t offset = maxInd - numSamplesUpsampled/2;
  value = crossCorrelationsUpsampled[pol][combo][maxInd];
  time = offset*correlationDeltaT;
}



void CrossCorrelator::getMaxUpsampledCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
							  Double_t& time, Double_t& value){
  Int_t combo = comboIndices[ant1][ant2];
  getMaxUpsampledCorrelationTimeValue(pol, combo, time, value);
}





Bool_t CrossCorrelator::useCombo(Int_t ant1, Int_t ant2, Int_t phiSector, Int_t deltaPhiSect){

  // I want to be able to choose whether or not require one of the antennas to be the phi-sector
  // of interest or just to have both in range of deltaPhiSect.
  // I'm going to use the sign of deltaPhiSect to do this.
  // If deltaPhiSect < 0 this implies one of the antenna pairs must be in an l3Triggered phi-sector

  Int_t absDeltaPhiSect = TMath::Abs(deltaPhiSect);

  // Require that the difference in phi-sectors be <= absDeltaPhiSect
  Int_t phiSectorOfAnt1 = ant1%NUM_PHI;
  Bool_t ant1InRange = TMath::Abs(phiSector - (phiSectorOfAnt1))<=absDeltaPhiSect;

  // Takes account of wrapping around payload (e.g. antennas in phi-sectors 1 and 16 are neighbours)
  ant1InRange = ant1InRange || TMath::Abs(phiSector - (phiSectorOfAnt1))>=(NUM_PHI-absDeltaPhiSect);

  Int_t phiSectorOfAnt2 = ant2%NUM_PHI;
  Bool_t ant2InRange = TMath::Abs(phiSector - (phiSectorOfAnt2))<=absDeltaPhiSect;
  ant2InRange = ant2InRange || TMath::Abs(phiSector - (phiSectorOfAnt2))>=(NUM_PHI-absDeltaPhiSect);

  // See rather rant above.
  Bool_t extraCondition = true;
  if(deltaPhiSect < 0){
    extraCondition = (phiSectorOfAnt1==phiSector || phiSectorOfAnt2==phiSector);
  }

  return (ant1InRange && ant2InRange && extraCondition);
}


void CrossCorrelator::fillCombosToUse(){

  for(Int_t phiSector = 0; phiSector<NUM_PHI; phiSector++){
    if(combosToUseGlobal[phiSector].size() == 0){
      for(Int_t combo=0; combo<numCombos; combo++){
	Int_t ant1 = comboToAnt1s.at(combo);
	Int_t ant2 = comboToAnt2s.at(combo);
	if(useCombo(ant1, ant2, phiSector, kDeltaPhiSect)){
	  combosToUseGlobal[phiSector].push_back(combo);
	}
      }
    }
  }
}




Double_t CrossCorrelator::getInterpolatedUpsampledCorrelationValue(AnitaPol::AnitaPol_t pol,
								   Int_t combo, Double_t deltaT){

  Int_t offsetLow = floor(deltaT/correlationDeltaT);

  Double_t dt1 = offsetLow*correlationDeltaT;
  // Double_t dt2 = offsetHigh*ptr->correlationDeltaT;

  const Int_t offset = numSamplesUpsampled/2;
  offsetLow += offset;
  Int_t offsetHigh = offsetLow+1;

  // offsetLow = offsetLow < 0 ? offsetLow + numSamplesUpsampled : offsetLow;
  // offsetHigh = offsetHigh < 0 ? offsetHigh + numSamplesUpsampled : offsetHigh;

  Double_t c1 = crossCorrelationsUpsampled[pol][combo][offsetLow];
  Double_t c2 = crossCorrelationsUpsampled[pol][combo][offsetHigh];

  Double_t cInterp = (deltaT - dt1)*(c2 - c1)/(correlationDeltaT) + c1;

  return cInterp;

}





TGraph* CrossCorrelator::getCrossCorrelationGraphWorker(Int_t numSamps, AnitaPol::AnitaPol_t pol,
							Int_t ant1, Int_t ant2){
  // Primarily for debugging, put cross correlations in a TGraph

  Int_t combo = comboIndices[ant1][ant2];
  if(combo < 0){
    return NULL;
  }

  Double_t graphDt = correlationDeltaT;
  Double_t* corrPtr = crossCorrelationsUpsampled[pol][combo];
  if(numSamps != numSamplesUpsampled){
    numSamps = numSamples;
    graphDt = nominalSamplingDeltaT;
    corrPtr = crossCorrelations[pol][combo];
  }

  // Could actually perform the calculation here... but I'll leave it NULL for now
  if(corrPtr==NULL){
    return NULL;
  }

  std::vector<Double_t> offsets = std::vector<Double_t>(numSamps, 0);
  std::vector<Double_t> corrs = std::vector<Double_t>(numSamps, 0);

  // put the correlations in, in time order in memory

  for(Int_t i=0; i<numSamps; i++){
    Int_t offset = (i - numSamps/2);
    offsets.at(i) = offset*graphDt;
    corrs.at(i) = corrPtr[i];
  }


  TGraph* gr = new TGraph(numSamps,  &offsets[0], &corrs[0]);
  if(numSamps != numSamplesUpsampled){
    gr->SetName(TString::Format("grCorr_%d_%d", ant1, ant2));
    gr->SetTitle(TString::Format("Cross Correlation ant1 = %d, ant2 = %d", ant1, ant2));
  }
  else{
    gr->SetName(TString::Format("grCorrUpsampled_%d_%d", ant1, ant2));
    gr->SetTitle(TString::Format("Upsampled Cross Correlation ant1 = %d, ant2 = %d", ant1, ant2));
  }

  return gr;
}

TGraph* CrossCorrelator::getCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2){
  // Primarily for debugging, put cross correlations in a TGraph
  return getCrossCorrelationGraphWorker(numSamples, pol, ant1, ant2);
}

TGraph* CrossCorrelator::getUpsampledCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2){
  return getCrossCorrelationGraphWorker(numSamplesUpsampled, pol, ant1, ant2);
}


