#include "CrossCorrelator.h"

Acclaim::CrossCorrelator::CrossCorrelator(){
  initializeVariables();
}


Acclaim::CrossCorrelator::~CrossCorrelator(){
}



void Acclaim::CrossCorrelator::initializeVariables(){

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
  numSamples = PAD_FACTOR*NUM_SAMP; // Factor of two for padding. Turns circular xcor into linear xcor.
  numSamplesUpsampled = numSamples*UPSAMPLE_FACTOR; // For upsampling

  nominalSamplingDeltaT = NOMINAL_SAMPLING_DELTAT;
  correlationDeltaT = nominalSamplingDeltaT/UPSAMPLE_FACTOR;

  do5PhiSectorCombinatorics();
}


void Acclaim::CrossCorrelator::do5PhiSectorCombinatorics(){
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






Double_t Acclaim::CrossCorrelator::getCrossCorrelation(AnitaPol::AnitaPol_t pol, Int_t combo, Double_t deltaT){

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



void Acclaim::CrossCorrelator::getNormalizedInterpolatedTGraphs(const FilteredAnitaEvent* fEv,
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




void Acclaim::CrossCorrelator::doFFTs(AnitaPol::AnitaPol_t pol){

  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
    FancyFFTs::doFFT(numSamples, fVolts[pol][ant], ffts[pol][ant]);
    renormalizeFourierDomain(pol, ant); 
  }
}



void Acclaim::CrossCorrelator::renormalizeFourierDomain(AnitaPol::AnitaPol_t pol, Int_t ant){

  const int numFreqs = FancyFFTs::getNumFreqs(numSamples);
  Double_t sumOfVSquared = 0;
  for(int freqInd=0; freqInd < numFreqs; freqInd++){
    Double_t re = ffts[pol][ant][freqInd].real();
    Double_t im = ffts[pol][ant][freqInd].imag();
    Double_t factor = freqInd == numFreqs - 1 ? 0.5 : 1;
    sumOfVSquared += factor*(re*re + im*im);
  }

  Double_t timeDomainVSquared = sumOfVSquared/(numSamples*numSamples/2);
  Double_t norm = TMath::Sqrt(timeDomainVSquared);

  for(int freqInd=0; freqInd < numFreqs; freqInd++){
    Double_t re = ffts[pol][ant][freqInd].real();
    Double_t im = ffts[pol][ant][freqInd].imag();

    ffts[pol][ant][freqInd].real(re/norm);
    ffts[pol][ant][freqInd].imag(im/norm);
  }
}







void Acclaim::CrossCorrelator::correlateEvent(const FilteredAnitaEvent* fEv){

  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    correlateEvent(fEv, (AnitaPol::AnitaPol_t)pol);
  }
}


void Acclaim::CrossCorrelator::correlateEvent(const FilteredAnitaEvent* fEv, AnitaPol::AnitaPol_t pol){

  // Read TGraphs from events into memory (also deletes old TGraphs)
  // getFftsAndStartTimes(fEv, pol);
  eventNumber[pol] = fEv->eventNumber;
  
  getNormalizedInterpolatedTGraphs(fEv, pol);
  doFFTs(pol);  
  
  // std::cout << "here" << std::endl;
  doCrossCorrelations(pol);

}

void Acclaim::CrossCorrelator::doCrossCorrelations(AnitaPol::AnitaPol_t pol){

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


void Acclaim::CrossCorrelator::doUpsampledCrossCorrelations(AnitaPol::AnitaPol_t pol, Int_t phiSector){


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





Bool_t Acclaim::CrossCorrelator::useCombo(Int_t ant1, Int_t ant2, Int_t phiSector, Int_t deltaPhiSect){

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


void Acclaim::CrossCorrelator::fillCombosToUse(){

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




Double_t Acclaim::CrossCorrelator::getInterpolatedUpsampledCorrelationValue(AnitaPol::AnitaPol_t pol,
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





TGraph* Acclaim::CrossCorrelator::getCrossCorrelationGraphWorker(Int_t numSamps, AnitaPol::AnitaPol_t pol,
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

TGraph* Acclaim::CrossCorrelator::getCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2){
  // Primarily for debugging, put cross correlations in a TGraph
  return getCrossCorrelationGraphWorker(numSamples, pol, ant1, ant2);
}

TGraph* Acclaim::CrossCorrelator::getUpsampledCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2){
  return getCrossCorrelationGraphWorker(numSamplesUpsampled, pol, ant1, ant2);
}


