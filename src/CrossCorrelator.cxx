#include "CrossCorrelator.h"
#include "InterferometricMap.h"
#include "AnitaDataset.h"
#include "AcclaimCorrelationSummary.h"
#include "FilterStrategy.h"
#include "TFile.h" ///@todo remove after debugging complete...

Acclaim::CrossCorrelator::CrossCorrelator(){
  initializeVariables();
}


Acclaim::CrossCorrelator::~CrossCorrelator(){
}


double Acclaim::CrossCorrelator::correlationIndexToTime(bool upsampled, int corrIndex, AnitaPol::AnitaPol_t pol, int combo) const { 
  return correlationIndexToTime(upsampled, corrIndex, pol, comboToAnt1s.at(combo), comboToAnt2s.at(combo));
}

double Acclaim::CrossCorrelator::correlationIndexToTime(bool upsampled, int corrIndex, AnitaPol::AnitaPol_t pol, int ant1, int ant2) const {

  const int n     = upsampled ? numSamplesUpsampled : numSamples;
  const double dt = upsampled ? correlationDeltaT : nominalSamplingDeltaT;
  
  double deltaT = dt*(corrIndex - n/2);
  
  deltaT -= startTimes[pol][ant1];
  deltaT += startTimes[pol][ant2];
 
  return deltaT;
}


std::shared_ptr<const Acclaim::CorrelationSummary> Acclaim::CrossCorrelator::makeSummary(AnitaPol::AnitaPol_t pol, const FilteredAnitaEvent* event, double waisPhi, double waisTheta,  Adu5Pat* pat){
  
  int phiSector = InterferometricMap::getPhiSectorFromPhiRough(waisPhi);
  phiSector = phiSector % NUM_PHI; // enforce wrapping
  // std::cout << phiSector << std::endl;
  correlateEvent(event, pol);
  doUpsampledCrossCorrelations(pol, phiSector);

  auto summary = std::make_shared<Acclaim::CorrelationSummary>(pol, eventNumber[pol], waisPhi, waisTheta, phiSector, pat);

  for(int combo : combosToUseGlobal[phiSector]){
    int ant1 = comboToAnt1s[combo];
    int ant2 = comboToAnt2s[combo];

    const double* cc = crossCorrelationsUpsampled[pol][combo];
    int maxIndex = RootTools::getIndexOfMaximum(numSamplesUpsampled, cc);

    double maxCC = cc[maxIndex]/(numSamples*numSamples);
    
    double deltaT = correlationIndexToTime(true, maxIndex, pol, ant1, ant2);
    
    // std::cout << ant1 << "\t" << ant2 << "\t" << maxIndex << " of " << numSamplesUpsampled << "\t" << dt << "\t" << cc << std::endl;
    summary->add(ant1, ant2, deltaT, maxCC);
  }
  
  return summary;
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






Double_t Acclaim::CrossCorrelator::getCrossCorrelation(AnitaPol::AnitaPol_t pol, Int_t combo, Double_t deltaT) const{

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
                                                                AnitaPol::AnitaPol_t pol, bool raw){

  for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){
    const AnalysisWaveform* wf = raw ? fEv->getRawGraph(ant, pol) : fEv->getFilteredGraph(ant, pol);
    // const AnalysisWaveform* wf = raw ? fEv->getRawGraph(0, pol) : fEv->getFilteredGraph(0, pol);    
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
  FancyFFTs::normalizeFourierDomain(numSamples, ffts[pol][ant]);
}




void Acclaim::CrossCorrelator::correlateEvent(const FilteredAnitaEvent* fEv){

  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    correlateEvent(fEv, (AnitaPol::AnitaPol_t)pol);
  }
}


void Acclaim::CrossCorrelator::correlateEvent(const FilteredAnitaEvent* fEv, AnitaPol::AnitaPol_t pol){

  // Read TGraphs from events into memory (also deletes old TGraphs)
  eventNumber[pol] = fEv->getHeader()->eventNumber;
  
  getNormalizedInterpolatedTGraphs(fEv, pol);
  doFFTs(pol);  
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
    for(Int_t samp=0; samp < offset; samp++){
      crossCorrelations[pol][combo][samp+offset] = ccInternalArray[samp];
    }
    // copies second half of original array (times < 0) into first half of internal storage
    for(Int_t samp=0; samp < numSamples; samp++){
      crossCorrelations[pol][combo][samp-offset] = ccInternalArray[samp];
    }
  }
  
}


void Acclaim::CrossCorrelator::doUpsampledCrossCorrelations(AnitaPol::AnitaPol_t pol, Int_t phiSector){

  Double_t* ccInternalArray = FancyFFTs::getRealArray(std::pair<int, int>(numSamplesUpsampled, 0));

  for(int combo : combosToUseGlobal[phiSector]){
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
    for(Int_t samp=0; samp < offset; samp++){
      crossCorrelationsUpsampled[pol][combo][samp+offset] = ccInternalArray[samp];
    }
    // copies second half of original array (times < 0) into first half of internal storage
    for(Int_t samp=offset; samp < numSamplesUpsampled; samp++){
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

  bool usePair = (ant1InRange && ant2InRange && extraCondition);
  // if(usePair){
  //   std::cout << "Phi-sector " << phiSector << "\t" << ant1 << "\t" << ant2  << "\t(" << comboIndices[ant1][ant2] << ")" << std::endl;
  // }
  return usePair;
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
      // std::cout << "combosToUse[" << phiSector << "].size()==" << combosToUseGlobal[phiSector].size() << std::endl;
    }
  }
}







TGraph* Acclaim::CrossCorrelator::getCrossCorrelationGraphWorker(Int_t numSamps, AnitaPol::AnitaPol_t pol,
								 Int_t ant1, Int_t ant2) const {
  // Primarily for debugging, put cross correlations in a TGraph

  Int_t combo = comboIndices[ant1][ant2];
  if(combo < 0){
    return NULL;
  }

  Double_t graphDt = correlationDeltaT;
  const Double_t* corrPtr = crossCorrelationsUpsampled[pol][combo];
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

TGraph* Acclaim::CrossCorrelator::getCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2) const {
  // Primarily for debugging, put cross correlations in a TGraph
  return getCrossCorrelationGraphWorker(numSamples, pol, ant1, ant2);
}

TGraph* Acclaim::CrossCorrelator::getUpsampledCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2) const {
  return getCrossCorrelationGraphWorker(numSamplesUpsampled, pol, ant1, ant2);
}

