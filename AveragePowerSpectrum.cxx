#include "AveragePowerSpectrum.h"

AveragePowerSpectrum::AveragePowerSpectrum(Double_t dt, Int_t nSamp){
  deltaT = dt;
  numSamples = nSamp;
  numFreqs = FancyFFTs::getNumFreqs(nSamp);
  freqArray = FancyFFTs::getFreqArray(nSamp, dt);
}

AveragePowerSpectrum::~AveragePowerSpectrum(){
  delete [] freqArray;
  this->empty();
}


size_t AveragePowerSpectrum::add(TGraph* gr){

  Double_t* ps = FancyFFTs::getPowerSpectrum(numSamples, gr->GetY(), deltaT, PowSpecNorm::kPowSpecDensity);
  storedPowSpecs.push_back(ps);
  
  return storedPowSpecs.size();;
}

void AveragePowerSpectrum::empty(){
  while(!storedPowSpecs.empty()){
    delete [] storedPowSpecs.back();
    storedPowSpecs.pop_back();    
  }
}


TGraph* AveragePowerSpectrum::get(TString name, TString title){

  std::vector<Double_t> avePowSpec(numFreqs);

  UInt_t numEvents = storedPowSpecs.size();
  for(UInt_t eventInd=0; eventInd < numEvents; eventInd++){
    for(Int_t freqInd=0; freqInd < numFreqs; freqInd++){
      avePowSpec.at(freqInd) += storedPowSpecs.at(eventInd)[freqInd];
    }
  }
  for(Int_t freqInd=0; freqInd < numFreqs; freqInd++){
    avePowSpec.at(freqInd)/=numEvents;
  }
  
  TGraph* gr = new TGraph(numFreqs, freqArray, &avePowSpec[0]);
  gr->SetName(name);
  gr->SetTitle(title);  
  return gr;
}


TGraph* AveragePowerSpectrum::getScaled(TString name, TString title){

  TGraph* gr = get(name, title);
  for(Int_t freqInd=0; freqInd < numFreqs; freqInd++){
    Double_t y = gr->GetY()[freqInd];
    gr->GetY()[freqInd] = 10*TMath::Log10(y);
    gr->GetX()[freqInd]*= 1e3;
  }
  return gr;
}

