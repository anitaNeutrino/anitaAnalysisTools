#include "AnitaAveragePowerSpectrum.h"


AnitaAveragePowerSpectrum::AnitaAveragePowerSpectrum(){

  fName = "AnitaAveragePowerSpectrum";
  fTitle = "Anita Average Power Spectrum";
  mode = AveragePowerSpectrum::kRolling;
  numSamps = 256;
  deltaT = 1./2.6;
  
  // Default constructor, just zero internals
  for(Int_t polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      avePowSpecs[polInd][ant] = NULL;
    }
  }
  
}


AnitaAveragePowerSpectrum::AnitaAveragePowerSpectrum(TString name, TString title, Double_t dt, Int_t n,
						     AveragePowerSpectrum::mode_t powSpecMode){

  // Record initialization options for contained AveragePowerSpectra
  fName = name;
  fTitle = title;
  mode = powSpecMode;
  numSamps = n;
  deltaT = dt;

  // Initialize them.
  initAllAvePowSpecs();
}


AnitaAveragePowerSpectrum::~AnitaAveragePowerSpectrum(){
  deleteAllAvePowSpecs();
}




void AnitaAveragePowerSpectrum::deleteAllAvePowSpecs(){
  // Delete all non-null internals
  for(Int_t polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){

      if(avePowSpecs[polInd][ant]!=NULL){
	delete avePowSpecs[polInd][ant];
	avePowSpecs[polInd][ant] = NULL;
      }
    }
  }
}


void AnitaAveragePowerSpectrum::initAllAvePowSpecs(){
  // Initialize the power spectra from options recorded in constructor
  
  for(Int_t polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){

      TString name = fName + TString::Format("_%d_%d", polInd, ant);

      TString title = fTitle;
      title += polInd == 0 ? " HPOL " : " VPOL";
      title += TString::Format(" antenna %d", ant);
      
      avePowSpecs[polInd][ant] = new AveragePowerSpectrum(name, title, deltaT, numSamps, mode);
    }
  }
}


AveragePowerSpectrum* AnitaAveragePowerSpectrum::get(AnitaPol::AnitaPol_t pol, Int_t ant){
  // Getter function with array size checks
  AveragePowerSpectrum* avePowSpecPtr = NULL;
  if(ant >= 0 && ant < NUM_SEAVEYS){
    avePowSpecPtr = avePowSpecs[pol][ant];
  }
  return avePowSpecPtr;
}


void AnitaAveragePowerSpectrum::add(AnitaPol::AnitaPol_t pol, Int_t ant, TGraph* gr){
  // Wrapper for adding waveform to particular channels power spectrum
  AveragePowerSpectrum* avePowSpecPtr = get(pol, ant);
  avePowSpecPtr->add(gr);
}


void AnitaAveragePowerSpectrum::reset(){
  // Function to delete and reinitialize all contained AveragePowerSpectra
  deleteAllAvePowSpecs();
  initAllAvePowSpecs();
}
