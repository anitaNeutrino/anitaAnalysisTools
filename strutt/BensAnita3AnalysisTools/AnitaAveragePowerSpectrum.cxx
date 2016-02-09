#include "AnitaAveragePowerSpectrum.h"


AnitaAveragePowerSpectrum::AnitaAveragePowerSpectrum(){

  fName = "AnitaAveragePowerSpectrum";
  fTitle = "Anita Average Power Spectrum";
  
  // Default constructor, just zero internals
  for(Int_t polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      avePowSpecs[polInd][ant] = NULL;
    }
  }
  
}


AnitaAveragePowerSpectrum::AnitaAveragePowerSpectrum(TString name, TString title){

  // Record initialization options for contained AveragePowerSpectra
  fName = name;
  fTitle = title;

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

      TString title = RootTools::getAntName(AnitaPol::AnitaPol_t(polInd), ant);

      avePowSpecs[polInd][ant] = new AveragePowerSpectrum(name, title);
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


TMultiGraph* AnitaAveragePowerSpectrum::drawSpectralSummary(AnitaPol::AnitaPol_t pol, AnitaRing::AnitaRing_t ring){
  TMultiGraph* mg = new TMultiGraph();

  TString grTitle = TString::Format("%s PSDs", fTitle.Data());

  switch (ring){
  case AnitaRing::kTopRing:
    grTitle += " Top ring";
    break;
  case AnitaRing::kMiddleRing:
    grTitle += " Middle ring";
    break;
  case AnitaRing::kBottomRing:
    grTitle += " Bottom ring";
    break;
  case AnitaRing::kNotARing:
    break;
  }

  grTitle += pol == AnitaPol::kHorizontal ? " HPol " : " VPol ";

  grTitle += ";Frequency (MHz); PSD (mV^{2}/MHz)";

  // std::cout << grTitle.Data() << std::endl;
  mg->SetTitle(grTitle);
  for(int phi=0; phi<NUM_PHI; phi++){
    Int_t ant = ring*NUM_PHI + phi;
    AveragePowerSpectrum* aps = get(pol, ant);
    TGraph* gr = aps->makeAvePowSpecTGraph_dB();
    gr->SetLineColor(gStyle->GetColorPalette(phi*Int_t(254./(NUM_PHI-1))));
    gr->SetFillColorAlpha(0, 0);
    mg->Add(gr);
  }

  return mg;
}
