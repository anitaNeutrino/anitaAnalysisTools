#include "AnitaAveragePowerSpectrum.h"




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Constructor
 *
 * Inserts default names for internal parameters and NULL pointers for avePowSpecs
 */
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




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Constructor
 *
 * @param name is the name of the AnitaAveragePowerSpectrum
 * @param title is the title of the AnitaAveragePowerSpectrum
 */
AnitaAveragePowerSpectrum::AnitaAveragePowerSpectrum(TString name, TString title){

  // Record initialization options for contained AveragePowerSpectra
  fName = name;
  fTitle = title;

  // Initialize them.
  initAllAvePowSpecs();
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Destructor
 */
AnitaAveragePowerSpectrum::~AnitaAveragePowerSpectrum(){
  deleteAllAvePowSpecs();
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Deletes all non-NULL interal pointers.
 */
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





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates empty AveragePowerSpectrums.
 *
 * Names each interal AveragePowerSpectrums with the fName plus polarization / antenna number suffix.
 */
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






//---------------------------------------------------------------------------------------------------------
/**
 * @brief Access internal AveragePowerSpectrum.
 *
 * @param pol is polarization of interest.
 * @param ant is the antenna of interest.
 */
AveragePowerSpectrum* AnitaAveragePowerSpectrum::get(AnitaPol::AnitaPol_t pol, Int_t ant){
  // Getter function with array size checks
  AveragePowerSpectrum* avePowSpecPtr = NULL;
  if(ant >= 0 && ant < NUM_SEAVEYS){
    avePowSpecPtr = avePowSpecs[pol][ant];
  }
  return avePowSpecPtr;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Add a waveform TGraph to an internal AveragePowerSpectrum.
 *
 * @param pol is polarization of interest.
 * @param ant is the antenna of interest.
 * @param gr is the waveform TGraph.
 */
void AnitaAveragePowerSpectrum::add(AnitaPol::AnitaPol_t pol, Int_t ant, TGraph* gr){
  // Wrapper for adding waveform to particular channels power spectrum
  AveragePowerSpectrum* avePowSpecPtr = get(pol, ant);
  avePowSpecPtr->add(gr);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Deletes and reinitializes all internal AveragePowerSpectrums.
 */
void AnitaAveragePowerSpectrum::reset(){
  // Function to delete and reinitialize all contained AveragePowerSpectra
  deleteAllAvePowSpecs();
  initAllAvePowSpecs();
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Get a ring of TGraphs of a single polarization as a TMultiGraph
 * 
 * @param pol is polarization of interest.
 * @param ring is the ring of interest.
 * @return TMultiGraph containing NUM_PHI TGraphs.
 *
 * NUM_PHI (16) seems like enough to examine at once.
 */
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
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
    std::cerr << "Your  version of ROOT does not support all the features of eventReaderRoot\n";
    gr->SetFillColorAlpha(0, 0);
#endif
    mg->Add(gr);
  }

  return mg;
}
