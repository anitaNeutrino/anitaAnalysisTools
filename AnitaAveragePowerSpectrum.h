/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             A holder for the whole payload's worth of AveragePowerSpectra. 
	     Because TTrees have a particularly inelegant interface in ROOT.
*************************************************************************************************************** */

#ifndef ANITAAVERAGEPOWERSPECTRUM_H
#define ANITAAVERAGEPOWERSPECTRUM_H

#include "AveragePowerSpectrum.h"
#include "RootTools.h"
#include "AnitaConventions.h"
#include "TMultiGraph.h"
#include "TStyle.h"


/**
 * @class AnitaAveragePowerSpectrum
 * @brief Essentially a wrapper around a whole payload of AveragePowerSpectrum.
 * 
 * I used to write these into TTrees, but the way the data was packed made access so slow that now I don't.
 * However, this class is still used in some spectral analysis code in some places.
*/
class AnitaAveragePowerSpectrum : public TNamed {

public:

  // Copied from AveragePowerSpectrum.h
  AnitaAveragePowerSpectrum();
  AnitaAveragePowerSpectrum(TString name, TString title);
  ~AnitaAveragePowerSpectrum();
  
  AveragePowerSpectrum* get(AnitaPol::AnitaPol_t pol, Int_t ant);
  void add(AnitaPol::AnitaPol_t pol, Int_t ant, TGraph* gr);
  void reset();

  // Produce summary information
  TMultiGraph* drawSpectralSummary(AnitaPol::AnitaPol_t pol, AnitaRing::AnitaRing_t ring);
  
private:
  AveragePowerSpectrum* avePowSpecs[AnitaPol::kNotAPol][NUM_SEAVEYS]; //!< A payload's worth of AveragePowerSpectrum
  void initAllAvePowSpecs();
  void deleteAllAvePowSpecs();
  
  ClassDef(AnitaAveragePowerSpectrum, 9);
};

  



#endif









