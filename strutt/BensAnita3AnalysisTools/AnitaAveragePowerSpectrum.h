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
  AveragePowerSpectrum* avePowSpecs[AnitaPol::kNotAPol][NUM_SEAVEYS];  
  void initAllAvePowSpecs();
  void deleteAllAvePowSpecs();
  
  ClassDef(AnitaAveragePowerSpectrum, 8);
};

  



#endif









