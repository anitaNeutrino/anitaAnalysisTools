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
#include "AnitaConventions.h"


class AnitaAveragePowerSpectrum : public TNamed {

public:

  // Copied from AveragePowerSpectrum.h
  AnitaAveragePowerSpectrum();
  AnitaAveragePowerSpectrum(TString name, TString title, Double_t dt, Int_t n,
			    AveragePowerSpectrum::mode_t powSpecMode=AveragePowerSpectrum::kRolling);

  ~AnitaAveragePowerSpectrum();
  
  AveragePowerSpectrum* get(AnitaPol::AnitaPol_t pol, Int_t ant);
  void add(AnitaPol::AnitaPol_t pol, Int_t ant, TGraph* gr);
  void reset();

  
private:
  AveragePowerSpectrum* avePowSpecs[AnitaPol::kNotAPol][NUM_SEAVEYS];  
  void initAllAvePowSpecs();
  void deleteAllAvePowSpecs();
  
  Double_t deltaT;
  Int_t numSamps;
  AveragePowerSpectrum::mode_t mode;
  
  ClassDef(AnitaAveragePowerSpectrum, 2);
};

  



#endif









