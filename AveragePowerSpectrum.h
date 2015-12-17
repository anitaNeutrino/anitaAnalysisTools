/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             A Class to manage power spectra, probably needs to do things like rolling averages.
*************************************************************************************************************** */

#ifndef AVERAGEPOWERSPECTRUM_H
#define AVERAGEPOWERSPECTRUM_H

#include <TGraph.h>
#include "FancyFFTs.h"

class AveragePowerSpectrum {

public:
  AveragePowerSpectrum(Double_t dt, Int_t n);
  ~AveragePowerSpectrum();

  size_t add(TGraph* gr);
  TGraph* get(TString name, TString title);
  TGraph* getScaled(TString name, TString title);  


private:
  Double_t deltaT;
  Double_t numSamples;
  Int_t numFreqs;

  Double_t* freqArray; ///< Array of frequency values
  std::vector<Double_t*> storedPowSpecs; ///< vector 

};


#endif
