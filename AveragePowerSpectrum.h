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
#include "TH1D.h"
#include "TF1.h"

class AveragePowerSpectrum {

public:

  enum mode_t{
    kRolling, /// < This mode keeps a pointer to each event's power spectrum
    kSummed /// < This mode sums each event into a total
  };
  
  AveragePowerSpectrum(TString histBaseNameIn, Double_t dt, Int_t n, mode_t powSpecMode=kRolling);
  ~AveragePowerSpectrum();

  size_t add(TGraph* gr);
  TGraph* get(TString name, TString title);
  TGraph* getScaled(TString name, TString title);
  void emptyRolling(); ///< Deletes all the power spectra and empties the vector
  void emptyRayleighs(); ///< Writes all the Rayleigh histograms and deletes them
  static TF1* makeRayleighFunction(TString name, Double_t xMin, Double_t xMax);

private:
  Double_t deltaT;
  Double_t numSamples;
  Int_t numFreqs;
  Int_t count;
  TString histBaseName;

  Double_t* freqArray; ///< Array of frequency values
  std::vector<Double_t*> storedPowSpecs; ///< vector
  AveragePowerSpectrum::mode_t mode; ///< How we store the average power spectrum in the rolling case
  std::vector<Double_t> summedPowSpec; ///< How we store the average power spectrum in the summed case
  std::vector<TH1D*> hRayleighs; ///< Histograms for Rayleigh distributions

  std::vector<std::vector<Double_t> > psdOutliers;///< Storage for outlier values, if get too many then add them to rayeligh histogram
  UInt_t maxNumOutliers; ///< Critical number of outliers before adding them, set in constructor.
  
};




#endif









