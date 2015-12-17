/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             A Class to manage power spectra, probably needs to do things like rolling averages.
*************************************************************************************************************** */

#ifndef AVERAGEPOWERSPECTRUM_H
#define AVERAGEPOWERSPECTRUM_H

#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TNamed.h"

#include "FancyFFTs.h"


#define NUM_AMPLITUDE_BINS 64
#define INITIAL_MAX_AMPLITUDE 4
#define INITIAL_MIN_AMPLITUDE 0

class AveragePowerSpectrum : public TNamed {

public:

  enum mode_t{
    kRolling, /// < This mode keeps a pointer to each event's power spectrum
    kSummed /// < This mode sums each event into a total
  };
  
  AveragePowerSpectrum();
  AveragePowerSpectrum(TString name, TString title, Double_t dt, Int_t n, mode_t powSpecMode=kRolling);
  ~AveragePowerSpectrum();

  size_t add(TGraph* gr);


  
  TGraph* makeAvePowSpecTGraph(); ///< Creates and returns a TGraph of the average power spectrum.
  TGraph* makeAvePowSpecTGraph_dB(); ///< Creates and returns a TGraph of the average power spectrum with dB scale and bins in MHz bins.


  
  void deleteRayleighDistributions(); ///< Deletes the Rayleigh Histograms

  void fitRayleighHistogram(Int_t freqInd);
  void fitAllRayleighHistograms();
  

  TH1D* getRayleighHistogram(Int_t freqInd);
  TH1D* getRayleighHistogramFromFrequencyMHz(Double_t freqMHz);
  TF1* getRayleighHistogramFit(Int_t freqInd);  

  TH2D* makeRayleigh2DHistogram();
  

  static TF1* makeRayleighFunction(TString name, Double_t xMin, Double_t xMax);
  static TString getRayleighFunctionText();


  Double_t numSamples;
  Double_t deltaT;
  Double_t deltaF;  
  Int_t numFreqs;
  Double_t* freqArray;//[numFreqs] ///< Array of frequency values  
  std::vector<Double_t> summedPowSpec; ///< How we store the average power spectrum in the summed case
  std::vector<TH1D*> hRayleighs; ///< Histograms for Rayleigh distributions
  std::vector<TF1*> hRayleighFits; ///< Fits to Rayleigh distributions histograms
  std::vector<std::vector<Double_t> > psdOutliers;///< Storage for outlier values, if get too many then add them to rayeligh histogram
  UInt_t maxNumOutliers; ///< Critical number of outliers before adding them, set in constructor.
  Int_t count;
  AveragePowerSpectrum::mode_t mode; ///< How we store the average power spectrum in the rolling case
  std::vector<std::vector<Double_t> > storedPowSpecs; ///< vector

  
  ClassDef(AveragePowerSpectrum, 1);
};

  



#endif









