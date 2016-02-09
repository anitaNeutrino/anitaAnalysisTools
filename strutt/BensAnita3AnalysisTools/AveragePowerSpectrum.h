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
#include "RootTools.h"
#include "CrossCorrelator.h"

#define NUM_AMPLITUDE_BINS 64
#define INITIAL_MAX_AMPLITUDE 4
#define INITIAL_MIN_AMPLITUDE 0

#define MAX_NUM_OUTLIERS 10

// NUM_SAMPLES is included from CrossCorrelator.h
#define NUM_FREQS ((NUM_SAMPLES/2)+1)

class AveragePowerSpectrum : public TNamed {

public:
  
  AveragePowerSpectrum();
  AveragePowerSpectrum(TString name, TString title);
  ~AveragePowerSpectrum();

  size_t add(TGraph* gr);


  
  TGraph* makeAvePowSpecTGraph(); ///< Creates and returns a TGraph of the average power spectrum.
  TGraph* makeAvePowSpecTGraph_dB(); ///< Creates and returns a TGraph of the average power spectrum with dB scale and bins in MHz bins.


  
  void deleteRayleighDistributions(); ///< Deletes the Rayleigh Histograms
  void rebinAllRayleighHistograms(Int_t rebinFactor); ///< Loops over them all and rebins

  
  void fitRayleighHistogramOverRange(Int_t freqInd, Double_t xLowVal, Double_t xHighVal,
				     Double_t* rAmplitudes,
				     Double_t* rChiSquares,
				     Int_t* rNdf,
				     Double_t* rChiSquaresFullRange,
				     Int_t* rNdfFullRange);


  void fitRayleighHistogram(Int_t freqInd);
  void fitRayleighHistogramRisingEdge(Int_t freqInd);
  void fitRayleighHistogramRisingEdgeAndHalfFallingEdge(Int_t freqInd);

  void fitAllRayleighHistograms();
  void fitAllRayleighHistogramsRisingEdge();
  void fitAllRayleighHistogramsRisingEdgeAndHalfFallingEdge();

  TH1D* getRayleighHistogram(Int_t freqInd);
  TH1D* getRayleighHistogramFromFrequencyMHz(Double_t freqMHz);
  TF1* getRayleighHistogramFit(Int_t freqInd);  

  TH2D* makeRayleigh2DHistogram();
  

  static TF1* makeRayleighFunction(TString name, Double_t xMin, Double_t xMax);
  static TString getRayleighFunctionText();
  TF1* constructFitFromAmplitude(Int_t freqInd, Double_t amplitude);
  TF1* constructFitFromRayleighAmplitude(Int_t freqInd);
  TF1* constructFitFromRayleighAmplitudeRisingEdge(Int_t freqInd);
  TF1* constructFitFromRayleighAmplitudeRisingEdgeAndHalfFalling(Int_t freqInd);  

  Double_t deltaFMHz;  
  Double_t summedPowSpec[NUM_FREQS]; ///< How we store the average power spectrum in the summed case
  TH1D* hRayleighs[NUM_FREQS]; ///< Histograms for Rayleigh distributions
  TF1* hRayleighFits[NUM_FREQS]; ///< Fits to Rayleigh distributions histograms
  std::vector<Double_t> psdOutliers[NUM_FREQS];///< Storage for outlier values, if get too many then add them to rayeligh histogram
  // UInt_t maxNumOutliers; ///< Critical number of outliers before adding them, set in constructor.
  Int_t count;

  
  Double_t rayleighFitChiSquares[NUM_FREQS];
  Double_t rayleighFitChiSquaresRisingEdge[NUM_FREQS];
  Double_t rayleighFitChiSquaresRisingEdgeAndHalfFalling[NUM_FREQS];

  Double_t rayleighAmplitudes[NUM_FREQS];
  Double_t rayleighAmplitudesRisingEdge[NUM_FREQS];
  Double_t rayleighAmplitudesRisingEdgeAndHalfFalling[NUM_FREQS];

  Int_t rayleighNdf[NUM_FREQS];
  Int_t rayleighNdfRisingEdge[NUM_FREQS];
  Int_t rayleighNdfRisingEdgeAndHalfFalling[NUM_FREQS];

  Double_t rayleighFitChiSquaresFullRange[NUM_FREQS];
  Double_t rayleighFitChiSquaresRisingEdgeFullRange[NUM_FREQS];
  Double_t rayleighFitChiSquaresRisingEdgeAndHalfFallingFullRange[NUM_FREQS];

  Int_t rayleighNdfFullRange[NUM_FREQS];
  Int_t rayleighNdfRisingEdgeFullRange[NUM_FREQS];
  Int_t rayleighNdfRisingEdgeAndHalfFallingFullRange[NUM_FREQS];

  Double_t xHigh[NUM_FREQS];
  Double_t xHighRisingEdge[NUM_FREQS];
  Double_t xHighRisingEdgeAndHalfFalling[NUM_FREQS];

  Double_t outliers[NUM_FREQS][MAX_NUM_OUTLIERS];
  Int_t numOutliers[NUM_FREQS];  
  
  ClassDef(AveragePowerSpectrum, 9);
};

  



#endif









