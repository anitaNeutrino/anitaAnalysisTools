/* -*- C++ -*-.*/
/**************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             Slightly complicated functions to manipulate ROOT objects should go in here.
	     I want to only ever write them once, so when I need to do something I'll add it here.
*************************************************************************************************************** */

#ifndef ROOTTOOLS_H
#define ROOTTOOLS_H

#include "TObjArray.h"
#include "TGraph.h"
#include "TPaveStats.h"
#include "TTree.h"
#include "TAxis.h"
#include "TMath.h"
#include "iostream"
#include "TH2.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "cfloat"
#include "TColor.h"
#include "TChain.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"

#include "RawAnitaHeader.h"
#include "Adu5Pat.h"

/** @namespace RootTools
 * @brief My commonly used, general functions to manipulate ROOT objects; so I only ever write them once.
 *
 * This lovingly curated namespace can be imported into CINT with gSystem->Load('libBensAnitaTools.so"). 
 * I find myself often writing little routines to modify things like TGraphs over and over again.
 * Functions that do little jobs, too simple to deserve a dedicated class of their own go in this namespace.
 * The idea that I only spend any time coding up a particular function once. 
 * Expect this namespace to be referenced a lot in my code.
*/

namespace RootTools{

  void writeTGraph(TGraph* gr, TString name);
  void printArray(int n, double* array, TString delimiter = ", ", TString start = "{" ,TString end = "}\n");
  void printYVals(TGraph* gr, TString delimiter = ", ", TString start = "{" ,TString end = "}\n");
  void printXVals(TGraph* gr, TString delimiter = ", ", TString start = "{" ,TString end = "}\n");
  void printTGraphInfo(TGraph* gr);


  /* Get info about input */
  Int_t getIndexOfMaximum(Int_t len, Double_t* arr);
  std::vector<Int_t> getIndicesOfNans(TGraph* gr);  
  Double_t getSumOfYVals(TGraph* gr);
  void getMaxMin(TGraph* gr, Double_t& max, Double_t& min);
  void getMeanAndRms(TGraph* gr, Double_t& mean, Double_t& rms);
  void getMaxMin(TGraph* gr, Double_t& maxY, Double_t& maxX, Double_t& minY, Double_t& minX);
  void getMaxMinWithinLimits(TGraph* gr, Double_t& maxY, Double_t& maxX, 
			     Double_t& minY, Double_t& minX, 
			     Double_t lowerLimit, Double_t upperLimit);
  void getLocalMaxToMin(TGraph* gr, 
			Double_t& maxY, Double_t& maxX, 
			Double_t& minY, Double_t& minX);
  void getLocalMaxToMinWithinLimits(TGraph* gr, 
				    Double_t& maxY, Double_t& maxX, 
				    Double_t& minY, Double_t& minX,
				    Double_t lowerLimit, Double_t upperLimit);

  Int_t getPeakBinOfHistogram(TH1D* h);
  Double_t getPeakBinOfHistogram(TH2D* hist, Int_t& binx, Int_t& biny);
  Double_t getLowBinEdgeOfHistogramPeak(TH1D* h);
  Double_t getFullWidthHalfMax(TH1D* h);
  Int_t getBit(UInt_t bitIndex, UInt_t bitMask);

  /* Do geometric things */
  Double_t getDeltaAngleDeg(Double_t angle1, Double_t angle2);
  

  /* Modify input */
  void subtractOffset(TGraph* gr, Double_t offset);
  void normalize(TGraph* gr, Double_t& mean, Double_t& rms);
  void normalize(TGraph* gr);
  void zeroPadTGraph(TGraph* gr, Int_t newLen, Double_t dt=0);
  void offsetTGraphXAxes(Int_t numGrs, TGraph* grs[], Double_t offsets[]);
  void multiplyTGraphYAxes(Int_t numGrs, TGraph* grs[], Double_t factors[]);


  /* Prettify */
  void makeZaxisScaleEqualAboutZero(TH2D* h);
  
  /* Make new output based on input */
  TGraph* makeNormalized(TGraph* gr); ///< Creates new TGraph (leaving original unchanged) with mean = 0 & RMS = 1
  TGraph* makeNormalized(TGraph* gr, Double_t& mean, Double_t& rms);
  TGraph* makeSortedTGraph(TTree* tree, TString drawText, TString cutString, Double_t wrapValue);
  TGraph* makeLinearlyInterpolatedGraph(TGraph* grIn, Double_t dt);
  TGraph* makeDerivativeTGraph(TGraph* gr);
  TGraph* makeUnwrappedCorrelationGraph(TGraph* gr);
  TGraph* interpolateWithStartTime(TGraph* grIn, Double_t startTime, Double_t dt, Int_t nSamp);

  TH1D* plotsZaxisDist(TH2* h2, TString hName, Int_t nBins, Double_t xMin, Double_t xMax);
  TCanvas* drawArrayOfHistosPrettily(TH1D* hs[], Int_t numHists, TCanvas* can=NULL, 
				     Double_t* colWeights = NULL);
  TCanvas* drawArrayOfTGraphsPrettily(TGraph* grs[], Int_t numGrs, 
				      TString drawOpt = "l", TCanvas* can=NULL,
				      Double_t* colWeights = NULL);
  TLegend* makeLegend(TGraph* grs[], Int_t numGrs, TString titles[], TString opt = "l", 
		      Double_t minX=0.8, Double_t minY=0.8,Double_t maxX=1, Double_t maxY=1);
  TLegend* makeLegend(TH1D* hs[], Int_t numHists, TString titles[], TString opt = "l", 
		      Double_t minX=0.8, Double_t minY=0.8,Double_t maxX=1, Double_t maxY=1);



  void saveCanvas(TCanvas* c1, TString fileName);
  void setWhiteZeroColorScale();
  void draw2D(TH2D* hist, TString opt);
  Int_t getColorFracThroughPalette(Int_t index, Int_t maxVal);
  TCanvas* drawHistsWithStatsBoxes(Int_t numHists, TH1D* hs[], TString drawOpt, TString statsOption);
  TString getAntName(AnitaPol::AnitaPol_t pol, Int_t antInd);
  

  /* Load ROOT data into chains quickly*/
  TChain* getHeadChain(Int_t firstRun, Int_t lastRun, RawAnitaHeader*& headPtr);
  TChain* getAdu5PatChain(Int_t firstRun, Int_t lastRun, Adu5Pat*& pat);  

}

#endif
