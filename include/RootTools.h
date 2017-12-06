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

#include "AnitaConventions.h"
#include <iostream>
#include "TString.h"

class TGraph;
class TH1D;
class TH2;
class TH2D;
class TTree;
class TCanvas;
class TLegend;
class TPad;
class TAxis;
class TChain;

class UsefulAnitaEvent;
class RawAnitaHeader;
class Adu5Pat;

/** @mainpage
 * Yes, the name ACCLAIM is kind of lame, but passed a certain point anitaAnalysisTools was a a bit generic, and it was the best I could do in a lunch break.
 * A lot has changed since I last made a concerted effort to document my code.
 * Significant restructuring makes the code fit in better with the AnalysisFramework originally designed back in Chicago back in 2015, and the code is much better modularised.
 */


/** @namespace Acclaim
 * @brief Namespace which wraps everything in the library
 *
 * 
 */

namespace Acclaim
{
  /** @namespace RootTools
   * @brief My commonly used, general functions to manipulate ROOT objects; so I only ever write them once.
   *
   * This lovingly curated namespace can be imported into cint/cling with gSystem->Load('anitaAnalysisTools.so"), replace .so with .dylib if you're on a Mac.
   * This namespace is a collection of little routines to modify things like TGraphs in a general way, for things that are too simple to deserve a dedicated class of their own.
   */

  namespace RootTools{

  /** 
   * Return the opposite polarisation, prints an error if there is non-sensical input
   * Useful for getting the cross-polarisation.
   * 
   * @param pol is the polarisation to swap
   * 
   * @return the other polarisation if valid input, kNotAPol otherwise
   */
    inline AnitaPol::AnitaPol_t swapPol(AnitaPol::AnitaPol_t pol){
      switch(pol){
        case AnitaPol::kHorizontal:
          return AnitaPol::kVertical;
        case AnitaPol::kVertical:
          return AnitaPol::kHorizontal;
        default:
          std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", unknown polarisation!\t" << pol << std::endl;        
          return AnitaPol::kNotAPol;
      }
    }

    // for histogramming the whole flight.
    static const UInt_t a3StartTime = 1418869215; //!< Time of first event in ANITA-3 run 352
    static const UInt_t a3EndTime = 1420782770; //!< Time of last event + 1 second
  
    void writeTGraph(TGraph* gr, TString name);
    void printArray(int n, double* array, TString delimiter = ", ", TString start = "{" ,TString end = "}\n");
    void printYVals(const TGraph* gr, TString delimiter = ", ", TString start = "{" ,TString end = "}\n");
    void printXVals(const TGraph* gr, TString delimiter = ", ", TString start = "{" ,TString end = "}\n");
    void printTGraphInfo(const TGraph* gr);


    /* Get info about input */
    Int_t getIndexOfMaximum(Int_t len, Double_t* arr);
    std::vector<Int_t> getIndicesOfNans(TGraph* gr);
    Double_t getSumOfYVals(const TGraph* gr);
    void getMaxMin(const TGraph* gr, Double_t& max, Double_t& min);
    void getMeanAndRms(TGraph* gr, Double_t& mean, Double_t& rms);
    void getMaxMin(const TGraph* gr, Double_t& maxY, Double_t& maxX, Double_t& minY, Double_t& minX);
    void getMaxMinWithinLimits(const TGraph* gr, Double_t& maxY, Double_t& maxX,
			       Double_t& minY, Double_t& minX,
			       Double_t lowerLimit, Double_t upperLimit);
    void getLocalMaxToMin(const TGraph* gr,
			  Double_t& maxY, Double_t& maxX,
			  Double_t& minY, Double_t& minX);
    void getLocalMaxToMinWithinLimits(const TGraph* gr,
				      Double_t& maxY, Double_t& maxX,
				      Double_t& minY, Double_t& minX,
				      Double_t lowerLimit, Double_t upperLimit);

    Int_t getPeakBinOfHistogram(TH1D* h);
    Double_t getPeakBinOfHistogram(TH2D* hist, Int_t& binx, Int_t& biny);
    Double_t getLowBinEdgeOfHistogramPeak(TH1D* h);
    Double_t getFullWidthHalfMax(TH1D* h);
    Int_t getBit(UInt_t bitIndex, UInt_t bitMask);
    Int_t getNumBitsSet(Int_t numBitsToCheck, UInt_t bitMask);

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
    void setBinLabelsFromLowEdge(TAxis* ax, const char* format);

    /* Make new output based on input */
    TGraph* makeNormalized(TGraph* gr); ///< Creates new TGraph (leaving original unchanged) with mean = 0 & RMS = 1
    TGraph* makeNormalized(TGraph* gr, Double_t& mean, Double_t& rms);
    TGraph* makeSortedTGraph(TTree* tree, TString drawText, TString cutString, Double_t wrapValue);
    TGraph* makeLinearlyInterpolatedGraph(TGraph* grIn, Double_t dt);
    TGraph* makeDerivativeTGraph(const TGraph* gr);
    TGraph* makeUnwrappedCorrelationGraph(const TGraph* gr);
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

    TPad* makeSubPad(TPad* parentPad, double xlow, double ylow, double xup, double yup, TString suffix);

    UsefulAnitaEvent* makeGaussianEvent(UInt_t eventNumber);


    /* Load ROOT data into chains quickly*/
    TChain* getHeadChain(Int_t firstRun, Int_t lastRun, RawAnitaHeader*& headPtr);
    TChain* getAdu5PatChain(Int_t firstRun, Int_t lastRun, Adu5Pat*& pat);


    Int_t isMinBiasSampleEvent(const RawAnitaHeader* header);

    enum {
      kUnknownPeakTime = -9999
    };
    std::pair<double, double> findSmallestWindowContainingFracOfPower(const TGraph* grPow, double fracOfPowerInWindow);
    double getTimeIntegratedPower(const TGraph* gr, int firstSamp=0, int lastSamp=-1);

    void tokenize(std::vector<TString>& tokenizedOutput, const char* inputString, const char* separator);
    void tokenize(std::vector<TString>& tokenizedOutput, const char* inputString, const std::vector<const char*>& separators);

    template <class T>
    inline bool vectorContainsValue(const std::vector<T>& vec, const T& value){
      return std::find(vec.begin(), vec.end(), value)!=vec.end();
    }
  }
}
#endif
