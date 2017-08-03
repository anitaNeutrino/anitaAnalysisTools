/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             Easily plot parameters falling into different cut bins
***********************************************************************************************************/

#ifndef ANALYSIS_PLOT_H
#define ANALYSIS_PLOT_H

#include "TNamed.h"

class AnitaEventSummary;
class TH1;

namespace Acclaim{

/** 
 * @class AnalysisPlot does the hard work of histogramming variables in the presence of cuts.
 * 
 * This class creates a set of histograms (either TH1D or TH2D) for a set of cuts it is told about.
 * The dimension of the histograms is determined in the constructor, if nBinsY = 0 it's a TH1D, TH2D otherwise.
 * The cuts take the form of functions, which take a const AnitaEventSummary* and return and int.
 * See the AnalysisCuts namespace for examples of functions to add with AnalysisPlot::addCut().
 * There are some STRICT rules the AnalysisCut functions must conform to for this class to work! See AnalysisCuts.h!
 * AnalysisPlot objects can be saved in ROOT files.
 * 
 * Please note that a histogram is allocated for every combination of return values from the AnalysisCut functions.
 * This means that memory requirements increase exponentially with added cuts so be careful!
 * 
 * For example, the memory required for typical 2D histogram with 1024 x 128 bins is
 * sizeof(double)*128*1024 = 1,048,576 ~ 1MB
 * For 6 cuts, I require 2^{6} = 64 of them, which is 64MB
 * For 10 cuts would be 2^{10} = 1024, so that would be 1GB...
 */

class AnalysisPlot : public TNamed {

 public:
  AnalysisPlot() {;}
  AnalysisPlot(const char* name, const char* title, int nBinsX, double xMin, double xMax, int nBinsY=0, double yMin=0, double yMax=0);
  virtual ~AnalysisPlot();

  // the first parameter is a pointer to a function that takes a const AnitaEventSummary* as an argument
  size_t addCut(int(*)(const AnitaEventSummary*), const char* nameStr, const char* titleStr);

  virtual int Fill(const AnitaEventSummary* sum, double xVal, double yVal=1, double zVal=1);
  void Draw(Option_t* opt="");
  void Draw(Option_t* opt, const TString& selection);
  void Print(Option_t* = "") const;

 protected:

  virtual int getIndexFromCuts(const AnitaEventSummary* sum);
  virtual TH1* makeHist(const char* name, const char* title) const;

  double fNx;
  double fMinX;
  double fMaxX;

  double fNy;
  double fMinY;
  double fMaxY;
  std::vector<TH1*> hs; /// Stores all the created histograms
  std::vector<Int_t> ns; /// Stores the nRetVals for the cuts
  std::vector<Int_t> indexMultipliers; /// Used to map cut to histogram index
  std::vector<TString> cutNames; /// the string to associate with the cut in histogram names
  std::vector<int(*)(const AnitaEventSummary*)> analysisCuts; /// pointers to the functions to apply

  // This object does not persist, and is only for ease of drawing histograms hence the !
  std::vector<TH1*> hDummies; //! Dummy histograms for drawing, don't persist
  


  ClassDef(AnalysisPlot, 1);
};






/** 
 * @class AnalysisProf is the same an analysis plot, but replaces histograms with profiles
 * i.e. TH1D -> TProfile, TH2D -> TProfile2D
 *  */

class AnalysisProf : public AnalysisPlot {

 public:
  AnalysisProf(){;}
  virtual ~AnalysisProf(){;}
  AnalysisProf(const char* name, const char* title, int nBinsX, double xMin, double xMax, int nBinsY=0, double yMin=0, double yMax=0);
  virtual int Fill(const AnitaEventSummary* sum, double xVal, double yVal=1, double zVal=1);
  
 protected:
  virtual TH1* makeHist(const char* name, const char* title) const;
  
  ClassDef(AnalysisProf, 1);
};


}




#endif








