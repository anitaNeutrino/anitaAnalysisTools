/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             Efficiently plot histograms/profiles for arbitrary subdivisions of the data.
***********************************************************************************************************/

#ifndef ANALYSIS_PLOT_H
#define ANALYSIS_PLOT_H

#include "TNamed.h"
#include "AnalysisCuts.h"

class AnitaEventSummary;
class TH1;

namespace Acclaim{

/** 
 * @class AnalysisPlot 
 * @brief Does the hard work of histogramming variables in the presence of cuts.
 * 
 * Purpose:
 * TChain::Draw is bloody slow even on the 10% ANITA data.
 * Suppose you want to make a histogram with a set of cuts, and quickly look at the opposite set of cuts.
 * Calling TChain::Draw N times, takes 5N mins (on my reletively fast laptop).
 * This class is designed to produce a plot with many possible sets of cuts applied looping through the data once.
 * The final desired histogram can be assembled from some combination of the produced subplots.
 * 
 * 
 * How it works:
 * This class creates a set of histograms (either TH1D or TH2D) for a set of cuts it is told about.
 * The dimension of the histograms is determined in the constructor, if nBinsY = 0 it's a TH1D, TH2D otherwise.
 * The cuts take the form of minimalistic objects with an apply(const AnitaEventSummary*) function, which retuns an int.
 * See AnalysisCuts.h/AnalysisCuts.cxx for examples.
 * The AnalysisCuts must conform to some rules for this class to work! See AnalysisCuts.h for more info!
 * 
 * Drawing:
 * This class then has it's own slightly different and hopefully useful implementation of Draw().
 * Draw(Option_t*, const TString& selection) creates a new histogram and adds all the contained histograms with names matching
 * the regexp contained in selection.
 * The histogram names are produced programatically form the combination of cuts applied.
 * Therefore by using an appropriate regexp one can draw the histogram matching your desired set of cuts. 
 * For example, one could draw a histogram matching the sum of all the member histograms with Draw("colz", "*")
 * (Actually this is the default behaviour of Draw("colz") if no argument is supplied).
 * Cut labels can be inspected with Print(), which dumps the contained cut labels to the terminal.
 * If that's not enough Print("v") prints the names of all the contained histograms too, although this list gets long for many cuts.
 * 
 * Behaviour:
 * AnalysisPlot objects can be saved in ROOT files.
 * It is designed to mimic the default histogram behaviour by automagically being appended to the current ROOT directory.
 * i.e. it should be saved to the current file without Write() being explicitly called.
 * 
 * Notes:
 * A histogram is allocated for every combination of cut function return values (from the AnalysisCut functions).
 * This means that memory requirements increase exponentially with added cuts so be careful!
 * For example, the memory required for just the bins of a TH2D histogram with 1024 x 128 bins is approximately
 * sizeof(double)*128*1024 = 1,048,576 ~ 1MB (yeah yeah over/underflow not accounted for I know)
 * For 5 cuts with just 2 possible return values, I would require 2^{5} = 32 of them, which is 32MB
 * For 10 cuts with just 2 possible return values, I would require 2^{10} = 1024, so that would be 1GB,
 * For 20 cuts with just 2 possible return values, I would require 2^{20} = 1,048,576 so that would be 1TB!!!.
 */

class AnalysisPlot : public TNamed {

  friend class SummarySelector;
  
 public:
  AnalysisPlot() {;}
  AnalysisPlot(const char* name, const char* title, int nBinsX, double xMin, double xMax, int nBinsY=0, double yMin=0, double yMax=0);
  virtual ~AnalysisPlot();

  size_t addCut(const AnalysisCuts::AnalysisCut* cut);

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
  std::vector<Int_t> ns; /// The maximum return values of the cuts
  std::vector<Int_t> indexMultipliers; /// Used to map cut to histogram index
  std::vector<TString> cutNames; /// the string to associate with the cut in histogram names

  // These don't persist (deliberately) due to the //! comment string
  std::vector<const AnalysisCuts::AnalysisCut*> analysisCuts; //! The cuts to apply, do not persist, the only data needed later are saved in cutNames
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








