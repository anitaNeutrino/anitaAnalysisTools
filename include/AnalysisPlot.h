#ifndef ANALYSIS_PLOT_H
#define ANALYSIS_PLOT_H

#include "TNamed.h"

class AnitaEventSummary;

class TH1;

namespace Acclaim{

class AnalysisPlot : public TNamed {

 public:
  AnalysisPlot() {;}
  AnalysisPlot(const char* name, const char* title, int nx, double xMin, double xMax, int ny=0, double yMin=0, double yMax=0);
  virtual ~AnalysisPlot();

  size_t addCut(int nRetVals, const char* nameStr, const char* titleStr, int(*AnalysisCut)(AnitaEventSummary*));

  int Fill(AnitaEventSummary* sum, double xVal, double yVal=1, double zVal=1);
  void Draw(Option_t* opt="");
  void Draw(Option_t* opt, const TString& selection);
  void Print(Option_t* = "") const;

 protected:

  TH1* makeHist(const char* name, const char* title);

  double fNx;
  double fMinX;
  double fMaxX;

  double fNy;
  double fMinY;
  double fMaxY;

  // typical 2D histogram size
  // sizeof(double)*128*1024 = 1,048,576 ~ 1MB
  // suppose I have 6 cuts, 2^{6} = 64. Which is 64MB
  // 10 cuts would be 2^{10} = 1024, so that would be 1GB...
  // I suppose I should be wary of storing more than that...

  std::vector<TH1*> hs; /// Stores all the created histograms
  std::vector<Int_t> ns; /// Stores the nRetVals for the cuts
  std::vector<Int_t> indexMultipliers; /// Used to map cut to histogram index
  std::vector<TString> cutNames; /// the string to associate with the cut in histogram names
  std::vector<int(*)(AnitaEventSummary*)> cutFuncs; /// pointers to the functions to apply

  // This object does not persist, and is only for ease of drawing histograms hence the !
  std::vector<TH1*> hDummies; //! Dummy histograms for drawing, don't persist
  


  ClassDef(AnalysisPlot, 1);
};
}


#endif








