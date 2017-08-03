#include "SummarySet.h"
#include "ProgressBar.h"
#include "OutputConvention.h"
#include "AnitaEventSummary.h"

#include "TH2D.h"
#include "AnalysisPlot.h"
#include "AnalysisCuts.h"

using namespace Acclaim;


/** 
 * A script initially developed for looking at the 10% data set.
 * Although, I suppose it can look at any set of AnitaEventSummary files
 * 
 * @param summaryTreeFileGlob is a globbing expression for trees containing AnitaEventSummaries
 * (make sure you add a single quotes at the start and end to prevent the shell expanding your wildcard!)
 * 
 * @return 0 on success, 1 if too many arguments
 */

int main(int argc, char* argv[]){

  if(argc!=2){
    std::cerr << argv[0] << ": [summaryTreeFileGlob]" << std::endl;
    return 1;
  }

  OutputConvention oc(1, argv);
  TFile* fOut = oc.makeFile();

  const char* summaryTreeFileGlob = argv[1];
  SummarySet ss(summaryTreeFileGlob);  

  const Long64_t N = ss.N();
  std::cout << "Processing " << summaryTreeFileGlob << " with " << N << " entries." << std::endl;

  AnalysisPlot* hPeakVsTime = ss.bookTimeAnalysisPlot("hPeakVsTime", "Higher map peak vs time; realTime; Map peak", 1024, 128, 0, 1);
  hPeakVsTime->addCut(AnalysisCuts::isAboveHorizontal, "isAboveHorizontal", "Above Horizontal");
  hPeakVsTime->addCut(AnalysisCuts::isTaggedAsWaisPulser, "isTaggedAsWaisPulser", "WAIS Pulser");
  hPeakVsTime->addCut(AnalysisCuts::higherPol, "higherPol", "Pol");

  AnalysisProf* hDeltaAngleSun = new AnalysisProf("hDeltaAngleSun", "Angle between peak and sun;#delta#phi (Degrees);#delta#theta (Degrees)", 128, -5, 5, 128, -5, 5);
  hDeltaAngleSun->addCut(AnalysisCuts::isAboveHorizontal, "isAboveHorizontal", "Above Horizontal");
  hDeltaAngleSun->addCut(AnalysisCuts::isTaggedAsWaisPulser, "isTaggedAsWaisPulser", "WAIS Pulser");
  hDeltaAngleSun->addCut(AnalysisCuts::higherPol, "higherPol", "Pol");

  AnalysisPlot* hPeakThetaVsTime = ss.bookTimeAnalysisPlot("hPeakThetaVsTime", "Higher map peak vs time; realTime; Map peak", 1024, 128, -90, 90);
  hPeakThetaVsTime->addCut(AnalysisCuts::isAboveHorizontal, "isAboveHorizontal", "Above Horizontal");
  hPeakThetaVsTime->addCut(AnalysisCuts::isTaggedAsWaisPulser, "isTaggedAsWaisPulser", "WAIS Pulser");
  hPeakThetaVsTime->addCut(AnalysisCuts::higherPol, "higherPol", "Pol");
  
  ProgressBar p(N);
  for(Long64_t entry=0; entry < N; entry++){

    ss.getEntry(entry);
    AnitaEventSummary* sum = ss.summary();

    hPeakVsTime->Fill(sum, sum->realTime, sum->higherPeak().value);
    hDeltaAngleSun->Fill(sum, sum->dPhiSun(), sum->dThetaSun(), sum->higherPeak().value);
    hPeakThetaVsTime->Fill(sum, sum->realTime, sum->higherPeak().theta);
    
    p.inc(entry, N);
  }

  hPeakVsTime->Write();
  delete hPeakVsTime;

  hDeltaAngleSun->Write();
  delete hDeltaAngleSun;

  hPeakThetaVsTime->Write();
  delete hPeakThetaVsTime;
  
  
  fOut->Write();
  fOut->Close();
  return 0;
}














