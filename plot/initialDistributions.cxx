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
  hPeakVsTime->addCut(&AnalysisCuts::isAboveHorizontal);
  hPeakVsTime->addCut(&AnalysisCuts::isTaggedAsWaisPulser);
  hPeakVsTime->addCut(&AnalysisCuts::higherPol);
  hPeakVsTime->addCut(&AnalysisCuts::hasSourceLocation);

  AnalysisProf* hDeltaAngleSun = new AnalysisProf("hDeltaAngleSun", "Angle between peak and sun;#delta#phi (Degrees);#delta#theta (Degrees)", 128, -5, 5, 128, -5, 5);
  hDeltaAngleSun->addCut(&AnalysisCuts::isAboveHorizontal);
  hDeltaAngleSun->addCut(&AnalysisCuts::higherPol);

  AnalysisPlot* hPeakThetaVsTime = ss.bookTimeAnalysisPlot("hPeakThetaVsTime", "Peak elevation vs time; realTime; Peak #theta (Degrees)", 1024, 128, -90, 90);
  hPeakThetaVsTime->addCut(&AnalysisCuts::isAboveHorizontal);
  hPeakThetaVsTime->addCut(&AnalysisCuts::isTaggedAsWaisPulser);
  hPeakThetaVsTime->addCut(&AnalysisCuts::higherPol);
  hPeakThetaVsTime->addCut(&AnalysisCuts::hasSourceLocation);
  
  AnalysisPlot* hPeakBearingVsTime = ss.bookTimeAnalysisPlot("hPeakBearingVsTime", "Bearing of map peak vs time; realTime; Peak #phi (Degrees)", 1024, 256, 0, 360);
  hPeakBearingVsTime->addCut(&AnalysisCuts::isAboveHorizontal);
  hPeakBearingVsTime->addCut(&AnalysisCuts::isTaggedAsWaisPulser);
  hPeakBearingVsTime->addCut(&AnalysisCuts::higherPol);
  hPeakBearingVsTime->addCut(&AnalysisCuts::hasSourceLocation);
  
  ProgressBar p(N);
  for(Long64_t entry=0; entry < N; entry++){

    ss.getEntry(entry);
    AnitaEventSummary* sum = ss.summary();

    hPeakVsTime->Fill(sum, sum->realTime, sum->higherPeak().value);
    hDeltaAngleSun->Fill(sum, sum->dPhiSun(), sum->dThetaSun(), sum->higherPeak().value);
    hPeakThetaVsTime->Fill(sum, sum->realTime, sum->higherPeak().theta);
    hPeakBearingVsTime->Fill(sum, sum->realTime, sum->peakBearing());
    
    p.inc(entry, N);
  }
  
  fOut->Write();
  fOut->Close();
  return 0;
}














