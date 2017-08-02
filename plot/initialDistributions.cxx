#include "SummarySet.h"
#include "ProgressBar.h"
#include "OutputConvention.h"
#include "AnitaEventSummary.h"

#include "TH2D.h"
#include "AnalysisPlot.h"
#include "AnalysisCuts.h"

using namespace Acclaim;

int main(int argc, char* argv[]){

  if(argc!=2){
    std::cerr << argv[0] << "summaryTreeFileGlob" << std::endl;
    return 1;
  }

  OutputConvention oc(1, argv);
  TFile* fOut = oc.makeFile();

  const char* summaryTreeFileGlob = argv[1];
  SummarySet ss(summaryTreeFileGlob);

  const Long64_t N = ss.N();
  std::cout << "Processing " << summaryTreeFileGlob << " with " << N << " entries." << std::endl;

  AnalysisPlot* hPeakVsTime = ss.bookTimeAnalysisPlot("hPeakVsTime", "Higher map peak vs time", 1024, 128, 0, 1);
  hPeakVsTime->addCut(AnalysisCuts::isAboveHorizontal, "isAboveHorizontal", "Above Horizontal");
  hPeakVsTime->addCut(AnalysisCuts::isTaggedAsWaisPulser, "isTaggedAsWaisPulser", "WAIS Pulser");
  
  ProgressBar p(N);
  for(Long64_t entry=0; entry < N; entry++){

    ss.getEntry(entry);
    AnitaEventSummary* sum = ss.summary();

    hPeakVsTime->Fill(sum, sum->realTime, sum->higherPeak().value);
    
    p.inc(entry, N);
  }

  hPeakVsTime->Write();
  delete hPeakVsTime;
  
  fOut->Write();
  fOut->Close();
  return 0;
}














