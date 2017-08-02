#include "SummarySet.h"
#include "ProgressBar.h"
#include "OutputConvention.h"
#include "AnitaEventSummary.h"

#include "TH2D.h"
#include "QualityCut.h"
#include "AnalysisPlot.h"

using namespace Acclaim;


int isAboveHorizontal(AnitaEventSummary* sum){
  return sum->higherPeak().theta > 0 ? 1 : 0;
}


int isTaggedAsWaisPulser(AnitaEventSummary* sum){
  return sum->flags.pulser == AnitaEventSummary::EventFlags::WAIS;
}


int main(int argc, char* argv[]){

  if(argc!=2){
    std::cerr << argv[0] << "summaryTreeFileGlob" << std::endl;
    return 1;
  }

  OutputConvention oc(1, argv);
  TFile* fOut = oc.makeFile();

  const char* summaryTreeFileGlob = argv[1];

  // make note of what we've done here...
  SummarySet ss(summaryTreeFileGlob);

  const Long64_t N = ss.N();
  std::cout << "Processing " << summaryTreeFileGlob << " with " << N << " entries." << std::endl;

  TString peakVsTimeName = TString::Format("hHigherPeakVsEventNumber");
  AnalysisPlot* hPeakVsTime = new AnalysisPlot("hPeakVsTime", "Higher map peak vs time", 1024, ss.getFirstTime(), ss.getLastTime(), 128, 0, 1);
  hPeakVsTime->addCut(2, "isAboveHorizontal", "Above Horizontal", isAboveHorizontal);
  hPeakVsTime->addCut(2, "isTaggedAsWaisPulser", "WAIS Pulser", isTaggedAsWaisPulser);  
    
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














