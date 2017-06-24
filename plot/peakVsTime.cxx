#include "SummarySet.h"
#include "ProgressBar.h"
#include "OutputConvention.h"
#include "AnitaEventSummary.h"

#include "TH2D.h"
#include "QualityCut.h"

using namespace Acclaim;

int main(int argc, char* argv[]){

  const int numSets = argc - 1;

  if(numSets==0){
    std::cerr << argv[0] << "summaryTreeFileGlob0 summaryTreeFileGlob1 summaryTreeFileGlob2..." << std::endl;
    return 1;
  }

  OutputConvention oc(1, argv);
  TFile* fOut = oc.makeFile();


  for(int set=0; set < numSets; set++){
    const char* summaryTreeFileGlob = argv[1+set];

    // make note of what we've done here...
    TString noteName = TString::Format("Set%d", set);
    TNamed note(noteName.Data(), summaryTreeFileGlob);
    note.Write();

    SummarySet ss(summaryTreeFileGlob);

    const Long64_t N = ss.N();
    std::cout << "Processing " << summaryTreeFileGlob << " with " << N << " entries." << std::endl;

    ss.first();
    AnitaEventSummary* sum = ss.summary();
    const UInt_t firstTime = sum->realTime;
    const UInt_t firstEventNumber = sum->eventNumber;
    ss.last();
    const UInt_t lastTime = sum->realTime;
    const UInt_t lastEventNumber = sum->eventNumber;

    TString peakVsTimeName = TString::Format("hHigherPeakVsEventNumber%d", set);
    // TH2D h(peakVsTimeName, "Image Peak vs. Time", 1024, firstTime, lastTime, 128, 0, 1);
    TH2D h(peakVsTimeName, "Image Peak vs. eventNumber", 1024, firstEventNumber, lastEventNumber, 128, 0, 1);

    ProgressBar p(N);
    for(Long64_t entry=0; entry < N; entry++){

      ss.getEntry(entry);
      AnitaEventSummary* sum = ss.summary();

      if(QualityCut::passedAll(sum, true)){

        // h.Fill(sum->realTime, sum->higherPeak().value);
      }
      h.Fill(sum->eventNumber, sum->higherPeak().value);

      p.inc(entry, N);
    }
    h.Write();
  }

  fOut->Write();
  fOut->Close();
  return 0;
}














