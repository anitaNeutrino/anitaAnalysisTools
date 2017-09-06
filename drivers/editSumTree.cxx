#include "AnitaEventSummary.h"
#include "TTree.h"
#include "TFile.h"
#include "SummarySet.h"
#include "AnitaDataset.h"
#include "ProgressBar.h"
#include "AnalysisReco.h"
#include "InterferometricMap.h"

using namespace Acclaim;


/** 
 * Script to make a copy of summary trees, and hackily edit a couple of values if possible
 */
int main(int argc, char* argv[]){

  if(argc!=4){
    std::cerr << "Usage: " <<  argv[0] << " inputDir inputGlob run" << std::endl;
    return 1;
  }

  TString inDir = argv[1];
  TString inFileBaseName = argv[2];
  int run = atoi(argv[3]);

  TString inGlob = inDir + "/" + inFileBaseName + TString::Format("_%d_*.root", run);
  SummarySet ss(inGlob);

  TString outFileName = inFileBaseName + TString::Format("_%d.root", run);
  TFile* fOut = new TFile(outFileName, "recreate");
  TTree* sumTree = new TTree("sumTree", "sumTree");
  AnitaEventSummary* sum = NULL;
  sumTree->Branch("sum", &sum);

  AnitaDataset d(run);

  const Long64_t N = ss.N();
  ProgressBar p(N);
  for(Long64_t entry=0; entry < N; entry++){
    ss.getEntry(entry);
    AnitaEventSummary* inSum = ss.summary();

    d.getEvent(inSum->eventNumber);
    

    *sum = *inSum;


    const double phi0 = InterferometricMap::getBin0PhiDeg();
    for(int polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
      for(int peakInd=0; peakInd < sum->nPeaks[pol]; peakInd++){

        // hack to recover peakPhiSector...
        double phi_rough = sum->peak[pol][peakInd].phi - sum->peak[pol][peakInd].dphi_rough;
        int peakPhiSector = floor((phi_rough - phi0)/22.5);
        // std::cout << peakPhiSector << "\t" << phi_rough << std::endl;
        AnalysisReco::setTriggerInfoFromPeakPhi(d.header(), pol, peakPhiSector, sum->peak[pol][peakInd]);
      }
    }

    sumTree->Fill();
    p.inc(entry, N);
  }
  fOut->Write();
  fOut->Close();
  
  return 0;
}
