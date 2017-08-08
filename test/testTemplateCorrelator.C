#include "CrossCorrelator.h"
#include "AnitaDataset.h"
#include "ProgressBar.h"

using namespace Acclaim;

void testTemplateCorrelator(){

  int run = 352; //170;
  UInt_t eventNumber = 60831701; //14102511;

  TemplateCorrelator* tc = new TemplateCorrelator(run, eventNumber);

  // AnitaDataset d(run, true);
  AnitaDataset d(run);  
  Long64_t N = d.N();
  ProgressBar p(N);

  FilterStrategy empty;

  TFile* fOut = new TFile("/tmp/testTemplateCorrelator.root", "recreate");
  TString treeName = TString::Format("templateTree_%d_%u", run, eventNumber);
  TTree* t = new TTree(treeName, treeName);
  double blastiness[AnitaPol::kNotAPol];
  t->Branch("blastiness[2]", blastiness, "blastiness[2]/D");
  
  for(Long64_t entry=0; entry < N; entry++){
    d.getEntry(entry);
    FilteredAnitaEvent fEv(d.useful(), &empty, d.gps(), d.header());

    tc->correlateEvent(&fEv);

    for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
      blastiness[pol] = tc->getPeakCorrelation(pol);
      // std::cerr << pol << "\t" << blastiness[pol] << std::endl;
    }

    t->Fill();
    p.inc(entry, N);
    std::cerr << std::endl; break; 
  }

  AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;
  int ant = 1;
  auto gr = tc->getCrossCorrelationGraph(pol, ant);
  gr->Draw();

  fOut->Write();
  fOut->Close();
  
}
