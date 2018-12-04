#include "AnitaDataset.h"
#include "AcclaimCmdLineArgs.h"
#include "AnalysisFlow.h"
#include "FilterStrategy.h"
#include "CrossCorrelator.h"
#include "AcclaimCorrelationSummary.h"
#include "ProgressBar.h"


#include "TTree.h"

int main(){

  int run = 125;
  std::cout <<  AnitaVersion::get()  <<  std::endl;
  AnitaDataset d(run);
  auto cc = std::make_shared<Acclaim::CrossCorrelator>();
  
  FilterStrategy fs;


  TString fileName = TString::Format("AcclaimCorrelationSummary_%d.root", run);

  TFile* fOut = new TFile(fileName,  "recreate");
  TTree* t = new TTree("corrTree", "Tree of WAIS correlation pairs");
  Acclaim::CorrelationPair* corrPair = nullptr;
  t->Branch("corrPair", &corrPair);
  
  const Long64_t n = d.N();
  Acclaim::ProgressBar p(n);
  for(Long64_t entry=0; entry < n; entry++){
    d.getEntry(entry);

    FilteredAnitaEvent ev(d.useful(), &fs, d.gps(),  d.header());

    UsefulAdu5Pat pat(d.gps());

    Double_t theta, phi;
    pat.getThetaAndPhiWaveWaisDivide(theta, phi);

    theta*=-TMath::RadToDeg();
    phi*=TMath::RadToDeg();
  
    auto s = cc->makeSummary(AnitaPol::kVertical, &ev,  phi, theta);

    for(int pair = 0; pair < s->N(); pair++){
      auto p = s->get(pair);
      corrPair = const_cast<Acclaim::CorrelationPair*>(&p);

      t->Fill();
    }

    p.inc(entry, n);
  }

  fOut->Write();
  fOut->Close();
  
}
