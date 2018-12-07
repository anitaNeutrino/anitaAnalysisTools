#include "AnitaDataset.h"
#include "AcclaimCmdLineArgs.h"
#include "AnalysisFlow.h"
#include "FilterStrategy.h"
#include "CrossCorrelator.h"
#include "AcclaimCorrelationSummary.h"
#include "ProgressBar.h"
#include "UCUtil.h" //UCorrelator

#include "TTree.h"

int main(int argc, char* argv[]){

  if(argc!=2){
    std::cerr << "Usage: " << argv[0] << " run "<<std::endl;
    return 1;
  }
  
  int run = atoi(argv[1]);
  AnitaVersion::set(4);
  AnitaDataset d(run);
  auto cc = std::make_shared<Acclaim::CrossCorrelator>();
  
  FilterStrategy fs;

  auto calib = AnitaEventCalibrator::Instance();
  for(auto& channelsPerSurf : calib->relativePhaseCenterToAmpaDelays){
    for(auto& chanExtraDt : channelsPerSurf){
      // std::cout << "Setting " << chanExtraDt << " to ";
      chanExtraDt = 0;
      // std::cout << chanExtraDt << std::endl;
    }
  }


  TString fileName = TString::Format("AcclaimCorrelationSummary_%d.root", run);

  TFile* fOut = new TFile(fileName,  "recreate");
  TTree* t = new TTree("corrTree", "Tree of WAIS correlation pairs");
  Acclaim::CorrelationSummary* cs = nullptr;
  t->Branch("correlationSummary", &cs);
  
  const Long64_t n = d.N();
  Acclaim::ProgressBar p(n);
  for(Long64_t entry=0; entry < n; entry++){
    d.getEntry(entry);

    FilteredAnitaEvent ev(d.useful(), &fs, d.gps(), d.header());

    UsefulAdu5Pat pat(d.gps());

    Double_t theta, phi;
    pat.getThetaAndPhiWaveWaisDivide(theta, phi);

    theta*=-TMath::RadToDeg();
    phi*=TMath::RadToDeg();

    AnitaPol::AnitaPol_t pol;
    if(UCorrelator::isWAISVPol(&pat, d.header())){
      pol = AnitaPol::kVertical;
    }
    else if(UCorrelator::isWAISHPol(&pat, d.header())){
      pol = AnitaPol::kHorizontal;
    }
    else {
      pol = AnitaPol::kNotAPol;
      std::cerr << "Unknown WAIS pol" << std::endl;
    }

    if(pol != AnitaPol::kNotAPol){
      auto s = cc->makeSummary(pol, &ev,  phi, theta, d.gps());
      cs = const_cast<Acclaim::CorrelationSummary*>(s.get());
      t->Fill();
    }
    p.inc(entry, n);
  }

  fOut->Write();
  fOut->Close();
  
}
