#include "AnitaDataset.h"
#include "AcclaimCmdLineArgs.h"
#include "AnalysisFlow.h"
#include "FilterStrategy.h"
#include "CrossCorrelator.h"
#include "AcclaimCorrelationSummary.h"
#include "ProgressBar.h"
#include "UCUtil.h" //UCorrelator

#include "TTree.h"
#include "TChain.h"


std::map<std::pair<int,  int>, TGraph* > fGrs;
std::map<std::pair<int,  int>, TH1D* > fH1s;
std::map<std::pair<int,  int>, TH2D* > fH2s;

TGraph* getGraph(int ant1, int ant2){
  auto it = fGrs.find({ant1, ant2});
  if(it != fGrs.end()){
    return it->second;
  }
  else{
    TGraph* gr = new TGraph();
    auto name = TString::Format("gr_%d_%d", ant1, ant2);
    auto title = name;
    title += ";Measured #deltat (ns); Expected #deltat (ns)";
    gr->SetName(name);
    gr->SetTitle(title);
    fGrs[{ant1, ant2}] = gr;
    return gr;
  }
}

TH1D* getHist(int ant1, int ant2){
  auto it = fH1s.find({ant1, ant2});
  if(it != fH1s.end()){
    return it->second;
  }
  else{
    auto name = TString::Format("h_%d_%d", ant1, ant2);
    
    auto title = name;
    title += ";(Measured - Expected) #deltat (ns); Events per bin";
    TH1D* h = new TH1D(name, title, 128, -5, 5);    
    fH1s[{ant1, ant2}] = h;
    return h;
  }
}


TH2D* getHist2(int ant1, int ant2){
  auto it = fH2s.find({ant1, ant2});
  if(it != fH2s.end()){
    return it->second;
  }
  else{
    auto name = TString::Format("h2_%d_%d", ant1, ant2);
    
    auto title = name;
    title += ";Measured #deltat (ns); Expected #deltat (ns); Events per bin";
    TH2D* h = new TH2D(name, title, 128, -20, 20, 128, -20, 20);
    fH2s[{ant1, ant2}] = h;
    return h;
  }
}


int main(int argc, char* argv[]){

  if(argc!=2){
    std::cerr << "Usage: " << argv[0] << " run "<<std::endl;
    return 1;
  }
  
  TChain* corrChain = new TChain("corrTree");
  corrChain->Add("AcclaimCorrelationSummary_125.root");

  Acclaim::CorrelationSummary* cs = nullptr;
  corrChain->SetBranchAddress("correlationSummary", &cs);

  auto geom = AnitaGeomTool::Instance();
  geom->usePhotogrammetryNumbers(true);
  for(int ant=0; ant < 48; ant++){
    for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
      Double_t phi1=geom->getAntPhiPositionRelToAftFore(ant, pol);
      Double_t r1=geom->getAntR(ant, pol);
      Double_t z1=geom->getAntZ(ant, pol);
      std::cout << "ant[pol] "  << ant << "[" << pol << "] " << r1 << ", " << phi1*TMath::RadToDeg() << ", " << z1 << "\n";
    }
  }
  return 0;

  const Long64_t n = corrChain->GetEntries();
  Acclaim::ProgressBar p(n);
  for(Long64_t entry=0; entry < n; entry++){
    corrChain->GetEntry(entry);

    double phiWave = cs->fPhiDeg*TMath::DegToRad();
    double thetaWave = -cs->fThetaDeg*TMath::DegToRad();
    UsefulAdu5Pat usefulPat(&cs->fPat);

    for(const auto& corrPair : cs->fPairs){

      if(corrPair.correlation > 0.4){
	double dtExpected = -usefulPat.getDeltaTExpected(corrPair.ant1, corrPair.ant2, phiWave, thetaWave);
	TGraph* gr = getGraph(corrPair.ant1, corrPair.ant2);      
	gr->SetPoint(gr->GetN(), corrPair.dt, dtExpected);

	TH1D* h = getHist(corrPair.ant1, corrPair.ant2);      
	h->Fill(corrPair.dt - dtExpected);

	TH2D* h2 = getHist2(corrPair.ant1, corrPair.ant2);      
	h2->Fill(corrPair.dt,  dtExpected);
	
	// std::cout << cs->fEventNumber << "\t" << corrPair.dt << "\t" << dtExpected << std::endl;
      }
    }    
    p.inc(entry,  n);
  }

  TFile* fOut = new TFile("measured_vs_expected.root", "recreate");
  for(auto it_pair_gr : fGrs){
    it_pair_gr.second->Write();
    delete it_pair_gr.second;
  }
  for(auto it_pair_h : fH1s){
    it_pair_h.second->Write();
    delete it_pair_h.second;
  }
  for(auto it_pair_h : fH2s){
    it_pair_h.second->Write();
    delete it_pair_h.second;
  }
  
  fOut->Write();
  fOut->Close();
}
