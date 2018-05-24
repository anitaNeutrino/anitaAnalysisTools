#include "ThermalChain.h"
#include "DrawStrings.h"
#include "ProgressBar.h"

using namespace Acclaim;

void plotWaisEfficiency(){

  ThermalChain tc("data/makeThermalTree_3*.root");

  tc.setCut(ThermalTree::isTaggedAsWaisPulser);
  std::cout << tc.N() << std::endl;

  double lastRealTime = 0;
  double lastCoherentFilteredSnr = 0;
  ProgressBar p(tc.N());

  const int nx = 50;
  const double xLow = 0;
  const double xUp = 10;  
  TH1D* hPassed = new TH1D("hPassed", "", nx, xLow, xUp);
  TH1D* hTotal = new TH1D("hTotal", "", nx, xLow, xUp);
  for(Long64_t entry=0; entry < tc.N(); entry++){
    tc.getEntry(entry);
    if(lastRealTime > 0){
      double deltaTime = tc.realTime - lastRealTime;
      // std::cout << deltaTime << "\t" << tc.coherent_filtered_snr << std::endl;
      // tc.coherent_filtered_snr;

      
      hPassed->Fill(tc.coherent_filtered_snr);
      hTotal->Fill(tc.coherent_filtered_snr);
      
      for(UInt_t dt=1; dt < deltaTime; dt++){
	hTotal->Fill(tc.coherent_filtered_snr);
      }
    }
    lastRealTime = tc.realTime;
    lastCoherentFilteredSnr = tc.coherent_filtered_snr;
    p.inc(entry);
  }

  TEfficiency* hEff = new TEfficiency(*hPassed, *hTotal);
  hEff->SetTitle("ANITA-3 WAIS pulser hardware efficiency; SNR (no units); Efficiency");
  hEff->Draw();
}
