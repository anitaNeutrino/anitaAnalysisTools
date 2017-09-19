#include "NoiseMonitor.h"
#include "TStyle.h"

void plotRmsPlots(int doSave=0, UInt_t hash=271896405){

  gStyle->SetPalette(kRainBow);
  gStyle->SetNumberContours(100);
  
  AnitaVersion::set(3);
  NoiseMonitor* nm = new NoiseMonitor(hash);

  

  // for(int run=130; run <= 131; run++){
  for(int run=130; run <= 439; run++){    
    if(run >= 257 && run <= 263){
      continue;
    }
    TString canName = TString::Format("can%d", run);    
    auto c1 = new TCanvas(canName, canName, 1600, 800);

    c1->Divide(AnitaPol::kNotAPol);
    
    for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
      c1->cd(polInd+1);
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
      TProfile2D* p = (TProfile2D* ) nm->getProfile(pol, run);
      p->GetYaxis()->SetTimeDisplay(1);
      p->SetMaximum(30);
      p->SetMinimum(0);
      if(p){
        p->Draw("colz");
      }
      c1->Modified();
      c1->Update();
    }
    if(doSave){
      c1->SaveAs(TString::Format("run%d_anitaRMS.gif", run));
      delete c1;      
    }
  }
// run 130 no
// run 131 no
// run 132 no
// run 142 no
// run 143 no
// run 144 no
// run 146 no
// run 346 no

// 130 131 132 142 143 144 146 346

}
