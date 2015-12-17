#include "AntarcticaMapPlotter.h"

void testAntarcticaMapPlotterInRoot(){

  AntarcticaMapPlotter* amp = new AntarcticaMapPlotter("hTest", "Testing", 128, 128);
  

  amp->Fill(90, 0);

  amp->Fill(85, 90);
  amp->Fill(80, 90);
  amp->Fill(75, 90);
  amp->Fill(70, 90);
  amp->Fill(65, 90);
  amp->Fill(60, 90);

  amp->Fill(85, 180);
  amp->Fill(80, 180);
  amp->Fill(75, 180);
  amp->Fill(70, 180);
  amp->Fill(65, 180);
  amp->Fill(60, 180);

  amp->Fill(85, 0);
  amp->Fill(80, 0);
  amp->Fill(75, 0);
  amp->Fill(70, 0);
  amp->Fill(65, 0);
  amp->Fill(60, 0);


  amp->Fill(85, 270);
  amp->Fill(80, 270);
  amp->Fill(75, 270);
  amp->Fill(70, 270);
  amp->Fill(65, 270);
  amp->Fill(60, 270);
  
  amp->DrawHist("colz");

  std::vector<Double_t> lats;
  std::vector<Double_t> longs;
  Double_t latVal = -90;
  Double_t longVal=0;
  
  Int_t numPoints = 0;
  while(latVal < -65){

    lats.push_back(latVal);
    longs.push_back(longVal);

    latVal += 5./360;
    longVal += 1;
    
    numPoints++;
  }
  
  amp->addTGraph("grTest", "More Testing", numPoints, &lats[0], &longs[0]);
  amp->getCurrentTGraph()->SetMarkerStyle(8);
  amp->getCurrentTGraph()->SetMarkerSize(0.2);
  amp->getCurrentTGraph()->SetMarkerColor(2);    
  amp->DrawTGraph("psame");
  
  
  
  
}
