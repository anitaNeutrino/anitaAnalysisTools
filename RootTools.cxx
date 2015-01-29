#include "RootTools.h"

ClassImp(RootTools)

void RootTools::normalize(TGraph* gr, Double_t& mean, Double_t& rms){
  RootTools::getMeanAndRms(gr, mean, rms);
  for(int i=0; i<gr->GetN(); i++){    
    gr->GetY()[i] -= mean;
    gr->GetY()[i] /= rms;
  }
}

void RootTools::getMeanAndRms(TGraph* gr, Double_t& mean, Double_t& rms){
  Double_t sum = 0;
  Double_t square = 0;
  for(int i=0; i<gr->GetN(); i++){
    sum += gr->GetY()[i];
    square += gr->GetY()[i]*gr->GetY()[i];
  }
  mean = sum/gr->GetN();
  rms = TMath::Sqrt(square/gr->GetN() - mean*mean);
}
