#include "RootTools.h"

ClassImp(RootTools)

RootTools::RootTools(){

}

RootTools::~RootTools(){

}


void RootTools::getMaxMin(TGraph* gr, Double_t& maxY, Double_t& minY){
  Double_t maxX=0;
  Double_t minX=0;
  RootTools::getMaxMin(gr, maxY, maxX, minY, minX);
}



void RootTools::getMaxMin(TGraph* gr, Double_t& maxY, Double_t& maxX, Double_t& minY, Double_t& minX){

  maxY=gr->GetY()[0];
  maxX=gr->GetX()[0];
  minY=gr->GetY()[0];
  minX=gr->GetX()[0];

  for(int i=0; i<gr->GetN(); i++){
    if(gr->GetY()[i] > maxY){
      maxY = gr->GetY()[i];
      maxX = gr->GetX()[i];
    }
    if(gr->GetY()[i] < minY){
      minY = gr->GetY()[i];
      minX = gr->GetX()[i];
    }
  } 
}


void RootTools::subtractOffset(TGraph* gr, Double_t offset){
  for(int i=0; i<gr->GetN(); i++){
    gr->GetY()[i] -= offset;
  }
}



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
