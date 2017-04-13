#include "RayleighHist.h"
#include "TMath.h"

ClassImp(Acclaim::RayleighHist);


Acclaim::RayleighHist::RayleighHist(const char* name, const char* title) : TH1D(name, title, 1, -1001, -1000){ // deliberately stupid initial binning
  init();
}


void Acclaim::RayleighHist::init(){
  
  fracOfEventsWanted = 0.99;
  risingEdgeBins = 4;
  maxOverFlowThresh = 0.75;
}



bool Acclaim::RayleighHist::axisRangeOK(double meanAmp) const{

  // assume I'm cutting all the serious crap before here (SURF saturation + self triggered blasts)
  // that means if we've got a lot of overflow we need to rebin.
  double overFlow = GetBinContent(GetNbinsX()+1);
  double integralWithOverUnderFlow = Integral() + GetBinContent(0) + overFlow;

  // or it could be that the axis limits are too large...
  // I guess, let's rebin if all the samples are below some fraction of the total
  // but I'll sort that later
  
  if(overFlow/integralWithOverUnderFlow > maxOverFlowThresh){
    return false;
  }
  else{
    return true;
  }
}


void Acclaim::RayleighHist::rebinAndEmptyHist(double meanAmp){

  const double fracOfEventsWanted = 0.99;
  const double desiredMaxAmp = meanAmp*TMath::Sqrt(-(4./TMath::Pi())*TMath::Log(1.0 - fracOfEventsWanted));
  const double sigmaGuess = TMath::Sqrt(2./TMath::Pi())*meanAmp;

  // get 4 bins from the rising edge to the peak
  const Int_t risingEdgeBins = 4;
  const double binWidth = sigmaGuess/risingEdgeBins;

  const int nBins = TMath::Nint(desiredMaxAmp/binWidth);
    
  SetBins(nBins, 0, desiredMaxAmp);

  // empty bins inc. overflow and underflow
  for(int bx = 0; bx <= GetNbinsX() + 1; bx++){
    SetBinContent(bx, 0);
    SetBinError(bx, 0);
  }
}


int Acclaim::RayleighHist::Fill(double amp, double sign){

  // Here I force poisson errors for bin content *even if removing events*
  // this allows this histogram to be used for a rolling average

  sign = sign >= 0 ? 1 : -1;

  int bx = TH1D::Fill(amp, sign);
  double n = GetBinContent(bx);
  SetBinError(bx, TMath::Sqrt(n));

  SetTitle(TString::Format("%d events", int(Integral() + GetBinContent(GetNbinsX()+1))));

  return bx;
}
