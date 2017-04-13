#include "RayleighHist.h"
#include "TMath.h"
#include "TF1.h"
#include "TList.h"
#include "FourierBuffer.h"

ClassImp(Acclaim::RayleighHist);


Acclaim::RayleighHist::RayleighHist(FourierBuffer* fb,
				    const char* name, const char* title) : TH1D(name, title, 1, -1001, -1000){ // deliberately stupid initial binning
  fParent = fb;
  freqMHz = -9999;
  fracOfEventsWanted = 0.99;
  risingEdgeBins = 4;
  maxOverFlowThresh = 0.9;
  fRay = (TF1*) fParent->fRay->Clone(TString::Format("%s_fit", name));
  // fRay = fParent->fRay;
}

Acclaim::RayleighHist::~RayleighHist(){

  if(fRay){
    delete fRay;
    fRay = NULL;
  }
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

  SetTitle(TString::Format("%4.1lf MHz - %d events in distribution", freqMHz, int(Integral() + GetBinContent(GetNbinsX()+1))));

  return bx;
}



void Acclaim::RayleighHist::Fit(){



  double integralWithOverFlow = Integral(0, GetNbinsX()+1);

  if(integralWithOverFlow >= 50){

    double sum=0;
    for(int bx=0; bx <= GetNbinsX(); bx++){
      sum += GetXaxis()->GetBinCenter(bx)*GetBinContent(bx);
    }
    double mean = sum/Integral();
    
    double binWidth = GetBinWidth(1);
    double rayFitNorm = integralWithOverFlow*binWidth;
    double rayGuessAmp = mean*TMath::Sqrt(2./TMath::Pi());

    Double_t xMin = GetBinLowEdge(1);
    Double_t xMax = GetBinLowEdge(GetNbinsX()+1);
    fRay->SetRange(xMin, xMax);
    // Set fit params
    fRay->FixParameter(0, rayFitNorm);
    // fRay->FixParameter(1, rayGuessAmp);
    fRay->SetParameter(1, rayGuessAmp);
    // fRay->SetParLimits(1, 0, 1e9); // essentially infinite  
    fRay->SetParLimits(1, rayGuessAmp, rayGuessAmp); // fixed width
    
    TH1::Fit(fRay, "Q0");

    if(fName.Contains(TString::Format("%d", fParent->fDrawFreqBin))){
      std::cout << fRay->GetName() << ":\t" << "rayGuessAmp = " << rayGuessAmp << "\tmean = " << mean << "\tintegral = " << Integral() << "\tintegral+of = " << integralWithOverFlow << "\tbinWidth = " << binWidth << "\t func eval at mean = " << fRay->Eval(mean) << std::endl;
    }
    
    // if(GetListOfFunctions()->GetEntries() > 0){      
    //   std::cout << fRay << ":\t" << "\t" << GetListOfFunctions()->At(0) << "\t" << GetListOfFunctions()->At(0)->GetName() << std::endl;
    //   GetListOfFunctions()->Print();
    // }      
  }  
}



void Acclaim::RayleighHist::Draw(Option_t* opt){

  TH1D::Draw(opt);  
  if(fRay){
    fRay->Draw("lsame");
  }
  
}

void Acclaim::RayleighHist::SetFreqBinToDraw(Int_t freqBin){
  fParent->fDrawFreqBin = freqBin;
}
