#include "RayleighHist.h"
#include "TMath.h"
#include "TF1.h"
#include "TList.h"
#include "FourierBuffer.h"

ClassImp(Acclaim::RayleighHist);


Acclaim::RayleighHist::RayleighHist(FourierBuffer* fb,
				    const char* name, const char* title) : TH1D(name, title, 1, -1001, -1000), // deliberately stupid initial binning
									   amplitudes(fb ? fb->bufferSize : 1000){ 
  fParent = fb;
  freqMHz = -9999;
  fracOfEventsWanted = 0.99;
  risingEdgeBins = 4;
  maxOverFlowThresh = 0.9;
  
  
  
  fRay = fParent ? (TF1*) fParent->fRay->Clone(TString::Format("%s_fit", name)) : NULL;
}

Acclaim::RayleighHist::~RayleighHist(){

  if(fRay){
    delete fRay;
    fRay = NULL;
  }
}


void Acclaim::RayleighHist::guessMaxBinLimitAndSigmaFromMean(double meanAmp, double& upperLimit, double& sigmaGuess, double fracOfEventsInsideMaxAmp){
  upperLimit = meanAmp*TMath::Sqrt(-(4./TMath::Pi())*TMath::Log(1.0 - fracOfEventsInsideMaxAmp));
  sigmaGuess = TMath::Sqrt(2./TMath::Pi())*meanAmp;
}


bool Acclaim::RayleighHist::axisRangeOK() const{

  const double fracOfEventsWanted = 0.99;
  double mean = amplitudes.getMean();
  double maxAmp, sigmaGuess;
  guessMaxBinLimitAndSigmaFromMean(mean, maxAmp, sigmaGuess, fracOfEventsWanted);

  // let's say, within 10% is OK?
  double xup = GetBinLowEdge(GetNbinsX()+1);
  if(xup/maxAmp > 0.9 && xup/maxAmp < 1.1){
    return true;
  }
  else{
    return false;
  }
  
}


void Acclaim::RayleighHist::rebinAndRefill(double meanAmp){

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
  }

  // set errors once at the end
  for(RingBuffer::iterator it=amplitudes.begin(); it != amplitudes.end(); ++it){
    double amp = (*it);
    int bx = TH1D::FindBin(amp);
    SetBinContent(bx, GetBinContent(bx)+1);
  }
  // set errors once at the end
  for(int bx=1; bx <= GetNbinsX(); bx++){
    SetBinError(bx, TMath::Sqrt(GetBinContent(bx)));
  }
  
}


int Acclaim::RayleighHist::Fill(double amp, double sign){

  // Here I force poisson errors for bin content *even if removing events*
  // this allows this histogram to be used for a rolling average

  sign = sign >= 0 ? 1 : -1;
  int bx = TH1D::Fill(amp, sign);
  double n = GetBinContent(bx);
  SetBinError(bx, TMath::Sqrt(n));

  return bx;
}


int Acclaim::RayleighHist::add(double newAmp){

  // first we remove the old value, should be zero if unused
  
  bool needRemoveOld = amplitudes.numElements() == amplitudes.size();
  
  double oldAmp = amplitudes.insert(newAmp);

  if(needRemoveOld){
    Fill(oldAmp, -1); // remove old event from hist
  }
  Fill(newAmp);

  numEvents = amplitudes.numElements();  

  double mean = amplitudes.getMean();
  if(!axisRangeOK()){
    rebinAndRefill(mean);
  }

  return numEvents;
}




void Acclaim::RayleighHist::getTF1Params(Double_t& normalization, Double_t& rayAmp){
  if(fRay){
    normalization = fRay->GetParameter(0);
    rayAmp = fRay->GetParameter(1);
  }
  else{
    normalization = -999;
    rayAmp = -999;
  }
}


  
void Acclaim::RayleighHist::Eval(Double_t& chiSquare, Int_t& ndf){  

  // double integralWithOverFlow = Integral(0, GetNbinsX()+1);

    
  double binWidth = GetBinWidth(1);
  double rayFitNorm = binWidth*numEvents;
  double rayGuessAmp = amplitudes.getMean()*TMath::Sqrt(2./TMath::Pi());

  Double_t xMin = GetBinLowEdge(1);
  Double_t xMax = GetBinLowEdge(GetNbinsX()+1);
  fRay->SetRange(xMin, xMax);
  fRay->FixParameter(0, rayFitNorm);
  fRay->FixParameter(1, rayGuessAmp);
  // fRay->SetParLimits(1, rayGuessAmp, rayGuessAmp); // fixed width
  // fRay->SetParLimits(1, 0, 3*rayGuessAmp); // fixed width
  // fRay->SetParLimits(1, 0, 3*rayGuessAmp); // fixed width        

  // TH1::Fit(fRay, "Q0");

  
  chiSquare=0;
  ndf = 0;
  for(int bx=1; bx <= GetNbinsX(); bx++){
    double y = GetBinContent(bx);

    if(y > 0){
      double x = GetXaxis()->GetBinCenter(bx);
      double yR = fRay->Eval(x);
      double yErr = GetBinError(bx);
      double deltaY = (y - yR)/yErr;
      chiSquare += deltaY*deltaY;
      ndf++;
    }
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

