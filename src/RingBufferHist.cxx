#include "RingBufferHist.h"
#include "TMath.h"
#include "TF1.h"
#include "TList.h"
#include "FourierBuffer.h"

#include "TMinuit.h"
#include "TPad.h"

ClassImp(Acclaim::RingBufferHist);

Acclaim::RingBufferHist::RingBufferHist(const char* name, const char* title, int nBins, double xMin, double xMax, int ringBufferSize) : 
  TH1D(name, title, 10, 0, 1), amplitudes(ringBufferSize), fNx(nBins), fNumEvents(0),  fNumNonEmptyBins(0) {

  binCentres.resize(fNx, 0);
  squaredBinCentres.resize(fNx, 0);
  binValues.resize(fNx, 0);
  squaredBinErrors.resize(fNx, 0);


  
}

Acclaim::RingBufferHist::~RingBufferHist(){
}



int Acclaim::RingBufferHist::Fill(double amp, double sign){

  // Here I force poisson errors for bin content *even if removing events*
  // this allows this histogram to be used for a rolling average

  sign = sign >= 0 ? 1 : -1;
  int bx = TH1D::Fill(amp, sign);
  double n = GetBinContent(bx);

  // if it currently equals 1 and we just filled it, 
  // then it must previously have been empty
  fNumNonEmptyBins += n == 1 ? 1 : 0;

  SetBinError(bx, TMath::Sqrt(n));
  if(bx > 0 && bx <= fNx){
    binValues[bx-1] = n;
    squaredBinErrors[bx-1] = binValues[bx-1];
  }
  fNumEvents += sign;
  return bx;
}


bool Acclaim::RingBufferHist::add(double newAmp){

  // first we remove the old value, should be zero if unused
  
  bool needRemoveOld = amplitudes.numElements() == amplitudes.size();

  double oldAmp = amplitudes.insert(newAmp);
  bool removedOld = false;
  if(needRemoveOld){
    Fill(oldAmp, -1); // remove old event from hist
    removedOld = true;    
  }
  Fill(newAmp);

  return removedOld;
}


