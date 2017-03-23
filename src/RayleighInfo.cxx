#include "RayleighInfo.h"

ClassImp(RayleighInfo)

RayleighInfo::RayleighInfo(){

  run = 0;
  eventNumber = 0;
  realTimeNs = 0;
  firstRun = 0;
  firstEventNumber = 0;
  firstRealTimeNs = 0;
  eventAmp = 0;
  rayFitNorm = 0;
  rayFitAmp = 0;
  rayFitAmpError = 0;
  rayGuessAmp = 0;
  rayChiSquare = 0;
  rayNdf = 0;
  nBins = 0;
  binWidth = 0;
  integralPlusOverUnderFlow = 0;

}

TH1D* RayleighInfo::makeHistogram(){
  return NULL;
}
