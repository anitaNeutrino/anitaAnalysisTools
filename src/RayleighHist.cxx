#include "RayleighHist.h"
#include "TMath.h"
#include "TF1.h"
#include "TList.h"
#include "FourierBuffer.h"

#include "TMinuit.h"
#include "TPad.h"

ClassImp(Acclaim::RayleighHist);

const double Acclaim::RayleighHist::dExponent = 0.001; // maybe
const double Acclaim::RayleighHist::minExponent = 10; // min is zero
std::vector<double> Acclaim::RayleighHist::exponentialCache;


Acclaim::RayleighHist::RayleighHist(FourierBuffer* fb,
				    const char* name, const char* title) 
  : TH1D(name, title, 1, -1001, -1000), // deliberately stupid initial binning
    fBinWidth(1), fFitEveryNAdds(10), fNumAddsMod10(0), fRayleighNorm(0), fNumEvents(0),
    amplitudes(fb ? fb->bufferSize : 1000), fNumFitParams(1), 
    theFitParams(std::vector<double>(fNumFitParams, 0)),
    theFitParamsSteps(std::vector<double>(fNumFitParams, 1e-3)),
    fChiSquaredFunc(this, &Acclaim::RayleighHist::getRayleighChiSquare, fNumFitParams)
{ 
  fParent = fb;
  freqMHz = -9999;
  fracOfEventsWanted = 0.99;
  risingEdgeBins = 4;

  fitMethod = FitMethod::Adaptive;
  // fitMethod = FitMethod::Scan;
  // fitMethod = FitMethod::JustEvalGuess;
  
  
  fRay = fParent ? (TF1*) fParent->fRay->Clone(TString::Format("%s_fit", name)) : NULL;

  fNx = GetNbinsX();
  binCentres.resize(fNx, 0);
  squaredBinCentres.resize(fNx, 0);
  binValues.resize(fNx, 0);
  squaredBinErrors.resize(fNx, 0);

  fNumNonEmptyBins = 0;

  fParamsTF1.resize(2, 0);

  fMinimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  // set tolerance , etc...
  fMinimizer->SetMaxFunctionCalls(1000); // for Minuit/Minuit2 
  fMinimizer->SetMaxIterations(1000);  // for GSL  
  fMinimizer->SetTolerance(0.001);
  fMinimizer->SetPrintLevel(0);

  // create funciton wrapper for minmizer
  // a IMultiGenFunction type
  
  fNumFitParams = 1;
  fMinimizer->SetFunction(fChiSquaredFunc);

  
  // someone has to initialize 
  if(exponentialCache.size()==0){
    unsigned size = TMath::Nint(TMath::Abs(minExponent)/dExponent);
    exponentialCache.reserve(size);
    double minExponentCalculated = 0;
    while(minExponentCalculated > minExponent){
      exponentialCache.push_back(exp(minExponentCalculated));
      minExponentCalculated -= dExponent;
    }
  }
}









Acclaim::RayleighHist::~RayleighHist(){

  if(fRay){
    delete fRay;
    fRay = NULL;
  }

  if(fMinimizer){
    delete fMinimizer;
    fMinimizer = NULL;
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


void Acclaim::RayleighHist::rebinAndRefill(double meanAmp, double sigmaGuess){

  const double desiredMaxAmp = meanAmp*TMath::Sqrt(-(4./TMath::Pi())*TMath::Log(1.0 - fracOfEventsWanted));

  // get 4 bins from the rising edge to the peak
  const Int_t risingEdgeBins = 4;
  fBinWidth = sigmaGuess/risingEdgeBins;

  fNx = TMath::Nint(desiredMaxAmp/fBinWidth);
  
  binCentres.resize(fNx, 0);
  squaredBinCentres.resize(fNx, 0);
  binValues.resize(fNx, 0);
  squaredBinErrors.resize(fNx, 0);

  SetBins(fNx, 0, desiredMaxAmp);
  
  // empty bins inc. overflow and underflow
  for(int bx = 0; bx <= fNx + 1; bx++){
    SetBinContent(bx, 0);
  }

  // set bin content by hand...
  for(RingBuffer::iterator it=amplitudes.begin(); it != amplitudes.end(); ++it){
    double amp = (*it);
    int bx = TH1D::FindBin(amp);
    SetBinContent(bx, GetBinContent(bx)+1);
  }

  // ... set errors once at the end
  fNumNonEmptyBins = 0;
  for(int bx=1; bx <= GetNbinsX(); bx++){
    double y = GetBinContent(bx);
    SetBinError(bx, TMath::Sqrt(y));

    // update cached values
    binCentres[bx-1] = fBinWidth*0.5 + (bx-1)*fBinWidth;
    squaredBinCentres[bx-1] = binCentres[bx-1]*binCentres[bx-1];
    binValues[bx-1] = y;
    squaredBinErrors[bx-1] = binValues[bx-1];
    fNumNonEmptyBins += y > 0 ? 1 : 0;

  }
  fRayleighNorm = amplitudes.size()*fBinWidth;
  
}


int Acclaim::RayleighHist::Fill(double amp, double sign){

  // Here I force poisson errors for bin content *even if removing events*
  // this allows this histogram to be used for a rolling average

  sign = sign >= 0 ? 1 : -1;
  int bx = TH1D::Fill(amp, sign);
  double n = GetBinContent(bx);
  SetBinError(bx, TMath::Sqrt(n));
  if(bx > 0 && bx <= fNx){
    binValues[bx-1] = n;
    squaredBinErrors[bx-1] = binValues[bx-1];
  }
  fNumEvents += sign;
  return bx;
}


bool Acclaim::RayleighHist::add(double newAmp){

  // first we remove the old value, should be zero if unused
  
  bool needRemoveOld = amplitudes.numElements() == amplitudes.size();

  double oldAmp = amplitudes.insert(newAmp);
  if(needRemoveOld){
    Fill(oldAmp, -1); // remove old event from hist
  }
  Fill(newAmp);

  // fNumEvents = amplitudes.numElements();

  double mean = amplitudes.getMean();
  fRayleighAmpGuess = TMath::Sqrt(2./TMath::Pi())*mean;
  if(!axisRangeOK()){
    rebinAndRefill(mean, fRayleighAmpGuess);
  }
  fRayleighNorm = fBinWidth*fNumEvents; // fNumEvents is updated in Fill()

  bool updatedFit = false;
  fNumAddsMod10++;
  if(fFitEveryNAdds > 0 && fNumAddsMod10 == fFitEveryNAdds){
    fitRayleigh(false);
    fNumAddsMod10 = 0;
    updatedFit = true;
  }


  return updatedFit;
}




void Acclaim::RayleighHist::getRayleighFitParams(double& rayAmp, double& chiSquare, int& ndf){
  rayAmp = fRayleighAmplitude;
  chiSquare = fChiSquare;
  ndf = fNDF;
}



double Acclaim::RayleighHist::getRayleighChiSquare(const double* params){
  
  double amplitude = params[0];
  double chiSquare=0;
  if(amplitude <= 0){ // does this fuck the fitter?
    chiSquare = 9999 + fabs(amplitude);
  }
  else{
    double invAmpSquare = 1./(amplitude*amplitude);
    double minusHalfInvAmpSquare = -0.5*invAmpSquare;
    double exponentPreFactorPart = fRayleighNorm*invAmpSquare;
    // std::cout << amplitude << "\t" << invAmpSquare << "\t" << exponentPreFactorPart << "\t" << fRayleighNorm << "\t" << fBinWidth << "\t" << fNumEvents << std::endl;
    for(int i=0; i < fNx; i++){
      if(squaredBinErrors[i] > 0){
	double yR = exponentPreFactorPart*binCentres[i]*exp(minusHalfInvAmpSquare*squaredBinCentres[i]) - binValues[i];
        chiSquare += yR*yR/squaredBinErrors[i];
      }
    }
  }
  return chiSquare;
}





void Acclaim::RayleighHist::fitRayleighScan(bool forGuiUpdateTF1){
  fNDF = fNumNonEmptyBins - 1; // -1 since we are varying the amplitude parameter

  // const int nStep = 2*rayleighAmplitude/theFitParamsSteps.at(0);
  const int nStep = 100; //2*rayleighAmplitude/theFitParamsSteps.at(0);
  double ds = 2*fRayleighAmpGuess/nStep;
  double bestAmp = 0;
  double minChiSquare = 1e9;
  int bestS = 0;
  for(int s=1; s < nStep; s++){
    double amp = s*ds;
    double cs = getRayleighChiSquare(&amp);
    // std::cout << s << "\t" << amp << "\t" << cs << "\t" << nStep << "\t" << ds << std::endl;

    if(cs < minChiSquare){
      bestAmp = amp;
      minChiSquare = cs;
      bestS = s;
    }
  }
  // if(bestS==1){
  //   std::cout << fName << "\t" << bestAmp << "\t" << minChiSquare << "\t" << fRayleighAmpGuess << std::endl;
  // }
  // if(fNumEvents){
  //   std::cout << bestAmp << "\t" << minChiSquare << "\t" << nStep << "\t" << fNumEvents << std::endl << std::endl;
  // }
  fRayleighAmplitude = bestAmp;
  fChiSquare = minChiSquare;  

  if(forGuiUpdateTF1){
    updateParamsOfTF1();
  }
}


void Acclaim::RayleighHist::fitRayleighAdaptive(bool forGuiUpdateTF1){
  // prepareRayleighFitCache();
  fNDF = fNumNonEmptyBins; // no varying parameters
  fChiSquare = getRayleighChiSquare(&fRayleighAmplitude);
  
  if(fChiSquare / fNDF >= 2){
    fitRayleighScan(forGuiUpdateTF1);
  }
  else{
    if(forGuiUpdateTF1){
      updateParamsOfTF1();
    }
  }
}


void Acclaim::RayleighHist::fitRayleighJustEvalGuess(bool forGuiUpdateTF1){
  // prepareRayleighFitCache();
  fNDF = fNumNonEmptyBins; // no varying parameters
  fChiSquare = getRayleighChiSquare(&fRayleighAmplitude);
  if(forGuiUpdateTF1){
    updateParamsOfTF1();
  }

}

void Acclaim::RayleighHist::fitRayleighTF1(){
  fRay->FixParameter(0, fRayleighNorm);
  fRay->SetParameter(1, fRayleighAmpGuess);
  fRay->SetParLimits(1, 0, fRayleighAmpGuess*2);
  this->Fit(fRay, "Q0");
  fChiSquare = fRay->GetChisquare();
  fNDF = fRay->GetNDF();
  fRayleighAmplitude = fRay->GetParameter(1);
}

void Acclaim::RayleighHist::fitRayleighMinuit(bool forGuiUpdateTF1){
  theFitParams.at(0) = fRayleighAmpGuess;
  fMinimizer->SetVariable(0, "amplitude", fRayleighAmpGuess, theFitParamsSteps.at(0));
  bool successfulMinimization = fMinimizer->Minimize();
  if(!successfulMinimization){
    std::cerr << fTitle << " failed the minimization, amplitude = " << theFitParams.at(0) 
	      << "\t" << ", intial guess was " << fRayleighAmpGuess << std::endl << std::endl;
  }
  fRayleighAmplitude = theFitParams.at(0); // this isn't quite right, but I'm not going to fix it now...
  // std::cout << fRayleighAmplitude << "\t" << fMinimizer->MinValue() << std::endl;
  fChiSquare = getRayleighChiSquare(&fRayleighAmplitude);
  fNDF = fNumNonEmptyBins - 1;
  if(forGuiUpdateTF1){
    updateParamsOfTF1();
  }

}

void Acclaim::RayleighHist::fitRayleigh(bool forGuiUpdateTF1){

  // the internals that need to be updated by the fit
  fRayleighAmplitude = fRayleighAmpGuess; // sensible initial guess calculated elsewhere
  fNDF = 0;
  fChiSquare = 0;

  if(fitMethod==FitMethod::Adaptive){
    fitRayleighAdaptive();
  }
  else if(fitMethod==FitMethod::Scan){
    fitRayleighScan(forGuiUpdateTF1);
  }
  else if(fitMethod==FitMethod::JustEvalGuess){
    fitRayleighJustEvalGuess(forGuiUpdateTF1);
  }
  else if(fitMethod==FitMethod::TF1){
    fitRayleighTF1();
  }
  else if(fitMethod==FitMethod::Minuit){
    fitRayleighMinuit(forGuiUpdateTF1);
  }
  else{
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ 
	      << ", you asked for an unknown fit method so I'm not doing anything." << std::endl;
  }    
}


void Acclaim::RayleighHist::updateParamsOfTF1(){
  if(fRay){
    if(fParamsTF1.at(0) != fRayleighNorm){
      fParamsTF1.at(0) = fRayleighNorm;
      fRay->SetParameter(0, fParamsTF1.at(0));
    }
    if(fParamsTF1.at(1) != fRayleighAmplitude){
      fParamsTF1.at(1) = fRayleighAmplitude;
      fRay->SetParameter(1, fParamsTF1.at(1));
    }
  }
}


void Acclaim::RayleighHist::Draw(Option_t* opt){

  TH1D::Draw(opt);  

  if(fRay){
    updateParamsOfTF1();
    fRay->Draw("lsame");
  }  
}

