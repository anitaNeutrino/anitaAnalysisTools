#include "RayleighHist.h"
#include "TMath.h"
#include "TF1.h"
#include "TList.h"
#include "FourierBuffer.h"

#include "TMinuit.h"
#include "TPad.h"

ClassImp(Acclaim::RayleighHist);

#define NUM_BINS 10


Acclaim::RayleighHist::RayleighHist(FourierBuffer* fb,
				    const char* name, const char* title) : 
    TH1D(name, title, 10, 0, 1), amplitudes(fb ? fb->bufferSize : 1000), fNx(10), fNumEvents(0),  fNumNonEmptyBins(0),  
    fNumFitParams(1),
    theFitParams(std::vector<double>(fNumFitParams, 0)),
    theFitParamsSteps(std::vector<double>(fNumFitParams, 1e-3)),
    fChiSquaredFunc(this, &Acclaim::RayleighHist::getRayleighChiSquare, fNumFitParams),
    fBinWidth(1), fFitEveryNAdds(10), fNumAddsMod10(fFitEveryNAdds), fRayleighNorm(0){

  binValues.resize(fNx, 0);
  squaredBinErrors.resize(fNx, 0),
  
      fParent = fb;
  freqMHz = -9999;
  // fracOfEventsWanted = (1 - 1e-8);
  fracOfEventsWanted = 0.999;  

  fitMethod = kAdaptive;
  // chiSquareErrorMethod = kPearson;
  chiSquareErrorMethod = kPoisson;
  // fitMethod = Scan;
  // fitMethod = JustEvalGuess;
  
  
  fRay = fParent ? (TF1*) fParent->fRay->Clone(TString::Format("%s_fit", name)) : NULL;
  fParamsTF1.resize(2, 0);


  grLastAddedAmp = new TGraph(2);
  grLastAddedAmp->SetName("grLastAddedAmplitude");
  grLastAddedAmp->SetTitle("Last added amplitude");
  
  fMinimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  // set tolerance , etc...
  fMinimizer->SetMaxFunctionCalls(1000); // for Minuit/Minuit2 
  fMinimizer->SetMaxIterations(1000);  // for GSL  
  fMinimizer->SetTolerance(0.001);
  fMinimizer->SetPrintLevel(0);

  // create funciton wrapper for minmizer
  // a IMultiGenFunction type
  
  fMinimizer->SetFunction(fChiSquaredFunc);  
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

  if(grLastAddedAmp){
    delete grLastAddedAmp;
    grLastAddedAmp = NULL;
  }
}


void Acclaim::RayleighHist::guessMaxBinLimitAndSigmaFromMean(double meanAmp, double& upperLimit, double& sigmaGuess, double fracOfEventsInsideMaxAmp){
  upperLimit = meanAmp*TMath::Sqrt(-(4./TMath::Pi())*TMath::Log(1.0 - fracOfEventsInsideMaxAmp));
  sigmaGuess = TMath::Sqrt(2./TMath::Pi())*meanAmp;  
}


bool Acclaim::RayleighHist::axisRangeOK() const{

  double mean = amplitudes.getMean();
  double maxAmp, sigmaGuess;
  guessMaxBinLimitAndSigmaFromMean(mean, maxAmp, sigmaGuess, fracOfEventsWanted);

  // let's say, within 10% is OK?
  // double xup = GetBinLowEdge(NUM_BINS+1);
  double xup = fBinWidth*(NUM_BINS+1);
  if(xup/maxAmp > 0.9 && xup/maxAmp < 1.1){
    return true;
  }
  else{
    return false;
  }
}


void Acclaim::RayleighHist::rebinAndRefill(double meanAmp){

  const double desiredMaxAmp = meanAmp*TMath::Sqrt(-(4./TMath::Pi())*TMath::Log(1.0 - fracOfEventsWanted));

  // // get 4 bins from the rising edge to the peak
  // const Int_t risingEdgeBins = 4;
  // fBinWidth = sigmaGuess/risingEdgeBins;

  fBinWidth = desiredMaxAmp/NUM_BINS;
  // NUM_BINS = TMath::Nint(desiredMaxAmp/fBinWidth);
  // NUM_BINS = TMath::Nint(desiredMaxAmp/fBinWidth);
  
  binCentres.resize(NUM_BINS, 0);
  squaredBinCentres.resize(NUM_BINS, 0);
  binValues.resize(NUM_BINS, 0);
  squaredBinErrors.resize(NUM_BINS, 0);

  SetBins(NUM_BINS, 0, desiredMaxAmp);
  
  // empty bins inc. overflow and underflow
  for(int bx = 0; bx <= NUM_BINS + 1; bx++){
    SetBinContent(bx, 0);    
  }

  for(int bin=0; bin < NUM_BINS; bin++){
    binValues[bin] = 0;
  }

  // set bin content by hand...
  for(RingBuffer::iterator it=amplitudes.begin(); it != amplitudes.end(); ++it){
    double amp = (*it);
    if(amp < desiredMaxAmp){
      int bin = floor(amp/fBinWidth);
      binValues[bin]++;
    }
  }

  // // set bin content by hand...
  // for(RingBuffer::iterator it=amplitudes.begin(); it != amplitudes.end(); ++it){
  //   double amp = (*it);
  //   int bx = TH1D::FindBin(amp);
  //   SetBinContent(bx, GetBinContent(bx)+1);
  // }  

  // ... set errors once at the end
  fNumNonEmptyBins = 0;
  for(int bx=1; bx <= NUM_BINS; bx++){
    // double y = GetBinContent(bx);
    double y = binValues[bx-1];
    
    SetBinContent(bx, y);
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

  // if it currently equals 1 and we just filled it, 
  // then it must previously have been empty
  fNumNonEmptyBins += n == 1 ? 1 : 0;

  SetBinError(bx, TMath::Sqrt(n));
  if(bx > 0 && bx <= NUM_BINS){
    binValues[bx-1] = n;
    squaredBinErrors[bx-1] = binValues[bx-1];
  }
  fNumEvents += sign;
  return bx;
}


bool Acclaim::RayleighHist::add(double newAmp){

  // first we remove the old value, should be zero if unused
  
  // first we remove the old value, should be zero if unused
  
  bool needRemoveOld = amplitudes.numElements() == amplitudes.size();
  double oldAmp = amplitudes.insert(newAmp);
  bool removedOld = false;
  if(needRemoveOld){
    Fill(oldAmp, -1); // remove old event from hist
    removedOld = true;    
  }
  Fill(newAmp);

  // fNumEvents = amplitudes.numElements();

  double mean = amplitudes.getMean();
  fRayleighAmpGuess = TMath::Sqrt(2./TMath::Pi())*mean;
  if(!axisRangeOK()){
    rebinAndRefill(mean);
  }
  fRayleighNorm = fBinWidth*fNumEvents; // fNumEvents is updated in Fill()

  bool updatedFit = false;
  if(fNumAddsMod10 == fFitEveryNAdds){
    fitRayleigh(false);
    fNumAddsMod10 = 0;
    updatedFit = true;
  }
  fNumAddsMod10++;  

  grLastAddedAmp->SetPoint(0, newAmp, 0);
  grLastAddedAmp->SetPoint(1, newAmp, 100000000);

  
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
    for(int i=0; i < NUM_BINS; i++){
      if(squaredBinErrors[i] > 0){
	double yTheoretical = exponentPreFactorPart*binCentres[i]*exp(minusHalfInvAmpSquare*squaredBinCentres[i]);
        double yMeasured = binValues[i];
        double dy = yTheoretical - yMeasured;
        if(chiSquareErrorMethod==kPoisson){
          chiSquare += dy*dy/squaredBinErrors[i];
        }
        else if(chiSquareErrorMethod==kPearson){
          chiSquare += dy*dy/yTheoretical;
        }
        else{
          std::cerr << "Error in " << __PRETTY_FUNCTION__
                    << ". Unknown chiSquareErrorMethod" << std::endl;
        }
      }
    }
  }
  return chiSquare;
}




void Acclaim::RayleighHist::fitRayleighAdaptive(bool forGuiUpdateTF1){
  // prepareRayleighFitCache();
  fitRayleighJustEvalGuess(forGuiUpdateTF1);
  
  if(fChiSquare / fNDF >= 2){
    fitRayleighScan(forGuiUpdateTF1);
  }
  else{
    if(forGuiUpdateTF1){
      updateParamsOfTF1();
      std::cout << __PRETTY_FUNCTION__ << "\t" << fChiSquare << "\t" << fNDF << "\t" << fNumNonEmptyBins << "\t" << std::endl;
    }
  }
}



void Acclaim::RayleighHist::fitRayleighScan(bool forGuiUpdateTF1){
  fNDF = fNumNonEmptyBins - 1; // -1 since we are varying the amplitude parameter

  // const int nStep = 2*rayleighAmplitude/theFitParamsSteps.at(0);
  const int nStep = 100; //2*rayleighAmplitude/theFitParamsSteps.at(0);
  double ds = 2*fRayleighAmpGuess/nStep;
  double bestAmp = 0;
  double minChiSquare = 1e9;
  for(int s=1; s < nStep; s++){
    double amp = s*ds;
    double cs = getRayleighChiSquare(&amp);
    // std::cout << s << "\t" << amp << "\t" << cs << "\t" << nStep << "\t" << ds << std::endl;

    if(cs < minChiSquare){
      bestAmp = amp;
      minChiSquare = cs;
    }
  }
  fRayleighAmplitude = bestAmp;
  fChiSquare = minChiSquare;  

  if(forGuiUpdateTF1){
    updateParamsOfTF1();
    std::cout << __PRETTY_FUNCTION__ << "\t" << fChiSquare << "\t" << fNDF << "\t" << fNumNonEmptyBins << "\t" << std::endl;
  }
}


void Acclaim::RayleighHist::fitRayleighJustEvalGuess(bool forGuiUpdateTF1){
  // prepareRayleighFitCache();
  fNDF = fNumNonEmptyBins; // no varying parameters
  fChiSquare = getRayleighChiSquare(&fRayleighAmplitude);
  if(forGuiUpdateTF1){
    updateParamsOfTF1();
    std::cout << __PRETTY_FUNCTION__ << "\t" << fChiSquare << "\t" << fNDF << "\t" << fNumNonEmptyBins << "\t" << std::endl;
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
  std::cout << __PRETTY_FUNCTION__ << "\t" << fChiSquare << "\t" << fNDF << "\t" << fNumNonEmptyBins << "\t" << std::endl;
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
    std::cout << __PRETTY_FUNCTION__ << "\t" << fChiSquare << "\t" << fNDF << "\t" << fNumNonEmptyBins << "\t" << std::endl;
  }
}



void Acclaim::RayleighHist::fitRayleigh(bool forGuiUpdateTF1){

  // the internals that need to be updated by the fit
  fRayleighAmplitude = fRayleighAmpGuess; // sensible initial guess calculated elsewhere
  fNDF = 0;
  fChiSquare = 0;

  if(fitMethod==kAdaptive){
    fitRayleighAdaptive(forGuiUpdateTF1);
  }
  else if(fitMethod==kScan){
    fitRayleighScan(forGuiUpdateTF1);
  }
  else if(fitMethod==kJustEvalGuess){
    fitRayleighJustEvalGuess(forGuiUpdateTF1);
  }
  else if(fitMethod==kTF1){
    fitRayleighTF1();
  }
  else if(fitMethod==kMinuit){
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
  if(grLastAddedAmp){
    grLastAddedAmp->SetLineColor(kMagenta);
    grLastAddedAmp->Draw("lsame");
    // std::cout << grLastAddedAmp->GetX()[0] << "\t" 
    // 	      << fRayleighAmplitude << "\t" 
    // 	      << getCDF(grLastAddedAmp->GetX()[0]) << std::endl;
    // for(int i=0; i <= 10; i++){
    //   std::cout << getCDF(double(i)*fRayleighAmplitude) << std::endl;
    // }
  }

}

