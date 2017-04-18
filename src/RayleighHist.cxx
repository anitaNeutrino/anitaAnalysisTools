#include "RayleighHist.h"
#include "TMath.h"
#include "TF1.h"
#include "TList.h"
#include "FourierBuffer.h"

#include "TMinuit.h"

ClassImp(Acclaim::RayleighHist);


Acclaim::RayleighHist::RayleighHist(FourierBuffer* fb,
				    const char* name, const char* title) : TH1D(name, title, 1, -1001, -1000), // deliberately stupid initial binning
									   fBinWidth(1),
									   fRayleighNorm(0),
									   amplitudes(fb ? fb->bufferSize : 1000),
									   fNumFitParams(1),
									   theFitParams(std::vector<double>(fNumFitParams, 0)),
									   theFitParamsSteps(std::vector<double>(fNumFitParams, 1e-3)),
									   fChiSquaredFunc(this, &Acclaim::RayleighHist::getRayleighChiSquare, fNumFitParams)
{ 
  fParent = fb;
  freqMHz = -9999;
  fracOfEventsWanted = 0.99;
  risingEdgeBins = 4;
  fitMethod = FitMethod::Default;
  
  
  fRay = fParent ? (TF1*) fParent->fRay->Clone(TString::Format("%s_fit", name)) : NULL;




  fMinimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  // set tolerance , etc...
  fMinimizer->SetMaxFunctionCalls(1000); // for Minuit/Minuit2 
  fMinimizer->SetMaxIterations(1000);  // for GSL
  // fMinimizer->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  // fMinimizer->SetMaxIterations(10000);  // for GSL
  
  fMinimizer->SetTolerance(0.001);
  fMinimizer->SetPrintLevel(0);

  // create funciton wrapper for minmizer
  // a IMultiGenFunction type
  
  fNumFitParams = 1;
  fMinimizer->SetFunction(fChiSquaredFunc);
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
  fBinWidth = sigmaGuess/risingEdgeBins;

  const int nBins = TMath::Nint(desiredMaxAmp/fBinWidth);
    
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
  fRayleighNorm = amplitudes.size()*fBinWidth;
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
  fRayleighNorm = fBinWidth*numEvents;
  return numEvents;
}



double Acclaim::RayleighHist::getRayleighChiSquare(const double* params){
  
  double amplitude = params[0];
  if(amplitude <= 0){ // does this fuck the fitter?
    return 9999;
  }
  double chiSquare=0;
  for(int i=0; i < 3*fNx; i+=3){
    double yR = EvalRayleigh(fRayleighNorm, amplitude, rayleighChiSquareCache[i]) - rayleighChiSquareCache[i+1];    
    chiSquare += rayleighChiSquareCache[i+2] > 0 ? yR*yR/rayleighChiSquareCache[i+2] : 0;
  }
  return chiSquare;
}


void Acclaim::RayleighHist::Fit(Double_t& rayleighAmplitude, Double_t& chiSquare, Int_t& ndf){


  rayleighAmplitude = amplitudes.getMean()*TMath::Sqrt(2./TMath::Pi());  
  ndf = 0;
  chiSquare=0;
  ndf = 0;
  
  // the number of degrees of freedom are the number of non-empty bins - the number of fit parameters (either 0 or 1)  
  switch(fitMethod){

  case FitMethod::TF1:
    {
      std::cerr << "Not currently implemented" << std::endl;
      return;
    }
  case FitMethod::Minuit:
    {
      double rayInitialGuess = rayleighAmplitude;
      prepareRayleighFitCache();
      fMinimizer->SetVariable(0, "amplitude", rayleighAmplitude, theFitParamsSteps.at(0));
      bool successfulMinimization = fMinimizer->Minimize();
      if(!successfulMinimization){
	std::cerr << fTitle << " failed the minimization, amplitude = " << theFitParams.at(0) << "\t" << ", intial guess was " << rayInitialGuess << std::endl << std::endl;
      }
      return;
    }
  case FitMethod::JustEvalGuess:
    {
      prepareRayleighFitCache();
      chiSquare = getRayleighChiSquare(&rayleighAmplitude);
      return;
    }
  case FitMethod::Scan:
  default:
    {
      prepareRayleighFitCache();

      // const int nStep = 2*rayleighAmplitude/theFitParamsSteps.at(0);
      const int nStep = 100; //2*rayleighAmplitude/theFitParamsSteps.at(0);
      const double ds = 2*rayleighAmplitude/nStep;
      double bestAmp = 0;
      double minChiSquare = 1e9;
      for(int s=1; s < nStep; s++){
	double amp = s*ds;
	double cs = getRayleighChiSquare(&amp);
	// std::cout << s << "\t" << amp << "\t" << cs << "\t" << nStep << "\t" << std::endl;

	if(cs < minChiSquare){
	  bestAmp = amp;
	  minChiSquare = cs;
	}	
      }
      // std::cout << bestAmp << "\t" << minChiSquare << "\t" << nStep << "\t" << std::endl << std::endl;
      rayleighAmplitude = bestAmp;
      chiSquare = minChiSquare;
      return;
    }
  }  
}


void Acclaim::RayleighHist::prepareRayleighFitCache(){

  // aggressive caching to try and speed up fitting function  
  const int nCacheVals = 3*GetNbinsX();
  rayleighChiSquareCache.resize(0);
  rayleighChiSquareCache.reserve(nCacheVals);
  fNx = GetNbinsX();
  for(int bx=1; bx <= fNx; bx++){
    double y = GetBinContent(bx);
    double x = GetXaxis()->GetBinCenter(bx);
    double yErr = GetBinError(bx);
    rayleighChiSquareCache.push_back(x);
    rayleighChiSquareCache.push_back(y);
    rayleighChiSquareCache.push_back(yErr*yErr);
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

