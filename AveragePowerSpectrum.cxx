#include "AveragePowerSpectrum.h"

ClassImp(AveragePowerSpectrum);



AveragePowerSpectrum::AveragePowerSpectrum(){
  // Don't use this.
  for(Int_t freqInd=0; freqInd<NUM_FREQS; freqInd++){  
    hRayleighs[freqInd] = NULL;
    hRayleighFits[freqInd] = NULL;
  }

}


AveragePowerSpectrum::AveragePowerSpectrum(TString name, TString title){


  deltaFMHz = 1e3/(NOMINAL_SAMPLING_DELTAT*NUM_SAMPLES);

  count=0;
  SetName(name);
  SetTitle(title);  

  // Default power spec histogram values.
  // maxNumOutliers = 5;
  // maxNumOutliers = 0;  

  for(Int_t freqInd=0; freqInd<NUM_FREQS; freqInd++){
    // psdOutliers[freqInd] = std::vector<Double_t>(0, 0);
    summedPowSpec[freqInd] = 0;

    TString histName = name + TString::Format("_%d", freqInd);
    TString histTitle = title + TString::Format(" Rayleigh Distribution for %4.2lf MHz bin",
						deltaFMHz*freqInd);
    histTitle += "; Amplitude (mV/MHz); Events/bin";
    TH1D* hTemp = new TH1D(histName, histTitle,
			   NUM_AMPLITUDE_BINS,
			   INITIAL_MIN_AMPLITUDE,
			   INITIAL_MAX_AMPLITUDE);
    hTemp->SetDirectory(0);
    hTemp->Sumw2();

    // This prepocessor variable is defined in the Makefile and is (I hope)
    // querying the ROOT major version number. At the moment it asks whether the
    // ROOT major version number is >= 6, hence the name.
    // It works for me, for now.
#ifdef IS_ROOT_6
    hTemp->SetCanExtend(TH1::kAllAxes);
#else
    hTemp->SetBit(TH1::kCanRebin);
#endif

    hRayleighs[freqInd] = hTemp;
    hRayleighFits[freqInd] = NULL;

    rayleighFitChiSquares[freqInd] = -1;
    rayleighFitChiSquaresRisingEdge[freqInd] = -1;
    rayleighFitChiSquaresRisingEdgeAndHalfFalling[freqInd] = -1;

    rayleighAmplitudes[freqInd] = -1;
    rayleighAmplitudesRisingEdge[freqInd] = -1;
    rayleighAmplitudesRisingEdgeAndHalfFalling[freqInd] = -1;

    rayleighNdf[freqInd] = -1;
    rayleighNdfRisingEdge[freqInd] = -1;
    rayleighNdfRisingEdgeAndHalfFalling[freqInd] = -1;

    for(int i=0; i<MAX_NUM_OUTLIERS; i++){
      outliers[freqInd][i] = -1;
    }
    numOutliers[freqInd] = 0;
  }
}

AveragePowerSpectrum::~AveragePowerSpectrum(){
  deleteRayleighDistributions();
}


TString AveragePowerSpectrum::getRayleighFunctionText(){
  return "([0]*x/([1]*[1]))*exp(-x*x/(2*[1]*[1]))";
}

TF1* AveragePowerSpectrum::makeRayleighFunction(TString name, Double_t xMin, Double_t xMax){
  TF1* fRay = new TF1(name, getRayleighFunctionText(), xMin, xMax);
  return fRay;
}

size_t AveragePowerSpectrum::add(TGraph* gr){

  // Double_t* ps = FancyFFTs::getPowerSpectrum(numSamples, gr->GetY(), 1e-3*deltaT, PowSpecNorm::kPowSpecDensity);
  // Double_t* ps = FancyFFTs::getPowerSpectrum(numSamples, gr->GetY(), deltaT, PowSpecNorm::kPowSpecDensity);

  // Returns the sum over V[i]*V[i], does not normalize bin width by frequency.
  Double_t* ps = FancyFFTs::getPowerSpectrum(NUM_SAMPLES, gr->GetY(),
					     NOMINAL_SAMPLING_DELTAT, PowSpecNorm::kSum);
  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    // Double_t sqrtPSD = TMath::Sqrt(ps[freqInd]);
    Double_t sqrtPSD = TMath::Sqrt(ps[freqInd])/(deltaFMHz);
    TH1D* h = hRayleighs[freqInd];
    Double_t histMaxVal = h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()+1);
    if(sqrtPSD < histMaxVal){
      h->Fill(sqrtPSD);
    }
    else{
      outliers[freqInd][numOutliers[freqInd]] = sqrtPSD;
      numOutliers[freqInd]++;


      if(numOutliers[freqInd]>=MAX_NUM_OUTLIERS){

	while(numOutliers[freqInd] >= MAX_NUM_OUTLIERS){
	  histMaxVal*=2;
	  
	  Double_t outliersTemp[MAX_NUM_OUTLIERS] = {0};
	  Int_t numOutliersTemp = 0;
	
	  for(int i=0; i<MAX_NUM_OUTLIERS; i++){
	    if(outliers[freqInd][i] < histMaxVal){
	      h->Fill(outliers[freqInd][i]);
	    }
	    else{
	      outliersTemp[numOutliersTemp] = outliers[freqInd][i];
	      numOutliersTemp++;
	    }
	  }

	  numOutliers[freqInd] = 0;
	  for(int i=0; i < numOutliersTemp; i++){
	    outliers[freqInd][i] = outliersTemp[i];
	    numOutliers[freqInd]++;
	  }
	  for(int i=numOutliersTemp; i < MAX_NUM_OUTLIERS; i++){
	    outliers[freqInd][i] = -1;
	  }
	}
      }
    }
  }

  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    summedPowSpec[freqInd] += ps[freqInd];
  }
  count++;
  
  delete [] ps;
  return count;
}


void AveragePowerSpectrum::rebinAllRayleighHistograms(Int_t rebinFactor){
  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    TH1D* h = hRayleighs[freqInd];
    h->Rebin(rebinFactor);
  }
}



TH1D* AveragePowerSpectrum::getRayleighHistogramFromFrequencyMHz(Double_t freqMHz){

  Int_t bestFreqInd=0;
  Double_t bestFreqDiff = DBL_MAX;
  for(Int_t freqInd=0; freqInd < NUM_FREQS-1; freqInd++){
    Double_t freqDiff = TMath::Abs(freqInd*deltaFMHz - freqMHz);
    if(freqDiff < bestFreqDiff){
      bestFreqInd = freqInd;
      bestFreqDiff = freqDiff;
    }
  }
  return hRayleighs[bestFreqInd];
}


TH1D* AveragePowerSpectrum::getRayleighHistogram(Int_t freqInd){
  return hRayleighs[freqInd];
}

TF1* AveragePowerSpectrum::getRayleighHistogramFit(Int_t freqInd){
  TH1D* h = hRayleighs[freqInd];
  TString funcName = TString::Format("fit_%d", freqInd);
  TF1* f = (TF1*) h->FindObject(funcName);
  return f;
}


TH2D* AveragePowerSpectrum::makeRayleigh2DHistogram(){

  Double_t maxVal = 0;
  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    TH1D* h = getRayleighHistogram(freqInd);
    Double_t xMax = h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()+1);
    if(xMax > maxVal){
      maxVal = xMax;
    }
  }


  TString name = TString::Format("h2D_%s", GetName());
  TString title = TString::Format("%s Rayleigh Distribution Summary", GetTitle());
  
  TH2D* h2 = new TH2D(name, title, NUM_FREQS, 0, NUM_FREQS*deltaFMHz,
		      NUM_AMPLITUDE_BINS, 0, maxVal);

  h2->GetXaxis()->SetTitle("Frequency (MHz)");
  h2->GetXaxis()->SetNoExponent(1);
  h2->GetYaxis()->SetTitle("Amplitude (mV/MHz)");  
  h2->GetYaxis()->SetNoExponent(1);
  h2->Sumw2();
  
  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    TH1D* h = getRayleighHistogram(freqInd);
    for(Int_t binx=1; binx<=h->GetNbinsX(); binx++){
      Double_t freqMHz = freqInd*deltaFMHz;
      Double_t amplitude = h->GetBinCenter(binx);
      Double_t weight = h->GetBinContent(binx);      
      
      h2->Fill(freqMHz, amplitude, weight);
    }
  }
  
  return h2;
  
}

void AveragePowerSpectrum::fitAllRayleighHistograms(){
  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    fitRayleighHistogram(freqInd);
  }
}


void AveragePowerSpectrum::fitAllRayleighHistogramsRisingEdge(){
  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    fitRayleighHistogramRisingEdge(freqInd);
  }
}

void AveragePowerSpectrum::fitAllRayleighHistogramsRisingEdgeAndHalfFallingEdge(){
  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    fitRayleighHistogramRisingEdgeAndHalfFallingEdge(freqInd);
  }
}


TF1* AveragePowerSpectrum::constructFitFromAmplitude(Int_t freqInd, Double_t amplitude){

  TH1D* h = getRayleighHistogram(freqInd);
  TString fitName = TString::Format("fit_%d_%lf", freqInd, amplitude);
  TF1* fit = makeRayleighFunction(fitName, 0, h->GetBinLowEdge(h->GetNbinsX()+1));
  fit->SetParameter(1, amplitude);
  fit->SetParameter(0, h->GetBinLowEdge(2)*h->Integral());  
  return fit;  
}


TF1* AveragePowerSpectrum::constructFitFromRayleighAmplitude(Int_t freqInd){
  return constructFitFromAmplitude(freqInd, rayleighAmplitudes[freqInd]);
}

TF1* AveragePowerSpectrum::constructFitFromRayleighAmplitudeRisingEdge(Int_t freqInd){
  return constructFitFromAmplitude(freqInd, rayleighAmplitudesRisingEdge[freqInd]);
}

TF1* AveragePowerSpectrum::constructFitFromRayleighAmplitudeRisingEdgeAndHalfFalling(Int_t freqInd){
  return constructFitFromAmplitude(freqInd, rayleighAmplitudesRisingEdgeAndHalfFalling[freqInd]);
}



void AveragePowerSpectrum::fitRayleighHistogram(Int_t freqInd){
  // std::cout << __PRETTY_FUNCTION__ << std::endl;
  
  TH1D* h = getRayleighHistogram(freqInd);
  xHigh[freqInd] = h->GetBinLowEdge(h->GetNbinsX()+1);

  fitRayleighHistogramOverRange(freqInd, h->GetBinLowEdge(1), xHigh[freqInd],
				rayleighAmplitudes,
				rayleighFitChiSquares,
				rayleighNdf,			
				rayleighFitChiSquaresFullRange,
				rayleighNdfFullRange);
}

void AveragePowerSpectrum::fitRayleighHistogramRisingEdge(Int_t freqInd){
  // std::cout << __PRETTY_FUNCTION__ << std::endl;

  TH1D* h = getRayleighHistogram(freqInd);
  Int_t peakBin = RootTools::getPeakBinOfHistogram(h);

  xHighRisingEdge[freqInd] = h->GetBinCenter(peakBin);
  fitRayleighHistogramOverRange(freqInd, h->GetBinLowEdge(1), xHighRisingEdge[freqInd],
				rayleighAmplitudesRisingEdge,
				rayleighFitChiSquaresRisingEdge,
				rayleighNdfRisingEdge,
				rayleighFitChiSquaresRisingEdgeFullRange,
				rayleighNdfRisingEdgeFullRange);
}


void AveragePowerSpectrum::fitRayleighHistogramRisingEdgeAndHalfFallingEdge(Int_t freqInd){
  // std::cout << __PRETTY_FUNCTION__ << std::endl;

  TH1D* h = getRayleighHistogram(freqInd);
  Double_t xHighTemp = -1;
  Int_t peakBin = RootTools::getPeakBinOfHistogram(h);
  Double_t peakVal = h->GetBinContent(peakBin);

  for(Int_t binx=peakBin; binx<=h->GetNbinsX(); binx++){
    Double_t binVal = h->GetBinContent(binx);
    if(binVal < peakVal*0.5){
      xHighTemp = h->GetBinCenter(binx);
      // std::cout << fName << "\t" << freqInd << "\t" << binx << "\t" << xHigh << std::endl;
      break;
    }
  }

  xHighRisingEdgeAndHalfFalling[freqInd] = xHighTemp;
  fitRayleighHistogramOverRange(freqInd, h->GetBinLowEdge(1), xHighTemp,
				rayleighAmplitudesRisingEdgeAndHalfFalling,
				rayleighFitChiSquaresRisingEdgeAndHalfFalling,
				rayleighNdfRisingEdgeAndHalfFalling,
				rayleighFitChiSquaresRisingEdgeAndHalfFallingFullRange,
				rayleighNdfRisingEdgeAndHalfFallingFullRange);
  // fitRayleighHistogramOverRange(freqInd, h->GetBinLowEdge(1), h->GetBinLowEdge(3));  
}



void AveragePowerSpectrum::fitRayleighHistogramOverRange(Int_t freqInd, Double_t xLowVal, Double_t xHighVal,
							 Double_t* rAmplitudes,
							 Double_t* rChiSquares,
							 Int_t* rNdf,
							 Double_t* rChiSquaresFullRange,
							 Int_t* rNdfFullRange){


  TH1D* h = getRayleighHistogram(freqInd);

  if(h->Integral() > 0){ // Fit will fail with empty histogram

    TString fitName = TString::Format("fit_%d", freqInd);
    TF1* fit = hRayleighFits[freqInd];

    if(fit==NULL){
      fit = makeRayleighFunction(fitName, xLowVal, xHighVal);
    }

    Double_t mean = h->GetMean();
    // Mean of rayleigh distribution = sigma* (pi/2)^{0.5}
    Double_t sigGuess = mean / TMath::Sqrt(0.5*TMath::Pi());
    
    fit->SetParameter(1, sigGuess);

    // Can't not set upper bound if setting lower bound..
    // Therefore set the upper bound to be insanely high.
    // fit->SetParLimits(1, 0, h->GetMean()*1e9); 

    Double_t histArea = h->Integral() * h->GetBinLowEdge(2);
    fit->FixParameter(0, histArea);

    h->Fit(fit, "RQ0");
    
    // h->GetFunction(fit->GetName())->ResetBit(TF1::kNotDraw);
    rAmplitudes[freqInd] = h->GetFunction(fit->GetName())->GetParameter(1);
    rChiSquares[freqInd] = h->GetFunction(fit->GetName())->GetChisquare();
    rNdf[freqInd] = h->GetFunction(fit->GetName())->GetNDF(); 

    rChiSquaresFullRange[freqInd] = 0;
    rNdfFullRange[freqInd] = 0;
    for(int binx=1; binx<=h->GetNbinsX(); binx++){
      Double_t binVal = h->GetBinContent(binx);
      Double_t binError = h->GetBinError(binx);
      Double_t fitVal = fit->Eval(h->GetBinCenter(binx));
      if(binError > 0){
	Double_t chi = (binVal - fitVal)/binError;
	rChiSquaresFullRange[freqInd] += chi*chi;
	rNdfFullRange[freqInd]++;
      }
    }
    
    delete fit;
  }
}


void AveragePowerSpectrum::deleteRayleighDistributions(){
  for(int freqInd=0; freqInd < NUM_FREQS; freqInd++){
    if(hRayleighs[freqInd]!=NULL){
      delete hRayleighs[freqInd];
      hRayleighs[freqInd] = NULL;
    }
  }
}

TGraph* AveragePowerSpectrum::makeAvePowSpecTGraph(){

  TString name = TString::Format("gr_%s", GetName());
  TString title = TString::Format("%s", GetTitle());
  
  std::vector<Double_t> avePowSpec(NUM_FREQS);
  
  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    if(count > 0){
      avePowSpec[freqInd] = summedPowSpec[freqInd]/count;
    }
  }
  
  // Double_t termResOhms = 50;
  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    avePowSpec[freqInd]/=(deltaFMHz);
  }
  
  Double_t freqArray[NUM_FREQS];
  for(int freqInd=0; freqInd < NUM_FREQS; freqInd++){
    freqArray[freqInd] = deltaFMHz*freqInd;
  }
  
  TGraph* gr = new TGraph(NUM_FREQS, freqArray, &avePowSpec[0]);
  gr->SetName(name);
  gr->SetTitle(title);
  return gr;
}


TGraph* AveragePowerSpectrum::makeAvePowSpecTGraph_dB(){
  
  TGraph* gr = makeAvePowSpecTGraph();
  TString name = TString::Format("%s_dB", gr->GetName());

  for(Int_t freqInd=0; freqInd < NUM_FREQS; freqInd++){
    Double_t y = gr->GetY()[freqInd];
    // gr->GetY()[freqInd] = 10*TMath::Log10(y);
    gr->GetY()[freqInd] = 10*TMath::Log10(y);
  }
  return gr;
}

  
