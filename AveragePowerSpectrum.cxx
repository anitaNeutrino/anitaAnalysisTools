#include "AveragePowerSpectrum.h"

ClassImp(AveragePowerSpectrum);



AveragePowerSpectrum::AveragePowerSpectrum(){
  // Don't use this.
  freqArray = NULL;
}


AveragePowerSpectrum::AveragePowerSpectrum(TString name, TString title, Double_t dt,
					   Int_t nSamp, AveragePowerSpectrum::mode_t psMode){


  deltaT = dt;
  numSamples = nSamp;
  numFreqs = FancyFFTs::getNumFreqs(nSamp);
  freqArray = FancyFFTs::getFreqArray(nSamp, dt);
  for(Int_t freqInd=0; freqInd<numFreqs; freqInd++){
    freqArray[freqInd] *= 1e3; // GHz -> MHz;
  }  
  deltaF = freqArray[1]-freqArray[0];

  mode = psMode;
  count=0;
  SetName(name);
  SetTitle(title);  

  // Default power spec histogram values.
  maxNumOutliers = 5;

  for(Int_t freqInd=0; freqInd<numFreqs; freqInd++){
    psdOutliers.push_back(std::vector<Double_t>(0, 0));
    summedPowSpec.push_back(0);

    TString histName = name + TString::Format("_%d", freqInd);
    TString histTitle = title + TString::Format(" Rayleigh Distribution for %4.2lf MHz bin",
					       freqArray[freqInd]);
    histTitle += "; Amplitude (mV/MHz); Number of events";
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

    hRayleighs.push_back(hTemp);
    hRayleighFits.push_back(NULL);    
  }
}

AveragePowerSpectrum::~AveragePowerSpectrum(){
  deleteRayleighDistributions();
  if(freqArray){
    delete [] freqArray;
    freqArray = NULL;
  }
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
  Double_t* ps = FancyFFTs::getPowerSpectrum(numSamples, gr->GetY(), deltaT, PowSpecNorm::kSum);
  for(Int_t freqInd=0; freqInd < numFreqs; freqInd++){
    // Double_t sqrtPSD = TMath::Sqrt(ps[freqInd]);
    Double_t sqrtPSD = TMath::Sqrt(ps[freqInd])/(deltaF);
    TH1D* h = hRayleighs.at(freqInd);
    Double_t histMaxVal = h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()+1);
    if(sqrtPSD < histMaxVal){
      h->Fill(sqrtPSD);
    }
    else{
      psdOutliers.at(freqInd).push_back(sqrtPSD);
      if(psdOutliers.at(freqInd).size() > maxNumOutliers){
	std::vector<Double_t> newOutliers;
	for(UInt_t index = 0; index < psdOutliers.at(freqInd).size(); index++){	    
	  Double_t outlier = psdOutliers.at(freqInd).at(index);
	  if(outlier < histMaxVal*2){
	    h->Fill(outlier);
	  }
	  else{
	    newOutliers.push_back(outlier);
	  }
	}
	psdOutliers.at(freqInd).clear();
	for(UInt_t newIndex=0; newIndex < newOutliers.size(); newIndex++){
	  psdOutliers.at(freqInd).push_back(newOutliers.at(newIndex));
	}
      }
    }
  }

  size_t retVal = 0;
  if(mode==AveragePowerSpectrum::kRolling){
    storedPowSpecs.push_back(std::vector<Double_t> (ps, ps+numFreqs));
    retVal = storedPowSpecs.size();
  }
  else{
    for(Int_t freqInd=0; freqInd < numFreqs; freqInd++){
      summedPowSpec.at(freqInd) += ps[freqInd];
    }
    count++;
    retVal = count;
  }

  
  delete [] ps;
  return retVal;
}


TH1D* AveragePowerSpectrum::getRayleighHistogramFromFrequencyMHz(Double_t freqMHz){

  Int_t bestFreqInd=0;
  Double_t bestFreqDiff = DBL_MAX;
  for(Int_t freqInd=0; freqInd < numFreqs-1; freqInd++){
    Double_t freqDiff = TMath::Abs(freqArray[freqInd] - freqMHz);
    if(freqDiff < bestFreqDiff){
      bestFreqInd = freqInd;
      bestFreqDiff = freqDiff;
    }
  }
  return hRayleighs.at(bestFreqInd);
}


TH1D* AveragePowerSpectrum::getRayleighHistogram(Int_t freqInd){
  return hRayleighs.at(freqInd);
}

TF1* AveragePowerSpectrum::getRayleighHistogramFit(Int_t freqInd){
  TH1D* h = hRayleighs.at(freqInd);
  TString funcName = TString::Format("fit_%d", freqInd);
  TF1* f = (TF1*) h->FindObject(funcName);
  return f;
}


TH2D* AveragePowerSpectrum::makeRayleigh2DHistogram(){

  Double_t maxVal = 0;
  for(Int_t freqInd=0; freqInd < numFreqs; freqInd++){
    TH1D* h = getRayleighHistogram(freqInd);
    Double_t xMax = h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()+1);
    if(xMax > maxVal){
      maxVal = xMax;
    }
  }


  TString name = TString::Format("h2D_%s", GetName());
  TString title = TString::Format("%s Rayleigh Distribution Summary", GetTitle());
  
  TH2D* h2 = new TH2D(name, title, numFreqs, 0, freqArray[numFreqs-1],
		      NUM_AMPLITUDE_BINS, 0, maxVal);

  h2->GetXaxis()->SetTitle("Frequency (MHz)");
  h2->GetXaxis()->SetNoExponent(1);
  h2->GetYaxis()->SetTitle("Amplitude (mV/MHz)");  
  h2->GetYaxis()->SetNoExponent(1);
  h2->Sumw2();
  
  for(Int_t freqInd=0; freqInd < numFreqs; freqInd++){
    TH1D* h = getRayleighHistogram(freqInd);
    for(Int_t binx=1; binx<=h->GetNbinsX(); binx++){
      Double_t freqMHz = freqArray[freqInd];
      Double_t amplitude = h->GetBinCenter(binx);
      Double_t weight = h->GetBinContent(binx);      
      
      h2->Fill(freqMHz, amplitude, weight);
    }
  }
  
  return h2;
  
}

void AveragePowerSpectrum::fitAllRayleighHistograms(){
  for(Int_t freqInd=0; freqInd < numFreqs; freqInd++){
    fitRayleighHistogram(freqInd);
  }
}

void AveragePowerSpectrum::fitRayleighHistogram(Int_t freqInd){

  TH1D* h = getRayleighHistogram(freqInd);

  if(h->Integral() > 0){ // Fit will fail with empty histogram

    TString fitName = TString::Format("fit_%d", freqInd);
    TF1* fit = hRayleighFits.at(freqInd);

    if(fit==NULL){
      fit = makeRayleighFunction(fitName,
				 h->GetXaxis()->GetBinLowEdge(1),
				 h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()+1));
    }
    fit->SetParameter(1, h->GetMean());
    h->Fit(fit, "Q0");
    
    h->GetFunction(fit->GetName())->ResetBit(TF1::kNotDraw);

    delete fit;
  }
}

void AveragePowerSpectrum::deleteRayleighDistributions(){
  while(!hRayleighs.empty()){
    TH1D* hTemp = hRayleighs.back();
    delete hTemp;
    hRayleighs.pop_back();    
  }
}

TGraph* AveragePowerSpectrum::makeAvePowSpecTGraph(){

  TString name = TString::Format("gr_%s", GetName());
  TString title = TString::Format("%s Average Power Spectrum", GetTitle());
  
  std::vector<Double_t> avePowSpec(numFreqs);
  
  if(mode==AveragePowerSpectrum::kRolling){

    UInt_t numEvents = storedPowSpecs.size();

    for(UInt_t eventInd=0; eventInd < numEvents; eventInd++){
      for(Int_t freqInd=0; freqInd < numFreqs; freqInd++){
	avePowSpec.at(freqInd) += storedPowSpecs.at(eventInd)[freqInd];
      }
    }
    for(Int_t freqInd=0; freqInd < numFreqs; freqInd++){
      avePowSpec.at(freqInd)/=numEvents;
    }
  }
  else{

    for(Int_t freqInd=0; freqInd < numFreqs; freqInd++){
      if(count > 0){
	avePowSpec.at(freqInd) = summedPowSpec.at(freqInd)/count;
      }
    }
  }

  // Double_t ohms = 50;
  for(Int_t freqInd=0; freqInd < numFreqs; freqInd++){
    avePowSpec.at(freqInd)/=(deltaF);
  }
  
  TGraph* gr = new TGraph(numFreqs, freqArray, &avePowSpec[0]);
  gr->SetName(name);
  gr->SetTitle(title);
  return gr;
}


TGraph* AveragePowerSpectrum::makeAvePowSpecTGraph_dB(){
  
  TGraph* gr = makeAvePowSpecTGraph();
  TString name = TString::Format("%s_dB", gr->GetName());

  for(Int_t freqInd=0; freqInd < numFreqs; freqInd++){
    Double_t y = gr->GetY()[freqInd];
    // gr->GetY()[freqInd] = 10*TMath::Log10(y);
    gr->GetY()[freqInd] = 10*TMath::Log10(y);
  }
  return gr;
}

  
