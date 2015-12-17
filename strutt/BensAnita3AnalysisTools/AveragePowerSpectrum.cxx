#include "AveragePowerSpectrum.h"

AveragePowerSpectrum::AveragePowerSpectrum(TString histBaseNameIn, Double_t dt, Int_t nSamp, AveragePowerSpectrum::mode_t psMode){


  deltaT = dt;
  numSamples = nSamp;
  numFreqs = FancyFFTs::getNumFreqs(nSamp);
  freqArray = FancyFFTs::getFreqArray(nSamp, dt);
  mode = psMode;
  count=0;
  histBaseName = histBaseNameIn;
  
  const Int_t numPowBins = 64;
  const Double_t maxPow = 128;
  const Double_t minPow = 0;
  maxNumOutliers = 10;
  
  for(Int_t freqInd=0; freqInd<numFreqs; freqInd++){
    psdOutliers.push_back(std::vector<Double_t>(0, 0));
    summedPowSpec.push_back(0);

    TString name = histBaseName + TString::Format("_%d_%d", freqInd, numFreqs);
    TH1D* hTemp = new TH1D(name, name, numPowBins, minPow, maxPow);
    hTemp->Sumw2(2);

    // This prepocessor variable should be defined in the Makefile and (I hope)
    // is querying the ROOT version. At the moment it asks whether the
    // ROOT major version number is >= 6, hence the name.
    // It works for me, for now.
#ifdef IS_ROOT_6
    hTemp->SetCanExtend(TH1::kAllAxes);
#else
    hTemp->SetBit(TH1::kCanRebin);
#endif


  }
}

AveragePowerSpectrum::~AveragePowerSpectrum(){
  emptyRolling();
  emptyRayleighs();  
  delete [] freqArray;
}


TF1* AveragePowerSpectrum::makeRayleighFunction(TString name, Double_t xMin, Double_t xMax){
  TF1* fRay = new TF1(name, "([0]*x/([1]*[1]))*exp(-x*x/(2*[1]*[1]))", xMin, xMax);
  return fRay;
}

size_t AveragePowerSpectrum::add(TGraph* gr){

  Double_t* ps = FancyFFTs::getPowerSpectrum(numSamples, gr->GetY(), deltaT, PowSpecNorm::kPowSpecDensity);
  for(Int_t freqInd=0; freqInd < numFreqs; freqInd++){
    Double_t sqrtPSD = TMath::Sqrt(ps[freqInd]);
    TH1D* h = hRayleighs.at(freqInd);
    Double_t histMaxVal = h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()+1);
    if(sqrtPSD < histMaxVal){
      h->Fill(sqrtPSD);
    }
    else{
      psdOutliers.at(freqInd).push_back(sqrtPSD);
      if(psdOutliers.at(freqInd).size() > maxNumOutliers){
	std::vector<Double_t> newOutliers;

	// for(int index = psdOutliers.at(freqInd).size()-1; index >= 0; index--){
	for(UInt_t index = 0; index < psdOutliers.at(freqInd).size(); index++){	    
	  Double_t outlier = psdOutliers.at(freqInd).at(index);
	  if(outlier < histMaxVal*2){


#ifndef IS_ROOT_6
	    // Need to rebin the axis ourselves for ROOT versions lower than 6?
	    if(h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()+1) <= histMaxVal){
	      h->RebinAxis(histMaxVal*2, h->GetXaxis());
	    }
#endif


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
  
  if(mode==AveragePowerSpectrum::kRolling){
    storedPowSpecs.push_back(ps);
    return storedPowSpecs.size();
  }
  else{
    for(Int_t freqInd=0; freqInd < numFreqs; freqInd++){
      // if(histBaseName.Contains("North")){
      // 	std::cout << histBaseName.Data() << "\t" << freqInd << ps[freqInd] << "\t" << summedPowSpec.at(freqInd) << "\t" << count << std::endl;
      // 	}
      summedPowSpec.at(freqInd) += ps[freqInd];
    }
    delete [] ps;
    count++;
    return count;
  }
}

void AveragePowerSpectrum::emptyRolling(){
  while(!storedPowSpecs.empty()){
    delete [] storedPowSpecs.back();
    storedPowSpecs.pop_back();    
  }
}

void AveragePowerSpectrum::emptyRayleighs(){
  while(!hRayleighs.empty()){
    TH1D* hTemp = hRayleighs.back();
    hTemp->Write();
    delete hTemp;
    hRayleighs.pop_back();    
  }
}

TGraph* AveragePowerSpectrum::get(TString name, TString title){

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
      // std::cerr << freqInd << "\t" << avePowSpec.at(freqInd) << "\t" << summedPowSpec.at(freqInd) << "\t" << count << std::endl;
      Double_t norm = count > 0 ? count : 1;
      avePowSpec.at(freqInd) = summedPowSpec.at(freqInd)/norm;
    }
  }

  TGraph* gr = new TGraph(numFreqs, freqArray, &avePowSpec[0]);
  gr->SetName(name);
  gr->SetTitle(title);
  return gr;
}


TGraph* AveragePowerSpectrum::getScaled(TString name, TString title){

  TGraph* gr = get(name, title);
  for(Int_t freqInd=0; freqInd < numFreqs; freqInd++){
    Double_t y = gr->GetY()[freqInd];
    y*=1e-3;
    // gr->GetY()[freqInd] = 10*TMath::Log10(y);
    gr->GetY()[freqInd] = 10*TMath::Log10(y);
    gr->GetX()[freqInd]*= 1e3;
  }
  return gr;
}

  
