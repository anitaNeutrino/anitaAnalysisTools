#include "FourierBuffer.h"
#include "TMath.h"
#include "FFTWComplex.h"

const char* rayleighFuncText = "([0]*x/([1]*[1]))*exp(-x*x/(2*[1]*[1]))";
const char* riceFuncText = "([0]*x/([1]*[1]))*exp(-(x*x+[2]*[2])/(2*[1]*[1]))*TMath::BesselI0([2]*x/([1]*[1]))";

FourierBuffer::FourierBuffer(Int_t timeScaleSeconds){
  timeScale = timeScaleSeconds;
}





size_t FourierBuffer::add(const RawAnitaHeader* header, const AnalysisWaveform& wave){

  // std::vector<FFTWComplex> tempFreqs(wave.freq(), wave.freq()+NUM_FREQS);


  // for old compilers, push back compy of emptry vector for speed.
  freqVecs.push_back(std::vector<FFTWComplex>());

  // then get reference to vector in the list
  std::vector<FFTWComplex>& freqVec = freqVecs.back();

  // copy frequencies into vector
  freqVec.assign(wave.freq(), wave.freq()+wave.Nfreq());


  // (wave.freq(), wave.freq()+wave.Nfreq());
  // std::cout << wave.Nfreq() << "\t" << tempFreqs.size() << std::endl;
  // std::cout << wave.Nfreq() << "\t" << tempFreqs.size() << std::endl;
  // std::cout << "\t" << freqVec.size() << "\t" << wave.Nfreq() << std::endl;
  eventNumbers.push_back(header->eventNumber);
  runs.push_back(header->run);

  Double_t realTime= header->realTime;
  realTime += 1e-9*header->triggerTimeNs;
  realTimesNs.push_back(realTime);

  
  // std::cout << freqVecs.size() << "\t" << eventNumbers.size() << "\t" << runs.size() << "\t" << realTimesNs.size() << std::endl;
  removeOld();
  // std::cout << freqVecs.size() << "\t" << eventNumbers.size() << "\t" << runs.size() << "\t" << realTimesNs.size() << std::endl;

  // std::cout << std::endl;
  
  return eventNumbers.size();
}



Int_t FourierBuffer::removeOld(){

  Int_t nPopped = 0;
  if(realTimesNs.size() > 0){

    Double_t mostRecentTime = realTimesNs.back();

    while(mostRecentTime - realTimesNs.front() > timeScale){
      // std::cout << realTimes.size() << "\t" << mostRecentTime <<  "\t" << realTimes.front() << "\t" << mostRecentTime - realTimes.front() << "\t" << timeScale << std::endl;
      realTimesNs.pop_front();
      eventNumbers.pop_front();
      runs.pop_front();
      freqVecs.pop_front();
      nPopped++;
    }
  }

  return nPopped;
}






TH1D* FourierBuffer::fillRayleighInfo(Int_t freqBin, RayleighInfo* info){

  TH1D* hRay = getRayleighDistribution(freqBin);
  TString fName = TString::Format("fRay%d", freqBin);
  Double_t xMin= hRay->GetXaxis()->GetBinLowEdge(1);
  Double_t xMax= hRay->GetXaxis()->GetBinLowEdge(hRay->GetNbinsX()+1);
  TF1* fRay = new TF1(fName, rayleighFuncText, xMin, xMax);

  // Fill event level info
  info->eventNumber = eventNumbers.back();
  info->run = runs.back();
  info->realTimeNs = realTimesNs.back();
  
  info->firstRun = runs.front();
  info->firstEventNumber = eventNumbers.front();  
  info->firstRealTimeNs = realTimesNs.front();  


  FFTWComplex val = freqVecs.back().at(freqBin);
  info->eventAmp = TMath::Sqrt(val.re*val.re + val.im*val.im);
  
  // Fill histogram related histogram info
  info->integralPlusOverUnderFlow = hRay->Integral() + hRay->GetBinContent(0) + hRay->GetBinContent(hRay->GetNbinsX()+1);
  info->binWidth = hRay->GetXaxis()->GetBinLowEdge(2) - hRay->GetXaxis()->GetBinLowEdge(1);
  // std::cout << info->binWidth << std::endl;
  info->nBins = hRay->GetNbinsX();
  info->rayFitNorm = info->integralPlusOverUnderFlow*info->binWidth;
  info->rayGuessAmp = hRay->GetMean()*TMath::Sqrt(2./TMath::Pi());  

  // Set fit params
  fRay->FixParameter(0, info->rayFitNorm);
  // fRay->SetParameter(1, info->rayGuessAmp);
  fRay->SetParameter(1, info->rayFitAmp);  
  fRay->SetParLimits(1, 0, 1e9); // essentially infinite

  // Do fit
  hRay->Fit(fRay, "Q0"); //, "", "", xMin, xMax);

  // Fill fit result params
  info->rayFitAmp = fRay->GetParameter(1);
  info->rayFitAmpError = fRay->GetParError(1);
  info->rayChiSquare = fRay->GetChisquare();
  info->rayNdf = fRay->GetNDF();


  // think copy lives with the histogram
  delete fRay;

  return hRay;
}



TH1D* FourierBuffer::fillRiceInfo(Int_t freqBin, RiceInfo* info){

  TH1D* hRay = fillRayleighInfo(freqBin, info);


  if(info->rayChiSquare/info->rayNdf > 2){

    TString fName = TString::Format("fRice%d", freqBin);
    Double_t xMin = hRay->GetXaxis()->GetBinLowEdge(1);
    Double_t xMax = hRay->GetXaxis()->GetBinLowEdge(hRay->GetNbinsX()+1);
    TF1* fRice = new TF1(fName, riceFuncText, xMin, xMax);

    fRice->FixParameter(0, info->rayFitNorm);
    fRice->SetParameter(1, info->rayGuessAmp);
    fRice->SetParLimits(1, 0, 1e9); // essentially infinite
    fRice->SetParameter(2, 0);
    fRice->SetParLimits(2, 0, 1e9); // essentially infinite

    hRay->Fit(fRice, "Q0");


    info->riceFitNorm = fRice->GetParameter(0); // fixed for the fit  
    info->riceFitAmp = fRice->GetParameter(1);
    info->riceFitAmpError = fRice->GetParError(1);  
    info->riceFitSignal = fRice->GetParameter(2);
    info->riceFitSignalError = fRice->GetParError(2);  

    info->riceChiSquare = fRice->GetChisquare();
    info->riceNdf = fRice->GetNDF();

    delete fRice;
  }
  
  return hRay;
}






TH1D* FourierBuffer::getRayleighDistribution(Int_t freqBin){





  double meanAmp = 0;
  int count = 0;
  std::list<std::vector<FFTWComplex> >::iterator it;
  for(it = freqVecs.begin(); it!=freqVecs.end(); ++it){

    std::vector<FFTWComplex>& freqVec = (*it);

    Double_t absSq = freqVec.at(freqBin).re*freqVec.at(freqBin).re + freqVec.at(freqBin).im*freqVec.at(freqBin).im;
    Double_t abs = TMath::Sqrt(absSq);

    // could add some logic to exclude outliers?
    meanAmp += abs;
    count++;

    // std::cout << abs << std::endl;
  }

  // now we have the mean
  meanAmp/=count;


  // Now we can make a guess at the axis limits we want
  const double fracOfEventsWanted = 0.99;

  const double maxAmp = meanAmp*TMath::Sqrt(-(4./TMath::Pi())*TMath::Log(1.0 - fracOfEventsWanted));

  // get 5 bins on the rising edge
  const Int_t risingEdgeBins = 4;

  // peak is at sigma =
  const double sigmaGuess = TMath::Sqrt(2./TMath::Pi())*meanAmp;
  const double binWidth = sigmaGuess/risingEdgeBins;

  const int nBins = TMath::Nint(maxAmp/binWidth);

  // put it all together...

  // std::cout << sigmaGuess << "\t" << maxAmp << "\t" << binWidth << "\t" << nBins << std::endl;



  TString name = TString::Format("hRayleigh%u", eventNumbers.back());
  TH1D* hRayleigh = new TH1D(name, name, nBins, 0, maxAmp);
  hRayleigh->Sumw2();


  // std::list<std::vector<FFTWComplex> >::iterator it;
  for(it = freqVecs.begin(); it!=freqVecs.end(); ++it){

    std::vector<FFTWComplex>& freqVec = (*it);

    Double_t absSq = freqVec.at(freqBin).re*freqVec.at(freqBin).re + freqVec.at(freqBin).im*freqVec.at(freqBin).im;
    Double_t abs = TMath::Sqrt(absSq);

    hRayleigh->Fill(abs);

    // std::cout << abs << std::endl;
  }

  return hRayleigh;

}
