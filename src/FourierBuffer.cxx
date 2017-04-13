#include "FourierBuffer.h"
#include "TMath.h"
#include "FFTWComplex.h"
#include "TSpectrum.h"

const char* rayleighFuncText = "([0]*x/([1]*[1]))*exp(-x*x/(2*[1]*[1]))";
const char* riceFuncText = "([0]*x/([1]*[1]))*exp(-(x*x+[2]*[2])/(2*[1]*[1]))*TMath::BesselI0([2]*x/([1]*[1]))";



Acclaim::FourierBuffer::FourierBuffer(double timeScaleSeconds, Int_t theAnt, AnitaPol::AnitaPol_t thePol){
  timeScale = timeScaleSeconds;
  ant = theAnt;
  pol = thePol;
  df = -1;

  // will initialize this dynamically to get around this no-copy-constructor bullshit
  spectrum = NULL; 
}


Acclaim::FourierBuffer::~FourierBuffer(){
  if(spectrum){
    delete spectrum;
    spectrum = NULL;
  }
}



size_t Acclaim::FourierBuffer::add(const RawAnitaHeader* header, const AnalysisWaveform& wave){

  // for old compilers, push back compy of emptry vector for speed.
  powerRingBuffer.push_back(std::vector<double>(0));

  // then get reference that vector in the list
  std::vector<double>& freqVec = powerRingBuffer.back();

  // copy frequencies into vector
  const TGraphAligned* grPower = wave.power();
  freqVec.assign(grPower->GetY(), grPower->GetY()+grPower->GetN());

  // dynamically assign deltaF variable
  if(df <= 0){
    df = grPower->GetX()[1] - grPower->GetX()[0];
  }

  if(sumPower.size()==0){
    sumPower.assign(grPower->GetY(), grPower->GetY()+grPower->GetN());
  }
  else{
    for(int freqInd=0; freqInd < grPower->GetN(); freqInd++){
      sumPower.at(freqInd) += grPower->GetY()[freqInd];
    }
  }

  if(grPower->GetN() != (int)sumPower.size()){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unexpected waveform size" << std::endl;
  }

  eventNumbers.push_back(header->eventNumber);
  runs.push_back(header->run);

  Double_t realTime= header->realTime;
  realTime += 1e-9*header->triggerTimeNs;
  realTimesNs.push_back(realTime);

  
  removeOld();
  
  return eventNumbers.size();
}



Int_t Acclaim::FourierBuffer::removeOld(){

  Int_t nPopped = 0;
  if(realTimesNs.size() > 0){

    Double_t mostRecentTime = realTimesNs.back();

    while(mostRecentTime - realTimesNs.front() > timeScale){
      realTimesNs.pop_front();
      eventNumbers.pop_front();
      runs.pop_front();
      std::vector<double>& removeThisPower = powerRingBuffer.front();
      for(unsigned int freqInd=0; freqInd < removeThisPower.size(); freqInd++){
	sumPower.at(freqInd) -= removeThisPower.at(freqInd);
      }
      powerRingBuffer.pop_front();
      nPopped++;
    }
  }

  return nPopped;
}






TH1D* Acclaim::FourierBuffer::fillRayleighInfo(Int_t freqBin, RayleighInfo* info) const{

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

  // FFTWComplex val = freqVecs.back().at(freqBin);  
  // info->eventAmp = TMath::Sqrt(val.re*val.re + val.im*val.im);
  double powerBin = powerRingBuffer.back().at(freqBin);
  info->eventAmp = TMath::Sqrt(powerBin);  
  
  // Fill histogram related histogram info
  info->nBins = hRay->GetNbinsX();  
  info->integralPlusOverUnderFlow = hRay->Integral(0, info->nBins+1);
  info->binWidth = hRay->GetBinWidth(1);
  // std::cout << info->binWidth << std::endl;
  info->rayFitNorm = info->integralPlusOverUnderFlow*info->binWidth;
  info->rayGuessAmp = hRay->GetMean()*TMath::Sqrt(2./TMath::Pi());  

  // Set fit params
  fRay->FixParameter(0, info->rayFitNorm);
  fRay->SetParameter(1, info->rayGuessAmp);  
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



TH1D* Acclaim::FourierBuffer::fillRiceInfo(Int_t freqBin, RiceInfo* info) const{

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






TH1D* Acclaim::FourierBuffer::getRayleighDistribution(Int_t freqBin) const{





  double meanAmp = 0;
  int count = 0;
  // std::list<std::vector<FFTWComplex> >::iterator it;
  // for(it = freqVecs.begin(); it!=freqVecs.end(); ++it){  
  std::list<std::vector<double> >::const_iterator it;  
  for(it = powerRingBuffer.begin(); it!=powerRingBuffer.end(); ++it){

    // std::vector<FFTWComplex>& freqVec = (*it);
    const std::vector<double>& power = (*it);    
    
    // Double_t absSq = freqVec.at(freqBin).re*freqVec.at(freqBin).re + freqVec.at(freqBin).im*freqVec.at(freqBin).im;
    Double_t absSq = power.at(freqBin);
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


  for(it = powerRingBuffer.begin(); it!=powerRingBuffer.end(); ++it){

    // std::vector<FFTWComplex>& freqVec = (*it);
    const std::vector<double>& power = (*it);    

    // Double_t absSq = freqVec.at(freqBin).re*freqVec.at(freqBin).re + freqVec.at(freqBin).im*freqVec.at(freqBin).im;
    Double_t absSq = power.at(freqBin);
    Double_t abs = TMath::Sqrt(absSq);

    hRayleigh->Fill(abs);

    // std::cout << abs << std::endl;
  }
  // for(it = freqVecs.begin(); it!=freqVecs.end(); ++it){

  //   std::vector<FFTWComplex>& freqVec = (*it);

  //   Double_t absSq = freqVec.at(freqBin).re*freqVec.at(freqBin).re + freqVec.at(freqBin).im*freqVec.at(freqBin).im;
  //   Double_t abs = TMath::Sqrt(absSq);

  //   hRayleigh->Fill(abs);

  //   // std::cout << abs << std::endl;
  // }

  return hRayleigh;

}


TGraphAligned* Acclaim::FourierBuffer::getAvePowSpec_dB(double thisTimeRange) const{
  
  TGraphAligned* gr = getAvePowSpec(thisTimeRange);
  gr->dBize();
  return gr;
}


TGraphAligned* Acclaim::FourierBuffer::getAvePowSpec(double thisTimeRange) const{

  // set default value (which is the whole range of the fourier buffer)
  thisTimeRange = thisTimeRange < 0 ? timeScale : thisTimeRange;

  int nEvents = eventNumbers.size();

  TGraphAligned* gr = new TGraphAligned(sumPower.size());
  
  TString name;
  TString title;  
  if(ant < 0){
    name = TString::Format("grAvePowSpec_%u_%u", eventNumbers.front(), eventNumbers.back());
    title = TString::Format("Average Power Spectrum event %u - %u (unknown channel)", eventNumbers.front(), eventNumbers.back());    
  }
  else{
    const char* polName = AnitaPol::kHorizontal ? "HPol" : "VPol";
    name = TString::Format("grAvePowSpec_%d_%s_%u_%u", ant, polName, eventNumbers.front(), eventNumbers.back());
    title = TString::Format("Average Power Spectrum %d %s event %u - %u", ant, polName, eventNumbers.front(), eventNumbers.back());        
  }
  gr->SetName(name);
  gr->SetTitle(title);
  

  for(unsigned freqInd=0; freqInd < sumPower.size(); freqInd++){
    gr->SetPoint(freqInd, freqInd*df, sumPower.at(freqInd)/nEvents);
  }

  return gr;
}



TGraphAligned* Acclaim::FourierBuffer::getBackground_dB(double thisTimeRange) const{
  
  TGraphAligned* gr = getBackground(thisTimeRange);
  gr->dBize();
  return gr;
}


TGraphAligned* Acclaim::FourierBuffer::getBackground(double thisTimeRange) const{

  TGraphAligned* gr = getAvePowSpec(thisTimeRange);

  if(!spectrum){
    spectrum = new TSpectrum();
  }

  // at some point this changed from doubles to floats, not added a new method, changed...
  // git blames Lorenzo Moneta... 
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)

  spectrum->Background(gr->GetY(),gr->GetN(),
                      6,TSpectrum::kBackDecreasingWindow,
                      TSpectrum::kBackOrder2,kFALSE,
                      TSpectrum::kBackSmoothing3,kFALSE);
#else
  std::vector<float> tempFloats(gr->GetN(), 0);
  for(int i=0; i < gr->GetN(); i++){
    tempFloats[i] = gr->GetY()[i];
  }
  spectrum->Background(&tempFloats[0],gr->GetN(),
                      6,TSpectrum::kBackDecreasingWindow,
                      TSpectrum::kBackOrder2,kFALSE,
                      TSpectrum::kBackSmoothing3,kFALSE);
  for(int i=0; i < gr->GetN(); i++){
    gr->GetY()[i] = tempFloats[i];
  }
#endif
  

  return gr;
}
