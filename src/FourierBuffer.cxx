#include "FourierBuffer.h"
#include "TMath.h"
#include "FFTWComplex.h"
#include "TSpectrum.h"
#include "RayleighHist.h"
#include "TF1.h"
 

Acclaim::FourierBuffer::FourierBuffer(Int_t theBufferSize, Int_t theAnt, AnitaPol::AnitaPol_t thePol){
  
  // timeScale = timeScaleSeconds;
  bufferSize = theBufferSize <= 0 ? 1000 : theBufferSize;
  ant = theAnt;
  pol = thePol;
  df = -1;
  fDrawFreqBin = 100; //26;

  // will initialize this dynamically to get around this no-copy-constructor bullshit
  spectrum = NULL;


  const char* rayleighFuncText = "([0]*x/([1]*[1]))*exp(-x*x/(2*[1]*[1]))";
  // const char* riceFuncText = "([0]*x/([1]*[1]))*exp(-(x*x+[2]*[2])/(2*[1]*[1]))*TMath::BesselI0([2]*x/([1]*[1]))";
  
  TString funcName = TString::Format("fRay_%d_%d", ant, pol);
  fRay = new TF1(funcName, rayleighFuncText, 0, 1);
  

  doneVectorInit = false;
}



void Acclaim::FourierBuffer::initVectors(int n){
  sumPower.resize(n, 0);
  hRays.resize(n, NULL);
  chiSquares.resize(n, 0);
  ndfs.resize(n, 0);  
  for(int freqBin=0; freqBin < n; freqBin++){
    TString name = TString::Format("hRayleigh");
    name += ant >= 0 && pol < AnitaPol::kNotAPol ? TString::Format("_%d_%d", ant, pol) : "";
    name += TString::Format("_%d", freqBin);
    hRays.at(freqBin) = new Acclaim::RayleighHist(this, name, name);
  }
  doneVectorInit = true;
}



Acclaim::FourierBuffer::~FourierBuffer(){
  if(spectrum){
    delete spectrum;
    spectrum = NULL;
  }
  
  for(unsigned i=0; i < hRays.size(); i++){
    delete hRays.at(i);
    hRays.at(i) = NULL;
  }
}


TGraphAligned* Acclaim::FourierBuffer::getReducedChiSquaresOfRayelighDistributions() const{
  TGraphAligned* gr = NULL;

  if(chiSquares.size() > 2){
    gr = new TGraphAligned(chiSquares.size() - 1);   
    for(unsigned i=1; i < chiSquares.size(); i++){
      double val = (i < chiSquares.size() - 1 && ndfs[i] > 0) ? chiSquares[i]/ndfs[i] : 0;
      gr->SetPoint(i, df*i, val);
    }
  }
  return gr;
}



size_t Acclaim::FourierBuffer::add(const RawAnitaHeader* header, const AnalysisWaveform* wave){

  // get the power
  const TGraphAligned* grPower = wave->power();


  // do dynamic initialization and sanity checking
  if(!doneVectorInit){
    initVectors(grPower->GetN());
  }
  if(grPower->GetN() != (int)sumPower.size()){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unexpected waveform size" << std::endl;
  }
  if(df <= 0){
    df = grPower->GetX()[1] - grPower->GetX()[0];
    for(size_t freqBin=0; freqBin < hRays.size(); freqBin++){
      hRays.at(freqBin)->freqMHz = df*1e3*freqBin;
    }
  }
  
  eventNumbers.push_back(header->eventNumber);
  runs.push_back(header->run);

  // for old compilers, push back copy of empty vector for speed.
  // then get reference that vector in the list
  powerRingBuffer.push_back(std::vector<double>(0));
  std::vector<double>& freqVec = powerRingBuffer.back();
  freqVec.assign(grPower->GetY(), grPower->GetY()+grPower->GetN());

  
  // update sum of power
  for(int freqInd=0; freqInd < grPower->GetN(); freqInd++){
    sumPower.at(freqInd) += grPower->GetY()[freqInd];

    double amp = TMath::Sqrt(grPower->GetY()[freqInd]);
    hRays.at(freqInd)->add(amp);
  }
  
  removeOld();

  for(int freqInd=0; freqInd < grPower->GetN(); freqInd++){
    hRays.at(freqInd)->Eval(chiSquares[freqInd], ndfs[freqInd]);
  }
  
  return eventNumbers.size();
}




Int_t Acclaim::FourierBuffer::removeOld(){

  Int_t nPopped = 0;
  // if(realTimesNs.size() > 0){
  if(eventNumbers.size() > 0){    

    while((int)realTimesNs.size() > bufferSize){
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







TGraphAligned* Acclaim::FourierBuffer::getAvePowSpec_dB(int lastNEvents) const{

  TGraphAligned* gr = getAvePowSpec(lastNEvents);
  gr->dBize();
  return gr;

}


TGraphAligned* Acclaim::FourierBuffer::getAvePowSpec(int lastNEvents) const{

  // set default value (which is the whole range of the fourier buffer)

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



TGraphAligned* Acclaim::FourierBuffer::getBackground_dB(int lastNEvents) const{
  
  TGraphAligned* gr = getBackground(lastNEvents);
  gr->dBize();
  return gr;
}


TGraphAligned* Acclaim::FourierBuffer::getBackground(int lastNEvents) const{

  TGraphAligned* gr = getAvePowSpec(lastNEvents);

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

