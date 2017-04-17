#include "FourierBuffer.h"
#include "TMath.h"
#include "FFTWComplex.h"
#include "TSpectrum.h"
#include "RayleighHist.h"
#include "TF1.h"
#include "FilteredAnitaEvent.h"
#include "TVirtualPad.h"
#include "TMath.h"
#include "TCanvas.h"
#include "RawAnitaHeader.h"


Acclaim::FourierBuffer::FourierBuffer(Int_t theBufferSize){
  
  // timeScale = timeScaleSeconds;
  bufferSize = theBufferSize <= 0 ? 1000 : theBufferSize;
  df = -1;
  fDrawFreqBin = 100; //26;

  // will initialize this dynamically to get around this no-copy-constructor bullshit

  for(int pol=0; pol < AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      spectrums[pol][ant] = NULL;
    }
  }

  const char* rayleighFuncText = "([0]*x/([1]*[1]))*exp(-x*x/(2*[1]*[1]))";
  // const char* riceFuncText = "([0]*x/([1]*[1]))*exp(-(x*x+[2]*[2])/(2*[1]*[1]))*TMath::BesselI0([2]*x/([1]*[1]))";
  
  TString funcName = TString::Format("fRay");
  fRay = new TF1(funcName, rayleighFuncText, 0, 1);
  
  doneVectorInit = false;
}



void Acclaim::FourierBuffer::initVectors(int n){

  
  for(int pol=0; pol < AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      sumPowers[pol][ant].resize(n, 0);
      hRays[pol][ant].resize(n, NULL);
      chiSquares[pol][ant].resize(n, 0);
      ndfs[pol][ant].resize(n, 0);
      for(int freqBin=0; freqBin < n; freqBin++){
	TString name = TString::Format("hRayleigh_%d_%d_%d", ant, pol, freqBin);
	hRays[pol][ant].at(freqBin) = new Acclaim::RayleighHist(this, name, name);
      }
    }
  }
  doneVectorInit = true;  
}



Acclaim::FourierBuffer::~FourierBuffer(){


  for(int pol=0; pol < AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      if(spectrums[pol][ant]){
	delete spectrums[pol][ant];
	spectrums[pol][ant] = NULL;
      }
  
      for(unsigned i=0; i < hRays[pol][ant].size(); i++){
	delete hRays[pol][ant].at(i);
	hRays[pol][ant].at(i) = NULL;
      }
    }
  }
}


Acclaim::TGraphFB* Acclaim::FourierBuffer::getReducedChiSquaresOfRayelighDistributions(Int_t ant, AnitaPol::AnitaPol_t pol) const{
  TGraphFB* gr = NULL;
  
  if(chiSquares[pol][ant].size() > 2){
    gr = new TGraphFB(this, ant, pol, chiSquares[pol][ant].size() - 1);   
    for(unsigned i=1; i < chiSquares[pol][ant].size(); i++){
      double val = (i < chiSquares[pol][ant].size() - 1 && ndfs[i] > 0) ? chiSquares[pol][ant][i]/ndfs[pol][ant][i] : 0;
      gr->SetPoint(i, df*i, val);
    }
  }
  return gr;
}



// size_t Acclaim::FourierBuffer::add(const RawAnitaHeader* header, const AnalysisWaveform* wave){
size_t Acclaim::FourierBuffer::add(const FilteredAnitaEvent* fEv){  

  const RawAnitaHeader* header = fEv->getHeader();
  eventNumbers.push_back(header->eventNumber);
  runs.push_back(header->run);

  // get the power


  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(int ant=0; ant < NUM_SEAVEYS; ant++){

      const AnalysisWaveform* wave = fEv->getFilteredGraph(ant, pol);
      const TGraphAligned* grPower = wave->power();


      // do dynamic initialization and sanity checking
      if(!doneVectorInit){
	initVectors(grPower->GetN());
      }
      if(grPower->GetN() != (int)sumPowers[pol][ant].size()){
	std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unexpected waveform size" << std::endl;
      }
      if(df <= 0){
	df = grPower->GetX()[1] - grPower->GetX()[0];
	for(size_t freqBin=0; freqBin < hRays[pol][ant].size(); freqBin++){
	  hRays[pol][ant].at(freqBin)->freqMHz = df*1e3*freqBin;
	}
      }
  

      // for old compilers, push back copy of empty vector for speed.
      // then get reference that vector in the list
      powerRingBuffers[pol][ant].push_back(std::vector<double>(0));
      std::vector<double>& freqVec = powerRingBuffers[pol][ant].back();
      freqVec.assign(grPower->GetY(), grPower->GetY()+grPower->GetN());

  
      // update sum of power
      for(int freqInd=0; freqInd < grPower->GetN(); freqInd++){
	sumPowers[pol][ant].at(freqInd) += grPower->GetY()[freqInd];

	double amp = TMath::Sqrt(grPower->GetY()[freqInd]);
	hRays[pol][ant].at(freqInd)->add(amp);
      }
    }
  }
  
  removeOld();

  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(int ant=0; ant < NUM_SEAVEYS; ant++){      
      for(int freqInd=0; freqInd < (int)sumPowers[pol][ant].size(); freqInd++){
	hRays[pol][ant].at(freqInd)->Eval(chiSquares[pol][ant][freqInd], ndfs[pol][ant][freqInd]);
      }
    }
  }
  return eventNumbers.size();
}




Int_t Acclaim::FourierBuffer::removeOld(){

  Int_t nPopped = 0;
  while((int)eventNumbers.size() > bufferSize){
    eventNumbers.pop_front();
    runs.pop_front();

    for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
      for(int ant=0; ant < NUM_SEAVEYS; ant++){      
      
	std::vector<double>& removeThisPower = powerRingBuffers[pol][ant].front();
	for(unsigned int freqInd=0; freqInd < removeThisPower.size(); freqInd++){
	  sumPowers[pol][ant].at(freqInd) -= removeThisPower.at(freqInd);
	}
	powerRingBuffers[pol][ant].pop_front();
      }
    }
    nPopped++;
  }
  return nPopped;
}







Acclaim::TGraphFB* Acclaim::FourierBuffer::getAvePowSpec_dB(Int_t ant, AnitaPol::AnitaPol_t pol, int lastNEvents) const{

  TGraphFB* gr = getAvePowSpec(ant, pol, lastNEvents);
  gr->dBize();
  return gr;

}


Acclaim::TGraphFB* Acclaim::FourierBuffer::getAvePowSpec(Int_t ant, AnitaPol::AnitaPol_t pol, int lastNEvents) const{

  // set default value (which is the whole range of the fourier buffer)

  int nEvents = eventNumbers.size();

  TGraphFB* gr = new TGraphFB(this, ant, pol, sumPowers[pol][ant].size());
  
  const char* polName = AnitaPol::kHorizontal ? "HPol" : "VPol";
  TString name = TString::Format("grAvePowSpec_%d_%s_%u_%u", ant, polName, eventNumbers.front(), eventNumbers.back());
  TString title = TString::Format("Average Power Spectrum %d %s event %u - %u", ant, polName, eventNumbers.front(), eventNumbers.back());        
  gr->SetNameTitle(name, title);
  

  for(unsigned freqInd=0; freqInd < sumPowers[pol][ant].size(); freqInd++){
    gr->SetPoint(freqInd, freqInd*df, sumPowers[pol][ant].at(freqInd)/nEvents);
  }

  return gr;
}



Acclaim::TGraphFB* Acclaim::FourierBuffer::getBackground_dB(Int_t ant, AnitaPol::AnitaPol_t pol, int lastNEvents) const{
  
  TGraphFB* gr = getBackground(ant, pol, lastNEvents);
  gr->dBize();
  return gr;
}


Acclaim::TGraphFB* Acclaim::FourierBuffer::getBackground(Int_t ant, AnitaPol::AnitaPol_t pol, int lastNEvents) const{

  TGraphFB* gr = getAvePowSpec(ant, pol, lastNEvents);

  if(!spectrums[pol][ant]){
    spectrums[pol][ant] = new TSpectrum();
  }

  // at some point this changed from doubles to floats, not added a new method, changed...
  // git blames Lorenzo Moneta... 
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)

  spectrums[pol][ant]->Background(gr->GetY(),gr->GetN(),
				  6,TSpectrum::kBackDecreasingWindow,
				  TSpectrum::kBackOrder2,kFALSE,
				  TSpectrum::kBackSmoothing3,kFALSE);
#else
  std::vector<float> tempFloats(gr->GetN(), 0);
  for(int i=0; i < gr->GetN(); i++){
    tempFloats[i] = gr->GetY()[i];
  }
  spectrums[pol][ant]->Background(&tempFloats[0],gr->GetN(),
				  6,TSpectrum::kBackDecreasingWindow,
				  TSpectrum::kBackOrder2,kFALSE,
				  TSpectrum::kBackSmoothing3,kFALSE);
  for(int i=0; i < gr->GetN(); i++){
    gr->GetY()[i] = tempFloats[i];
  }
#endif
  

  return gr;
}














void Acclaim::TGraphFB::ExecuteEvent(int event, int x, int y){
  if(event == kButton1Double){
    Double_t freq = gPad->AbsPixeltoX(x);
    double df = GetX()[1] - GetX()[0];
    if(df > 0){
      int freqBin = TMath::Nint(freq/df);
      TCanvas* c1 = new TCanvas();
      c1->cd();
      // need to pretend it's non-const
      // don't delete this!
      RayleighHist* h = (RayleighHist*) fb->getRayleighDistribution(ant, pol, freqBin);
      h->Draw();
      c1->Modified();
      c1->Update();
    }
  }
  TGraphAligned::ExecuteEvent(event, x, y);
}
