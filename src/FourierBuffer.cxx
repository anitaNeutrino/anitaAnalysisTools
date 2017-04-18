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
#include "RootTools.h"
#include "TROOT.h"

Acclaim::FourierBuffer::FourierBuffer(Int_t theBufferSize) : doneVectorInit(false), eventsInBuffer(0), fMinFitFreq(0.15), fMaxFitFreq(1.3)
{
  
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
  for(int ant=0; ant < NUM_SEAVEYS; ant++){
    summaryPads[ant] = NULL;
  }


  const char* rayleighFuncText = "([0]*x/([1]*[1]))*exp(-x*x/(2*[1]*[1]))";
  // const char* riceFuncText = "([0]*x/([1]*[1]))*exp(-(x*x+[2]*[2])/(2*[1]*[1]))*TMath::BesselI0([2]*x/([1]*[1]))";
  
  fRay = new TF1("fRay", rayleighFuncText, 0, 1, TF1::EAddToList::kNo);
}



void Acclaim::FourierBuffer::initVectors(int n, double df){

  
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      sumPowers[pol][ant].resize(n, 0);

      chiSquares[pol][ant].resize(n, 0);
      ndfs[pol][ant].resize(n, 0);
      fitAmplitudes[pol][ant].resize(n, 0);


      hRays[pol][ant].resize(n, NULL);      
      for(int freqBin=0; freqBin < n; freqBin++){
	double f = df*1e3*freqBin;
	TString name = TString::Format("hRayleigh_%d_%d_%d", ant, pol, freqBin);
	const char* polChar = pol == AnitaPol::kHorizontal ? "H" : "V";
	int phiName = (ant%NUM_PHI)+1;
	int ring = ant / NUM_PHI;
	const char* ringChar = ring == 0 ? "T" : ring == 1 ? "M" : "V";
	TString title = TString::Format("Distribution of amplitudes %d%s%s %4.1lf MHz", phiName, ringChar, polChar, f);
	title += ";Amplitude (mV/MHz?); Events per bin";
	hRays[pol][ant].at(freqBin) = new Acclaim::RayleighHist(this, name, title);
	hRays[pol][ant].at(freqBin)->SetDirectory(0); // hide from when .ls is done in MagicDisplay
	hRays[pol][ant].at(freqBin)->freqMHz = df*1e3*freqBin;	
      }

      // summary graphs for drawSummary
      grChiSquares[pol].push_back(TGraphFB(this, ant, pol, n));
      grNDFs[pol].push_back(TGraphFB(this, ant, pol, n));
      grAmplitudes[pol].push_back(TGraphFB(this, ant, pol, n));
      grReducedChiSquares[pol].push_back(TGraphFB(this, ant, pol, n));

      for(int freqBin=0; freqBin < n; freqBin++){
	double f = df*freqBin;
	grChiSquares[pol][ant].GetX()[freqBin] = f;
	grNDFs[pol][ant].GetX()[freqBin] = f;
	grAmplitudes[pol][ant].GetX()[freqBin] = f;
	grReducedChiSquares[pol][ant].GetX()[freqBin] = f;

	grChiSquares[pol][ant].GetY()[freqBin] = 0;
	grNDFs[pol][ant].GetY()[freqBin] = 0;
	grAmplitudes[pol][ant].GetY()[freqBin] = 0;
	grReducedChiSquares[pol][ant].GetY()[freqBin] = 0;
	
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

  for(int ant=0; ant < NUM_SEAVEYS; ant++){
    if(summaryPads[ant]){
      delete summaryPads[ant];
      summaryPads[ant] = NULL;
    }
  }
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
	df = grPower->GetX()[1] - grPower->GetX()[0];
	initVectors(grPower->GetN(), df);
      }
      if(grPower->GetN() != (int)sumPowers[pol][ant].size()){
	std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unexpected waveform size" << std::endl;
      }
  

      // for old compilers, push back copy of empty vector for speed.
      // then get reference that vector in the list
      powerRingBuffers[pol][ant].push_back(std::vector<double>(0));
      std::vector<double>& freqVec = powerRingBuffers[pol][ant].back();
      freqVec.assign(grPower->GetY(), grPower->GetY()+grPower->GetN());
  
      // update sum of power
      for(int freqInd=0; freqInd < grPower->GetN(); freqInd++){
	sumPowers[pol][ant].at(freqInd) += grPower->GetY()[freqInd];
	double f = df*freqInd;
	if(f >= fMinFitFreq && f < fMaxFitFreq){
	
	  double amp = TMath::Sqrt(grPower->GetY()[freqInd]);
	
	  bool updated = hRays[pol][ant].at(freqInd)->add(amp);
	  if(updated){
	    hRays[pol][ant].at(freqInd)->getRayleighFitParams(fitAmplitudes[pol][ant][freqInd], chiSquares[pol][ant][freqInd], ndfs[pol][ant][freqInd]);
	    // std::cout << fitAmplitudes[pol][ant][freqInd] << "\t" << chiSquares[pol][ant][freqInd] << "\t" << ndfs[pol][ant][freqInd] << std::endl;
	    grChiSquares[pol][ant].GetY()[freqInd] = chiSquares[pol][ant][freqInd];
	    if(ndfs[pol][ant][freqInd] > 0){
	      grReducedChiSquares[pol][ant].GetY()[freqInd] = chiSquares[pol][ant][freqInd]/ndfs[pol][ant][freqInd];
	    }
	    grNDFs[pol][ant].GetY()[freqInd] = ndfs[pol][ant][freqInd];
	    grAmplitudes[pol][ant].GetY()[freqInd] = fitAmplitudes[pol][ant][freqInd];
	  }
	}

	// // is there a more elegant way to do this?
	// TSeqCollection* cans = gROOT->GetListOfCanvases();
	// for(int i=0; i < cans->GetEntries(); i++){
	//   TCanvas* can = (TCanvas*) cans->At(i);
	//   TObject* objInCan = can->FindObject(hRays[pol][ant][freqInd]->GetName());
	//   if(objInCan){
	//     can->Modified();
	//     can->Update();
	//   }
	// }
      }
    }
  }
  eventsInBuffer++;
  
  removeOld();
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
    eventsInBuffer--;
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





void Acclaim::FourierBuffer::drawSummary(TPad* pad) const{

  pad->cd();

  // find the maximum and minimum values before plotting
  Double_t yMax = -DBL_MAX;
  Double_t yMin = DBL_MAX;  
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;  
    for(int ant=0; ant < NUM_SEAVEYS; ant++){

      if(ant == 4 && pol == AnitaPol::kHorizontal){ // skip alfa
	continue;
      }
      
      TGraphFB& gr = (TGraphFB&) grReducedChiSquares[pol][ant];

      for(int i=0; i < gr.GetN(); i++){
	if(gr.GetY()[i] > yMax){
	  yMax = gr.GetY()[i];
	}
	if(gr.GetY()[i] < yMin){
	  yMin = gr.GetY()[i];
	}	
      }
    }
  }

  // now do the plotting (setting the max/min)
  
  for(int ant=0; ant < NUM_SEAVEYS; ant++){


    // dynamically initialize the summary pads
    if(!summaryPads[ant]){
      int ring = ant/NUM_PHI;
      int phi = ant % NUM_PHI;

      // instead of 3 by 16
      // do 6 by 8
      int xBin = phi/2;
      int yBin = (2*ring + phi%2);
    
      summaryPads[ant] = RootTools::makeSubPad(pad, double(xBin)/8, double(5 - yBin)/6,
					       double(xBin+1)/8, double(6 - yBin)/6,
					       TString::Format("FourierBuffer%d", ant));
      summaryPads[ant]->SetRightMargin(0);
      summaryPads[ant]->SetLeftMargin(0.1 * (xBin == 0));
      summaryPads[ant]->SetTopMargin(0.1 * (yBin == 0));
      summaryPads[ant]->SetBottomMargin(0.1 * (yBin == 5));

      // avoid deletion?
      summaryPads[ant]->SetBit(kCanDelete);
      summaryPads[ant]->InvertBit(kCanDelete);      
    }
    else{
      // already have the pad, so append it to the boss pad
      pad->cd();
      summaryPads[ant]->AppendPad();
    }


    // do the actual plotting
    summaryPads[ant]->cd();
    for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
      
      // TGraphFB* gr = getReducedChiSquaresOfRayelighDistributions(ant, pol);
      TGraphFB& gr = (TGraphFB&) grReducedChiSquares[pol][ant];
      gr.SetEditable(kFALSE);
      
      const char* opt = pol == AnitaPol::kVertical ? "lsame" : "al";
      EColor lineCol = pol == AnitaPol::kHorizontal ? kBlue : kBlack;
      gr.SetLineColor(lineCol);
      gr.SetMaximum(yMax);
      gr.SetMinimum(yMin);
      
      gr.Draw(opt);
    }
  }
  
}












void Acclaim::TGraphFB::drawCopy() const{
  TCanvas* c1 = new TCanvas();
  c1->cd();
  TGraphFB* gr2 = (TGraphFB*) Clone();
  gr2->pol = pol;
  gr2->ant = ant;
  gr2->fb = fb;
  gr2->SetBit(kCanDelete);
  gr2->fDoubleClickOption = kDrawRayleigh;
  gr2->Draw();
  c1->Modified();
  c1->Update();

}
void Acclaim::TGraphFB::ExecuteEvent(int event, int x, int y){

  if(event == kButton1Double){    
    if(fDoubleClickOption==kDrawRayleigh){
      drawRayleighHistNearMouse(x, y);
    }
    else if(fDoubleClickOption==kDrawCopy){
      drawCopy();
    }
  }  
  TGraphAligned::ExecuteEvent(event, x, y);
}


void Acclaim::TGraphFB::drawRayleighHistNearMouse(int x, int y) const{

  int selectedPoint = -1;
  double minD2 = 1e9;
  for(int i=0; i < GetN(); i++){
    Int_t pixelX = gPad->XtoAbsPixel(GetX()[i]);
    Int_t pixelY = gPad->YtoAbsPixel(GetY()[i]);

    const int closeDist = 4; // number of pixels that means you're close to a point
    if(TMath::Abs(pixelX - x) < closeDist && TMath::Abs(pixelY - y) < closeDist){
      double d2 = (pixelX - x)*(pixelX - x) + (pixelY - y)*(pixelY - y);
      
      if(d2 < minD2){
	minD2 = d2;
	selectedPoint = i;
      }
    }
  }
      
  double df = GetX()[1] - GetX()[0];
  if(selectedPoint > -1 && df > 0){
    TCanvas* c1 = new TCanvas();
    c1->cd();
    // need to pretend it's non-const
    // don't delete this!
    RayleighHist* h = (RayleighHist*) fb->getRayleighDistribution(ant, pol, selectedPoint);
    h->SetLineColor(GetLineColor());
    h->Draw();
    c1->Modified();
    c1->Update();
  }

}
