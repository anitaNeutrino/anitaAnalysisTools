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
#include "AnitaDataset.h"
#include "FilterStrategy.h"
#include "ProgressBar.h"
#include "QualityCut.h"
#include "AcclaimFilters.h"
#include "FFTtools.h"
#include "AcclaimFilters.h"
#define IS_ROOT_6 (ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0))
#define IS_ROOT_6_06_08 (ROOT_VERSION_CODE >= ROOT_VERSION(6,6,8))


Acclaim::FourierBuffer::FourierBuffer(Int_t theBufferSize ) :
    doneVectorInit(false), fCurrentlyLoadingHistory(false), fForceLoadHistory(false), eventsInBuffer(0), fMinFitFreq(Filters::Bands::anitaHighPassGHz), fMaxFitFreq(Filters::Bands::anitaLowPassGHz), fNumSkipped(0){
  
  bufferSize = theBufferSize <= 0 ? 1000 : theBufferSize;  
  df = -1;
  
  // will initialize this dynamically to get around this no-copy-constructor bullshit
  fSpectrum = NULL;
  const char* rayleighFuncText = "([0]*x/([1]*[1]))*exp(-x*x/(2*[1]*[1]))";
  for(int ant=0; ant < NUM_SEAVEYS; ant++){
    summaryPads[ant] = NULL;
  }

#if IS_ROOT_6_06_08
  fRay = new TF1("fRay", rayleighFuncText, 0, 1e4, TF1::EAddToList::kNo);
#else
  fRay = new TF1("fRay", rayleighFuncText, 0, 1e4);
#endif
}


void Acclaim::FourierBuffer::initGraphAndVector(std::vector<double> vec[][NUM_SEAVEYS],
						std::vector<TGraphFB>* gr,
						int n, double df, double defaultVal){
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;

    if(gr){
      gr[pol].reserve(NUM_SEAVEYS);
    }
    
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      if(vec){
	vec[pol][ant].resize(n, defaultVal);
      }

      if(gr){
	gr[pol].push_back(TGraphFB(this, ant, pol, n));
	for(int freqBin=0; freqBin < n; freqBin++){
	  double f = df*freqBin;
	  gr[pol][ant].GetX()[freqBin] = f;
	  gr[pol][ant].GetY()[freqBin] = defaultVal;
	}
      }
      
    }    
  }
}



void Acclaim::FourierBuffer::initVectors(int n, double df){

  initGraphAndVector(spectrumAmplitudes, grSpectrumAmplitudes, n, df, 0);
  initGraphAndVector(chiSquares, grChiSquares, n, df, 0);
  initGraphAndVector(chiSquaresRelativeToSpectrum, grChiSquaresRelativeToSpectrum, n, df, 0);
  initGraphAndVector(NULL, grReducedChiSquares, n, df, 0);
  initGraphAndVector(NULL, grReducedChiSquaresRelativeToSpectrum, n, df, 0);
  initGraphAndVector(probs, grProbs, n, df, 0);
  initGraphAndVector(fitAmplitudes, grAmplitudes, n, df, 0);  
  initGraphAndVector(sumPowers, NULL, n, df, 0);
  initGraphAndVector(lastAmps, grLastAmps, n, df, 0);
  initGraphAndVector(NULL, grNDFs, n, df, 0);
  initGraphAndVector(chanChiSquares, NULL, n, df, 0);  
  
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      ndfs[pol][ant].resize(n, 0);
          
      hRays[pol][ant].resize(n, NULL);
      for(int freqBin=0; freqBin < n; freqBin++){
	double f = df*1e3*freqBin;
	TString name = TString::Format("hRayleigh_%d_%d_%d", ant, pol, freqBin);
	const char* polChar = pol == AnitaPol::kHorizontal ? "H" : "V";
	int phiName = (ant%NUM_PHI)+1;
	int ring = ant / NUM_PHI;
	const char* ringChar = ring == 0 ? "T" : ring == 1 ? "M" : "B";
	TString title = TString::Format("Distribution of amplitudes %d%s%s %4.1lf MHz", phiName, ringChar, polChar, f);
	title += ";Amplitude (mV/MHz?); Events per bin";
	hRays[pol][ant].at(freqBin) = new Acclaim::RayleighHist(this, name, title);
	hRays[pol][ant].at(freqBin)->SetDirectory(0); // hide from when .ls is done in MagicDisplay
	hRays[pol][ant].at(freqBin)->freqMHz = df*1e3*freqBin;	
      }
    }
  }

  // need to do this at the end so the vector doesn't reallocate the memory (also reserved it too now so doubly fixed)
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(int ant=0; ant < NUM_SEAVEYS; ant++){

      // get them to know about eachother
      
      
      std::vector<Acclaim::TGraphFB*> grAmpVec;
      grAmpVec.push_back(&grLastAmps[pol][ant]);
      grAmpVec.push_back(&grAmplitudes[pol][ant]);
      grAmpVec.push_back(&grSpectrumAmplitudes[pol][ant]);
      Acclaim::TGraphFB::setDrawingDependencies(grAmpVec);      

      std::vector<Acclaim::TGraphFB*> grChiSqVec;
      grChiSqVec.push_back(&grReducedChiSquares[pol][ant]);
      grChiSqVec.push_back(&grReducedChiSquaresRelativeToSpectrum[pol][ant]);
      
      Acclaim::TGraphFB::setDrawingDependencies(grChiSqVec);      
    }
  }

  doneVectorInit = true;  
}



Acclaim::FourierBuffer::~FourierBuffer(){

  if(fSpectrum){
    delete fSpectrum;
    fSpectrum = NULL;
  }

  for(int pol=0; pol < AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){

      for(unsigned i=0; i < hRays[pol][ant].size(); i++){
        if(hRays[pol][ant].at(i)){
          delete hRays[pol][ant].at(i);          
          hRays[pol][ant].at(i) = NULL;
        }
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




size_t Acclaim::FourierBuffer::add(const FilteredAnitaEvent* fEv){



  // do dynamic initialization, need this first in case we return from this function before the end
  if(!doneVectorInit){
    // arbitrarily choose 1TH    
    const AnalysisWaveform* wave = fEv->getFilteredGraph(0, AnitaPol::kHorizontal);
    const TGraphAligned* grPower = wave->power();
    df = grPower->GetX()[1] - grPower->GetX()[0];
    initVectors(grPower->GetN(), df);
  }


  
  // don't add nasty events that have crazy amounts of power or other pathologies
  Bool_t isGoodEvent = QualityCut::applyAll(fEv->getUsefulAnitaEvent());
  if(!isGoodEvent){
    fNumSkipped++;
    // then do nothing
    return eventNumbers.size();
  }

  // handle use in MagicDisplay where refreshing event display on same event...
  // things still won't be exactly right if you do anything other than repeat the same event
  // or move forwards sequentially...
  if(eventNumbers.size() > 0 && fEv->getHeader()->eventNumber == eventNumbers.back()){
    // However, if we've asked to load the history of an event, we've probably just jumped to that event in MagicDisplay
    // I'll create an exception to this rule if fForceLoadHistory is true
    if(!fForceLoadHistory){
      return eventNumbers.size();
    }
  }  

   
  if(!fCurrentlyLoadingHistory && fForceLoadHistory){
    automagicallyLoadHistory(fEv);
    fForceLoadHistory = false;
  }
  
  const RawAnitaHeader* header = fEv->getHeader();
  eventNumbers.push_back(header->eventNumber);
  runs.push_back(header->run);

  bool anyUpdated = false; // should all get updated at the same time, but whatever...

  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(int ant=0; ant < NUM_SEAVEYS; ant++){

      const AnalysisWaveform* wave = fEv->getFilteredGraph(ant, pol);
      const TGraphAligned* grPower = wave->power();

      // quick sanity check
      if(grPower->GetN() != (int)sumPowers[pol][ant].size()){
	std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unexpected waveform size" << std::endl;
      }
  
      // for old compilers, push back copy of empty vector for speed.
      // then get reference that vector in the list
      powerRingBuffers[pol][ant].push_back(std::vector<double>(0));
      std::vector<double>& freqVec = powerRingBuffers[pol][ant].back();
      freqVec.assign(grPower->GetY(), grPower->GetY()+grPower->GetN()); 

      // update sum of power
      //g_power.SetPoint(i, i * df, the_fft[i].getAbsSq()*2/fft_len/50/1000);	  

      // need *50*1000/2 = 25000
      const double cosminPowerConversionFactor = 25000*grPower->GetN();

      for(int freqInd=0; freqInd < grPower->GetN(); freqInd++){
	sumPowers[pol][ant].at(freqInd) += grPower->GetY()[freqInd];
	double f = df*freqInd;
	if(f >= fMinFitFreq && f < fMaxFitFreq){

	  double amp = TMath::Sqrt(grPower->GetY()[freqInd]*cosminPowerConversionFactor);
	  
	  bool updated = hRays[pol][ant].at(freqInd)->add(amp);
	  grLastAmps[pol][ant].GetY()[freqInd] = amp;
          lastAmps[pol][ant][freqInd] = amp;
	  
	  if(updated){
	    anyUpdated = true;
	    hRays[pol][ant].at(freqInd)->getRayleighFitParams(fitAmplitudes[pol][ant][freqInd],
							      chiSquares[pol][ant][freqInd],
							      ndfs[pol][ant][freqInd]);
	    
	    grChiSquares[pol][ant].GetY()[freqInd] = chiSquares[pol][ant][freqInd];
	    if(ndfs[pol][ant][freqInd] > 0){
	      
	      grReducedChiSquares[pol][ant].GetY()[freqInd] = chiSquares[pol][ant][freqInd]/ndfs[pol][ant][freqInd];

	      grNDFs[pol][ant].GetY()[freqInd] = ndfs[pol][ant][freqInd];
	      grAmplitudes[pol][ant].GetY()[freqInd] = fitAmplitudes[pol][ant][freqInd];	    
	    }
	  }
	  double distAmp = 0;
	  if(spectrumAmplitudes[pol][ant][freqInd] > 0){
	    distAmp = spectrumAmplitudes[pol][ant][freqInd];
	    chiSquaresRelativeToSpectrum[pol][ant].at(freqInd) = hRays[pol][ant].at(freqInd)->getRayleighChiSquare(&spectrumAmplitudes[pol][ant][freqInd]);
	    
	  }
	  else{
	    distAmp = hRays[pol][ant].at(freqInd)->getAmplitude();
	    chiSquaresRelativeToSpectrum[pol][ant].at(freqInd) = chiSquares[pol][ant][freqInd];
	  }
	  grChiSquaresRelativeToSpectrum[pol][ant].GetY()[freqInd] = chiSquaresRelativeToSpectrum[pol][ant].at(freqInd);

	  if(ndfs[pol][ant][freqInd] > 0){
	    grReducedChiSquaresRelativeToSpectrum[pol][ant].GetY()[freqInd] = chiSquaresRelativeToSpectrum[pol][ant].at(freqInd)/ndfs[pol][ant][freqInd];
	  }

          // if(distAmp > 0){
          //   chanChisquare[pol][ant] += (amp*amp)/(distAmp*distAmp);
          //   chanNdf[pol][ant] += 2;
          // }


	  double prob = hRays[pol][ant].at(freqInd)->getOneMinusCDF(amp, distAmp);
	  
	  double probVal = probVal = prob > 0 ? TMath::Log10(prob) : 0;

	  // if(!TMath::Finite(probVal)){
	  //   std::cout << probVal << "\t" << prob << std::endl;
	  // }
	  
	  probs[pol][ant].at(freqInd) = probVal;
	  grProbs[pol][ant].GetY()[freqInd] = probVal;
	}
      }

      double sumChanChiSquares = 0;
      int n=0;
      double sumSquaredChanChiSquares = 0;
      for(int freqInd=0; freqInd < grPower->GetN(); freqInd++){
        double amp = lastAmps[pol][ant][freqInd];
        if(amp > 0){
          double realAmp = fitAmplitudes[pol][ant][freqInd];
          chanChiSquares[pol][ant][freqInd] = (amp*amp)/(realAmp*realAmp);
          sumChanChiSquares += chanChiSquares[pol][ant][freqInd];
          sumSquaredChanChiSquares += chanChiSquares[pol][ant][freqInd]*chanChiSquares[pol][ant][freqInd];
          n++;

          // std::cout << "(" << amp << ", " << realAmp << ")" << std::endl;
          // std::cout << "(" << chanChiSquares[pol][ant][freqInd] << ")" << std::endl;          
        }
      }
      // std::cout << std::endl;

      if(n > 0){
        meanChanChiSquare[pol][ant] = sumChanChiSquares/n;
        varChanChiSquare[pol][ant] = sumSquaredChanChiSquares/n - meanChanChiSquare[pol][ant]*meanChanChiSquare[pol][ant];
      }
      else{
        meanChanChiSquare[pol][ant] = 0;
        varChanChiSquare[pol][ant] = 0;
      }

      if(!fCurrentlyLoadingHistory){
        // std::cout << pol << "\t" << ant << "\t" << sumChanChiSquares << "\t" << sumSquaredChanChiSquares << "\t" << meanChanChiSquare[pol][ant] << "\t" << varChanChiSquare[pol][ant] << "\t" << n << std::endl;
        // std::cout << pol << "\t" << ant << "\t" << meanChanChiSquare[pol][ant] << "\t" << varChanChiSquare[pol][ant] << "\t" << n << std::endl;
      }
      // if(!fCurrentlyLoadingHistory){
      //   std::cout << pol << "\t" << ant << "\t" << hChanChiSquares[pol][ant]->GetMean() << "\t" << hChanChiSquares[pol][ant]->GetRMS() << std::endl;
      // }
      // if(!fCurrentlyLoadingHistory){
      //   std::cout  << pol << "\t" << ant << "\t" << chanChisquare[pol][ant] << "\t" << chanNdf[pol][ant] << std::endl;
      // }
    }
  }


  if(anyUpdated){
    int firstSpecInd = 0; 

    for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
      for(int ant=0; ant < NUM_SEAVEYS; ant++){

	int lastSpecInd = int((isAlfaBandpassed(ant,pol) ? Filters::Bands::alfaLowPassGHz : fMaxFitFreq)/df);
	int nf = fitAmplitudes[pol][ant].size();

	// copy the fitted amplitudes into the background array
	// (before applying the TSpectrum background thing)
	for(int freqInd = 0; freqInd < nf; freqInd++){
	  spectrumAmplitudes[pol][ant][freqInd] = fitAmplitudes[pol][ant][freqInd];
	}

	// do the spectral analysis of the amplitudes
	getBackgroundSpectrum(&spectrumAmplitudes[pol][ant][firstSpecInd],
			      lastSpecInd - firstSpecInd);

	// update summary graphs	
	for(int freqInd = 0; freqInd < nf; freqInd++){
	  grSpectrumAmplitudes[pol][ant].GetY()[freqInd] = spectrumAmplitudes[pol][ant][freqInd];
	}
      }
    }
  }


  // if(!fCurrentlyLoadingHistory){
  //   for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
  //     AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
  //     for(int ant=0; ant < NUM_SEAVEYS; ant++){
  // 	std::cout << pol << "\t" << ant << "\t" << eventPower[pol][ant] << "\t" << expectedThermalPower[pol][ant] << std::endl;
  //     }
  //   }
  // }  
  eventsInBuffer++;
  
  removeOld();
  return eventNumbers.size();
}



TH1D* Acclaim::FourierBuffer::makeChanChisquareHist(AnitaPol::AnitaPol_t pol, Int_t ant, TPad* pad, const char* drawOpt) const{

  TString hName = TString::Format("hChanChiSquare_%d_%d", pol, ant);  
  TH1D* h = new TH1D(hName, hName, 8, 0, 24);
  // std::cout << "filling " << pol << "\t" << ant << ":";
  for(unsigned freqInd=0; freqInd < sumPowers[pol][ant].size(); freqInd++){
    h->Fill(chanChiSquares[pol][ant][freqInd]);
  }
  // std::cout << std::endl;

  if(!pad){
    new TCanvas();    
  }
  else{
    pad->cd();
  }
  h->Draw(drawOpt);
  return h;
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



void Acclaim::FourierBuffer::automagicallyLoadHistory(const FilteredAnitaEvent* fEv){

  fCurrentlyLoadingHistory = true;
  bool loadedHistory = false;

  const int safetyMargin = 100; // we don't add self triggered blasts or SURF saturated events, but I don't know how many there are before I try
  const int numEventsToLoopOver = bufferSize + safetyMargin; // so do a little more than the minimum

  Int_t run = fEv->getHeader()->run;
  UInt_t eventNumber = fEv->getHeader()->eventNumber;

  // if this FourierBuffer is inside an operation in the filter, then we can just process all the old events
  // otherwise we need to manually add them to this FourierBuffer
  // I imagine that needToAddManually will basically always be false...
  // bool needToAddManually = true;
  // const FilterStrategy* fs = fEv->getStrategy();
  // for(unsigned i=0; i < fs->nOperations(); i++){
  //   const Filters::RayleighMonitor* rm = dynamic_cast<const Filters::RayleighMonitor*>(fs->getOperation(i));
  //   if(rm){
  //     FourierBuffer fb = rm->getFourierBuffer();

  //     std::cout << fb.getAddress() << "\t" << this << std::endl;
      
  //     if(this == fb.getAddress()){
  // 	needToAddManually = false;
  // 	break;
  //     }
  //   }
  // }
  
  bool needToAddManually = false;

  std::cerr << "Info in " << __PRETTY_FUNCTION__ << ": Loading history for FourierBuffer, will loop over the last " << bufferSize << " events plus a few." << std::endl;
  // std::cerr << "In " << __PRETTY_FUNCTION__ << ", I think I am here " << this << std::endl;
  
  while(!loadedHistory){

    std::vector<AnitaDataset*> ds;
    std::vector<Int_t> firstEntries;
    std::vector<Int_t> lastEntries;
    Long64_t eventsInQueue = 0;
    int thisRun = run;
    while(eventsInQueue < numEventsToLoopOver){
      ds.push_back(new AnitaDataset(thisRun));
      
      AnitaDataset* d = ds.back();

      d->first();
      UInt_t firstEventThisRun = d->header()->eventNumber;

      d->last();
      UInt_t lastEventThisRun = d->header()->eventNumber;

      UInt_t maxPossEventNumber = run == thisRun ? eventNumber : lastEventThisRun;
      
      Int_t potentialEvents = maxPossEventNumber - firstEventThisRun;

      Int_t eventsThisRun = potentialEvents > numEventsToLoopOver ? numEventsToLoopOver : potentialEvents;

      firstEntries.push_back(potentialEvents - eventsThisRun);
      lastEntries.push_back(potentialEvents);
      eventsInQueue += lastEntries.back() - firstEntries.back();

      thisRun--;
    }

    // now we have the queue of events, iterate backwards (so run number is increasing), so they go in the buffer in the proper order...

    Long64_t eventsDone = 0;
    ProgressBar p(eventsInQueue);
    for(int i = (int)ds.size() - 1; i >= 0; i--){    

      AnitaDataset* d = ds.at(i);

      for(int entry = firstEntries.at(i); entry < lastEntries.at(i); entry++){

	int thisEntry = d->getEntry(entry);

	// std::cout << thisEntry << "\t" << entry << "\t" << d->header()->eventNumber << "\t" << std::endl;
      
	if(thisEntry > 0){
	
	  FilteredAnitaEvent fEv2(d->useful(), (FilterStrategy*) fEv->getStrategy(), d->gps(), d->header(), false);
	  if(needToAddManually){
	    add(&fEv2);
	  }
	}
	p.inc(eventsDone, eventsInQueue);
	eventsDone++;
      }

      delete d;
    }

    if(eventNumbers.size() != (unsigned) bufferSize){
      std::cerr << std::endl;
      std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", tried to load history, but only got "
		<< eventNumbers.size() << " events when I have a buffer size " << bufferSize << std::endl;
    }

    loadedHistory = true;
  }
  fCurrentlyLoadingHistory = false;
  
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

  getBackgroundSpectrum(gr->GetY(), gr->GetN());

  return gr;
}


// wrapper for spectrum since it's interface changes in a fucking stupid way between ROOT versions
void Acclaim::FourierBuffer::getBackgroundSpectrum(double* y, int n) const{

  if(!fSpectrum){
    fSpectrum = new TSpectrum();
  }

#if IS_ROOT_6
  double* tempPtr = y;
#else
  std::vector<float> tempFloats(n, 0);
  for(int i=0; i < n; i++){
    tempFloats[i] = y[i];
  }
  float* tempPtr = &tempFloats[0];
#endif

  const int numIter = 4;
  fSpectrum->Background(tempPtr,n,
			numIter,TSpectrum::kBackDecreasingWindow,
			TSpectrum::kBackOrder4,kFALSE,
			TSpectrum::kBackSmoothing3,kFALSE);
  
#if !(IS_ROOT_6)
  for(int i=0; i < n; i++){
    y[i] = tempFloats[i];
  }
#endif
}






void Acclaim::FourierBuffer::drawSummary(TPad* pad, SummaryOption_t summaryOpt) const{

  if(summaryOpt == None){
    return;
  }

  pad->cd();

  // find the maximum and minimum values before plotting
  Double_t yMax = -DBL_MAX;
  Double_t yMin = DBL_MAX;
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(int ant=0; ant < NUM_SEAVEYS; ant++){

      if(isAlfaBandpassed(ant, pol)){
	continue;
      }
      
      // TGraphFB& gr = (TGraphFB&) grReducedChiSquares[pol][ant]; 
      const TGraphFB* gr = getSelectedGraphForSummary(summaryOpt, ant, pol);
      // TGraphFB& gr = (TGraphFB&) grProbs[pol][ant];

      for(int i=0; i < gr->GetN(); i++){
	double x = gr->GetX()[i];
	if(x > fMinSpecFreq && x < fMaxFitFreq){
	  if(gr->GetY()[i] > yMax){
	    yMax = gr->GetY()[i];
	  }
	  if(gr->GetY()[i] < yMin){
	    yMin = gr->GetY()[i];
	  }
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
    summaryPads[ant]->Clear();
    
    summaryPads[ant]->cd();
    
    if(summaryOpt == FourierBuffer::ReducedChisquare){

      TH1D* hH = makeChanChisquareHist(AnitaPol::kHorizontal, ant, summaryPads[ant], "e");
      TH1D* hV = makeChanChisquareHist(AnitaPol::kVertical, ant, summaryPads[ant], "esame");
      hV->SetBit(kCanDelete);
      hH->SetBit(kCanDelete);
      hH->SetLineColor(kBlue);
      hV->SetLineColor(kBlack);
      gPad->SetLogy(1);
      
      continue;
      yMax = 10;
      yMin = 0;
      // gPad->SetLogy(1);
    }
    else{
      gPad->SetLogy(0);
    }

    summaryPads[ant]->cd();
    
    if(summaryOpt == FourierBuffer::Prob){
      yMax = 0; // otherwise some events get a little too crazy..      
      yMin = -10; // otherwise some events get a little too crazy..
    }
    
    for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
      
      TGraphFB* gr = getSelectedGraphForSummary(summaryOpt, ant, pol);      
      gr->SetEditable(kFALSE);
      
      const char* opt = pol == AnitaPol::kVertical ? "lsame" : "al";
      EColor lineCol = pol == AnitaPol::kHorizontal ? kBlue : kBlack;
      gr->SetLineColor(lineCol);
      gr->SetMaximum(yMax);
      gr->SetMinimum(yMin);
      
      gr->Draw(opt);

      for(unsigned i=0; i < gr->fDerivatives.size(); i++){
	TGraphFB* gr2 = gr->fDerivatives[i];
	gr2->Draw("lsame");
	gr2->SetLineColor(gr->GetLineColor());
	gr2->SetLineStyle(i+1);
      }      
    }
  }
}













void Acclaim::TGraphFB::drawCopy() const{


  if(fDerivedFrom){ // this setup should recursively call draw copy with the parent
    fDerivedFrom->drawCopy();
  }  
  else{
    int logx = gPad->GetLogx();
    int logy = gPad->GetLogy();
    int logz = gPad->GetLogz();  

    TCanvas* c1 = new TCanvas();
    c1->cd();
    c1->SetLogx(logx);
    c1->SetLogy(logy);
    c1->SetLogz(logz);

    TGraphFB* gr2 = new TGraphFB(this);    
    gr2->SetBit(kCanDelete);
    gr2->fDoubleClickOption = kDrawRayleigh;
    gr2->SetTitle(TString::Format(""));
  
    gr2->Draw();
    for(unsigned i=0; i < fDerivatives.size(); i++){
      TGraphFB* grD = new TGraphFB(fDerivatives.at(i));
      grD->SetBit(kCanDelete);
      grD->fDoubleClickOption = kDrawRayleigh;
      grD->SetTitle(TString::Format(""));
      grD->Draw("lsame");
    }
   
    gPad->Modified();
    gPad->Update();
  }
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

  int logx = gPad->GetLogx();
  int logy = gPad->GetLogy();
  int selectedPoint = -1;
  double minD2 = 1e99;
  for(int i=0; i < GetN(); i++){
    double xVal = logx ? TMath::Log10(GetX()[i]) : GetX()[i];
    double yVal = logy ? TMath::Log10(GetY()[i]) : GetY()[i];
    Int_t pixelX = gPad->XtoAbsPixel(xVal);
    Int_t pixelY = gPad->YtoAbsPixel(yVal);

    // std::cout << selectedPoint << "\t" << i << "\t" << x << "\t" << y << "\t" << pixelX << "\t" << pixelY << std::endl;
    

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


// Utility function to set fDerives from and fDerivatives (i.e. the drawing ownership)
void Acclaim::TGraphFB::setDrawingDependencies(const std::vector<TGraphFB*> grs){
  if(grs.size() > 0){
    for(unsigned i=1; i < grs.size(); i++){
      grs[0]->fDerivatives.push_back(grs[i]);
      grs[i]->fDerivedFrom = grs[0];
    }    
  }
}
