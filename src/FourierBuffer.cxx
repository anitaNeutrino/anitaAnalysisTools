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

Acclaim::FourierBuffer::FourierBuffer(Int_t theBufferSize) :
  doneVectorInit(false), fCurrentlyLoadingHistory(false), eventsInBuffer(0), fMinFitFreq(0.15), fMaxFitFreq(1.3), fMinSpecFreq(0.2), fMaxSpecFreq(1.3)
{
  
  // timeScale = timeScaleSeconds;
  bufferSize = theBufferSize <= 0 ? 1000 : theBufferSize;
  df = -1;
  fDrawFreqBin = 100; //26;

  // will initialize this dynamically to get around this no-copy-constructor bullshit
  fSpectrum = NULL;

  const char* rayleighFuncText = "([0]*x/([1]*[1]))*exp(-x*x/(2*[1]*[1]))";
  // const char* riceFuncText = "([0]*x/([1]*[1]))*exp(-(x*x+[2]*[2])/(2*[1]*[1]))*TMath::BesselI0([2]*x/([1]*[1]))";

  for(int ant=0; ant < NUM_SEAVEYS; ant++){
    summaryPads[ant] = NULL;
  }
  
  fRay = new TF1("fRay", rayleighFuncText, 0, 100, TF1::EAddToList::kNo);
}



void Acclaim::FourierBuffer::initVectors(int n, double df){

  
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      sumPowers[pol][ant].resize(n, 0);

      chiSquares[pol][ant].resize(n, 0);
      ndfs[pol][ant].resize(n, 0);

      probs[pol][ant].resize(n, 0);      

      fitAmplitudes[pol][ant].resize(n, 0);
      spectrumAmplitudes[pol][ant].resize(n, 0);


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

      // summary graphs for drawSummary
      grChiSquares[pol].push_back(TGraphFB(this, ant, pol, n));
      grNDFs[pol].push_back(TGraphFB(this, ant, pol, n));
      grReducedChiSquares[pol].push_back(TGraphFB(this, ant, pol, n));
      
      grProbs[pol].push_back(TGraphFB(this, ant, pol, n));

      grAmplitudes[pol].push_back(TGraphFB(this, ant, pol, n));
      grSpectrumAmplitudes[pol].push_back(TGraphFB(this, ant, pol, n));
      


      for(int freqBin=0; freqBin < n; freqBin++){
	double f = df*freqBin;
	grChiSquares[pol][ant].GetX()[freqBin] = f;
	grNDFs[pol][ant].GetX()[freqBin] = f;
	grReducedChiSquares[pol][ant].GetX()[freqBin] = f;

	grAmplitudes[pol][ant].GetX()[freqBin] = f;	
	grSpectrumAmplitudes[pol][ant].GetX()[freqBin] = f;

	grProbs[pol][ant].GetX()[freqBin] = f;

	grChiSquares[pol][ant].GetY()[freqBin] = 0;
	grNDFs[pol][ant].GetY()[freqBin] = 0;
	grReducedChiSquares[pol][ant].GetY()[freqBin] = 0;
	
	grAmplitudes[pol][ant].GetY()[freqBin] = 0;	
	grSpectrumAmplitudes[pol][ant].GetY()[freqBin] = 0;
	grProbs[pol][ant].GetY()[freqBin] = -0.1;
	
      }	
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



size_t Acclaim::FourierBuffer::add(const FilteredAnitaEvent* fEv){  


  if(!fCurrentlyLoadingHistory && eventNumbers.size() == 0){
    automagicallyLoadHistory(fEv);
  }
  
  const RawAnitaHeader* header = fEv->getHeader();
  eventNumbers.push_back(header->eventNumber);
  runs.push_back(header->run);


  bool anyUpdated = false; // should all get updated at the same time, but whatever...
      
  
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
	    anyUpdated = true;
	    hRays[pol][ant].at(freqInd)->getRayleighFitParams(fitAmplitudes[pol][ant][freqInd],
							      chiSquares[pol][ant][freqInd],
							      ndfs[pol][ant][freqInd]);

	    // if(ant == 47 && pol == AnitaPol::kHorizontal && freqInd == 44){
	    //   std::cout << pol << "\t" << ant << freqInd << "\t"
	    // 		<< fitAmplitudes[pol][ant][freqInd] << "\t"
	    // 		<< chiSquares[pol][ant][freqInd] << "\t"
	    // 		<< ndfs[pol][ant][freqInd] << std::endl;
	    // }
	    
	    grChiSquares[pol][ant].GetY()[freqInd] = chiSquares[pol][ant][freqInd];
	    if(ndfs[pol][ant][freqInd] > 0){
	      
	      grReducedChiSquares[pol][ant].GetY()[freqInd] = chiSquares[pol][ant][freqInd]/ndfs[pol][ant][freqInd];

	    // if(ant == 47 && pol == AnitaPol::kHorizontal && freqInd == 44){
	    //   std::cout << "inside ... " << std::endl;
	    //   std::cout << pol << "\t" << ant << "\t" << freqInd << "\t"
	    // 		<< fitAmplitudes[pol][ant][freqInd] << "\t"
	    // 		<< chiSquares[pol][ant][freqInd] << "\t"
	    // 		<< ndfs[pol][ant][freqInd] << std::endl;
	    // }
	      
	    }
	    grNDFs[pol][ant].GetY()[freqInd] = ndfs[pol][ant][freqInd];
	    grAmplitudes[pol][ant].GetY()[freqInd] = fitAmplitudes[pol][ant][freqInd];
	    
	  }

	  // this is N, if the probability of getting this amplitude or higher, given the rayeligh distribution is 1/N.
	  double probDenominator = 1./hRays[pol][ant].at(freqInd)->getOneMinusCDF(amp);
	  probs[pol][ant].at(freqInd) = probDenominator;
	  grProbs[pol][ant].GetY()[freqInd] = probDenominator;
	}

      }
    }
  }


  if(anyUpdated){
    int firstSpecInd = int(fMinSpecFreq/df);
    int lastSpecInd = int(fMaxSpecFreq/df);

    // double sumLowEdgeBins = 0;
    // double sumHighEdgeBins = 0;    
    // int numSum = 0;
    // for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    //   AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    //   for(int ant=0; ant < NUM_SEAVEYS; ant++){
    // 	sumLowEdgeBins += fitAmplitudes[pol][ant][firstSpecInd];
    // 	sumHighEdgeBins += fitAmplitudes[pol][ant][lastSpecInd];
    // 	numSum++;
    //   }
    // }

    // double meanSpecLowEdge = sumLowEdgeBins/numSum;
    // double meanSpecHighEdge = sumHighEdgeBins/numSum;    
    
    for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
      for(int ant=0; ant < NUM_SEAVEYS; ant++){

	int nf = fitAmplitudes[pol][ant].size();
	// std::cout << fitAmplitudes[pol][ant].size() << "\t" << sumPowers[ant][pol].size() << "\t" << nf << std::endl;
	for(int freqInd = 0; freqInd < nf; freqInd++){
	  // if(freqInd < firstSpecInd){
	  //   spectrumAmplitudes[pol][ant][freqInd] = meanSpecLowEdge;
	  // }
	  // else if(freqInd >= lastSpecInd){
	  //   spectrumAmplitudes[pol][ant][freqInd] = meanSpecHighEdge;
	  // }
	  // else{
	    spectrumAmplitudes[pol][ant][freqInd] = fitAmplitudes[pol][ant][freqInd];
	  // }
	}
	
	// getSpectrum(&spectrumAmplitudes[pol][ant][0], nf);
	getSpectrum(&spectrumAmplitudes[pol][ant][firstSpecInd], lastSpecInd - firstSpecInd);	
	for(int freqInd = 0; freqInd < nf; freqInd++){
	  grSpectrumAmplitudes[pol][ant].GetY()[freqInd] = spectrumAmplitudes[pol][ant][freqInd];
	}
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



void Acclaim::FourierBuffer::automagicallyLoadHistory(const FilteredAnitaEvent* fEv){

  fCurrentlyLoadingHistory = true;
  
  int fuck = 0;
  bool loadedHistory = false;

  UInt_t desiredFirstEvent = fEv->getHeader()->eventNumber;
  Int_t run = fEv->getHeader()->run;


  while(!loadedHistory){

    AnitaDataset d(run);

    d.first();
    UInt_t firstEventThisRun = d.header()->eventNumber;

    if(firstEventThisRun > desiredFirstEvent){
      run--;
    }
    else{

      // assume monotonically increasing eventNumber at one per event...
      int desiredFirstEntry = desiredFirstEvent - firstEventThisRun;

      for(int bi = 0; bi < bufferSize; bi++){
	int thisEntry = d.getEntry(bi + desiredFirstEntry);

	if(thisEntry > 0){
	
	  if(bi==0 && d.header()->eventNumber != desiredFirstEvent){
	    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", couldn't find the right first event, wanted "
		      << desiredFirstEvent << ", but got " << d.header()->eventNumber << ", while looking in run "
		      << run << std::endl;
	  }
	
	  FilteredAnitaEvent fEv2(d.useful(), (FilterStrategy*) fEv->getStrategy(), d.gps(), d.header(), false);
	  add(&fEv2);
	}	
      }

      if(eventNumbers.size() != (unsigned) bufferSize){
	std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", tried to load " << bufferSize
		  << " events of history but instead I got " <<  eventNumbers.size() << " events..." << std::endl;
      }
      
      loadedHistory = true;
    }
    fuck++;
    if(fuck==4){
      std::cerr << "quadruple fuck" << std::endl;
      exit(4);
    }
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


// wrapper for spectrum since it's interface changes in a fucking stupid way between ROOT versions
void Acclaim::FourierBuffer::getSpectrum(double* y, int n) const{

  if(!fSpectrum){
    fSpectrum = new TSpectrum();
  }

  const int numIter = 6;

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)  
  double* tempPtr = y;
#else
  std::vector<float> tempFloats(n, 0);
  for(int i=0; i < n; i++){
    tempFloats[i] = y[i];
  }
  float* tempPtr = &tempFloats[0];
#endif

  fSpectrum->Background(tempPtr,n,
			numIter,TSpectrum::kBackDecreasingWindow,
			TSpectrum::kBackOrder2,kFALSE,
			TSpectrum::kBackSmoothing3,kFALSE);
  
#if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
  for(int i=0; i < n; i++){
    y[i] = tempFloats[i];
  }
#endif
  
  
}

Acclaim::TGraphFB* Acclaim::FourierBuffer::getBackground(Int_t ant, AnitaPol::AnitaPol_t pol, int lastNEvents) const{

  TGraphFB* gr = getAvePowSpec(ant, pol, lastNEvents);

  getSpectrum(gr->GetY(), gr->GetN());

  return gr;
}





void Acclaim::FourierBuffer::drawSummary(TPad* pad, SummaryOption_t summaryOpt) const{

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
      
      // TGraphFB& gr = (TGraphFB&) grReducedChiSquares[pol][ant]; 
      const TGraphFB* gr = getSelectedGraphForSummary(summaryOpt, ant, pol);
      // TGraphFB& gr = (TGraphFB&) grProbs[pol][ant];

      for(int i=0; i < gr->GetN(); i++){
	if(gr->GetY()[i] > yMax){
	  yMax = gr->GetY()[i];
	}
	if(gr->GetY()[i] < yMin){
	  yMin = gr->GetY()[i];
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
    if(summaryOpt == FourierBuffer::Prob){
      gPad->SetLogy(1);
      yMin = 1;
    }
    else{
      gPad->SetLogy(0);
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

      if(summaryOpt == FourierBuffer::RayleighAmplitude){
	grSpectrumAmplitudes[pol][ant].SetLineColor(gr->GetLineColor());
	grSpectrumAmplitudes[pol][ant].SetLineStyle(2);
	grSpectrumAmplitudes[pol][ant].Draw("lsame");
      }           
    }
  }
  
}












void Acclaim::TGraphFB::drawCopy() const{

  int logx = gPad->GetLogx();
  int logy = gPad->GetLogy();
  int logz = gPad->GetLogz();  

  TCanvas* c1 = new TCanvas();
  c1->cd();
  c1->SetLogx(logx);
  c1->SetLogy(logy);
  c1->SetLogz(logz);    

  TGraphFB* gr2 = (TGraphFB*) Clone();
  gr2->pol = pol;
  gr2->ant = ant;
  gr2->fb = fb;
  gr2->SetBit(kCanDelete);
  gr2->fDoubleClickOption = kDrawRayleigh;
  gr2->SetTitle(TString::Format(""));
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
