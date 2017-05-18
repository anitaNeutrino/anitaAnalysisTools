#include "AcclaimFilters.h"
#include "AnalysisWaveform.h"
#include "FilteredAnitaEvent.h"
#include "BasicFilters.h" // Cosmin's example filters
#include "ResponseManager.h"
#include "FourierBuffer.h"
#include "FFTWComplex.h"
#include <iostream>
#include "TPad.h"
#include "RayleighHist.h"
#include "RootTools.h"

static std::map<TString, FilterStrategy*> acclaimDefaults;

/** 
 * Checks all operations in the Filter Strategy and forces the strategy to load history the before it processes the next event
 * This is useful for breaking up analysis runs into multiple jobs wthout disturbing the rolling averages.
 *
 * @param fs is the FilterStrategy
 */
void Acclaim::Filters::makeFourierBuffersLoadHistoryOnNextEvent(FilterStrategy* fs){
  // load history if needed
  for(unsigned i=0; i < fs->nOperations(); i++){
    const Acclaim::Filters::RayleighMonitor* rmConst = dynamic_cast<const Acclaim::Filters::RayleighMonitor*>(fs->getOperation(i));
    if(rmConst){
      const Acclaim::FourierBuffer* fb = rmConst->getFourierBuffer();
      fb->setForceLoadHistory(true);
    }
  }
}



/** 
 * This function adds my custom strategies to a map of TString to strategies...
 */
void Acclaim::Filters::appendFilterStrategies(std::map<TString, FilterStrategy*>& filterStrats, bool saveOutput){


  if(acclaimDefaults.size()==0){

    // first make the operations...

    // should always use some kind of alfa filter unless you have a good reason
    const double fpeSafety = 1e-10;    
    ALFAFilter* alfaFilter = new ALFAFilter(Bands::alfaLowPassGHz - fpeSafety);
    
    Notch* bandHighPass = new Notch(0, Bands::anitaHighPassGHz);
    Notch* bandLowPass = new Notch(Bands::anitaLowPassGHz, 2);

    // double log10ProbThresh = -100; //2.5;
    double reducedChiSquareThresh = 5;
    double channelChiSquareCdfThresh = 0.995;
    const int numEventsInRayleighDistributions = 1500;
    RayleighFilter* rf = new RayleighFilter(channelChiSquareCdfThresh, reducedChiSquareThresh, numEventsInRayleighDistributions);



    AnitaResponse::AllPassDeconvolution* allPass = new AnitaResponse::AllPassDeconvolution();    
    TString responseDir = TString::Format("%s/share/AnitaAnalysisFramework/responses/IndividualBRotter", getenv("ANITA_UTIL_INSTALL_DIR"));
    AnitaResponse::ResponseManager* responseManager = new AnitaResponse::ResponseManager(responseDir.Data(), 0, allPass);
    AnitaResponse::DeconvolveFilter* df = new AnitaResponse::DeconvolveFilter(responseManager, allPass);
    // then make the strategies

    // every operation is going to use these default strategies
    FilterStrategy* defaultOps = new FilterStrategy();
    defaultOps->addOperation(alfaFilter, saveOutput); // has internal check for ANITA version
    acclaimDefaults["Minimum"] = defaultOps;

    // every operation is going to use these default strategies
    FilterStrategy* defaultDeco = new FilterStrategy();
    (*defaultDeco) = (*defaultOps);
    defaultDeco->addOperation(df, saveOutput); // has internal check for ANITA version
    acclaimDefaults["Deconvolve"] = defaultDeco;
  
    FilterStrategy* fs = new FilterStrategy();
    (*fs) = (*defaultOps);
    fs->addOperation(bandHighPass, saveOutput);
    fs->addOperation(bandLowPass, saveOutput);
    fs->addOperation(rf, saveOutput);
    acclaimDefaults["RayleighFilter"] = fs;
  }

  std::map<TString, FilterStrategy*>::iterator it;
  for(it = acclaimDefaults.begin(); it != acclaimDefaults.end(); ++it){
    filterStrats[it->first] = it->second;
  }  
}

FilterStrategy* Acclaim::Filters::findDefaultStrategy(const TString& stratName){

  appendFilterStrategies(acclaimDefaults, false);
  return findStrategy(acclaimDefaults, stratName);
}

FilterStrategy* Acclaim::Filters::findStrategy(const std::map<TString, FilterStrategy*>& filterStrats, const TString& stratName){

  FilterStrategy* fs = NULL;
  std::map<TString, FilterStrategy*>::const_iterator it = filterStrats.find(stratName);

  if(it!=filterStrats.end()){
    fs = it->second;
  }
  else{
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unable to find filter strategy " << stratName << std::endl;
  }
  return fs;
}






Acclaim::Filters::Notch::Notch(Double_t lowEdgeGHz, Double_t highEdgeGHz){

  // Freq bins are currently in 0.1 GHz steps
  fTag = TString::Format("notch%.lfGHzto%.lfGHz", lowEdgeGHz, highEdgeGHz);
  fDescription = TString::Format("Notch filter from %.lf GHz to %.lf GHz", lowEdgeGHz, highEdgeGHz);
  fLowEdgeGHz = lowEdgeGHz;
  fHighEdgeGHz = highEdgeGHz;
  // fLowNotchIndex = -1;
  // fHighNotchIndex = -1;  

  for(int pol=0; pol < AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      fPowerRemovedByNotch[pol][ant] = 0;
    }
  }
}


void Acclaim::Filters::Notch::processOne(AnalysisWaveform * g){
  // add small offset to frequencies to avoid inconsistent notch edges due to
  // floating point error...
  const double floatPointError = 1e-10; // plently 
      
  const double deltaF_GHz = g->deltaF();
  const int nf = g->Nfreq();
  FFTWComplex* theFreqs = g->updateFreq();
  for(int freqInd=0; freqInd < nf; freqInd++){
    const double freqGHz = deltaF_GHz* freqInd + floatPointError;

    if(freqGHz >= fLowEdgeGHz && freqGHz < fHighEdgeGHz){
      theFreqs[freqInd] = 0;
    }
  }
}


void Acclaim::Filters::Notch::process(FilteredAnitaEvent * ev) 
{
  for(int pol=0; pol < AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      fPowerRemovedByNotch[pol][ant] = 0;
      AnalysisWaveform* wf = getWf(ev, ant, (AnitaPol::AnitaPol_t) pol);
      const TGraphAligned* grPower = wf->power();
      for(int i=0; i < grPower->GetN(); i++){
	fPowerRemovedByNotch[pol][ant] += grPower->GetY()[i];
      }
    }
  }
  
  UniformFilterOperation::process(ev);

  for(int pol=0; pol < AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      AnalysisWaveform* wf = getWf(ev, ant, (AnitaPol::AnitaPol_t) pol);
      const TGraphAligned* grPower = wf->power();
      for(int i=0; i < grPower->GetN(); i++){
	fPowerRemovedByNotch[pol][ant] -= grPower->GetY()[i];
      }
    }
  }
}


















Acclaim::Filters::SpikeSuppressor::SpikeSuppressor(double spikeThresh_dB, int numEvents) : fRandy(1234), fourierBuffer(numEvents) {
  fSpikeThresh_dB = spikeThresh_dB;
  fNumEvents = numEvents;
  fDescription = Form("Finds spikes (as I define them) greater than %4.2lf dB and removes them, replacing them with noise", fSpikeThresh_dB);
}

void Acclaim::Filters::SpikeSuppressor::processOne(AnalysisWaveform* wave){
  (void) wave;
  std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", I need more information, this filter is designed to be called with process(FilteredAnitaEvent*) instead!" << std::endl;
}



void Acclaim::Filters::SpikeSuppressor::process(FilteredAnitaEvent* fEv){

  setSeed(fEv->getHeader()->eventNumber); // for deterministic randomness

  fourierBuffer.add(fEv);
  
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(int ant=0; ant < NUM_SEAVEYS; ant++){      
      AnalysisWaveform* wave = getWf(fEv, ant, pol);

      TGraphAligned* grAvePowSpec_dB = fourierBuffer.getAvePowSpec_dB(ant, pol);

      const double cosminFactor = 1000*50*grAvePowSpec_dB->GetN(); // from the fft normlization option in AnalysisWaveform
      

      // TGraphAligned grFiltered_dB = suppressSpikes(grAvePowSpec_dB);


      TGraphAligned* grBackground_dB = fourierBuffer.getBackground_dB(ant, pol);
      TGraphAligned grFiltered_dB = suppressSpikes(grAvePowSpec_dB, grBackground_dB);      
      TGraphAligned grFiltered_not_dB = grFiltered_dB;
      grFiltered_not_dB.undBize();


      FFTWComplex* fft = wave->updateFreq();
      for(int i=0; i < grFiltered_not_dB.GetN(); i++){
	if(grFiltered_dB.GetY()[i] != grAvePowSpec_dB->GetY()[i]){

	  double mag = TMath::Sqrt(grFiltered_not_dB.GetY()[i]*cosminFactor);
	  double phase = TMath::TwoPi()*fRandy.Uniform();
	  fft[i].setMagPhase(mag, phase);
	}
      }

      delete grAvePowSpec_dB;
      grAvePowSpec_dB = NULL;
      
      delete grBackground_dB;
      grBackground_dB = NULL;
      
    }
  }
}




TGraphAligned Acclaim::Filters::SpikeSuppressor::suppressSpikes(const TGraphAligned* grPower, const TGraphAligned* grBackground){
  // expecting both in dB...  

  
  TGraphAligned grFilteredSpec = (*grPower);
  
  for(int i=0; i < grPower->GetN(); i++){
    double x = grPower->GetX()[i];

    double y = grPower->GetY()[i];    

    double yBack = grBackground->Eval(x);

    if(y - yBack > fSpikeThresh_dB){
      grFilteredSpec.GetY()[i] = yBack;
    }
  }

  return grFilteredSpec;
}







TGraphAligned Acclaim::Filters::SpikeSuppressor::suppressSpikes(const TGraphAligned* grPower){
  
  std::vector<int> extremaSamps;

  // this will leave out the edges, which is fine.
  for(Int_t samp=1; samp < grPower->GetN()-1; samp++){
    Double_t y0 = grPower->GetY()[samp-1];
    Double_t y1 = grPower->GetY()[samp];
    Double_t y2 = grPower->GetY()[samp+1];


    if(y0 > y1 && y1 < y2){
      // Is a local minimum
      extremaSamps.push_back(samp);
    }
    else if(y0 < y1 && y1 > y2){
      // Is a local maximum
      extremaSamps.push_back(samp);
    }
  }

  std::vector<double> x0s;
  std::vector<double> x2s;

  std::vector<double> y0s;
  std::vector<double> y2s;
    
  if(extremaSamps.size() <=1) {
    // this shouldn't be true but there might be a weird waveform somewhere...	
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << "You have a very strange waveform, this was unexpected..." << std::endl;
  }
  else{
	
    for(unsigned i=1; i < extremaSamps.size() - 1; i++){	
      int samp0 = extremaSamps.at(i-1);
      int samp1 = extremaSamps.at(i);
      int samp2 = extremaSamps.at(i+1);

      double y0 = grPower->GetY()[samp0];
      double y1 = grPower->GetY()[samp1];
      double y2 = grPower->GetY()[samp2];

      // it's a maxima
      if(y0 < y1 && y1 > y2){

	// now we figure out the spike width...
	// using the next local minima as the spike edges isn't good enough...
	// try defining the width as the distance from the from the closest minima
	// to the peak (extended to the other side)	    
	    
	double x0 = grPower->GetX()[samp0];
	double x1 = grPower->GetX()[samp1];
	double x2 = grPower->GetX()[samp2];

	double halfDx = TMath::Min(x1 - x0, x2 - x1);

	// redefine spike edges using halfDx
	x0 = x1 - halfDx;
	x2 = x1 + halfDx;

	y0 = grPower->Eval(x0);
	y2 = grPower->Eval(x2);

	double yInterp = interpolate_dB(x1, x0, x2, y0, y2);



	double spikeSize_dB = y1 - yInterp;


	// std::cout << fSpikeThresh_dB << "\t" << spikeSize_dB << "\t" << yInterp << "\t" << y1 << std::endl;

	if(spikeSize_dB > fSpikeThresh_dB){
	  x0s.push_back(x0);
	  x2s.push_back(x2);
	  y0s.push_back(y0);
	  y2s.push_back(y2);
	}
      }	  
    }
  }
   

  // std::cout << " I found " << y2s.size() << " spikes!" << std::endl;
  
  // std::cout << "I found spikes at: " << std::endl;
  // for(unsigned spike = 0; spike < x0s.size(); spike++){
  //   std::cout << 1e3*x0s.at(spike) << ", " << 1e3*x2s.at(spike) << std::endl;
  // }

  TGraphAligned grFilteredSpec = (*grPower);
  // std::vector<int> modified(grFilteredSpec.GetN(), 0);
  for(unsigned spike = 0; spike < x0s.size(); spike++){

    double x0 = x0s.at(spike);
    double x2 = x2s.at(spike);
    
    double y0 = y0s.at(spike);

    double y2 = y2s.at(spike);

    for(int i=0; i < grFilteredSpec.GetN(); i++){

      double x1 = grFilteredSpec.GetX()[i];

      if(x1 >= x0 && x1 < x2){
	double yInterp = interpolate_dB(x1, x0, x2, y0, y2);
	grFilteredSpec.GetY()[i] = yInterp;
	// modified.at(i) = 1;	
      }
      // std::cout << modified.at(i) << ", ";
    }
  }

  return grFilteredSpec;
}



double Acclaim::Filters::SpikeSuppressor::interpolate_dB(double x, double xLow, double xHigh, double yLow, double yHigh){

  double deltaX = xHigh - xLow;
  double deltaY = pow(10, yHigh) - pow(10, yLow);
  double gradient = deltaY/deltaX;
  // std::cout << deltaY << "\t" << deltaX << std::endl;
  // y = mx + c
  double yInterp = (x - xLow)*gradient + pow(10, yLow);
  yInterp = TMath::Log10(yInterp);

  return yInterp;
}





Acclaim::Filters::UniformMagnitude::UniformMagnitude(){
}

void Acclaim::Filters::UniformMagnitude::processOne(AnalysisWaveform* wf){

  FFTWComplex* fft = wf->updateFreq();  
  for(int i=0; i< wf->Nfreq(); i++){
    fft[i].setMagPhase(i > 0 ? 1 : 0, fft[i].getPhase());
  }
}


Acclaim::Filters::SpectrumMagnitude::SpectrumMagnitude(Int_t numEvents) : RayleighMonitor(numEvents){
  fDescription = "Sets the magnitude of each frequency bin equal to the fourier buffer spectrum";
}

void Acclaim::Filters::SpectrumMagnitude::process(FilteredAnitaEvent* fEv){
  RayleighMonitor::process(fEv);
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      AnalysisWaveform* wf = getWf(fEv, ant, pol);

      const int nf = wf->Nfreq();

      FFTWComplex* theFreqs = wf->updateFreq();
      for(int freqInd=0; freqInd < nf; freqInd++){
	double mag = fourierBuffer.getBackgroundSpectrumAmp(pol, ant, freqInd);
	double phase = theFreqs[freqInd].getPhase();
	theFreqs[freqInd].setMagPhase(mag, phase);
      }
    }
  }
}





























Acclaim::Filters::RayleighMonitor::RayleighMonitor(int numEvents) : fourierBuffer(numEvents) {
  fNumEvents = numEvents;
  fDescription = TString::Format("Decides whether or not to filter events based on characteristics of the event amplitude and Rayleigh distribution over %d events", fNumEvents);
  fNumOutputs = 6;
  fOutputAnt = 4;
  fOutputPol = AnitaPol::kHorizontal;
}


void Acclaim::Filters::RayleighMonitor::process(FilteredAnitaEvent* fEv){  
  fourierBuffer.add(fEv);
}




unsigned Acclaim::Filters::RayleighMonitor::outputLength(unsigned i) const{
  (void) i;
  if(i==0){
    return 1;
  }
  else{
    // return AnitaPol::kNotAPol*NUM_SEAVEYS*131; // TODO, check this at runtime maybe?
    return 131; // TODO, check this at runtime maybe?    
  }
}



const char* Acclaim::Filters::RayleighMonitor::outputName(unsigned i) const{
  // so what do I want to output?
  switch(i){
  case 0:
    return "numEvents";
  case 1:
    return "chiSquares";    
  case 2:
    return "ndfs";
  case 3:
    return "rayleighAmps";
  case 4:
    return "spectrumAmps";    
  case 5:
    return "prob";
  default:
    return NULL;
  };
}


void Acclaim::Filters::RayleighMonitor::fillOutput(unsigned i, double* v) const{

  if(i==0){
    v[0] = fourierBuffer.getNumEventsInBuffer();
    return;
  }
  
  int outInd = 0;
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      if(ant==fOutputAnt && pol == fOutputPol){
      
	const double* vals = 0;
	const int* valsInt = 0;
	unsigned n = 0;
	switch(i){
	case 0:
	  std::cerr << "uh oh!" << std::endl;
	  break;
	case 1:
	  vals = &fourierBuffer.getChiSquares(ant, pol)[0];
	  n = fourierBuffer.getChiSquares(ant, pol).size();
	  break;
	case 2:
	  valsInt = &fourierBuffer.getNDFs(ant, pol)[0];
	  n = fourierBuffer.getNDFs(ant, pol).size();
	  break;
	case 3:
	  vals = &fourierBuffer.getRayleighAmplitudes(ant, pol)[0];
	  n = fourierBuffer.getRayleighAmplitudes(ant, pol).size();
	  break;
	case 4:
	  vals = &fourierBuffer.getBackgroundSpectrumAmplitudes(ant, pol)[0];
	  n = fourierBuffer.getBackgroundSpectrumAmplitudes(ant, pol).size();
	  break;
	case 5:
	  vals = &fourierBuffer.getProbabilities(ant, pol)[0];
	  n = fourierBuffer.getProbabilities(ant, pol).size();
	  break;
	}

	// std::cout << vals << "\t" << valsInt << "\t" << n << std::endl;
	for(unsigned point=0; point < n; point++){
	  if(vals){
	    v[outInd] = vals[point];
	  }
	  else if(valsInt){
	    v[outInd] = valsInt[point];
	  }
	  outInd++;
	}
      }
    }
  }
}






















Acclaim::Filters::RayleighFilter::RayleighFilter(double channelChiSquareCdfThresh, double chiSquarePerDofThreshold, int numEvents) : RayleighMonitor(numEvents), fChiSquarePerDofThreshold(chiSquarePerDofThreshold), fChanChiSquareCdfThreshold(channelChiSquareCdfThresh)
{
  fRandy = new TRandom3(1234); // seed will be reset on a per event basis using the eventNumber
  fDescription = TString::Format("Tracks frequency bin amplitudes over %d events, the chi squared threshold is %lf", fNumEvents, fChiSquarePerDofThreshold);
  
  fChanChiSquareThreshold = -2*TMath::Log(1 - fChanChiSquareCdfThreshold);
}

Acclaim::Filters::RayleighFilter::~RayleighFilter()
{
  if(fRandy){
    delete fRandy;
    fRandy = NULL;
  }
}


void Acclaim::Filters::RayleighFilter::process(FilteredAnitaEvent* fEv){

  UInt_t eventNumber = fEv->getHeader()->eventNumber;  
  fRandy->SetSeed(eventNumber);
  RayleighMonitor::process(fEv);
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      AnalysisWaveform* wf = getWf(fEv, ant, pol);

      const int nf = wf->Nfreq();

      double meanChanChiSquare = 0;
      double varChanChiSquare = 0;
      fourierBuffer.getMeanVarChanChiSquares(ant, pol, meanChanChiSquare, varChanChiSquare);

      // thisEventChiSquareThreshold =
      double thisChannelExtraThreshold = meanChanChiSquare > 2 ? meanChanChiSquare - 2 : 0;
      double thisChannelChiChiSquareThreshold = fChanChiSquareThreshold + thisChannelExtraThreshold;
      // fChanChiSquareCdfThreshold

      FFTWComplex* theFreqs = wf->updateFreq();
      for(int freqInd=0; freqInd < nf; freqInd++){

	// double probVal = fourierBuffer.getProb(pol, ant, freqInd);
	int ndf = fourierBuffer.getNDFs(ant, pol)[freqInd];
	// double reducedChiSquare = fourierBuffer.getChiSquares(ant, pol)[freqInd]/ndf;
	double reducedChiSquareRelativeToSpectrum = ndf > 0 ? fourierBuffer.getChiSquaresRelativeToSpectrum(ant, pol)[freqInd]/ndf : -1;
        
        double freqBinChiSquare = fourierBuffer.getChanChiSquares(ant, pol)[freqInd];
	
	if(reducedChiSquareRelativeToSpectrum > fChiSquarePerDofThreshold || freqBinChiSquare > thisChannelChiChiSquareThreshold){
	  
	  // double specAmp = fourierBuffer.getBackgroundSpectrumAmp(pol, ant, freqInd);
	  // // std::cout << pol << "\t" << ant << "\t" << freqInd << "\t" << probVal << "\t" << specAmp << std::endl;
	  
	  // double x1 = fRandy->Gaus(0, specAmp);
	  // double x2 = fRandy->Gaus(0, specAmp);
	  // double newAmp = TMath::Sqrt(x1*x1 + x2*x2);
	  
	  // double phase = fRandy->Uniform(0, TMath::TwoPi());

	  // theFreqs[freqInd].setMagPhase(newAmp, phase);
	  theFreqs[freqInd].re = 0;
	  theFreqs[freqInd].im = 0;
	}        
      }
      
      int nZero=0;
      for(int freqInd=0; freqInd < nf; freqInd++){
        if(theFreqs[freqInd].re == 0 && theFreqs[freqInd].im == 0){
          nZero++;
        }
      }
      if(nZero==nf){
        std::cerr << "Event " << fEv->getHeader()->eventNumber << ", channel " << pol << ", " << ant << " has all zero frequency bins...";
        int nFailReducedChiSquare = 0;
        int nFailChannelChiSquare = 0;
        for(int freqInd=0; freqInd < nf; freqInd++){
          if(freqInd >= 26 && freqInd < 120){
            if(theFreqs[freqInd].re == 0 && theFreqs[freqInd].im == 0){
              int ndf = fourierBuffer.getNDFs(ant, pol)[freqInd];
              double reducedChiSquareRelativeToSpectrum = ndf > 0 ? fourierBuffer.getChiSquaresRelativeToSpectrum(ant, pol)[freqInd]/ndf : -1;
        
              double freqBinChiSquare = fourierBuffer.getChanChiSquares(ant, pol)[freqInd];
	
              if(reducedChiSquareRelativeToSpectrum > fChiSquarePerDofThreshold){
                nFailReducedChiSquare++;
              }
              if(freqBinChiSquare > thisChannelChiChiSquareThreshold){
                nFailChannelChiSquare++;                
              }
            }
          }
        }
        std::cerr << nFailReducedChiSquare << " failed the reducedChiSquare, ";
        std::cerr << nFailChannelChiSquare << " failed the channelChiSquare " << std::endl;
      }
      
      // std::cout << std::endl;
    }
  }  
}
