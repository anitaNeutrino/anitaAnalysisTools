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
#include "RawAnitaHeader.h"

static std::map<TString, FilterStrategy*> acclaimDefaults;


/**
 * Returns a c-style string containing the name of Cosmin's current favourite sine sub filter.
 * What this actually is, may change with time.
 * @return Cosmin's fave.
 */
const char* Acclaim::Filters::getCosminsFavouriteSineSubName()
{
  return "sinsub_10_3_ad_2";
}


//---------------------------------------------------------------------------------------------------------
/**
 * @brief Checks all operations and forces the first FourierBuffer to load history the before it processes the next event
 *
 * This is useful for breaking up analysis runs into multiple jobs wthout disturbing the rolling averages.
 * The FourierBuffer is a memeber of the RayleighMonitor class and its descendents.
 *
 * @param fs is the FilterStrategy in which to look for FourierBuffers
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


void Acclaim::Filters::generateFilterStrategies(bool saveOutput){

  if(acclaimDefaults.size()==0){
    // first make the operations...

    // should always use some kind of alfa filter unless you have a good reason
    ALFAFilter* alfaFilter = new ALFAFilter(Bands::alfaLowPassGHz);
    
    Notch* bandHighPass = new Notch(0, Bands::anitaHighPassGHz);
    Notch* bandLowPass = new Notch(Bands::anitaLowPassGHz, 2);

    Notch* brickWall260 = new Notch(0.23, 0.29);
    Notch* brickWall370 = new Notch(0.34, 0.40);    

    // double log10ProbThresh = -100; //2.5;
    // double reducedChiSquareThresh = 5;
    // double channelChiSquareCdfThresh = 0.995;
    // const int numEventsInRayleighDistributions = 1500;
    // RayleighFilter* rf = new RayleighFilter(channelChiSquareCdfThresh, reducedChiSquareThresh, numEventsInRayleighDistributions);



    AnitaResponse::AllPassDeconvolution* allPass = new AnitaResponse::AllPassDeconvolution();    
    TString responseDir = TString::Format("%s/share/AnitaAnalysisFramework/responses/IndividualBRotter", getenv("ANITA_UTIL_INSTALL_DIR"));
    AnitaResponse::ResponseManager* responseManager = new AnitaResponse::ResponseManager(responseDir.Data(), 0, allPass);
    AnitaResponse::DeconvolveFilter* df = new AnitaResponse::DeconvolveFilter(responseManager, allPass);
    // then make the strategies

    // every operation is going to use these default strategies
    FilterStrategy* defaultOps = new FilterStrategy();
    defaultOps->addOperation(alfaFilter, saveOutput); // has internal check for ANITA version
    acclaimDefaults["Minimum"] = defaultOps;

    FilterStrategy* brick = new FilterStrategy();
    (*brick) = (*defaultOps);
    brick->addOperation(brickWall260, saveOutput);
    brick->addOperation(brickWall370, saveOutput);
    acclaimDefaults["BrickWallSatellites"] = brick;

    FilterStrategy* defaultDeco = new FilterStrategy();
    (*defaultDeco) = (*defaultOps);
    defaultDeco->addOperation(df, saveOutput); // has internal check for ANITA version
    acclaimDefaults["Deconvolve"] = defaultDeco;
  
    // FilterStrategy* fs = new FilterStrategy();
    // (*fs) = (*defaultOps);
    // fs->addOperation(bandHighPass, saveOutput);
    // fs->addOperation(bandLowPass, saveOutput);
    // fs->addOperation(rf, saveOutput);
    // acclaimDefaults["RayleighFilter"] = fs;
  }

}


//---------------------------------------------------------------------------------------------------------
/**
 * @brief Adds my custom strategies to a map of TString to strategies
 *
 * My analysis will only use a handful of filters strategies.
 * This function appends the "defaultStrategies" to a TString/FilterStrategy* map like what is used in MagidDisplay.
 * This function allows the strategies to be loaded painlessly into MagicDisplay and analysis scripts
 * Note that the strategies are initialized once so adding the saveOutput bool will currently only work on the first call.
 *
 * @param filterStrats is a the TString/FilterStrategy* map.
 * @param saveOutput is passed to the filterStrategies when they are created.
 */
void Acclaim::Filters::appendFilterStrategies(std::map<TString, FilterStrategy*>& filterStrats, bool saveOutput){

  generateFilterStrategies(saveOutput);

  std::map<TString, FilterStrategy*>::iterator it;
  for(it = acclaimDefaults.begin(); it != acclaimDefaults.end(); ++it){
    filterStrats[it->first] = it->second;
  }  
}



//---------------------------------------------------------------------------------------------------------
/**
 * @brief Wrapper function to find a strategy from the defaults
 *
 * This function is pretty dumb and doesn't do any error checking and could fail with undefined behaviour.
 *
 * @param stratName is the name of the strategy (i.e. the key in the TString/FilterStrategy* map. 
 */

FilterStrategy* Acclaim::Filters::findDefaultStrategy(const TString& stratName){

  appendFilterStrategies(acclaimDefaults, false);
  return findStrategy(acclaimDefaults, stratName);
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Wrapper function to do a std::map::iterator lookup of a TString/filterStrategy* map
 *
 * @param filterStrats is the map of filterStrategies
 * @param stratName is the name of the strategy (i.e. the key in the TString/FilterStrategy* map. 
 */
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







/** 
 * @brief Constructor for notch
 *
 * @param lowEdgeGHz is the frequency in GHz of the low edge of the notch (i.e. the high pass frequency)
 * @param highEdgeGHz is the frequency in GHz of the high edge of the notch (i.e. the low pass frequency)
 */
Acclaim::Filters::Notch::Notch(Double_t lowEdgeGHz, Double_t highEdgeGHz)
    : fLowEdgeGHz(lowEdgeGHz), fHighEdgeGHz(highEdgeGHz) {

  // Freq bins are currently in 0.01 GHz steps
  fTag = TString::Format("notch%4.2lfGHzto%4.2lfGHz", lowEdgeGHz, highEdgeGHz);
  fDescription = TString::Format("Notch filter from %4.2lf GHz to %4.2lf GHz", lowEdgeGHz, highEdgeGHz);
  // fLowNotchIndex = -1;
  // fHighNotchIndex = -1;  

  for(int pol=0; pol < AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      fPowerRemovedByNotch[pol][ant] = 0;
    }
  }
}


/** 
 * @brief applies the notch to a single waveform
 * 
 * @param g is the waveform to which the notch filters should be applied
 */
void Acclaim::Filters::Notch::processOne(AnalysisWaveform * g){
  // add small offset to frequencies to avoid inconsistent notch edges due to
  // floating point error...
  const double floatPointError = 1e-10; // plenty 
      
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


/** 
 * @brief Applies the notch to all waveforms in an event
 *
 * Also stores the power removed from each channel
 * 
 * @param ev is the event to which the notch should be applied to all channels
 */
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























/** 
 * @brief Constructor for uniform magnitude
 */
Acclaim::Filters::UniformMagnitude::UniformMagnitude(){
}


/** 
 * @brief Sets each frequency bin to have identical magnitude, without disturbing the phase
 *
 * Pretty dumb really.
 *  
 * @param wf is the waveform to be filtered
 */
void Acclaim::Filters::UniformMagnitude::processOne(AnalysisWaveform* wf){

  FFTWComplex* fft = wf->updateFreq();  
  for(int i=0; i< wf->Nfreq(); i++){
    fft[i].setMagPhase(i > 0 ? 1 : 0, fft[i].getPhase());
  }
}





/** 
 * @brief Constructor for SoectrumMagnitude filter
 * 
 * @param numEvents is the number of events over which to derived the TSpectrum averaged amplitudes
 */
Acclaim::Filters::SpectrumMagnitude::SpectrumMagnitude(Int_t numEvents) : RayleighMonitor(numEvents){
  fDescription = TString::Format("Sets the magnitude of each frequency bin equal to the fourier buffer TSpectrum, averaged over %d events", numEvents);
}



/** 
 * @brief Set each frequency bin of each waveform to the FourierBuffer derived TSpectrum amplitude, without chaning the phase.
 * 
 * @param fEv is the event to be filtered
 */
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




























 /** 
 * @brief Constructor for theh RayleighMonitor
 *
 * This class is the mother of all filter classes which interface with the FourierBuffer.
 * @param numEvents is the number of events over which the FourierBuffer is to track frequency amplitudes
 */
Acclaim::Filters::RayleighMonitor::RayleighMonitor(int numEvents) : fourierBuffer(numEvents) {
  fNumEvents = numEvents;
  fDescription = TString::Format("Decides whether or not to filter events based on characteristics of the event amplitude and Rayleigh distribution over %d events", fNumEvents);
  fNumOutputs = 6;
  fOutputAnt = 4;
  fOutputPol = AnitaPol::kHorizontal;
}




/** 
 * @brief Wrapper for FourierBuffer::add(FilteredAnitaEvent*)
 * 
 * @param fEv is the event to pass to FourierBuffer
 */
void Acclaim::Filters::RayleighMonitor::process(FilteredAnitaEvent* fEv){  
  fourierBuffer.add(fEv);
}



/** 
 * @brief returns the number of doubles each output array 
 *
 * Since the FourierBuffer produces a high density of output, the standard FilterOperation i/o is not particularly useful.
 * Therefore the associated functions are not well tested.
 *
 * @param i is the element of the output array
 * 
 * @return the length of the output 
 */
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



/** 
 * @brief maps the output array index to a name
 * 
 * @param i is the output index.
 * 
 * @return a pointer to a c-style string containing the name of the output, or NULL if i is invalid.
 */
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




/** 
 * @brief puts the ith output buffer into the array v.
 * 
 * @param i us the output index
 * @param v points to the array into which the output should be written
 */
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






















/** 
 * @brief Constructor for the RayleighFilter
 *
 * This filter operation actually takes some action based on the values calculated by the FourierBuffer.
 * This is by far the most developed filter in this library (although, again, the hard work is all inside FourierBuffer).
 *
 * The parameters passed in the constructor here too complicated to be explained in a single comment line, so here goes:
 * @param channelChiSquareCdfThresh
 * If one divides each frequency bin amplitude squared, by the rayleigh amplitude squared and histograms them
 * (and the amplitudes are well behaved) you should find the distibution is a chi-square with 2 degrees of freedom.
 * (This is because underlying the Rayleigh distribution you have two gaussian distributed degrees of freedom)
 * By histogramming all the frequency amplitudes in a channel one can find outliers.
 
 * @param chiSquarePerDofThreshold:
 * The rayleigh amplitude chi square per degree of freedom is evaluated assuming the TSpectrum amplitude.
 * If the distribution of amplitudes is not well described by a Rayligh distribution with a amplitude equal
 * to the TSpectrum amplitude, that frequency bin is filtered.
 *
 * @param numEvents is the number of events over which to track frequency amplitudes in the FourierBuffer
 */
Acclaim::Filters::RayleighFilter::RayleighFilter(double channelChiSquareCdfThresh, double chiSquarePerDofThreshold, int numEvents) : RayleighMonitor(numEvents), fChiSquarePerDofThreshold(chiSquarePerDofThreshold), fChanChiSquareCdfThreshold(channelChiSquareCdfThresh)
{
  fRandy = new TRandom3(1234); // seed will be reset on a per event basis using the eventNumber
  fDescription = TString::Format("Tracks frequency bin amplitudes over %d events, the chi squared threshold is %lf", fNumEvents, fChiSquarePerDofThreshold);
  
  fChanChiSquareThreshold = -2*TMath::Log(1 - fChanChiSquareCdfThreshold);
}


/** 
 * @brief Destructor
 */
Acclaim::Filters::RayleighFilter::~RayleighFilter()
{
  if(fRandy){
    delete fRandy;
    fRandy = NULL;
  }
}


/** 
 * Applies the RayleighFilter operation
 * 
 * @param fEv is the event to filter
 */
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

      // some debugging statements
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
