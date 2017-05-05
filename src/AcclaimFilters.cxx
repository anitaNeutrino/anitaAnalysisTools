#include "AcclaimFilters.h"
#include "AnalysisWaveform.h"
#include "FilteredAnitaEvent.h"
#include "BasicFilters.h" // Cosmin's example filters
#include "FourierBuffer.h"
#include "FFTWComplex.h"
#include <iostream>
#include "TPad.h"
#include "RayleighHist.h"
#include "RootTools.h"



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

  // first make the operations

  double alfaLowPassFreq = 0.7;
  ALFASincFilter* alfaFilter = new ALFASincFilter(alfaLowPassFreq); // should always use some kind of alfa filter unless you have a good reason

  Double_t widthMHz = 26;
  Double_t centreMHz = 260;
  Notch* n260 = new Notch(centreMHz-widthMHz, centreMHz+widthMHz);

  // Double_t widthMHz = 26;
  centreMHz = 370;
  Notch* n370 = new Notch(centreMHz-widthMHz, centreMHz+widthMHz);

  double log10ProbThresh = -2.5;
  double reducedChiSquareThresh = 3;
  RayleighFilter* rf = new RayleighFilter(log10ProbThresh, reducedChiSquareThresh, 1500, alfaLowPassFreq);  

  // then make the strategies
  
  FilterStrategy* stupidNotchStrat = new FilterStrategy();
  stupidNotchStrat->addOperation(alfaFilter, saveOutput);
  stupidNotchStrat->addOperation(n260, saveOutput);
  stupidNotchStrat->addOperation(n370, saveOutput);

  filterStrats["BrickWallSatellites"] = stupidNotchStrat;

  FilterStrategy* fs = new FilterStrategy();
  fs->addOperation(alfaFilter, saveOutput);
  fs->addOperation(rf, saveOutput);
  filterStrats["RayleighFilter"] = fs;
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






Acclaim::Filters::Notch::Notch(Double_t lowEdgeMHz, Double_t highEdgeMHz){

  // Freq bins are currently in 10MHz steps
  fTag = TString::Format("notch%.lfMHzto%.lfMHz", lowEdgeMHz, highEdgeMHz);
  fDescription = TString::Format("Notch filter from %.lf MHz to %.lf MHz", lowEdgeMHz, highEdgeMHz);
  fLowEdgeMHz = lowEdgeMHz;
  fHighEdgeMHz = highEdgeMHz;

  for(int pol=0; pol < AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      fPowerRemovedByNotch[pol][ant] = 0;
    }
  }
}


void Acclaim::Filters::Notch::processOne(AnalysisWaveform * g){//, AnitaPol::AnitaPol_t pol, int ant){

  // std::cout << __PRETTY_FUNCTION__ << "\t" << g << std::endl;
  const double deltaFMHz = 1e3*g->deltaF();
  const int nf = g->Nfreq();

  // // std::cout << g << "\t" << g->deltaF() << "\t" << deltaFMHz << "\t" << nf << "\t" << std::endl;
  // Double_t removedPower = 0;
  
  FFTWComplex* theFreqs = g->updateFreq();
  for(int freqInd=0; freqInd < nf; freqInd++){
    const double freqMHz = deltaFMHz* freqInd;

    if(freqMHz >= fLowEdgeMHz && freqMHz < fHighEdgeMHz){
      // removedPower += theFreqs[freqInd].getAbsSq();
      theFreqs[freqInd] = 0;
    }
  }

  // if(pol!=AnitaPol::kNotAPol && ant >= 0 && ant < NUM_SEAVEYS){
  //   fPowerRemovedByNotch[pol][ant] = removedPower;
  // }
// return removedPower;
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
































Acclaim::Filters::RayleighMonitor::RayleighMonitor(int numEvents, double alfaLowPassFreqGHz) : fourierBuffer(numEvents, alfaLowPassFreqGHz) {
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
	  vals = &fourierBuffer.getSpectrumAmplitudes(ant, pol)[0];
	  n = fourierBuffer.getSpectrumAmplitudes(ant, pol).size();
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






















Acclaim::Filters::RayleighFilter::RayleighFilter(double log10ProbThreshold, double chiSquarePerDofThreshold, int numEvents, double alfaLowPassFreqGHz) : RayleighMonitor(numEvents, alfaLowPassFreqGHz), fLog10ProbThreshold(log10ProbThreshold), fChiSquarePerDofThreshold(chiSquarePerDofThreshold)
{
  fRandy = new TRandom3(1234); // seed will be reset on a per event basis using the eventNumber
  fDescription = TString::Format("Tracks frequency bin amplitudes over %d events", fNumEvents);
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

      FFTWComplex* theFreqs = wf->updateFreq();
      for(int freqInd=0; freqInd < nf; freqInd++){

	double probVal = fourierBuffer.getProb(pol, ant, freqInd);
	int ndf = fourierBuffer.getNDFs(ant, pol)[freqInd];
	double reducedChiSquare = fourierBuffer.getChiSquares(ant, pol)[freqInd]/ndf;

	if(reducedChiSquare > fChiSquarePerDofThreshold || probVal < fLog10ProbThreshold){
	  double specAmp = fourierBuffer.getSpectrumAmp(pol, ant, freqInd);
	  // std::cout << pol << "\t" << ant << "\t" << freqInd << "\t" << probVal << "\t" << specAmp << std::endl;
	  
	  double x1 = fRandy->Gaus(0, specAmp);
	  double x2 = fRandy->Gaus(0, specAmp);
	  double newAmp = TMath::Sqrt(x1*x1 + x2*x2);
	  
	  double phase = fRandy->Uniform(0, TMath::TwoPi());

	  // std::cout << probVal << "\t" << x1 << "\t" << x2 << "\t" << newAmp << "\t" << phase << "\t" << theFreqs[freqInd] << "\t";
	  
	  theFreqs[freqInd].setMagPhase(newAmp, phase);
	  

	  // std::cout << theFreqs[freqInd] << std::endl;
	}
	// else{
	//   std::cout << pol << "\t" << ant << "\t" << theFreqs[freqInd] << "\t" << probVal << "\t" << std::endl;
	// }
      }
    }
  }  
}
