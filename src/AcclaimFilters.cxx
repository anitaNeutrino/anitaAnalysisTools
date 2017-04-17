
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
 * This function adds my custom strategies to a map of TString to strategies...
 */
void Acclaim::Filters::appendFilterStrategies(std::map<TString, FilterStrategy*>& filterStrats, bool saveOutput){

  // first make the operations

  ALFAFilter* alfaFilter = new ALFAFilter(); // should always use this unless you have a good reaon...
  
  Double_t widthMHz = 26;
  Double_t centreMHz = 260;
  Notch* n260 = new Notch(centreMHz-widthMHz, centreMHz+widthMHz);

  // Double_t widthMHz = 26;
  centreMHz = 370;
  Notch* n370 = new Notch(centreMHz-widthMHz, centreMHz+widthMHz);


  // then make the strategies
  
  FilterStrategy* stupidNotchStrat = new FilterStrategy();
  stupidNotchStrat->addOperation(alfaFilter, saveOutput);  
  stupidNotchStrat->addOperation(n260, saveOutput);
  stupidNotchStrat->addOperation(n370, saveOutput);  

  filterStrats["BrickWallSatellites"] = stupidNotchStrat;

  // FilterStrategy* spikeKiller = new FilterStrategy();
  // spikeKiller->addOperation(alfaFilter);
  // SpikeSuppressor* ss = new SpikeSuppressor(3, 1000);
  // spikeKiller->addOperation(ss, saveOutput);
  // filterStrats["SpikeSuppressor"] = spikeKiller;

  FilterStrategy* justRm = new FilterStrategy();
  RayleighMonitor* rm = new RayleighMonitor(1500);
  justRm->addOperation(rm, saveOutput);
  filterStrats["RayleighMonitor"] = justRm;  
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
































Acclaim::Filters::RayleighMonitor::RayleighMonitor(int numEvents) : fourierBuffer(numEvents) {
  fNumEvents = numEvents;
}




void Acclaim::Filters::RayleighMonitor::process(FilteredAnitaEvent* fEv){
  
  fourierBuffer.add(fEv);

  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
    }
  }  
}

unsigned Acclaim::Filters::RayleighMonitor::outputLength(unsigned i) const{
  (void) i;
  if(i==2){
    return 1;
  }
  else{
    return AnitaPol::kNotAPol*NUM_SEAVEYS*131; // TODO, check this at runtime maybe?
  }
}





const char* Acclaim::Filters::RayleighMonitor::outputName(unsigned i) const{
  // so what do I want to output?
  switch(i){
  case 0:
    return "chiSquares";
  case 1:
    return "ndfs";    
  case 2:
    return "nEvents";
  default:
    return NULL;
  };
}


void Acclaim::Filters::RayleighMonitor::fillOutput(unsigned i, double* v) const{

  if(i==2){
    v[0] = fourierBuffer.getNumEventsInBuffer();
    return;
  }
  
  int outInd = 0;
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      const double* vals = 0;
      const int* valsInt = 0;
      unsigned n = 0;
      switch(i){
      case 0:
	vals = &fourierBuffer.getChiSquares(ant, pol)[0];
	n = fourierBuffer.getChiSquares(ant, pol).size();
      case 1:
	valsInt = &fourierBuffer.getNDFs(ant, pol)[0];
	n = fourierBuffer.getNDFs(ant, pol).size();
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

// /** The number of output variables (doubles or double arrays) defined by this operation */
// virtual unsigned nOutputs() const  { return 0; } 

// /** The name of the ith output variable */ 
// virtual const char *  outputName(unsigned i) const  { (void) i; return ""; } 

// /** The length of the ith output variable  (it's a double array of this size)*/ 
// virtual unsigned outputLength(unsigned i) const { (void) i; return 0; } 

// /** Fill the ith output */ 
// virtual void fillOutput(unsigned i, double * v) const{ (void) v; (void) i;  return; } 


