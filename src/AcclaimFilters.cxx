#include "AcclaimFilters.h"
#include "AnalysisWaveform.h"
#include "BasicFilters.h" // Cosmin's example filters
#include "FFTWComplex.h"
#include <iostream>


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
