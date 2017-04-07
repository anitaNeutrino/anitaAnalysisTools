/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             Potential filters to use in ANITA-3 analysis.
	     These all use the AnitaAnalysisFramework filter framework, which defines a bunch of virtual functions.
	     These are to be overloaded by very simple filter operations.
***********************************************************************************************************/

#ifndef ACCLAIM_FILTERS_H
#define ACCLAIM_FILTERS_H

#include "FilterStrategy.h"
#include "FilterOperation.h" // contains the filter operation syntax we must adhere to
#include "TString.h"
#include "TRandom3.h"
#include <map>

class TGraphAligned;

namespace Acclaim
{
  class FourierBuffer;

  namespace Filters
  {
    void appendFilterStrategies(std::map<TString, FilterStrategy*>& filterStrats, bool saveOutput = false); //!< Utility function for MagicDisplay
    FilterStrategy* findStrategy(const std::map<TString, FilterStrategy*>& filterStrats, const TString& stratName);


    // base notch class
    class Notch: public UniformFilterOperation
    {
    protected:
      TString fTag, fDescription, fOutputName;
      Double_t fLowEdgeMHz, fHighEdgeMHz;
      mutable Double_t fPowerRemovedByNotch[AnitaPol::kNotAPol][NUM_SEAVEYS];
      virtual void processOne(AnalysisWaveform * g);
  
    public:
      Notch(Double_t lowEdgeMHz, Double_t highEdgeMHz);
  
      virtual const char * tag () const {return fTag.Data();};
      virtual const char * description () const {return fDescription.Data();} 
      virtual unsigned nOutputs() const  { return AnitaPol::kNotAPol; } 
      virtual const char *  outputName(unsigned i) const {return i == AnitaPol::kHorizontal ? "powerRemovedByNotchHPol" : "powerRemovedByNotchVPol";}
      virtual unsigned outputLength(unsigned i) const {(void) i; return NUM_SEAVEYS;}
      virtual void fillOutput(unsigned i, double * v) const 
      {
	for(int ant=0; ant < NUM_SEAVEYS; ant++){
	  v[ant] = fPowerRemovedByNotch[i][ant];
	}
      }
      virtual void process(FilteredAnitaEvent* fe);
  
    };






    class SpikeSuppressor : public UniformFilterOperation {
    protected:
      double fSpikeThresh_dB;
      double fTimeScale;
      TRandom3 fRandy;
      TString fDescription;
      std::vector<std::vector<FourierBuffer> > fourierBuffers;

      TGraphAligned suppressSpikes(const TGraphAligned* grPower);
      TGraphAligned suppressSpikes(const TGraphAligned* grPower, const TGraphAligned* grBackground);      
      double interpolate_dB(double x, double xLow, double xHigh, double yLow, double yHigh);
      
    public:
      void setSeed(UInt_t seed){fRandy.SetSeed(seed);}
      
      SpikeSuppressor(double spikeThresh_dB, double timeScale);

      virtual const char * tag () const {return "SpikeSuppressor";};
      virtual const char * description () const {return fDescription.Data();}
      virtual void processOne(AnalysisWaveform* wave);
      virtual void process(FilteredAnitaEvent* fEv);


      
    };

    
    
  }
}

#endif
