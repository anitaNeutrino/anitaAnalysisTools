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
#include <iostream>
#include "FourierBuffer.h"

class TPad;
class TGraphAligned;

namespace Acclaim
{
  /** @namespace Filters
   * @brief Contains all my filter classes, some of which have less than obvious names
   *
   * All frequencies contained in here are in GHz unless otherwise stated in a variable name.
   * (Sometimes the GHz is stated anyway for good measure).
   *
   */
  namespace Filters
  {

  const char* getCosminsFavouriteSineSubName();

  /** @namespace Bands
   * @brief Hard coded frequencies in GHz which can get used in any given filter
   *
   */  
    namespace Bands {
      // zero everything outside of these
      const double anitaHighPassGHz = 0.2;
      const double anitaLowPassGHz = 1.2;
      const double alfaLowPassGHz = 0.7;
    }

    void generateFilterStrategies(bool saveOutput = false);
    void appendFilterStrategies(std::map<TString, FilterStrategy*>& filterStrats, bool saveOutput = false);
    FilterStrategy* findStrategy(const std::map<TString, FilterStrategy*>& filterStrats, const TString& stratName);
    FilterStrategy* findDefaultStrategy(const TString& stratName);

    void makeFourierBuffersLoadHistoryOnNextEvent(FilterStrategy* fs);



  /** 
   * @Class Notch
   * @brief Filters contiguous frequency bins in the FourierDomain
   * 
   */
    class Notch: public UniformFilterOperation
    {
    protected:
      TString fTag, fDescription, fOutputName;
      Double_t fLowEdgeGHz, fHighEdgeGHz;
      mutable Double_t fPowerRemovedByNotch[AnitaPol::kNotAPol][NUM_SEAVEYS];
      virtual void processOne(AnalysisWaveform * g);
  
    public:
      Notch(Double_t lowEdgeGHz, Double_t highEdgeGHz);
  
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




  /** 
   * @Class RayleighMonitor
   * @brief Tracks the amplitudes of frequencies but doesn't do anything else, to be inherited from
   *
   * (Actually all the hard work is done by the FourierBuffer class, which is a rather complicated beastie.)
   */  
    class RayleighMonitor : public UniformFilterOperation {
    protected:
      int fNumEvents;
      FourierBuffer fourierBuffer;
      TString fDescription;
      unsigned fNumOutputs;
      AnitaPol::AnitaPol_t fOutputPol;
      int fOutputAnt;
    public:
      explicit RayleighMonitor(int numEvents);
      virtual const char * tag () const {return "RayleighMonitor";};
      virtual const char * description () const {return fDescription.Data();}
      virtual void processOne(AnalysisWaveform* wave)
      {
	(void) wave;
	std::cerr << "Error in " << __PRETTY_FUNCTION__
		  << " function not implemented, use process(FilteredAnitaEvent*) instead" << std::endl;
      }
      virtual void process(FilteredAnitaEvent* fEv);
      virtual unsigned outputLength(unsigned i) const;
      virtual unsigned nOutputs() const{return fNumOutputs;}
      virtual const char* outputName(unsigned i) const;
      virtual void fillOutput(unsigned i, double* v) const;

      // const FourierBuffer& getFourierBuffer() const{return fourierBuffer;}
      const FourierBuffer* getFourierBuffer() const{return &fourierBuffer;}
    };


  
  /** 
   * @Class RayleigFilter
   * @brief Tracks the amplitudes of frequencies, how well they track a Rayleigh distribution and filters if they do not 
   *
   * As with RayleighMonitor, most all the hard work in tracking these amplitudes is done by the FourierBuffer class.
   * However, the process function extracts the numbers from FourierBuffer and applies the filtering logic there.  
   */    
    class RayleighFilter : public RayleighMonitor {
    public:
      explicit RayleighFilter(double channelChiSquareCdfThresh, double chiSquarePerDofThresh, Int_t numEvents);
      virtual ~RayleighFilter();
      virtual void process(FilteredAnitaEvent* fEv);
      // virtual unsigned nOutputs() const {return 0;}
      virtual const char * tag () const {return "RayleighFilter";};
      virtual const char * description () const {return fDescription.Data();}
    protected:
      TRandom3* fRandy;
      // double fLog10ProbThreshold; //!< What was the probability
      double fChiSquarePerDofThreshold; // Remove frequency if our fit of the rayleigh amplitude is bad
      double fChanChiSquareCdfThreshold; // threshold in chiSquareCdf
      double fChanChiSquareThreshold; // threshold in chiSquare (computed from threshold in chiSquareCdf)
    };


  /** 
   * @Class SpectrumMagnitude
   * @brief Silly filtering class, don't use it.
   *
   * Forces magnitude of each frequency bin to match the TSpectrum derived magnitudes inside FourierBuffer
   * (thereby deweighting CW frequency bins). 
   * Don't use this.
   */      
    class SpectrumMagnitude : public RayleighMonitor
    {
    public:
      explicit SpectrumMagnitude(Int_t numEvents);
      virtual ~SpectrumMagnitude() {;}
      virtual void process(FilteredAnitaEvent* fEv);
      // virtual unsigned nOutputs() const {return 0;}
      virtual const char * tag () const {return "SpectrumMagnitude";};
      virtual const char * description () const {return fDescription.Data();}
    protected:
    };
    
    

  /** 
   * @Class UniformMagnitude
   * @brief Silly filtering class, don't use it.
   *
   * Forces magnitude of each frequency bin to be equal (thereby deweighting CW frequency bins).
   * Don't use this.
   */      
    class UniformMagnitude : public UniformFilterOperation {
    public:
      explicit UniformMagnitude();
      virtual ~UniformMagnitude() { ;}
      virtual void processOne(AnalysisWaveform* wf);
      // virtual unsigned nOutputs() const {return 0;}
      virtual const char * tag () const {return "UniformMagnitude";};
      virtual const char * description () const {return "Gives every frequency bin the same magnitude, keeping the phase constant";}
    };


    

  }
}

#endif
