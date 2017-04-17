/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             A FIFO queue for frequency domain information
*************************************************************************************************************** */

#ifndef FOURIERBUFFER_H
#define FOURIERBUFFER_H

#include "AnalysisWaveform.h"

#include "RayleighInfo.h"
#include "RiceInfo.h"
#include "TH1D.h"
#include <complex>
#include <list>
#include "AnitaConventions.h"

class TSpectrum;
class FilteredAnitaEvent;

namespace Acclaim
{

  class RayleighHist;
  class TGraphFB;
  
  /**
   * @class FourierBuffer
   * @brief A class to hold frequency domain info in memory
   * I've removed the mag/phase info for now as it's probably a bit gratuitous
   */
  class FourierBuffer {

    friend class RayleighHist;
  public:

    virtual ~FourierBuffer();
    explicit FourierBuffer(Int_t theBufferSize=1000);

    size_t add(const FilteredAnitaEvent* fEv);    
    
    const RayleighHist* getRayleighDistribution(Int_t ant, AnitaPol::AnitaPol_t pol, Int_t freqBin=-1) const {return hRays[pol][ant].at(freqBin >= 0 ? freqBin : fDrawFreqBin);}
    TGraphFB* getAvePowSpec_dB(Int_t ant, AnitaPol::AnitaPol_t pol, int lastNEvents = -1) const;
    TGraphFB* getAvePowSpec(Int_t ant, AnitaPol::AnitaPol_t pol, int lastNEvents = -1) const;
    TGraphFB* getBackground_dB(Int_t ant, AnitaPol::AnitaPol_t pol, int lastNEvents = -1) const;
    TGraphFB* getBackground(Int_t ant, AnitaPol::AnitaPol_t pol, int lastNEvents = -1) const;
    TGraphFB* getReducedChiSquaresOfRayelighDistributions(Int_t ant, AnitaPol::AnitaPol_t pol) const;
    
  private:
    Int_t bufferSize;
    Int_t removeOld();
    void initVectors(int n);

    // list of events
    std::list<UInt_t> eventNumbers;
    std::list<Int_t> runs;    

    std::list<std::vector<double> > powerRingBuffers[AnitaPol::kNotAPol][NUM_SEAVEYS];

    // vectors of frequency bins
    std::vector<double> sumPowers[AnitaPol::kNotAPol][NUM_SEAVEYS];
    std::vector<RayleighHist*> hRays[AnitaPol::kNotAPol][NUM_SEAVEYS];

    std::vector<double> chiSquares[AnitaPol::kNotAPol][NUM_SEAVEYS];
    std::vector<int> ndfs[AnitaPol::kNotAPol][NUM_SEAVEYS];

    // it turns out that initialising a TF1 is very slow,
    // so I initialize a master here (owned by FourierBuffer) and clone others from this one.    
    TF1* fRay;
    double df;
    mutable TSpectrum* spectrums[AnitaPol::kNotAPol][NUM_SEAVEYS]; // to estimate the background
    bool doneVectorInit;
    int fDrawFreqBin;
  };






  // little class for some GUI i/o magic
  class TGraphFB : public TGraphAligned {
  public:
    TGraphFB(const FourierBuffer* theFb=NULL, Int_t theAnt=-1, AnitaPol::AnitaPol_t thePol=AnitaPol::kNotAPol,
	     int n=0) : TGraphAligned(n)
    {
      fb = theFb;
      ant = theAnt;
      pol = thePol;
    }
    virtual ~TGraphFB(){;}
    void ExecuteEvent(Int_t event, Int_t x, Int_t y);
  private:
    const FourierBuffer* fb; // pointer to parent
    Int_t ant;
    AnitaPol::AnitaPol_t pol;
  };
  
}



#endif //FOURIERBUFFER_H
