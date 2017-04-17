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
#include "TObject.h"

class TSpectrum;
class FilteredAnitaEvent;
class TPad;

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
    void drawSummary(TPad* pad) const;
    unsigned getN(int ant, AnitaPol::AnitaPol_t pol) const{return sumPowers[pol][ant].size();}
    unsigned getCurrentBufferSize();
    const std::vector<double>& getChiSquares(int ant, AnitaPol::AnitaPol_t pol) const {return chiSquares[pol][ant];};
    const std::vector<int>& getNDFs(int ant, AnitaPol::AnitaPol_t pol) const {return ndfs[pol][ant];};    
    int getNumEventsInBuffer() const {return eventsInBuffer;}
    void loadHistory();
    
  protected:
    Int_t bufferSize;
    Int_t removeOld();
    void initVectors(int n, double df);

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

    bool doneVectorInit;
    int fDrawFreqBin;
    int eventsInBuffer;
    
    double df; // frequency bin width (from AnalysisWaveform so probably in GHz)
    mutable TSpectrum* spectrums[AnitaPol::kNotAPol][NUM_SEAVEYS]; // to estimate the background
    mutable TPad* summaryPads[NUM_SEAVEYS]; // for drawSummary
    std::vector<TGraphFB> grReducedChiSquares[AnitaPol::kNotAPol]; // for drawSummary
    std::vector<TGraphFB> grChiSquares[AnitaPol::kNotAPol]; // for drawSummary
    std::vector<TGraphFB> grNDFs[AnitaPol::kNotAPol]; // for drawSummary
    std::vector<TGraphFB> grAmplitudes[AnitaPol::kNotAPol]; // for drawSummary
  };






  // little class for some GUI i/o magic
  class TGraphFB : public TGraphAligned {
  public:
    enum EDoubleClickOption{
      kDrawRayleigh,
      kDrawCopy
    };
    
    TGraphFB(const FourierBuffer* theFb=NULL, Int_t theAnt=-1, AnitaPol::AnitaPol_t thePol=AnitaPol::kNotAPol,
	     int n=0) : TGraphAligned(n), fDoubleClickOption(kDrawCopy)
    {
      fb = theFb;
      ant = theAnt;
      pol = thePol;
    }
    virtual ~TGraphFB(){;}
    virtual void ExecuteEvent(Int_t event, Int_t x, Int_t y);
    void drawCopy() const;
    void drawRayleighHistNearMouse(int x, int y) const;// *MENU*
  private:
    const FourierBuffer* fb; // pointer to parent, don't delete
    Int_t ant;
    AnitaPol::AnitaPol_t pol;
    EDoubleClickOption fDoubleClickOption;
    ClassDef(Acclaim::TGraphFB, 0);
  };

}



#endif //FOURIERBUFFER_H
