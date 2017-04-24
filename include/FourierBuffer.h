/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             A FIFO queue for frequency domain information
*************************************************************************************************************** */

#ifndef FOURIERBUFFER_H
#define FOURIERBUFFER_H

#include "AnalysisWaveform.h"

#include "TH1D.h"
#include <complex>
#include <deque>
#include "AnitaConventions.h"
#include "TObject.h"

class TSpectrum;
class FilteredAnitaEvent;
class TPad;
class UsefulAnitaEvent;

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

    enum SummaryOption_t{
      None,
      Chisquare,
      ReducedChisquare,
      NDF,
      RayleighAmplitude,
      Prob
    };

    
    virtual ~FourierBuffer();
    explicit FourierBuffer(Int_t theBufferSize=1000);

    size_t add(const FilteredAnitaEvent* fEv);

    Double_t getProb(AnitaPol::AnitaPol_t pol, Int_t ant, Int_t freqBin) const{
      return probs[pol][ant][freqBin];
    }
    Double_t getSpectrumAmp(AnitaPol::AnitaPol_t pol, Int_t ant, Int_t freqBin) const{
      return spectrumAmplitudes[pol][ant][freqBin];
    }

    
    const RayleighHist* getRayleighDistribution(Int_t ant, AnitaPol::AnitaPol_t pol, Int_t freqBin) const {return hRays[pol][ant].at(freqBin);}
    TGraphFB* getAvePowSpec_dB(Int_t ant, AnitaPol::AnitaPol_t pol, int lastNEvents = -1) const;
    TGraphFB* getAvePowSpec(Int_t ant, AnitaPol::AnitaPol_t pol, int lastNEvents = -1) const;
    TGraphFB* getBackground_dB(Int_t ant, AnitaPol::AnitaPol_t pol, int lastNEvents = -1) const;
    TGraphFB* getBackground(Int_t ant, AnitaPol::AnitaPol_t pol, int lastNEvents = -1) const;
    TGraphFB* getReducedChiSquaresOfRayelighDistributions(Int_t ant, AnitaPol::AnitaPol_t pol) const;
    void drawSummary(TPad* pad, SummaryOption_t) const;
    unsigned getN(int ant, AnitaPol::AnitaPol_t pol) const{return sumPowers[pol][ant].size();}
    unsigned getCurrentBufferSize();
    const std::vector<double>& getChiSquares(int ant, AnitaPol::AnitaPol_t pol) const {return chiSquares[pol][ant];};
    const std::vector<int>& getNDFs(int ant, AnitaPol::AnitaPol_t pol) const {return ndfs[pol][ant];};    
    int getNumEventsInBuffer() const {return eventsInBuffer;}
    void setForceLoadHistory(bool f) const {fForceLoadHistory=f;}
    bool isASelfTriggeredBlastOrHasSurfSaturation(const UsefulAnitaEvent* useful);

    const FourierBuffer* getAddress(){return this;}
  protected:
    Int_t bufferSize;
    Int_t removeOld();
    void initVectors(int n, double df);

    // list of events
    std::deque<UInt_t> eventNumbers;
    std::deque<Int_t> runs;

    std::deque<std::vector<double> > powerRingBuffers[AnitaPol::kNotAPol][NUM_SEAVEYS];

    // vectors of frequency bins
    std::vector<double> sumPowers[AnitaPol::kNotAPol][NUM_SEAVEYS];
    std::vector<RayleighHist*> hRays[AnitaPol::kNotAPol][NUM_SEAVEYS];

    std::vector<double> chiSquares[AnitaPol::kNotAPol][NUM_SEAVEYS];
    std::vector<int> ndfs[AnitaPol::kNotAPol][NUM_SEAVEYS];
    std::vector<double> fitAmplitudes[AnitaPol::kNotAPol][NUM_SEAVEYS];
    std::vector<double> spectrumAmplitudes[AnitaPol::kNotAPol][NUM_SEAVEYS];
    std::vector<double> probs[AnitaPol::kNotAPol][NUM_SEAVEYS];

    // it turns out that initialising a TF1 is very slow,
    // so I initialize a master here (owned by FourierBuffer) and clone others from this one.    
    TF1* fRay;

    bool doneVectorInit;


    void automagicallyLoadHistory(const FilteredAnitaEvent* fEv);
    bool fCurrentlyLoadingHistory;
    mutable bool fForceLoadHistory;
    
    int eventsInBuffer;
    
    double df; // frequency bin width (from AnalysisWaveform so probably in GHz)
    double fMinFitFreq;
    double fMaxFitFreq;
    mutable TSpectrum* fSpectrum; // to estimate the background    
    double fMinSpecFreq;
    double fMaxSpecFreq;
    void getSpectrum(double* y, int n) const;


    mutable TPad* summaryPads[NUM_SEAVEYS]; // for drawSummary
    mutable std::vector<TGraphFB> grReducedChiSquares[AnitaPol::kNotAPol]; // for drawSummary
    mutable std::vector<TGraphFB> grChiSquares[AnitaPol::kNotAPol]; // for drawSummary
    mutable std::vector<TGraphFB> grNDFs[AnitaPol::kNotAPol]; // for drawSummary
    mutable std::vector<TGraphFB> grSpectrumAmplitudes[AnitaPol::kNotAPol]; // for drawSummary
    mutable std::vector<TGraphFB> grAmplitudes[AnitaPol::kNotAPol]; // for drawSummary
    mutable std::vector<TGraphFB> grProbs[AnitaPol::kNotAPol]; // for drawSummary

    
    TGraphFB* getSelectedGraphForSummary(SummaryOption_t choice, int ant, AnitaPol::AnitaPol_t pol) const{
      switch(choice){
      case None:
	return NULL;
      case Chisquare:
	return &grChiSquares[pol][ant];
      case ReducedChisquare:
	return &grReducedChiSquares[pol][ant];
      case NDF:
	return &grNDFs[pol][ant];
      case RayleighAmplitude:
	return &grAmplitudes[pol][ant];
      case Prob:
	return &grProbs[pol][ant];
      }
    }    
  };


  



  // little class for some GUI i/o magic
  class TGraphFB : public TGraphAligned {
    friend class FourierBuffer;
  public:
    enum EDoubleClickOption{
      kDrawRayleigh,
      kDrawCopy
    };
    
    TGraphFB(const FourierBuffer* theFb=NULL, Int_t theAnt=-1, AnitaPol::AnitaPol_t thePol=AnitaPol::kNotAPol,
	     int n=0) : TGraphAligned(n), fDoubleClickOption(kDrawCopy), fDerivedFrom(NULL)
    {
      SetTitle("");
      fb = theFb;
      ant = theAnt;
      pol = thePol;
    }
    explicit TGraphFB(const TGraphFB* gr) : TGraphAligned(gr->GetN(), gr->GetX(), gr->GetY()), fDoubleClickOption(gr->fDoubleClickOption), fDerivatives(gr->fDerivatives)
    {
      fb = gr->fb;
      ant = gr->ant;
      pol = gr->pol;
      SetLineColor(gr->GetLineColor());
      SetLineStyle(gr->GetLineStyle());      
      fDerivedFrom = gr->fDerivedFrom;
    }
    virtual ~TGraphFB(){;}
    virtual void ExecuteEvent(Int_t event, Int_t x, Int_t y);
    void drawCopy() const;
    void drawRayleighHistNearMouse(int x, int y) const;
  private:
    const FourierBuffer* fb; // pointer to parent, don't delete
    Int_t ant;
    AnitaPol::AnitaPol_t pol;
    EDoubleClickOption fDoubleClickOption;

    // just raw pointers, no one "owns" anyone, rely on FourierBuffer and ROOT's garbage collection as appropriate
    const TGraphFB* fDerivedFrom; // e.g. spectrum is derived from amplitudes, want to draw them together
    std::vector<TGraphFB*> fDerivatives; // amplitudes point to spectrum
    ClassDef(Acclaim::TGraphFB, 0);
  };

}



#endif //FOURIERBUFFER_H
