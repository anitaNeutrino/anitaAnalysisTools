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
   * @brief A a glorified ring buffer of frequency amplitudes with a TSpectrum to look for CW spikes
   *
   * More specifically, it holds an array of specialised histograms (RayleighHist), one for each channel.
   * Each of those histograms has an internal RingBuffer, the histogram is updated as events are added/removed.
   * The RayleighHists can try to fit a Rayleigh distribution, or can quantify how well the distribution matches a given frequency,
   *
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
    Double_t getBackgroundSpectrumAmp(AnitaPol::AnitaPol_t pol, Int_t ant, Int_t freqBin) const{
      return spectrumAmplitudes[pol][ant][freqBin];
    }
    
    // void getChanChiSquareAndNDF(AnitaPol::AnitaPol_t pol, Int_t ant,
    //     			double& chiSquare, int& ndf) const{
    //   chiSquare = chanChisquare[pol][ant];
    //   ndf = chanNdf[pol][ant];      
    // }
    const std::vector<double>& getPowerRingBufferBack(AnitaPol::AnitaPol_t pol, int ant){
      return powerRingBuffers[pol][ant].back();
    }

    const RayleighHist* getRayleighDistribution(Int_t ant, AnitaPol::AnitaPol_t pol,
						Int_t freqBin) const {
      return hRays[pol][ant].at(freqBin);
    }

    unsigned getN(int ant, AnitaPol::AnitaPol_t pol) const {
      return sumPowers[pol][ant].size();
    }
    
    TGraphFB* getAvePowSpec_dB(Int_t ant, AnitaPol::AnitaPol_t pol, int lastNEvents = -1) const;
    TGraphFB* getAvePowSpec(Int_t ant, AnitaPol::AnitaPol_t pol, int lastNEvents = -1) const;
    TGraphFB* getBackground_dB(Int_t ant, AnitaPol::AnitaPol_t pol, int lastNEvents = -1) const;
    TGraphFB* getBackground(Int_t ant, AnitaPol::AnitaPol_t pol, int lastNEvents = -1) const;
    TGraphFB* getReducedChiSquaresOfRayelighDistributions(Int_t ant, AnitaPol::AnitaPol_t pol) const;
    void drawSummary(TPad* pad, SummaryOption_t) const;
    unsigned getCurrentBufferSize() const;

    const std::vector<double>& getChiSquares(int ant,
					     AnitaPol::AnitaPol_t pol) const {
      return chiSquares[pol][ant];
    }
    const std::vector<double>& getChiSquaresRelativeToSpectrum(int ant,
							       AnitaPol::AnitaPol_t pol) const {
      return chiSquaresRelativeToSpectrum[pol][ant];
    }
    const std::vector<int>& getNDFs(int ant,
				    AnitaPol::AnitaPol_t pol) const{
      return ndfs[pol][ant];
    }
    const std::vector<double>& getRayleighAmplitudes(int ant,
						     AnitaPol::AnitaPol_t pol) const {
      return fitAmplitudes[pol][ant];
    }
    const std::vector<double>& getBackgroundSpectrumAmplitudes(int ant,
							       AnitaPol::AnitaPol_t pol) const {
      return spectrumAmplitudes[pol][ant];
    }
    const std::vector<double>& getProbabilities(int ant, AnitaPol::AnitaPol_t pol) const {
      return probs[pol][ant];
    }

    const std::vector<double>& getChanChiSquares(int ant, AnitaPol::AnitaPol_t pol) const{
      return chanChiSquares[pol][ant];
    }

    TH1D* makeChanChisquareHist(AnitaPol::AnitaPol_t pol, Int_t ant, TPad* pad, const char* drawOpt) const; // for summary
    
    void getMeanVarChanChiSquares(int ant, AnitaPol::AnitaPol_t pol, double mean, double var) const{
      mean = meanChanChiSquare[pol][ant];
      var = varChanChiSquare[pol][ant];
    }
    
    int getNumEventsInBuffer() const {
      return eventsInBuffer;
    }
    void setForceLoadHistory(bool f) const {
      fForceLoadHistory=f;
    }

    bool isAlfaBandpassed(int ant, AnitaPol::AnitaPol_t pol) const {
      return (ant == 4 && pol == AnitaPol::kHorizontal) || (ant == 12 && pol == AnitaPol::kHorizontal);
    }

  protected:
    Int_t bufferSize;
    Int_t removeOld();
    void initVectors(int n, double df);

    // list of events
    std::deque<UInt_t> eventNumbers;
    std::deque<Int_t> runs;

    std::deque<std::vector<double> > powerRingBuffers[AnitaPol::kNotAPol][NUM_SEAVEYS];

    // utility function to sanitize vector/graph initialization
    void initGraphAndVector(std::vector<double> vec[][NUM_SEAVEYS],
			    std::vector<TGraphFB>* gr,
			    int n, double df, double defaultVal);

    TGraphFB* getSelectedGraphForSummary(SummaryOption_t choice, int ant, AnitaPol::AnitaPol_t pol) const;

    
    // vectors of frequency bins
    std::vector<double> sumPowers[AnitaPol::kNotAPol][NUM_SEAVEYS];
    std::vector<RayleighHist*> hRays[AnitaPol::kNotAPol][NUM_SEAVEYS];

    std::vector<double> chiSquares[AnitaPol::kNotAPol][NUM_SEAVEYS];    
    std::vector<double> chiSquaresRelativeToSpectrum[AnitaPol::kNotAPol][NUM_SEAVEYS];    
    std::vector<int> ndfs[AnitaPol::kNotAPol][NUM_SEAVEYS];
    std::vector<double> fitAmplitudes[AnitaPol::kNotAPol][NUM_SEAVEYS];
    std::vector<double> spectrumAmplitudes[AnitaPol::kNotAPol][NUM_SEAVEYS];
    std::vector<double> lastAmps[AnitaPol::kNotAPol][NUM_SEAVEYS];    
    std::vector<double> probs[AnitaPol::kNotAPol][NUM_SEAVEYS];


    std::vector<double> chanChiSquares[AnitaPol::kNotAPol][NUM_SEAVEYS];
    double meanChanChiSquare[AnitaPol::kNotAPol][NUM_SEAVEYS];
    double varChanChiSquare[AnitaPol::kNotAPol][NUM_SEAVEYS];


    // it turns out that initialising a TF1 is very slow,
    // so I initialize a master here (owned by FourierBuffer) and clone others from this one.    
    TF1* fRay;

    Int_t fNumSkipped; //!< Incremented when skipping a payload blast or SURF saturation event, currently only for debugging
    bool doneVectorInit; //!< Do we need to allocate a bunch of memory? (Do this dynamically to avoid MagicDisplay being slow)


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
    void getBackgroundSpectrum(double* y, int n) const;

    // double eventPower[AnitaPol::kNotAPol][NUM_SEAVEYS];
    // double expectedThermalPower[AnitaPol::kNotAPol][NUM_SEAVEYS];

    mutable TPad* summaryPads[NUM_SEAVEYS]; // for drawSummary
    mutable std::vector<TGraphFB> grReducedChiSquares[AnitaPol::kNotAPol]; // for drawSummary
    mutable std::vector<TGraphFB> grChiSquares[AnitaPol::kNotAPol]; // for drawSummary
    mutable std::vector<TGraphFB> grChiSquaresRelativeToSpectrum[AnitaPol::kNotAPol]; // for drawSummary        
    mutable std::vector<TGraphFB> grReducedChiSquaresRelativeToSpectrum[AnitaPol::kNotAPol]; // for drawSummary    
    mutable std::vector<TGraphFB> grNDFs[AnitaPol::kNotAPol]; // for drawSummary
    mutable std::vector<TGraphFB> grSpectrumAmplitudes[AnitaPol::kNotAPol]; // for drawSummary
    mutable std::vector<TGraphFB> grAmplitudes[AnitaPol::kNotAPol]; // for drawSummary
    mutable std::vector<TGraphFB> grLastAmps[AnitaPol::kNotAPol]; // for drawSummary    
    mutable std::vector<TGraphFB> grProbs[AnitaPol::kNotAPol]; // for drawSummary
    
    
  };


  /**
   * @class TGraphFB 
   * @brief A little class for some GUI magic (for MagicDisplay)
   */
  class TGraphFB : public TGraphAligned {
    friend class FourierBuffer;
  public:

    // Utility function to set fDerives from and fDerivatives (i.e. the drawing ownership)
    static void setDrawingDependencies(const std::vector<Acclaim::TGraphFB*> grs);
    
    enum EDoubleClickOption{
      kDrawRayleigh,
      kDrawCopy
    };

    /** 
     * @brief Constructor
     * 
     * @param theFb the fourier buffer, the parent of this TGraphFB
     * @param theAnt is the antenna number (0-47)
     * @param thePol is the polarisation
     * @param n is the number of points to construct the TGraphFB with
     */
    TGraphFB(const FourierBuffer* theFb=NULL, Int_t theAnt=-1, AnitaPol::AnitaPol_t thePol=AnitaPol::kNotAPol,
	     int n=0) : TGraphAligned(n), fDoubleClickOption(kDrawCopy), fDerivedFrom(NULL)
    {
      SetTitle("");
      fb = theFb;
      ant = theAnt;
      pol = thePol;
    }

    /** 
     * @brief Creates a copy of the passed TGraphFB
     * 
     * @param gr is the TGraphFB to copy 
     */    
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
