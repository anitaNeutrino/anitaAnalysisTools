#ifndef RAYLEIGH_HIST_H
#define RAYLEIGH_HIST_H

#include "TH1D.h"
#include "RingBuffer.h" // dumb class I put in eventReaderRoot ages ago, looks like I can actually use it again!
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

class TF1;

namespace Acclaim{
  class FourierBuffer;

  class RayleighHist : public TH1D {

    friend class FourierBuffer;

  public:

    enum class FitMethod{
      TF1, // slow 
      Minuit, // less slow
      Scan, // probably the fastest useful option, won't do errors
      JustEvalGuess, // fastest but obviously the least accurate
      Default = Scan
      // Default = JustEvalGuess
    };
    
    RayleighHist(FourierBuffer* fb=NULL, const char* name = "", const char* title = "");
    virtual ~RayleighHist();

    virtual void Draw(Option_t* opt="");
    int add(double newAmp); // Main interaction method
    void Fit(Double_t& rayleighAmplitude, Double_t& chiSquare, Int_t& ndf);
    
    

    void SetFreqBinToDraw(Int_t freqBin); // *MENU*
    

    static void guessMaxBinLimitAndSigmaFromMean(double meanAmp, double& maxAmp, double& sigmaGuess, double fracOfEventsInsideMaxAmp);

    inline static double EvalRayleigh(double normalization, double amplitude, double x){
      return (normalization*x/(amplitude*amplitude))*exp(-x*x/(2*amplitude*amplitude));
    }
    
  protected:
    RingBuffer amplitudes; //!< Tracks all the amplitudes

    virtual int Fill(double amp, double sign=1); //!< Fill the histogram, this is called by add(double)
    bool axisRangeOK() const; //!< Checks current axis range is reasonable
    void rebinAndRefill(double meanAmp); //!< Dynamically rebin and refill histogram with contents of RingBuffer of amplitudes
    
    FourierBuffer* fParent; //!< Daddy
    TF1* fRay; //!< Pointer to the Rayeligh TF1 cloned from parent FourierBuffer

    Int_t fNumFitParams; //!< Will be equal to one as we only try and fit the amplitude (normalization is fixed by fNumEvents and bin width)
    ROOT::Math::Minimizer* fMinimizer; //!< The minuit minimizer object
    ROOT::Math::Functor fChiSquaredFunc; //!< For minuit interface, will point to getRayelighChiSquare(const double*)
    
    double getRayleighChiSquare(const double* params); // for the minuit fitter
    std::vector<double> theFitParams; // for minuit fitter interface
    std::vector<double> theFitParamsSteps; // for minuit fitter interface (I think minuit ignores this, but oh well)
    
    double fracOfEventsWanted; //!< Fraction of events to be in the histogram bin limits using the guessed amplitude (don't set to 1 as this requires an infinite axis range)
    Int_t risingEdgeBins; //!< Number of bins between 0 and where we guess the histogram peak is, for dynamic rebinning
    Int_t fDrawFreqBin;
    double freqMHz; //!< The frequency (MHz) of this Rayleigh distribution
    Int_t numEvents; //!< Tracks the number of events in the RingBuffer/histogram (faster than integral)




    FitMethod fitMethod;

    // caching for fit functions
    double fBinWidth;
    double fRayleighNorm;
    void prepareRayleighFitCache();
    Int_t fNx;
    std::vector<double> rayleighChiSquareCache;
    
    ClassDef(RayleighHist, 0);
    
  };
}
#endif
