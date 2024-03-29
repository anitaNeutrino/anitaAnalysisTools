#ifndef RAYLEIGH_HIST_H
#define RAYLEIGH_HIST_H

#include "TH1D.h"
#include "RingBuffer.h"
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include "TF1.h"


class TGraph;

namespace Acclaim{
  class FourierBuffer;

  /*!
    A class to histogram frequency amplitudes on the fly and "fit" a Rayleigh distribution to the result
  */
  class RayleighHist : public TH1D {

    friend class FourierBuffer;

  public:

    typedef enum {
      kTF1, // slow 
      kMinuit, // less slow
      kScan, // probably the fastest useful option, won't do errors
      kJustEvalGuess, // fastest but obviously the least accurate
      kAdaptive, // Tries just evaluating the guess, if that's good enough (chiSquare < 2), stops there, otherwise does a scan
      kDefault = kAdaptive
      // Default = JustEvalGuess
    } FitMethod;

    typedef enum {
      kPoisson,
      kPearson
    } ChiSquareErrorMethod;
    
    RayleighHist(FourierBuffer* fb=NULL, const char* name = "", const char* title = "");
    virtual ~RayleighHist();

    virtual void Draw(Option_t* opt="");
    virtual bool add(double newAmp); //!< Input amplitudes events
    void getRayleighFitParams(double& rayAmp, double& chiSquare, int& ndf); //!< Output Rayleigh distribution parameters

    void fitRayleigh(bool forGuiUpdateTF1=true); // *MENU* Fit the Rayleigh distribution using the selected fit method

    // Mostly for validating fit sanity with a GUI
    void fitRayleighAdaptive(bool forGuiUpdateTF1=true); // *MENU* Tries to evaluate the guess, if it's good enough move on, otherwise scan.
    void fitRayleighScan(bool forGuiUpdateTF1=true); // *MENU* Fit the Rayleigh distribution scanning through amplitude values and choosing the lowest chiSquare residual
    void fitRayleighJustEvalGuess(bool forGuiUpdateTF1=true); // *MENU* "Fit" the rayleigh distribution just using the guess from the mean/integral of the histogram
    void fitRayleighTF1(); // *MENU* Fit using the TF1
    void fitRayleighMinuit(bool forGuiUpdateTF1=true); // *MENU* Fit using Minuit
    
    static void guessMaxBinLimitAndSigmaFromMean(double meanAmp, double& maxAmp, double& sigmaGuess, double fracOfEventsInsideMaxAmp);
    
    
    inline double getOneMinusCDF(double amp, double distAmp = -1) const { //!< This is the probability of getting this amplitude (amp) or higher
      distAmp = distAmp < 0 ? fRayleighAmplitude : distAmp; // use this histograms rayleigh distribution amplitude if one wasn't specified      
      return exp((-0.5*amp*amp)/(distAmp*distAmp));
    }
    
    inline double getCDF(double amp, double distAmp = -1) const{ // This is the fraction of amplitudes lower than amp 
      distAmp = distAmp < 0 ? fRayleighAmplitude : distAmp; // Use this histograms rayleigh distribution amplitude if one wasn't specified      
      return 1 - getOneMinusCDF(amp, distAmp);
    }
    inline double getAmplitude() const {return fRayleighAmplitude;}
    
  protected:
    RingBuffer amplitudes; //!< Tracks all the amplitudes
    
    virtual int Fill(double amp, double sign=1); //!< Fill the histogram, this is called by add(double)
    Int_t fNumNonEmptyBins; //!< Cache number of non empty bins
    Int_t fNx; //!< The number of bins (faster than GetNbinsX())
    Int_t fNumEvents; //!< Tracks the number of events in the RingBuffer/histogram (faster than integral)    
    std::vector<double> binCentres;
    std::vector<double> squaredBinCentres;    
    std::vector<int> binValues; // cache histogram bin contents, should be integers
    std::vector<double> squaredBinErrors;    
    



    // virtual int Fill(double amp, double sign=1); //!< Fill the histogram, this is called by add(double)
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
    double freqMHz; //!< The frequency (MHz) of this Rayleigh distribution

    TGraph* grLastAddedAmp; //!< A pretty visual representation of the last added amplitude

    Int_t fNDF;
    Double_t fChiSquare;
    Double_t fRayleighAmplitude;

    void updateParamsOfTF1();    
    std::vector<Double_t> fParamsTF1;
    
    FitMethod fitMethod;
    ChiSquareErrorMethod chiSquareErrorMethod;
    const Int_t fFitEveryNAdds;
    Int_t fNumAddsMod10;

    // caching for fit functions
    double fBinWidth;
    double fRayleighAmpGuess;
    double fRayleighNorm;
    // Int_t fNx;
    // std::vector<double> binCentres;
    // std::vector<double> squaredBinCentres;    
    // std::vector<int> binValues; // cache histogram bin contents, should be integers
    // std::vector<double> squaredBinErrors;    

    ClassDef(RayleighHist, 0);
  };


  
}
#endif
