#ifndef RAYLEIGH_HIST_H
#define RAYLEIGH_HIST_H

#include "TH1D.h"
#include "RingBuffer.h" // dumb class I put in eventReaderRoot ages ago, looks like I can actually use it again!
class TF1;

namespace Acclaim{
  class FourierBuffer;

  class RayleighHist : public TH1D {

    friend class FourierBuffer;

  public:
    RayleighHist(FourierBuffer* fb=NULL, const char* name = "", const char* title = "");
    virtual ~RayleighHist();

    virtual void Draw(Option_t* opt="");

    void Eval(Double_t& chiSquare, Int_t& ndf);
    int add(double newAmp);

    void SetFreqBinToDraw(Int_t freqBin); // *MENU*
    
    RingBuffer amplitudes;

    static void guessMaxBinLimitAndSigmaFromMean(double meanAmp, double& maxAmp, double& sigmaGuess, double fracOfEventsInsideMaxAmp);

        
  protected:

    virtual int Fill(double amp, double sign=1);
    bool axisRangeOK() const;
    void rebinAndRefill(double meanAmp);
    
    FourierBuffer* fParent;    
    TF1* fRay;    
    double fracOfEventsWanted;
    double maxOverFlowThresh;    
    Int_t risingEdgeBins;    
    Int_t fDrawFreqBin;
    double freqMHz;
    Int_t numEvents;

    ClassDef(RayleighHist, 0);
    
  };
}
#endif
