#ifndef RAYLEIGH_HIST_H
#define RAYLEIGH_HIST_H

#include "TH1D.h"


class TF1;

namespace Acclaim{
  class FourierBuffer;

  class RayleighHist : public TH1D {

    friend class FourierBuffer;

  public:
    RayleighHist(FourierBuffer* fb=NULL, const char* name = "", const char* title = "");
    virtual ~RayleighHist();

    virtual int Fill(double amp, double sign=1);
    virtual void Draw(Option_t* opt="");

    bool axisRangeOK(double meanAmp) const;
    void rebinAndEmptyHist(double meanAmp);
    void Fit();

    void SetFreqBinToDraw(Int_t freqBin); // *MENU*
    Int_t numEvents;
    
  protected:

    FourierBuffer* fParent;
    TF1* fRay;    
    double fracOfEventsWanted;
    double maxOverFlowThresh;    
    Int_t risingEdgeBins;    
    Int_t fDrawFreqBin;
    double freqMHz;
    ClassDef(RayleighHist, 0);
    
  };
}
#endif
