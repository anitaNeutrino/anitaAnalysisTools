#ifndef RAYLEIGH_HIST_H
#define RAYLEIGH_HIST_H

#include "TH1D.h"


namespace Acclaim{
  class RayleighHist : public TH1D {

    friend class FourierBuffer;
  public:
    RayleighHist(const char* name = "", const char* title = "");
    virtual ~RayleighHist(){;}

    virtual int Fill(double amp, double sign=1);

    bool axisRangeOK(double meanAmp) const;
    void rebinAndEmptyHist(double meanAmp);
    
  protected:

    void init();
    double fracOfEventsWanted;
    double maxOverFlowThresh;    
    Int_t risingEdgeBins;    

    
    ClassDef(RayleighHist, 1);
  };
}
#endif
