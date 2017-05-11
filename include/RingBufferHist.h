#ifndef RING_BUFFER_HIST_H
#define RING_BUFFER_HIST_H

#include "TH1D.h"
#include "RingBuffer.h"
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include "TF1.h"

class TGraph;

namespace Acclaim{

  class RingBufferHist : public TH1D {

    friend class FourierBuffer;

  public:

    RingBufferHist(const char* name = "", const char* title = "", int nBins=1, double xMin=0, double xMax=1, int ringBufferSize=1000);
    virtual ~RingBufferHist();

    virtual bool add(double newAmp); //!< Input amplitudes events
        
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

    ClassDef(RingBufferHist, 0);
  };


  
}
#endif
