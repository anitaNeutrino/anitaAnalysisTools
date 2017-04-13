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
#include <vector>
#include "RawAnitaHeader.h"

class TSpectrum;

namespace Acclaim
{

  class RayleighHist;
  
  /**
   * @class FourierBuffer
   * @brief A class to hold frequency domain info in memory
   * I've removed the mag/phase info for now as it's probably a bit gratuitous
   */
  class FourierBuffer {

    friend class RayleighHist;
  public:

    virtual ~FourierBuffer();
    explicit FourierBuffer(Double_t timeScaleSeconds=10, Int_t theAnt=-1, AnitaPol::AnitaPol_t thePol = AnitaPol::kNotAPol);

    size_t add(const RawAnitaHeader* header, const AnalysisWaveform* wave);
    
    const RayleighHist* getRayleighDistribution(Int_t freqBin=-1) const {return hRays.at(freqBin >= 0 ? freqBin : fDrawFreqBin);}
    TGraphAligned* getAvePowSpec_dB(double timeRange = -1) const;
    TGraphAligned* getAvePowSpec(double timeRange = -1) const;
    TGraphAligned* getBackground_dB(double timeRange = -1) const;
    TGraphAligned* getBackground(double timeRange = -1) const;

    void setAntPol(Int_t theAnt, AnitaPol::AnitaPol_t thePol){
      ant = theAnt;
      pol = thePol;      
    }
    
  private:
    Int_t removeOld();
    void initVectors(int n);
    Int_t ant;
    AnitaPol::AnitaPol_t pol;

    std::list<std::vector<double> > powerRingBuffer;
    std::list<UInt_t> eventNumbers;
    std::list<Int_t> runs;
    std::list<Double_t> realTimesNs;

    std::vector<double> sumPower;
    std::vector<RayleighHist*> hRays;
    std::vector<double> sumAmps;

    Double_t timeScale;

    double df;
    mutable TSpectrum* spectrum; // to estimate the background
    bool doneVectorInit;
    int fDrawFreqBin;
    TF1* fRay;
  };
}



#endif //FOURIERBUFFER_H
