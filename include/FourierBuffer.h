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

namespace Acclaim
{

  /**
   * @class FourierBuffer
   * @brief A class to hold frequency domain info in memory
   * I've removed the mag/phase info for now as it's probably a bit gratuitous
   */
  class FourierBuffer {

  public:

    explicit FourierBuffer(Double_t timeScaleSeconds=60, Int_t theAnt=-1, AnitaPol::AnitaPol_t thePol = AnitaPol::kNotAPol);

    size_t add(const RawAnitaHeader* header, const AnalysisWaveform& wave);

    TH1D* getRayleighDistribution(Int_t freqBin) const;
    TH1D* fillRayleighInfo(Int_t freqBin, RayleighInfo* info) const;
    TH1D* fillRiceInfo(Int_t freqBin, RiceInfo* info) const;
    TGraphAligned* getAvePowSpec_dB(double timeRange = -1) const;
    TGraphAligned* getAvePowSpec(double timeRange = -1) const;

    void setAntPol(Int_t theAnt, AnitaPol::AnitaPol_t thePol){
      ant = theAnt;
      pol = thePol;
    }

  private:

    Int_t removeOld();

    Int_t ant;
    AnitaPol::AnitaPol_t pol;
    
    // std::list<std::vector<FFTWComplex> > freqVecs;
    std::list<std::vector<double> > powerRingBuffer;    
    std::list<UInt_t> eventNumbers;
    std::list<Int_t> runs;  
    std::list<Double_t> realTimesNs;

    std::vector<double> sumPower;
    Double_t timeScale;


    double df;

  };

}



#endif //FOURIERBUFFER_H
