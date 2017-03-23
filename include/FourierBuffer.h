/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             A FIFO queue for frequency domain information
*************************************************************************************************************** */

#ifndef FOURIERBUFFER_H
#define FOURIERBUFFER_H

#include "RTypes.h"
#include "AnalysisWaveform.h"

#include "RayleighInfo.h"
#include "RiceInfo.h"

#include "TH1D.h"
#include <complex>
#include <list>
#include <vector>
#include "RawAnitaHeader.h"

//#define NUM_FREQS ((NUM_SAMPLES/2) + 1)


/**
 * @class FourierBuffer
 * @brief A class to hold fourier domain representations in memory
*/
class FourierBuffer {

public:

  explicit FourierBuffer(Int_t timeScaleSeconds=60);

  size_t add(const RawAnitaHeader* header, const AnalysisWaveform& wave);

  TH1D* getRayleighDistribution(Int_t freqBin);
  TH1D* fillRayleighInfo(Int_t freqBin, RayleighInfo* info);
  TH1D* fillRiceInfo(Int_t freqBin, RiceInfo* info);  

private:


  Int_t removeOld();

  std::list<std::vector<FFTWComplex> > freqVecs;
  std::list<UInt_t> eventNumbers;
  std::list<Int_t> runs;  
  std::list<Double_t> realTimesNs;

  Int_t timeScale;





};





#endif //FOURIERBUFFER_H
