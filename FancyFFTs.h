/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             C++ ROOT friendly class to do FFTs faster than Ryan. 
	     Will probably be pretty bare bones intially.
	     I only really want this for doing Cross Correlations.
*************************************************************************************************************** */

#ifndef FANCYFFTS_H
#define FANCYFFTS_H

#include <iostream>
#include <map>
#include <algorithm>
#include <TObject.h>
#include <TSystem.h>
#include <TMath.h>


/* 
   Will use std::complex<double> for i/o as should be bit-to-bit identical to typdef fftw_complex double[2].
   So long as <complex> is included in front of fftw3.h.
*/
#include <complex>
#include <fftw3.h>

namespace PowSpecNorm {
  enum conventionFlag {
    kSum = 0,
    kAverage = 1,
    kTimeIntegral = 2,
    kPowSpecDensity = 3
  };
}

class FancyFFTs : public TObject{

public:
  FancyFFTs(); /* Not used */
  ~FancyFFTs(); /* Not used */


  /* Real-to-complex and complex-to real functions */
  static std::complex<double>* doFFT(Int_t len, double* input, bool copyOutputToNewArray);
  static double* doInvFFT(int len, std::complex<double>* input, bool copyOutputToNewArray);

  static double* getPowerSpectrum(int len, double* input, double dt, PowSpecNorm::conventionFlag normFlag);
  static double* getFreqArray(int len, double dt);
  static int getNumFreqs(int len);
  static int printListOfKeys();
  static double* crossCorrelate(int len, double* v1, double* v2);

private:

  static int extendToPowerOfTwo(int len);
  static bool makeNewPlanIfNeeded(int len); /* Takes care of checking whether a plan exists or not */
  /*
     std::maps which hold all the fftw goodies.
     The length is the key so we have an easy way to check if a plan already exists.
     The values are the input/ouput arrays and the plans, real-to-complex and complex-to-real.
  */
  static std::map<int, fftw_plan> fRealToComplex;
  static std::map<int, fftw_plan> fComplexToReal;
  static std::map<int, double*> fReals;
  // static std::map<int, fftw_complex*> fComplex;
  static std::map<int, std::complex<double>*> fComplex;

  ClassDef(FancyFFTs, 0);
};


#endif //FANCYFFTS_H
