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
#include "TObject.h"
#include "TSystem.h"
#include "TMath.h"
#include "TGraph.h"
#include "TThread.h"

/* 
   Will use std::complex<double> for i/o as should be bit-to-bit identical to typdef fftw_complex double[2].
   So long as <complex> is included in front of fftw3.h.
*/
#include <complex>
#include <fftw3.h>


#ifdef __CINT__
#ifdef FFTW_64_BIT // Hack for Hawaii install of FFTW. Is there ever a reason to now have this defed? 
typedef struct {char a[16];} __float128; /* 16 chars have the same size as one __float128 */
#endif
#endif 


namespace PowSpecNorm {
  enum conventionFlag {
    kSum = 0,
    kAverage = 1,
    kTimeIntegral = 2,
    kPowSpecDensity = 3,
    kPowSpecDensity_dBm = 4
  };
}

class FancyFFTs{

  // This class needs to a friend to CrossCorrelator so it can assign
  // multiple plans of the same length, one for each thread
  // and I've chosen to make the assign plan thing private.
  friend class CrossCorrelator;
  
public:
  FancyFFTs(); /* Not used */
  ~FancyFFTs(); /* Not used */


  /* Real-to-complex and complex-to real functions */
  static std::complex<double>* doFFT(Int_t len, double* input, bool copyOutputToNewArray, int threadInd=0);
  static std::complex<double>* doFFT(Int_t len, double* input, std::complex<double>* output, int threadInd=0);
  static std::complex<double>* doFFT(Int_t len, double* input, std::complex<double>* output,
				     bool copyOutputToNewArray, int threadInd=0);  
  static double* doInvFFT(int len, std::complex<double>* input, bool copyOutputToNewArray, int threadInd=0);
  static double* doInvFFT(int len, std::complex<double>* input, double* output, int threadInd=0);  
  static double* doInvFFT(int len, std::complex<double>* input, double* output,
			  bool copyOutputToNewArray, int threadInd=0);

  static double* getPowerSpectrum(int len, double* input, double dt, 
				  PowSpecNorm::conventionFlag normFlag,
				  int threadInd=0);

  static double* getPowerSpectrum(int len, double* input, double dt, 
				  PowSpecNorm::conventionFlag normFlag,
				  double* outputPtr, int threadInd=0);

  static TGraph* getPowerSpectrumTGraph(int len, double* input, double dt, 
					PowSpecNorm::conventionFlag normFlag,
					bool dBScale, int threadInd=0);
  static double* getFreqArray(int len, double dt);
  static int getNumFreqs(int len);
  static int printListOfKeys();
  static double* crossCorrelate(int len, double* v1, double* v2, int threadInd=0);
  static double* crossCorrelate(int len, std::complex<double>* fft1, std::complex<double>* fft2, int threadInd=0);
  static double* crossCorrelate(int len, std::complex<double>* fft1, std::complex<double>* fft2,
				double* output, int threadInd=0);
  static int extendToPowerOfTwo(int len);
  static std::complex<double>* zeroPadFFT(std::complex<double>* fft, int numFreqs, int numFreqsPadded);
  static std::complex<double>* zeroPadFFT(std::complex<double>* fft, std::complex<double>* output,
					  int numFreqs, int numFreqsPadded);  


private:

 /* Takes care of checking whether a plan exists or not */
  static bool makeNewPlanIfNeeded(int len,int threadInd=0);

    /*
     std::maps which hold all the fftw goodies.
     The length is the key so we have an easy way to check if a plan already exists.
     The values are the input/ouput arrays and the plans, real-to-complex and complex-to-real.
  */

  static std::map<std::pair<int, int>, fftw_plan> fRealToComplex;
  static std::map<std::pair<int, int>, fftw_plan> fComplexToReal;
  static std::map<std::pair<int, int>, double*> fReals;
  // static std::map<std::pair<int, int>, fftw_complex*> fComplex;
  static std::map<std::pair<int, int>, std::complex<double>*> fComplex;
};


#endif //FANCYFFTS_H
