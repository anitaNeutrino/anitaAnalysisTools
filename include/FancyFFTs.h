/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             C++ ROOT friendly namespace to do FFTs faster than Ryan. 
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



using std::complex;

/** @namespace FancyFFTs
 * @brief My implementation of a wrapper for FFTW for use with ROOT thing.
 *
 * Most functionality has been merged into FFTtools, except some of the threading options.
 * Hopefully that will get done at some point in the future and this namespace can be done away with.
 */
namespace FancyFFTs {

  /**
   * @brief Flag to pass to FancyFFT functions specifying normalization of power spectra.
   */  
  enum conventionFlag {
    kSum = 0,
    kAverage = 1,
    kTimeIntegral = 2,
    kPowSpecDensity = 3,
    kPowSpecDensity_dBm = 4
  };

  complex<double>* doFFT(int len, double* input, bool copyOutputToNewArray, int threadInd=0);
  complex<double>* doFFT(int len, double* input, complex<double>* output, int threadInd=0);
  complex<double>* doFFT(int len, double* input, complex<double>* output,
				     bool copyOutputToNewArray, int threadInd=0);  
  double* doInvFFT(int len, complex<double>* input, bool copyOutputToNewArray, int threadInd=0);
  double* doInvFFT(int len, complex<double>* input, double* output, int threadInd=0);  
  double* doInvFFT(int len, complex<double>* input, double* output,
			  bool copyOutputToNewArray, int threadInd=0);

  double* getPowerSpectrum(int len, double* input, double dt, 
				  conventionFlag normFlag,
				  int threadInd=0);

  double* getPowerSpectrum(int len, double* input, double dt, 
				  conventionFlag normFlag,
				  double* outputPtr, int threadInd=0);

  TGraph* getPowerSpectrumTGraph(int len, double* input, double dt, 
					conventionFlag normFlag,
					bool dBScale, int threadInd=0);
  double* getFreqArray(int len, double dt);
  int getNumFreqs(int len);
  int printListOfKeys();
  double* crossCorrelate(int len, double* v1, double* v2, int threadInd=0);
  double* crossCorrelate(int len, complex<double>* fft1, complex<double>* fft2, int threadInd=0);
  double* crossCorrelate(int len, complex<double>* fft1, complex<double>* fft2,
				double* output, int threadInd=0);
  int extendToPowerOfTwo(int len);
  complex<double>* zeroPadFFT(complex<double>* fft, int numSamples, int numSamplesUpsampled);
  complex<double>* zeroPadFFT(complex<double>* fft, complex<double>* output,
			      int numSamples, int numSamplesUpsampled);  


 /* Takes care of checking whether a plan exists or not */
  bool makeNewPlanIfNeeded(int len,int threadInd=0);

  double* getRealArray(std::pair<Int_t, Int_t> key);
  

};


#endif //FANCYFFTS_H
