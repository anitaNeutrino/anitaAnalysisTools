#include "FancyFFTs.h"
#include <fftw3.h>

std::map<std::pair<int, int>, fftw_plan> fRealToComplex;
std::map<std::pair<int, int>, fftw_plan> fComplexToReal;
std::map<std::pair<int, int>, double*> fReals;
std::map<std::pair<int, int>, complex<double>*> fComplex;





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates new fftw plan to handle 1D transforms of length len if they haven't been allocated already.
 *
 * @param len is the length of the input.
 * @param threadInd is the index of the thread that needs an fftw plan.
*/
bool FancyFFTs::makeNewPlanIfNeeded(int len, int threadInd){

  std::pair<int, int> key(len, threadInd);
  std::map<std::pair<int, int>,fftw_plan>::iterator it = fRealToComplex.find(key);
  if(it==fRealToComplex.end()){
    // std::cout << len << "\t" << threadInd << std::endl;
    fReals[key] = (double*) fftw_malloc(sizeof(double)*len);
    fComplex[key] = (complex<double>*) fftw_malloc(sizeof(fftw_complex)*len);
    fRealToComplex[key] = fftw_plan_dft_r2c_1d(len,fReals[key],(fftw_complex*)fComplex[key],FFTW_MEASURE);
    fComplexToReal[key] = fftw_plan_dft_c2r_1d(len,(fftw_complex*)fComplex[key],fReals[key],FFTW_MEASURE);
    return true;
  }
  else{
    return false;
  }
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates a power spectrum TGraph from an array of doubles.
 *
 * @param len is the length of the input.
 * @param input is a pointer to an array of doubles of length len.
 * @param dt is the time step between samples.
 * @param normFlag is a flag specifying how to normalize the output TGraph.
 * @param dBScale specifies whether to return the TGraph in dB or linear scale.
 * @param threadInd specifies the thread index of the fftw plans to use.
 * @returns the power spectrum TGraph.
*/
TGraph* FancyFFTs::getPowerSpectrumTGraph(int len, double* input, double dt, FancyFFTs::conventionFlag normFlag, bool dBScale, int threadInd){
  double* powSpec = getPowerSpectrum(len, input, dt, normFlag, threadInd);
  int numFreqs = getNumFreqs(len);
  if(dBScale==true){
    for(int freqInd=0; freqInd < numFreqs; freqInd++){
      powSpec[freqInd] = 10*TMath::Log10(powSpec[freqInd]);
    }
  }
  return new TGraph(numFreqs, getFreqArray(len, dt), powSpec);
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates the power spectrum in an array of doubles of length (len/2 + 1).
 *
 * @param len is the length of the input.
 * @param input is a pointer to an array of doubles of length len.
 * @param dt is the time step between samples.
 * @param normFlag is a flag specifying how to normalize the output TGraph.
 * @param threadInd uses a particular threads plans to do the ffts.
 * @returns a pointer to the array containing the power spectrum.
*/
double* FancyFFTs::getPowerSpectrum(int len, double* input, double dt, FancyFFTs::conventionFlag normFlag, int threadInd){
  
  return FancyFFTs::getPowerSpectrum(len, input, dt, normFlag, NULL, threadInd);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates the power spectrum in an array of doubles of length (len/2 + 1).
 *
 * @param len is the length of the input.
 * @param input is a pointer to an array of doubles of length len.
 * @param dt is the time step between samples.
 * @param normFlag is a flag specifying how to normalize the output TGraph.
 * @param threadInd uses a particular threads plans to do the ffts.
 * @param outputPtr can be used to direct FancyFFTs to copy the output into pre-allocatred memory.
 * @returns a pointer to the array containing the power spectrum.
*/
double* FancyFFTs::getPowerSpectrum(int len, double* input, double dt, FancyFFTs::conventionFlag normFlag, double* outputPtr, int threadInd){

  /* 
     FancyFFTs::conventionFlag determines (you guessed it) the normalization of the power spectrum.
     For a thorough explanation see Ryan's presentation:
     http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf
  */
  
  const double ohms = 50; // Assuming 50 Ohm termination
  double conventionNorm = 1;
  switch (normFlag){
  case kSum:
    /* Sum of output powSpec == sum over V[i]*V[i] for each V[i] in input */
    conventionNorm = 1./len;
    break;

  case kAverage:
    /* 
       Sum of output powSpec == sum over V[i]*V[i]/len for each V[i] in input
       Tells you power per unit sample? ANITAns probably don't want this one.
    */
    conventionNorm = 1./(len*len);
    break;

  case kTimeIntegral:
    /*
      Sum of output powSpec == sum over dt*V[i]*V[i] for each V[i] in input 
      Tell you total time integrated power. 
    */
    conventionNorm = dt/len;
    break;


  case kPowSpecDensity:
    /* 
       Sum of df*psd[i] for each psd[i] in output == sum over dt*V[i]*V[i]/df for each V[i] in input
                                                  == sum over N*dt*dt*V[i]*V[i] for each V[i] in input
       Makes this function return the "power spectral density".
       i.e. frequency bins are normalized such that zero padding the waveform won't affect bin content.
    */

    conventionNorm = dt*dt; /* since df = 1./(len*dt) */
    break;


  case kPowSpecDensity_dBm:
    /* 
       Sum of df*psd[i] for each psd[i] in output == sum over dt*V[i]*V[i]/df/ohms for each V[i] in input
                                                  == sum over N*dt*dt*V[i]*V[i]/ohms for each V[i] in input
       Makes this function return the "power spectral density" (assuming 50 Ohm termination).
       i.e. frequency bins are normalized such that zero padding the waveform won't affect bin content.
    */

    conventionNorm = dt*dt/ohms; /* since df = 1./(len*dt) */
    break;
    
  default:
    /* You shouldn't get here and now I'm going to tell you that. */
    std::cerr << "Invalid FancyFFTs::conventionFlag in FancyFFTs::getPowerSpectrum(int len, double* input, double dt, FancyFFTs::conventionFlag normFlag) " << std::endl;
  }
  
  
  /* 
     Do FFT without putting the output in a new array.
     we need to normalize the output so lets do that when we move it.
  */
  doFFT(len, input, false, threadInd);

  
  std::pair<int, int> key(len, threadInd);
  
  /* Get the fftw_malloc'd array that the plan uses */
  complex<double>* rawFftOutputPtr = (complex<double>*) fComplex[key];

  const int powSpecLen = getNumFreqs(len);
  
  double* powSpec = NULL;
  if(outputPtr==NULL){
    powSpec = new double[powSpecLen];
  }
  else{
    powSpec = outputPtr;
  }

  /* Pedantry with the normalization (only half the nyquist bin if there's an even number of bins) */
  bool halfPowerNquist = (len%2)==0 ? true : false;

  for(int i=0; i<powSpecLen; i++){

    /* Take account of symetry in normalization of power spectrum */
    double symmetryFactor = 2;
    if(i==0 || (halfPowerNquist && i==powSpecLen-1)){
      symmetryFactor = 1;
    }

    powSpec[i] = symmetryFactor*conventionNorm*(pow(std::abs(rawFftOutputPtr[i]), 2));
  }
  return powSpec;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Does a forward fast fourier transform on a real input array.
 *
 * @param len is the length of the input.
 * @param input is a pointer to an array of doubles of length len.
 * @param copyOutputToNewArray leave output in internal memory (false) or allocate new memory and copy (true).
 * @param threadInd uses a particular threads plans to do the ffts.
 * @returns a pointer to an array of complex<double>s containing the fft.
 */    
complex<double>* FancyFFTs::doFFT(int len, double* input, bool copyOutputToNewArray, int threadInd){
  return doFFT(len, input, NULL, copyOutputToNewArray, threadInd);
}



//---------------------------------------------------------------------------------------------------------
/**
 * @brief Does a forward fast fourier transform on a real input array.
 *
 * @param len is the length of the input.
 * @param input is a pointer to an array of doubles of length len.
 * @param output is a pointer to an array of complex<double> to copy the output to.
 * @param threadInd uses a particular threads plans to do the ffts.
 * @returns a pointer to an array of complex<double>s containing the fft.
*/
complex<double>* FancyFFTs::doFFT(int len, double* input, complex<double>* output, int threadInd){
  return doFFT(len, input, output, true, threadInd);
}



//---------------------------------------------------------------------------------------------------------
/**
 * @brief Does a forward fast fourier transform on a real input array.
 *
 * @param len is the length of the input.
 * @param input is a pointer to an array of doubles of length len.
 * @param output is a pointer to an array of complex<double> to copy the output to.
 * @param copyOutputToNewArray leave output in internal memory (false) or allocate new memory and copy (true).
 * @param threadInd uses a particular threads plans to do the ffts.
 * @returns a pointer to an array of complex<double>s containing the fft.
*/
complex<double>* FancyFFTs::doFFT(int len, double* input, complex<double>* output, bool copyOutputToNewArray, int threadInd){
  /* 
     Using complex<double> instead of the typdef fftw_complex double[2] 
     because CINT has a better time with it, even though it's (apparently) exactly the same.
  */

  std::pair<int, int> key(len, threadInd);
  makeNewPlanIfNeeded(len, threadInd);

  memcpy(fReals[key], input, sizeof(double)*len);

  fftw_execute(fRealToComplex[key]);

  if(copyOutputToNewArray==true){
    int numFreqs = getNumFreqs(len);
    complex<double>* theOutput = output;
    if(theOutput==NULL){
      theOutput = new complex<double>[numFreqs];
    }
    
    /* Seems to work, see http://www.fftw.org/doc/Complex-numbers.html */
    memcpy(theOutput, fComplex[key], sizeof(fftw_complex)*numFreqs);
    return theOutput;
  }
  else{
    return NULL;
  }
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Does an inverse fast fourier transform on an array of complex<double>s.
 *
 * @param len is the length of the real output array.
 * @param input is a pointer to an array of complex<double>s of length (len/2 +1).
 * @param copyOutputToNewArray leave output in internal memory (false) or allocate new memory and copy (true).
 * @param threadInd uses a particular threads plans to do the ffts.
 * @returns a pointer to an array of doubles containing the (real) inverse FFT.
*/
double* FancyFFTs::doInvFFT(int len, complex<double>* input, bool copyOutputToNewArray, int threadInd){
  return doInvFFT(len, input, NULL, copyOutputToNewArray, threadInd);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Does an inverse fast fourier transform on an array of complex<double>s.
 *
 * @param len is the length of the real output array.
 * @param input is a pointer to an array of complex<double>s of length (len/2 +1).
 * @param output is a pointer to an array of doubles to copy the output to.
 * @param threadInd uses a particular threads plans to do the ffts.
 * @returns a pointer to an array of doubles containing the (real) inverse FFT.
*/
double* FancyFFTs::doInvFFT(int len, complex<double>* input, double* output, int threadInd){
  return doInvFFT(len, input, output, true, threadInd);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Does an inverse fast fourier transform on an array of complex<double>s.
 *
 * @param len is the length of the real output array.
 * @param input is a pointer to an array of complex<double>s of length (len/2 +1).
 * @param output is a pointer to an array of doubles to copy the output to.
 * @param copyOutputToNewArray leave output in internal memory (false) or allocate new memory and copy (true).
 * @param threadInd uses a particular threads plans to do the ffts.
 * @returns a pointer to an array of doubles containing the (real) inverse FFT.
*/
double* FancyFFTs::doInvFFT(int len, complex<double>* input, double* output, bool copyOutputToNewArray, int threadInd){
  
  /* 
     Normalization of 1/N done in this function. 
     Note: fftw_plan_c2r_1d USES AND MESSES UP THE INPUT ARRAY when executed.
  */

  std::pair<int, int> key(len, threadInd);

  makeNewPlanIfNeeded(len, threadInd);
  int numFreqs = getNumFreqs(len);

  complex<double>* tempVals = (complex<double>*) fComplex[key];
  
  // In the case that we cleverly left the fft in the internal array skip the copy
  if(tempVals!=input){ 
    // In fact this is undefined behaviour!
    memcpy(tempVals, input, sizeof(fftw_complex)*numFreqs);
  }

  // Do inverse FFT.
  fftw_execute(fComplexToReal[key]);
  
  /* Normalization needed on the inverse transform */
  double* invFftOutPtr = fReals[key];
  for(int i=0; i<len; i++){
    invFftOutPtr[i]/=len;
  }

  double* theOutput = output;
  if(copyOutputToNewArray==true){
    if(theOutput==NULL){
      theOutput = new double[len];
    }  
    memcpy(theOutput, invFftOutPtr, sizeof(double)*len);
    return theOutput;    
  }
  else{
    return NULL;
  }
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Find the next highest power of two.
 *
 * @param len is the number to find the next highest power of two.
 * @returns the next highest power of two.
*/
int FancyFFTs::extendToPowerOfTwo(int len){
  return pow(2, TMath::CeilNint(TMath::Log2(len)));
}






//---------------------------------------------------------------------------------------------------------
/**
 * @brief Get the array of frequencies to go with FFT output, determined by the tranform length and dt.
 *
 * @param len is the length of the real fft input.
 * @param dt is the time between real input samples.
 * @returns a pointer to an array of doubles containing the frequencies.
*/
double* FancyFFTs::getFreqArray(int len, double dt){
  /* Once and only once... */  
  int numFreq = getNumFreqs(len);
  double df = 1./(dt*len);
  double* freqArray = new double[numFreq];
  for(int i=0; i<numFreq; i++){
    freqArray[i] = df * i;
  }
  return freqArray;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Get the number of frequencies from the number of samples.
 *
 * @param len is the length of the real fft input.
 * @returns the number of frequencies.
*/
int FancyFFTs::getNumFreqs(int len){
  return (len/2 + 1);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Zero padds FFTts, 
 *
 * @param fft is the unpadded fft.
 * @param numSamples is the length of the time domain of the fft (i.e. NOT the number of frequency bins)
 * @param numSamplesUpsampled is the length of the time domain of the padded fft (i.e. NOT the number of frequency bins)
 * @returns the a pointer to an array of complex<doubles> containing the padded fft.
 *
 * This is for doing interpolation in the time domain.
 * Zero padding in the time/frequency domain is equivelent to interpolation in the time/frequency domain.
*/
complex<double>* FancyFFTs::zeroPadFFT(complex<double>* fft, int numSamples, int numSamplesUpsampled){
  return zeroPadFFT(fft, NULL, numSamples, numSamplesUpsampled);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Zero padds FFTts, 
 *
 * @param fft is the unpadded fft.
 * @param numSamples is the length of the time domain of the fft (i.e. NOT the number of frequency bins)
 * @param numSamplesUpsampled is the length of the time domain of the padded fft (i.e. NOT the number of frequency bins)
 * @param output is where the padded fft should be copied to.
 * @returns the a pointer to an array of complex<doubles> containing the padded fft.
 *
 * This is for doing interpolation in the time domain.
 * Zero padding in the time/frequency domain is equivelent to interpolation in the time/frequency domain.
*/
complex<double>* FancyFFTs::zeroPadFFT(complex<double>* fft, complex<double>* output, int numSamples, int numSamplesUpsampled){

  const int numFreqs = getNumFreqs(numSamples);  
  const int numFreqsPadded = getNumFreqs(numSamplesUpsampled);
  
  complex<double>* fftPadded = NULL;
  if(output==NULL){
    fftPadded = new complex<double>[numFreqsPadded];
  }
  else{
    fftPadded = output;
  }

  // this is the mistake!!!!
  // Double_t scale = numFreqsPadded/numFreqs;

  // Here I scale the padded FFT so that it is as if I fourier transformed a longer waveform.
  // (There is a scale factor of length picked up from a forward FFT.)
  Double_t scale = numSamplesUpsampled/numSamples;

  for(int freqInd=0; freqInd<numFreqs; freqInd++){
    fftPadded[freqInd].real(fft[freqInd].real()*scale);
    fftPadded[freqInd].imag(fft[freqInd].imag()*scale);
  }
  for(int freqInd=numFreqs; freqInd<numFreqsPadded; freqInd++){
    fftPadded[freqInd] = 0;
  }

  // this factor undoes the effect of the the nyquist frequency containing half the power relative
  // to the other bins (due to them containing the negative frequency components)
  const Double_t sqrtHalf = TMath::Sqrt(0.5);
  fftPadded[numFreqs-1].real(sqrtHalf*fftPadded[numFreqs-1].real());
  fftPadded[numFreqs-1].imag(sqrtHalf*fftPadded[numFreqs-1].imag());  
  
  return fftPadded;
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Prints the list of keys of plans to stdout.
 */
int FancyFFTs::printListOfKeys(){

  std::vector<std::pair<int, int> > keys;

  std::map<std::pair<int, int>,fftw_plan>::iterator it;
  for(it = fRealToComplex.begin(); it != fRealToComplex.end(); it++){
    keys.push_back(it->first);
  }
  
  /* Pretty sure in advance of testing that this list is not guarenteed to be sorted. */
  // std::vector<std::pair<int, int>> sortedIndices(keys.size());
  // TMath::Sort(int(keys.size()), &keys[0], &sortedIndices[0], kFALSE);

  /* Print to terminal */
  std::cout << "Plan lengths in memory = [";
  for(int i=0; i<int(keys.size()); i++){
    // int j = sortedIndices.at(i);
    // std::cout << keys.at(j);
    std::cout << "(" << keys.at(i).first << ", " << keys.at(i).second << ")";
    if(i < int(keys.size())-1 ) std::cout << ", ";
  }
  std::cout << "]" << std::endl;

  return int(keys.size());
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Cross correlates the two input arrays.
 * 
 * @param len is the length of both arrays, v1 and 2.
 * @param v1 is the first input array.
 * @param v2 is the second input array.
 * @param threadInd uses a particular threads plans to do the ffts.
 * @returns a pointer to an array of doubles containing the cross correlations.
 */
double* FancyFFTs::crossCorrelate(int len, double* v1, double* v2, int threadInd){
  /* 
     Cross correlation is the same as bin-by-bin multiplication in the frequency domain (with a conjugation).
     Will assume lengths are the same for now.
  */

  /* Store output of FFT2 in tempVals */
  complex<double>* tempVals2 = doFFT(len, v2, true, threadInd);

  /* Leave output of FFT1 in internal arrays to avoid an unnecessary copy */
  doFFT(len, v1, false, threadInd);

  std::pair<int, int> key(len, threadInd);
  
  /* Get pointer to internal array */
  complex<double>* tempVals1 = (complex<double>*) fComplex[key];

  /* Take the product */
  int numFreqs = getNumFreqs(len);
  for(int i=0; i<numFreqs; i++){
    tempVals1[i] *= std::conj(tempVals2[i]);
  }

  delete [] tempVals2;
  
  /* Product back to time domain */
  double* crossCorr = doInvFFT(len, tempVals1, true, threadInd);

  /* 
     Picked up two factors of len when doing forward FFT, only removed one doing invFFT.
     This takes out the second factor.
  */
  for(int i=0; i<len; i++){
    crossCorr[i] /= len;
  }

  return crossCorr;
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Cross correlates the two input arrays.
 * 
 * @param len is the length of both arrays, v1 and 2.
 * @param fft1 is the fft of the first input array.
 * @param fft2 is the fft of the second input array.
 * @param threadInd uses a particular threads plans to do the ffts.
 * @returns a pointer to an array of doubles containing the cross correlations.
 */
double* FancyFFTs::crossCorrelate(int len, complex<double>* fft1, complex<double>* fft2, int threadInd){
  return crossCorrelate(len, fft1, fft2, NULL, threadInd);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Cross correlates the two input arrays.
 * 
 * @param len is the length of the time domain input (not the length of the ffts).
 * @param fft1 is the fft of the first input array.
 * @param fft2 is the fft of the second input array.
 * @param output is a pointer to copy the outputted cross correlations to.
 * @param threadInd uses a particular threads plans to do the ffts.
 * @returns a pointer to an array of doubles containing the cross correlations.
 */
double* FancyFFTs::crossCorrelate(int len, complex<double>* fft1, complex<double>* fft2,
				  double* output, int threadInd){
  /* 
     Cross correlation is the same as bin-by-bin multiplication (and some conjugation) in the frequency domain.
     Will assume lengths are the same for now.
  */


  /* Stops tempVals returning NULL, normally done in doFFT step. 
     But the nice thing about this class is it means that I won't be duplicating work. */
  makeNewPlanIfNeeded(len, threadInd);

  std::pair<int, int> key(len, threadInd);
  
  /* Grab array associated with plan from internal memory */
  complex<double>* tempVals = (complex<double>*) fComplex[key];

  // TThread::Lock();
  // std::cout << "threadInd = " << threadInd << "\tfComplex[key] = " << tempVals << std::endl << std::endl;
  // TThread::UnLock();
  
  
  /* Take the product */
  int numFreqs = getNumFreqs(len);

  for(int i=0; i<numFreqs; i++){
    tempVals[i] = fft1[i]*std::conj(fft2[i]);
  }

  /* Product back to time domain */  
  double* crossCorr = output;
  if(crossCorr==NULL){
    /* Allocates new memory */
    crossCorr = doInvFFT(len, tempVals, true, threadInd);
  }
  else{
    /* Does not allocate new memory */
    crossCorr = doInvFFT(len, tempVals, crossCorr, true, threadInd);    
  }


  /* 
     Picked up two factors of len when doing forward FFT, only removed one doing invFFT.
     This takes out the second factor.
  */
  for(int i=0; i<len; i++){
    // std::cout << crossCorr[i] << std::endl;
    crossCorr[i] /= len;
  }

  return crossCorr;
}


//---------------------------------------------------------------------------------------------------------
/**
 * @brief Get the real array of the internal fftw memory. Do not delete this!
 * 
 * @param key is an std::pair of the real array length and the thread index of the plan.
 * @returns a pointer to fftw's internal real array. DO NOT DELETE THIS!
 */
double* FancyFFTs::getRealArray(std::pair<Int_t, Int_t> key){
  return fReals[key];
}
