/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             C++ ROOT friendly class to do FFTs faster than Ryan, but is mostly a shameless copy.
	     Will probably be pretty bare bones intially.
	     I only really want this for doing Cross Correlations.
*************************************************************************************************************** */

#include "FancyFFTs.h"

ClassImp(FancyFFTs)

/* Define static members */
/* https://stackoverflow.com/questions/18433752/c-access-private-static-member-from-public-static-method */
std::map<std::pair<int, int>, fftw_plan> FancyFFTs::fRealToComplex;
std::map<std::pair<int, int>, fftw_plan> FancyFFTs::fComplexToReal;
std::map<std::pair<int, int>, double*> FancyFFTs::fReals;
std::map<std::pair<int, int>, std::complex<double>*> FancyFFTs::fComplex;
FancyFFTsWisdomManager FancyFFTs::myWisdom;

FancyFFTs::FancyFFTs(){
  std::cout << "FancyFFTs::FancyFFTs()" << std::endl;
}

FancyFFTs::~FancyFFTs(){
  std::cout << "FancyFFTs::~FancyFFTs()" << std::endl;
}


bool FancyFFTs::makeNewPlanIfNeeded(int len, int threadInd){
  /* 
     Function which checks whether we've encountered a request to do an FFT of this length before.
     If we haven't then we need a new plan!
  */

  std::pair<int, int> key(len, threadInd);
  std::map<std::pair<int, int>,fftw_plan>::iterator it = fRealToComplex.find(key);
  if(it==fRealToComplex.end()){
    // std::cout << len << "\t" << threadInd << std::endl;
    fReals[key] = (double*) fftw_malloc(sizeof(double)*len);
    fComplex[key] = (std::complex<double>*) fftw_malloc(sizeof(fftw_complex)*len);
    fRealToComplex[key] = fftw_plan_dft_r2c_1d(len,fReals[key],(fftw_complex*)fComplex[key],FFTW_MEASURE);
    fComplexToReal[key] = fftw_plan_dft_c2r_1d(len,(fftw_complex*)fComplex[key],fReals[key],FFTW_MEASURE);
    return true;
  }
  else{
    return false;
  }
}

TGraph* FancyFFTs::getPowerSpectrumTGraph(int len, double* input, double dt, PowSpecNorm::conventionFlag normFlag, bool dBScale, int threadInd){
  double* powSpec = getPowerSpectrum(len, input, dt, normFlag, threadInd);
  int numFreqs = getNumFreqs(len);
  if(dBScale==true){
    for(int freqInd=0; freqInd < numFreqs; freqInd++){
      powSpec[freqInd] = 10*TMath::Log10(powSpec[freqInd]);
    }
  }
  return new TGraph(numFreqs, getFreqArray(len, dt), powSpec);
}

double* FancyFFTs::getPowerSpectrum(int len, double* input, double dt, PowSpecNorm::conventionFlag normFlag, int threadInd){

  /* 
     PowSpecNorm::conventionFlag determines (you guessed it) the normalization of the power spectrum.
     For a thorough explanation see Ryan's presentation:
     http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf
  */

  double conventionNorm = 1;
  switch (normFlag){
  case PowSpecNorm::kSum:
    /* Sum of output powSpec == sum over V[i]*V[i] for each V[i] in input */
    conventionNorm = 1./len;
    break;

  case PowSpecNorm::kAverage:
    /* 
       Sum of output powSpec == sum over V[i]*V[i]/len for each V[i] in input
       Tells you power per unit sample? ANITAns probably don't want this one.
    */
    conventionNorm = 1./(len*len);
    break;

  case PowSpecNorm::kTimeIntegral:
    /*
      Sum of output powSpec == sum over dt*V[i]*V[i] for each V[i] in input 
      Tell you total time integrated power. 
    */
    conventionNorm = dt/len;
    break;


  case PowSpecNorm::kPowSpecDensity:
    /* 
       Sum of df*psd[i] for each psd[i] in output == sum over dt*V[i]*V[i]/df for each V[i] in input
                                                  == sum over N*dt*dt*V[i]*V[i] for each V[i] in input
       Makes this function return the "power spectral density".
       i.e. frequency bins are normalized such that zero padding the waveform won't affect bin content.
    */

    conventionNorm = dt*dt; /* since df = 1./(len*dt) */
    break;

  default:
    /* You shouldn't get here and now I'm going to tell you that. */
    std::cerr << "Invalid PowSpecNorm::conventionFlag in FancyFFTs::getPowerSpectrum(int len, double* input, double dt, PowSpecNorm::conventionFlag normFlag) " << std::endl;
  }
  
  
  /* 
     Do FFT without putting the output in a new array.
     we need to normalize the output so lets do that when we move it.
  */
  doFFT(len, input, false, threadInd);

  
  std::pair<int, int> key(len, threadInd);
  
  /* Get the fftw_malloc'd array that the plan uses */
  std::complex<double>* rawFftOutputPtr = (std::complex<double>*) fComplex[key];

  const int powSpecLen = getNumFreqs(len);
  double* powSpec = new double[powSpecLen];

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

std::complex<double>* FancyFFTs::doFFT(int len, double* input, bool copyOutputToNewArray, int threadInd){
  /* 
     Using std::complex<double> instead of the typdef fftw_complex double[2] 
     because CINT has a better time with it, even though it's (apparently) exactly the same.
  */

  std::pair<int, int> key(len, threadInd);
  makeNewPlanIfNeeded(len, threadInd);

  memcpy(fReals[key], input, sizeof(double)*len);

  fftw_execute(fRealToComplex[key]);

  if(copyOutputToNewArray==true){
    int numFreqs = getNumFreqs(len);
    std::complex<double>* output = new std::complex<double>[numFreqs];

    /* Seems to work, see http://www.fftw.org/doc/Complex-numbers.html */
    memcpy(output, fComplex[key], sizeof(fftw_complex)*numFreqs);
    return output;
  }
  else{
    return NULL;
  }
}

double* FancyFFTs::doInvFFT(int len, std::complex<double>* input, bool copyOutputToNewArray, int threadInd){
  
  /* 
     Normalization of 1/N done in this function. 
     Note: fftw_plan_c2r_1d DESTROYS THE INPUT ARRAY when executed.
  */

  std::pair<int, int> key(len, threadInd);

  makeNewPlanIfNeeded(len, threadInd);
  int numFreqs = getNumFreqs(len);

  std::complex<double>* tempVals = (std::complex<double>*) fComplex[key];
  
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

  if(copyOutputToNewArray==true){
    double* output = new double[len];
    memcpy(output, invFftOutPtr, sizeof(double)*len);

    // std::cout << output << "\t" << invFftOutPtr << "\t" << len << "\t" << sizeof(double)*len << std::endl;
    
    return output;
  }
  else{
    return NULL;
  }
}


int FancyFFTs::extendToPowerOfTwo(int len){
  return pow(2, TMath::CeilNint(TMath::Log2(len)));
}


/* Once and only once... */
double* FancyFFTs::getFreqArray(int len, double dt){
  int numFreq = getNumFreqs(len);
  double df = 1./(dt*len);
  double* freqArray = new double[numFreq];
  for(int i=0; i<numFreq; i++){
    freqArray[i] = df * i;
  }
  return freqArray;
}

int FancyFFTs::getNumFreqs(int len){
  return (len/2 + 1);
}

std::complex<double>* FancyFFTs::zeroPadFFT(std::complex<double>* fft, int numFreqs, int numFreqsPadded){

  std::complex<double>* fftPadded = new std::complex<double>[numFreqsPadded];

  // Here I scale the padded FFT so that it is as if I fourier transformed a longer waveform.
  // (There is a scale factor of length picked up from a forward FFT.)
  Double_t scale = numFreqsPadded/numFreqs;
  for(int freqInd=0; freqInd<numFreqs; freqInd++){
    fftPadded[freqInd].real(fft[freqInd].real()*scale);
    fftPadded[freqInd].imag(fft[freqInd].imag()*scale);
  }
  for(int freqInd=numFreqs; freqInd<numFreqsPadded; freqInd++){
    fftPadded[freqInd] = 0;
  }
  
  return fftPadded;
}



int FancyFFTs::printListOfKeys(){
  /* Returns number of plans, prints the lengths of the plans to screen.*/

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


double* FancyFFTs::crossCorrelate(int len, double* v1, double* v2, int threadInd){
  /* 
     Cross correlation is the same as bin-by-bin multiplication in the frequency domain.
     Will assume lengths are the same for now.
  */

  /* Store output of FFT2 in tempVals */
  std::complex<double>* tempVals2 = doFFT(len, v2, true, threadInd);

  /* Leave output of FFT1 in internal arrays to avoid an unnecessary copy */
  doFFT(len, v1, false, threadInd);

  std::pair<int, int> key(len, threadInd);
  
  /* Get pointer to internal array */
  std::complex<double>* tempVals1 = (std::complex<double>*) fComplex[key];

  /* Take the product */
  int numFreqs = getNumFreqs(len);
  for(int i=0; i<numFreqs; i++){
    tempVals1[i] *= std::conj(tempVals2[i]);
  }
  
  /* Product back to time domain */
  double* crossCorr = doInvFFT(len, tempVals1, true, threadInd);
  delete [] tempVals1;

  /* 
     Picked up two factors of len when doing forward FFT, only removed one doing invFFT.
     This takes out the second factor.
  */
  for(int i=0; i<len; i++){
    crossCorr[i] /= len;
  }

  return crossCorr;
}


double* FancyFFTs::crossCorrelate(int len, std::complex<double>* fft1, std::complex<double>* fft2,
				  int threadInd){
  /* 
     Cross correlation is the same as bin-by-bin multiplication in the frequency domain.
     Will assume lengths are the same for now.
  */


  /* Stops tempVals returning NULL, normally done in doFFT step. 
     But the nice thing about this class is it means that I won't be duplicating work. */
  makeNewPlanIfNeeded(len, threadInd);

  std::pair<int, int> key(len, threadInd);
  
  /* Grab array associated with plan from internal memory */
  std::complex<double>* tempVals = (std::complex<double>*) fComplex[key];

  // TThread::Lock();
  // std::cout << "threadInd = " << threadInd << "\tfComplex[key] = " << tempVals << std::endl << std::endl;
  // TThread::UnLock();
  
  
  /* Take the product */
  int numFreqs = getNumFreqs(len);

  for(int i=0; i<numFreqs; i++){
    tempVals[i] = fft1[i]*std::conj(fft2[i]);
  }
  
  /* Product back to time domain */
  double* crossCorr = doInvFFT(len, tempVals, true, threadInd);

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
