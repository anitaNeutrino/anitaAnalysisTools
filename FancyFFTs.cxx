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
std::map<int, fftw_plan> FancyFFTs::fRealToComplex;
std::map<int, fftw_plan> FancyFFTs::fComplexToReal;
std::map<int, double*> FancyFFTs::fReals;
std::map<int, std::complex<double>*> FancyFFTs::fComplex;

FancyFFTs::FancyFFTs(){
  std::cout << "FancyFFTs::FancyFFTs()" << std::endl;
}

FancyFFTs::~FancyFFTs(){
  std::cout << "FancyFFTs::~FancyFFTs()" << std::endl;
}


bool FancyFFTs::makeNewPlanIfNeeded(int len){
  /* 
     Function which checks whether we've encountered a request to do an FFT of this length before.
     If we haven't then we need a new plan!
  */
  std::map<int,fftw_plan>::iterator it = fRealToComplex.find(len);
  if(it==fRealToComplex.end()){
    fReals[len] = (double*) fftw_malloc(sizeof(double)*len);
    fComplex[len] = (std::complex<double>*) fftw_malloc(sizeof(fftw_complex)*len);
    fRealToComplex[len] = fftw_plan_dft_r2c_1d(len,fReals[len],(fftw_complex*)fComplex[len],FFTW_MEASURE);
    fComplexToReal[len] = fftw_plan_dft_c2r_1d(len,(fftw_complex*)fComplex[len],fReals[len],FFTW_MEASURE);
    return true;
  }
  else{
    return false;
  }
}

double* FancyFFTs::getPowerSpectrum(int len, double* input, double dt, PowSpecNorm::conventionFlag normFlag){

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
  doFFT(len, input, false);

  /* Get the fftw_malloc'd array that the plan uses */
  std::complex<double>* rawFftOutputPtr = (std::complex<double>*) fComplex[len];

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

std::complex<double>* FancyFFTs::doFFT(int len, double* input, bool copyOutputToNewArray){
  /* 
     Using std::complex<double> instead of the typdef fftw_complex double[2] 
     because CINT has a better time with it, even though it's (apparently) exactly the same.
  */

  makeNewPlanIfNeeded(len);

  std::memcpy(fReals[len], input, sizeof(double)*len);

  fftw_execute(fRealToComplex[len]);

  if(copyOutputToNewArray==true){
    int numFreqs = getNumFreqs(len);
    std::complex<double>* output = new std::complex<double>[numFreqs];

    /* Seems to work, see http://www.fftw.org/doc/Complex-numbers.html */
    std::memcpy(output, fComplex[len], sizeof(fftw_complex)*numFreqs);
    return output;
  }
  else{
    return NULL;
  }
}


double* FancyFFTs::doInvFFT(int len, std::complex<double>* input, bool copyOutputToNewArray){
  
  /* 
     Normalization of 1/N done in this function. 
     Note: fftw_plan_c2r_1d DESTROYS THE INPUT ARRAY when executed.
  */

  makeNewPlanIfNeeded(len);
  int numFreqs = getNumFreqs(len);

  std::memcpy(fComplex[len], input, sizeof(fftw_complex)*numFreqs);
  fftw_execute(fComplexToReal[len]);
  
  /* Normalization needed on the inverse transform */
  double* invFftOutPtr = fReals[len];
  for(int i=0; i<len; i++){
    invFftOutPtr[i]/=len;
  }

  if(copyOutputToNewArray==true){
    double* output = new double[len];
    std::memcpy(output, fReals[len], sizeof(double)*len);
    return output;
  }
  else{
    return NULL;
  }
}


int FancyFFTs::extendToPowerOfTwo(int len){
  return TMath::CeilNint(TMath::Log2(len));
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


int FancyFFTs::printListOfKeys(){
  /* Returns number of plans, prints the lengths of the plans to screen.*/

  std::vector<int> keys;

  for(std::map<int,fftw_plan>::iterator it = fRealToComplex.begin(); it != fRealToComplex.end(); it++){
    keys.push_back(it->first);
  }
  
  /* Pretty sure in advance of testing that this list is not guarenteed to be sorted. */
  std::vector<int> sortedIndices(keys.size());
  TMath::Sort(int(keys.size()), &keys[0], &sortedIndices[0], kFALSE);

  /* Print to terminal */
  std::cout << "Plan lengths in memory = [";
  for(int i=0; i<int(keys.size()); i++){
    int j = sortedIndices.at(i);
    std::cout << keys.at(j);
    if(i < int(keys.size())-1 ) std::cout << ", ";
  }
  std::cout << "]" << std::endl;

  return int(keys.size());
}


double* FancyFFTs::crossCorrelate(int len, double* v1, double* v2){
  /* 
     Cross correlation is the same as bin-by-bin multiplication in the frequency domain.
     Will assume lengths are the same for now.
  */

  /* Store output of FFT1 in tempVals */
  std::complex<double>* tempVals1 = doFFT(len, v1, true);

  /* Leave output of FFT2 in internal arrays to avoid an unnessary copy */
  doFFT(len, v2, false);

  /* Get pointer to internal array */
  std::complex<double>* tempVals2 = (std::complex<double>*) fComplex[len];

  /* Take the product */
  int numFreqs = getNumFreqs(len);
  for(int i=0; i<numFreqs; i++){
    // std::cout << tempVals1[i] <<"\t" << tempVals2[i] << std::endl;
    tempVals1[i] *= std::conj(tempVals2[i]);
  }
  
  /* Product back to time domain */
  double* crossCorr = doInvFFT(len, tempVals1, true);

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
