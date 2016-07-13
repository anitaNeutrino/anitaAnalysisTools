/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             C++ ROOT friendly class to do FFTs faster than FFTtools, but is mostly a shameless copy.
	     Will probably be pretty bare bones intially as I only really want this for doing Cross Correlations.
*************************************************************************************************************** */

#include "TMath.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TSystem.h"

#include <assert.h>


#include "FancyFFTs.h"
#include "RootTools.h"

#include <fftw3.h>

double sum(int n, double* x){
  double sum = 0;
  for(int i=0; i<n; i++){
    sum += x[i];
  }
  return sum;
}

int main(int argc, char *argv[]){

  TApplication* theApp = new TApplication("App", &argc, argv);

  const int n = 256;
  double omega1 = 0.7;
  double omega2 = 1.8;

  double ts[n];
  double dt = 0.1;
  double phi = TMath::Pi()/2;
  double sineWave1[n];
  double sineWave2[n];

  for(int i=0; i<n; i++){
    double t = dt*i;
    ts[i] = t;
    sineWave1[i] = 3.3*TMath::Sin(omega1*t + phi);
    sineWave2[i] = 1.6*TMath::Sin(omega2*t);
  }
  
  TGraph* gr1 = new TGraph(n, ts, sineWave1);
  gr1->SetName("grSineWave");
  gr1->SetTitle("A sine wave; t (s); Amplitude (V)");
  TGraph* gr2 = new TGraph(n, ts, sineWave2);
  gr2->SetLineColor(kRed);

  TCanvas* c1 = new TCanvas("c1", "Fancy FFT Frivolities", 1200, 1600);
  c1->Divide(2, 2);
  c1->cd(1);
  gr1->Draw();
  gr2->Draw("l same");

  double* fs = FancyFFTs::getFreqArray(n, dt);

  std::cout << "Before doing anything: ";
  FancyFFTs::printListOfKeys();

  double* ps1 = FancyFFTs::getPowerSpectrum(n, sineWave1, dt, FancyFFTs::kSum);
  double* ps2 = FancyFFTs::getPowerSpectrum(n, sineWave2, dt, FancyFFTs::kSum);

  std::cout << "After makeing two power spectrums of length " << n<< ": ";
  FancyFFTs::printListOfKeys();  

  TGraph* gr3 = new TGraph(FancyFFTs::getNumFreqs(n), fs, ps1);
  TGraph* gr4 = new TGraph(FancyFFTs::getNumFreqs(n), fs, ps2);
  delete [] ps1;
  delete [] ps2;

  gr4->SetLineColor(kRed);
  gr3->SetTitle("Power spectrum of a sine wave; Frequency (Hz); Power (V^{2})");
  c1->cd(2);
  gr3->Draw();
  gr4->Draw("same");
  gPad->SetLogy(1);









  /* Tests...*/
  std::cout << std::endl << std::endl << std::endl;
  std::cout << "Doing some tests of normalization a la Parsevals Theorem." << std::endl;
  std::cout << "If the program doesn't creash unceremoniously then we're all good :)" << std::endl;

  /* Parsevals theroem (used to normalize) */
  double powTime = 0;
  for(int i=0; i<gr1->GetN(); i++){
    powTime += gr1->GetY()[i]*gr1->GetY()[i];
  }
  std::cout << "powTime = " << powTime << std::endl;

  double powFreq = RootTools::getSumOfYVals(gr3);

  std::cout << "powFreq = " << powFreq << std::endl;
  std::cout << "powFreq/powTime = " << powFreq/powTime << std::endl;
  std::cout << "powFreq-powTime = " << powFreq-powTime << std::endl;
  std::cout << std::endl << std::endl << std::endl;

  /* Hardcore test of parsevals theorem: kSum normalization. */
  assert(TMath::Abs(powFreq-powTime) < 1e-10);






  double* ps1_ave = FancyFFTs::getPowerSpectrum(n, sineWave1, dt, FancyFFTs::kAverage);
  double powFreq_ave = sum(FancyFFTs::getNumFreqs(n), ps1_ave);
  std::cout << "powFreq_ave = " << powFreq_ave << std::endl;
  std::cout << "powFreq_ave/powTime_ave = " << powFreq_ave/(powTime/n) << std::endl;
  std::cout << "powFreq_ave-powTime_ave = " << powFreq_ave-(powTime/n) << std::endl;
  std::cout << std::endl << std::endl << std::endl;

  /* Hardcore test of parsevals theorem: kAverage normalization. */
  assert(TMath::Abs(powFreq_ave-(powTime/n)) < 1e-10);  

  delete [] ps1_ave;


  double* ps1_timeIntegral = FancyFFTs::getPowerSpectrum(n, sineWave1, dt, FancyFFTs::kTimeIntegral);
  double powFreq_timeIntegral = sum(FancyFFTs::getNumFreqs(n), ps1_timeIntegral);  
  std::cout << "powFreq_timeIntegral = " << powFreq_timeIntegral << std::endl;
  std::cout << "powFreq_timeIntegral/powTime_ave = " << powFreq_timeIntegral/(dt*powTime) << std::endl;
  std::cout << "powFreq_timeIntegral-powTime_ave = " << powFreq_timeIntegral-(dt*powTime) << std::endl;
  std::cout << std::endl << std::endl << std::endl;

  /* Hardcore test of parsevals theorem: kTimeIntegral normalization. */
  assert(TMath::Abs(powFreq_timeIntegral-(dt*powTime)) < 1e-10);  

  delete [] ps1_timeIntegral;









  c1->cd(3);
  std::complex<double>* testFFT = FancyFFTs::doFFT(gr1->GetN(), gr1->GetY(), true);
  double* testInvFFT = FancyFFTs::doInvFFT(gr1->GetN(), testFFT, true);
  TGraph* grTestInv = new TGraph(gr1->GetN(), ts, testInvFFT);
  grTestInv->SetTitle("Test of inverse fourier transform; Time (s); Amplitude (V)");
  grTestInv->Draw();
  delete [] testFFT;
  delete [] testInvFFT;



  std::cout << "Checking conversion between fftw_complex and std::complex<double>... " << std::endl; 
  std::cout << "sizeof(std::complex<double>) = " << sizeof(std::complex<double>) << std::endl;
  std::cout << "sizeof(fftw_complex) = " << sizeof(fftw_complex) << std::endl;
  std::cout << "sizeof(double) = " << sizeof(double) << std::endl;

  std::complex<double> x(1, 1);
  std::cout << "std::complex<double> x(1, 1) = " << x << std::endl;
  fftw_complex y;
  std::cout << "Doing memcpy(&y, &x, sizeof(fftw_complex))... " << std::endl;
  memcpy(&y, &x, sizeof(fftw_complex));  
  std::cout << "fftw_complex y has y[0] = " << y[0] << " and y[1] = " << y[1] << std::endl;
  std::cout << std::endl << std::endl;

  std::cout << "Now copying array of fftw_complex to std::complex<double> with mempcy..." << std::endl;
  fftw_complex y2[2];
  y2[0][0] = 1;
  y2[0][1] = 2;
  y2[1][0] = 3;
  y2[1][1] = 4;
  std::cout << "y2[2] = {{" << y2[0][0] << ", " << y2[0][1] << "}, {" << y2[1][0] << ", " << y2[1][1] << "}}"
  	    << std::endl;
  std::complex<double> x2[2];
  std::cout << "Doing memcpy(x2, y2, 2*sizeof(fftw_complex))" << std::endl;
  memcpy(x2, y2, 2*sizeof(fftw_complex));
  std::cout << "x[0] = " << x2[0] << ", x[1] = " << x2[1] << std::endl;
  std::cout << std::endl << std::endl;




  c1->cd(4);
  TGraph* gr1norm = RootTools::makeNormalized(gr1);
  TGraph* gr2norm = RootTools::makeNormalized(gr2);
  double* crossCorr1 = FancyFFTs::crossCorrelate(gr1norm->GetN(), gr1norm->GetY(), gr1norm->GetY());
  double* crossCorr2 = FancyFFTs::crossCorrelate(gr1norm->GetN(), gr1norm->GetY(), gr2norm->GetY());
  TGraph* grCor1 = new TGraph(gr1->GetN(), ts, crossCorr1);
  TGraph* grCor2 = new TGraph(gr1->GetN(), ts, crossCorr2);
  grCor1->SetTitle("Cross correlation");
  grCor1->Draw();
  grCor2->SetLineColor(kRed);
  grCor2->Draw("l same");
  gPad->Update();









  /* Draw pretty things */

  c1->Update();
  gSystem->ProcessEvents();
  std::cerr << "Select File->Quit to quit." << std::endl;
  theApp->Run();
  
  return 0;
}
