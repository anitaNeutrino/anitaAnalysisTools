/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             Hold some statistics about an event related to the rayliegh distributions.
***********************************************************************************************************/

#ifndef RAYLEIGH_INFO_H
#define RAYLEIGH_INFO_H

#include "TH1D.h"
#include "TF1.h"
#include "TObject.h"

class RayleighInfo : public TObject {


public:

  // Rayleigh distributions are over certain number of events...
  // so track the first and last event number.

  // the final event in the distribution  
  Int_t run;
  UInt_t eventNumber;
  Double_t realTimeNs;

 // the first event in the distribution  
  Int_t firstRun;  
  UInt_t firstEventNumber;
  Double_t firstRealTimeNs;
  
  // the event amplitude
  Double_t eventAmp;

  
  // tf1 parameters
  Double_t rayFitNorm; // fixed for the fit
  Double_t rayFitAmp; //
  Double_t rayFitAmpError; //  
  Double_t rayGuessAmp; //
  
  // fit statistics
  Double_t rayChiSquare;
  Int_t rayNdf;

  // histogram statistics
  Int_t nBins;
  Double_t binWidth;  
  Double_t integralPlusOverUnderFlow;
  


  RayleighInfo();
  
  // For viewing
  TH1D* makeHistogram();


  ClassDef(RayleighInfo, 1)

};


#endif
