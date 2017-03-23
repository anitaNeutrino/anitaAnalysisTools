/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             Hold some statistics about an event related to the rayliegh distributions.
***********************************************************************************************************/

#ifndef RICE_INFO_H
#define RICE_INFO_H

#include "RayleighInfo.h"

class RiceInfo : public RayleighInfo {


public:

  // Rice distributions are over certain number of events...
  // so track the first and last event number.

  // tf1 parameters
  Double_t riceFitNorm; // fixed for the fit
  Double_t riceFitAmp; //
  Double_t riceFitAmpError; //
  Double_t riceFitSignal; //
  Double_t riceFitSignalError; //
  
  // fit statistics
  Double_t riceChiSquare;
  Int_t riceNdf;

  RiceInfo();

  ClassDef(RiceInfo, 1)

};


#endif
