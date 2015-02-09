/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             Slightly complicated functions to manipulate ROOT objects should go in here.
	     I want to only ever write them once.
	     Let's make everything static.
*************************************************************************************************************** */

#ifndef ROOTTOOLS_H
#define ROOTTOOLS_H

#include <TGraph.h>
#include <TMath.h>

class RootTools: public TObject{

public:
  RootTools();
  ~RootTools();

  static void getMaxMin(TGraph* gr, Double_t& max, Double_t& min);
  static void subtractOffset(TGraph* gr, Double_t offset);
  static void normalize(TGraph*, Double_t& mean, Double_t& rms);
  static void getMeanAndRms(TGraph* gr, Double_t& mean, Double_t& rms);
  static void getMaxMin(TGraph* gr, Double_t& maxY, Double_t& maxX, Double_t& minY, Double_t& minX);
  
  ClassDef(RootTools, 1);


};

#endif
