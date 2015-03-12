/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             Slightly complicated functions to manipulate ROOT objects should go in here.
	     I want to only ever write them once, so when I need to do something I'll add it here.
*************************************************************************************************************** */

#ifndef ROOTTOOLS_H
#define ROOTTOOLS_H

#include <TObjArray.h>
#include <TGraph.h>
#include <TTree.h>
#include <TAxis.h>
#include <TMath.h>
#include <iostream>


namespace RootTools{

  void getMaxMin(TGraph* gr, Double_t& max, Double_t& min);
  void subtractOffset(TGraph* gr, Double_t offset);
  void normalize(TGraph* gr, Double_t& mean, Double_t& rms);
  void normalize(TGraph* gr);
  void getMeanAndRms(TGraph* gr, Double_t& mean, Double_t& rms);
  void getMaxMin(TGraph* gr, Double_t& maxY, Double_t& maxX, Double_t& minY, Double_t& minX);

  void printArray(int n, double* array, TString delimiter, TString start ,TString end);
  void printYVals(TGraph* gr, TString delimiter, TString start, TString end);
  void printXVals(TGraph* gr, TString delimiter, TString start, TString end);
  double getSumOfYVals(TGraph* gr);

  TGraph* makeNormalized(TGraph* gr);
  TGraph* makeNormalized(TGraph* gr, Double_t& mean, Double_t& rms);

  std::vector<Int_t> getIndicesOfNans(TGraph* gr);
  void printTGraphInfo(TGraph* gr);
  
  TGraph* makeSortedTGraph(TTree* tree, TString drawText, TString cutString, Double_t wrapValue);
};

#endif
