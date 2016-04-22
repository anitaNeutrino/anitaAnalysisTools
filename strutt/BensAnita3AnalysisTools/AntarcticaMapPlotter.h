/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             Antarctica canvas for plotting locations on the continent class
************************************************************************************************************/


#ifndef ANTARCTICAMAPPLOTTER_H
#define ANTARCTICAMAPPLOTTER_H

#include "TCanvas.h"
#include "TH2D.h"
#include "TMath.h"
#include "TROOT.h"
#include "TImage.h"
#include "TGraph.h"

#include <iostream>
#include <map>

/**
 * @class AntarcticaMapPlotter
 * @brief Class to plot reconstructed events onto a picture of Antarctica.
 *
 * Much of the hard geometry stuff was done by Ryan and I've copied it from some EventCorrelator macros.
*/
class AntarcticaMapPlotter : public TNamed{

public:

  AntarcticaMapPlotter();
  AntarcticaMapPlotter(TString name, TString title, Int_t nBinsX, Int_t nBinsY);  
  ~AntarcticaMapPlotter();

  void addHistogram(TString name, TString title, Int_t nBinsX, Int_t nBinsY);
  void addTGraph(TString name, TString title, Int_t n=0, Double_t* latitude=NULL, Double_t* longitude=NULL);

  Int_t setCurrentHistogram(TString name);  
  void DrawHist(TString opt);
  TH2D* getCurrentHistogram();


  // Int_t GetN();
  // void SetPoint(Int_t n, Double_t latitude, Double_t longitude);
  Int_t Fill(Double_t latitude, Double_t longitude, Double_t weight=1);
  
  Int_t setCurrentTGraph(TString name);
  void DrawTGraph(TString opt);
  TGraph* getCurrentTGraph();

  
private:

  // Image scaling factors
  Double_t TrueScaleLat; //!< ????
  // Double_t CentralMeridian; //!< ????
  Double_t RadiusOfEarth; //!< ????
  Double_t xOffest; //!< ????
  Double_t yOffset; //!< ????
  Double_t scale; //!< ????
  Double_t xSize; //!< ????
  Double_t ySize; //!< ????

  // PNG of antarctica
  TImage *img; //!< The png image of Antarctica.

  void getRelXYFromLatLong(Double_t latitude, Double_t longitude,Double_t &x, Double_t &y);
  void initializeInternals();

  TH2D* hCurrent; //!< Current histogram
  TGraph* grCurrent; //!< Current TGraph

  std::map<TString, TH2D*> hists; //!< Stored histograms by name
  std::map<TString, TGraph*> grs; //!< Stored graphs by name

  ClassDef(AntarcticaMapPlotter,2);
};






#endif //ANTARCTICAMAPPLOTTER_H
