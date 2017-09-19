#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "CutOptimizer.h"
#include <iostream>
#include "TF1.h"
#include "TGraph.h"
#include "TKey.h"
#include "RootTools.h"
#include "TH2D.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"

std::vector<TString> fisherComponents;
TTree* signalTree;
TTree* backgroundTree;

TString makeCommand(std::vector<TString>& components,  const double* weights){
  TString drawCommand;
  for(UInt_t i=0; i < components.size(); i++){
    if(i > 0){
      drawCommand += " + ";
    }
    drawCommand += TString::Format("(%lf*%s)", weights[i], components[i].Data());

  }
  return drawCommand;
}

Double_t overlap(const double* pars){

  TString drawCommand = makeCommand(fisherComponents, pars);

  TH1D hs("hs", "hs", 1024, -1, 1);
  hs.SetCanExtend(TH1::kAllAxes);
  TH1D hb("hb", "hb", 1024, -1, 1);
  hb.SetCanExtend(TH1::kAllAxes);
  
  signalTree->Draw(drawCommand + ">>hs", "weight", "goff");
  backgroundTree->Draw(drawCommand + ">>hb", "weight", "goff");
  


  int lastNonZeroBackgroundBin = hb.FindLastBinAbove(0);

  double signalBelowLastBackground = hs.Integral(1, lastNonZeroBackgroundBin);
  double signalAboveLastBackground = hs.Integral(lastNonZeroBackgroundBin, hs.GetNbinsX()+1);  

  // double sMean = hs.GetMean();
  // double bMean = hb.GetMean();

  // std::cout << drawCommand << std::endl;
  double fom = signalBelowLastBackground/signalAboveLastBackground;
  std::cout << signalBelowLastBackground << "\t" << signalAboveLastBackground << "\t" << fom << std::endl;

  return signalBelowLastBackground/signalAboveLastBackground;
}


void minimizeOverlap(const char* fileName){

  TFile* f = TFile::Open(fileName);
  signalTree = (TTree*) f->Get("signalTree");
  backgroundTree = (TTree*) f->Get("backgroundTree");
  signalTree->Show(0);
  backgroundTree->Show(0);  
  Acclaim::CutOptimizer::FisherResult* fr = (Acclaim::CutOptimizer::FisherResult*) f->Get("FisherResult");
  fr->getExpressions(fisherComponents);
  // TString fisherForm = fr->getFisherFormula();
  // Acclaim::RootTools::tokenize(fisherComponents,  fisherForm.Data(), std::vector<const char*>{"+(",  ")"});
  
  std::vector<double> vars(fisherComponents.size());
  for(UInt_t i=0; i < fisherComponents.size(); i++){
    vars[i] = fr->getWeight(fisherComponents[i]);
    std::cerr << fisherComponents[i] << "\t" << vars[i] << std::endl;
  }
  
  std::cout << "Initial overlap = " << overlap(&vars[0]) << std::endl;

  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2");

  // set tolerance , etc...
  min->SetMaxFunctionCalls(1000); // for Minuit/Minuit2 
  // min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);
   
  // create funciton wrapper for minmizer
  // a IMultiGenFunction type 
  ROOT::Math::Functor func(&overlap,fisherComponents.size());
   
  std::vector<double> steps(fisherComponents.size(), 1);
  // starting point
    
  min->SetFunction(func);
 
  // Set the free variables to be minimized!
  for(UInt_t i=0; i < fisherComponents.size(); i++){
    min->SetVariable(i,fisherComponents[i].Data(),vars[i], steps[i]);
  }
 
  // do the minimization
  min->Minimize(); 
 
  const double *xs = min->X();
  for(unsigned i=0; i < fisherComponents.size(); i++){
    std::cout << fisherComponents[i] << "\t" << xs[i] << std::endl;
  }
  
  TString drawCommand = makeCommand(fisherComponents, xs);
  TH1D* hs = new TH1D("hs", "hs", 1024, -1, 1);
  hs->SetCanExtend(TH1::kAllAxes);
  TH1D* hb = new TH1D("hb", "hb", 1024, -1, 1);
  hb->SetCanExtend(TH1::kAllAxes);

  signalTree->Draw(drawCommand + ">>hs", "weight");
  backgroundTree->Draw(drawCommand + ">>hb", "weight", "same");

  hs->SetLineColor(kRed);
  hb->SetLineColor(kBlue);
  
  // expected minimum is 0
  if ( min->MinValue()  < 1.E-4  && func(xs) < 1.E-4) {}
  else {
    std::cout << "Minimizer failed to converge !!!" << std::endl;
    Error("NumericalMinimization","fail to converge");
  }
 
  return;
}



