#include "DrawStrings.h"
#include "RootTools.h"
#include "ThermalChain.h"
#include "FFTtools.h"

using namespace Acclaim;

Double_t weightedNumPassing(TChain* c,  const TCut& cut){
  c->Draw("weight", cut, "goff");
  Int_t n = c->GetSelectedRows();
  Double_t* v1 = c->GetV1();
  return RootTools::sum(n, v1);
}

Int_t numPassing(TChain* c,  const TCut& cut){
  c->Draw("weight", cut, "goff");
  Int_t n = c->GetSelectedRows();
  return n;
}

inline TString result(double n, double total){
  return TString::Format("%.2lf (%.2lf%%)", n, 100*n/total);  
}

inline TString result(Int_t n, Int_t total){
  return TString::Format("%d (%.2lf%%)", n, 100*double(n)/total);  
}



void makeCutTables(const char* thermalTreeDataGlob="data/makeThermalTree*.root", const char*  thermalTreeMcGlob="mc/makeThermalTree_mc_*.root"){

  gROOT->ProcessLine("#include \"FFTtools.h\"");
  
  ThermalChain c(thermalTreeDataGlob);
  ThermalChain c2(thermalTreeMcGlob);

  // Double_t w = numPassing(c, "");
  // Double_t n = weightedNumPassing(c2, "");
  Int_t nData0 = numPassing(c.getChain(), ThermalTree::analysisSample);
  Int_t nWais0 = numPassing(c.getChain(), ThermalTree::isTaggedAsWaisPulser);
  Double_t nMc0 = weightedNumPassing(c2.getChain(), "");

  Int_t nData1 = numPassing(c.getChain(), ThermalTree::analysisSample       + ThermalTree::passAllQualityCuts);  
  Int_t nWais1 = numPassing(c.getChain(), ThermalTree::isTaggedAsWaisPulser + ThermalTree::passAllQualityCuts);
  Double_t nMc1 = weightedNumPassing(c2.getChain(),                           ThermalTree::passAllQualityCuts);


  Int_t nData2 = nData1; //numPassing(c.getChain(), ThermalTree::analysisSample + thermalTree::passAllQualityCuts);
  Int_t nWais2 = numPassing(c.getChain(), ThermalTree::isTaggedAsWaisPulser + ThermalTree::passAllQualityCuts + ThermalTree::closeToWais);
  Double_t nMc2 = weightedNumPassing(c2.getChain(),                           ThermalTree::passAllQualityCuts + ThermalTree::closeToMC);

  const TCut fisherCut("fisherCut", ThermalTree::fisherDiscriminant + " > 5.8007812");

  Int_t nData3 = numPassing(c.getChain(), ThermalTree::analysisSample       + ThermalTree::passAllQualityCuts                            + fisherCut);
  Int_t nWais3 = numPassing(c.getChain(), ThermalTree::isTaggedAsWaisPulser + ThermalTree::passAllQualityCuts + ThermalTree::closeToWais + fisherCut);
  Double_t nMc3 = weightedNumPassing(c2.getChain(),                           ThermalTree::passAllQualityCuts + ThermalTree::closeToMC   + fisherCut);
  
  // reconstuction efficiency
  std::cout << "Cut\t"                       << "MC (%)"           << "\t" << "Wais (%)"             << "\t" << "Analysis Sample"      << std::endl;
  std::cout << "All\t"                       << result(nMc0, nMc0) << "\t" << result(nWais0, nWais0) << "\t" << result(nData0, nData0) << std::endl;
  std::cout << "QualityCuts\t"               << result(nMc1, nMc0) << "\t" << result(nWais1, nWais0) << "\t" << result(nData1, nData0) << std::endl;
  std::cout << "Reconstruction efficiency\t" << result(nMc2, nMc0) << "\t" << result(nWais2, nWais0) << "\t" << result(nData2, nData0) << std::endl;
  std::cout << "Thermal Cut\t"               << result(nMc3, nMc0) << "\t" << result(nWais3, nWais0) << "\t" << result(nData3, nData0) << std::endl;

  
  
}
