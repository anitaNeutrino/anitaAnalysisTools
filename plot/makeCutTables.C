#include "DrawStrings.h"
#include "RootTools.h"
#include "ThermalChain.h"
#include "FFTtools.h"
#include "TProof.h"

using namespace Acclaim;




Double_t numPassing(TChain* c,  const TCut& cut, const char* hName = NULL){
  TString cutText = cut.GetTitle();
  TString histName = hName ? hName : "hCuts";
  c->Draw(cutText + ">>" + histName + "(2, 0, 2)", Acclaim::ThermalTree::weight(cut), "goff");  
  TH1F* hCut = (TH1F*) gROOT->FindObject(histName);
  bool needHandleProofNonsense = false;
  if(!hCut && gProof){
    hCut = (TH1F*) gProof->GetOutputList()->FindObject(histName);
    needHandleProofNonsense = true;
  }

  Double_t n = -1;
  if(hCut){
    n = hCut->GetBinContent(2); // should be weight passing
  }
  
  if(needHandleProofNonsense){
    gProof->GetOutputList()->Clear();
  }
  delete hCut;
  
  return n;
}


void getCutTableParams(TChain* c, bool mc, std::vector<TCut>& cuts, std::vector<Double_t>& inSeq, std::vector<Double_t>& ifOnly, std::vector<Double_t>& ifNot){
  // cut 0 is preselection
  std::cout << "cuts.size() == " << cuts.size() << std::endl;
  inSeq.reserve(cuts.size());
  ifOnly.reserve(cuts.size());
  ifNot.reserve(cuts.size());
  
  if(cuts.size() > 1){
    int numCuts = cuts.size();
    const int whichCut = 0;
    for(int cutInd=whichCut; cutInd < TMath::Max(whichCut + 1, numCuts); cutInd++){
    
    // for(int cutInd=0; cutInd < numCuts; cutInd++){

      TCut seqCut = cuts[0];
      for(int c=1; c <= cutInd; c++){
	seqCut += cuts[c];
      }
      // std::cout << seqCut.GetTitle() << std::endl;
      
      TCut onlyCut = cuts[0];
      if(cutInd > 0){
	onlyCut += cuts[cutInd];
      }
      
      // std::cout << onlyCut.GetTitle() << std::endl;

      TCut notCut = cuts[0];
      for(int c=1; c < numCuts; c++){
	if(c!=cutInd){
	  notCut += cuts[c];
	}
      }
      // std::cout << notCut.GetTitle() << std::endl;
      // std::cout << std::endl << std::endl << std::endl;
      double seq = 0;
      double only = 0;
      double no = 0;

      TString histNameBase = TString::Format("hCut_%d", cutInd);
      TString hNameSeq = histNameBase + "_seq";
      TString hNameOnly  = histNameBase + "_only";
      TString hNameNo = histNameBase + "_no";
      seq = numPassing(c, seqCut, hNameSeq.Data());
      only = numPassing(c, onlyCut, hNameOnly.Data());
      no = numPassing(c, notCut, hNameNo.Data());

      inSeq.push_back(seq);
      ifOnly.push_back(only);
      ifNot.push_back(no);
      std::cout << cutInd << "\t" << int(seq) << "\t" << int(only) << "\t" << int(no) << std::endl;
      // std::cout << cutInd << "\t" << seq << "\t" << only << "\t" << no << std::endl;      
    }
  }
}



inline TString result(double n, double total){
  return TString::Format("%.2lf (%.2lf%%)", n, 100*n/total);
}

inline TString result(Int_t n, Int_t total){
  return TString::Format("%d (%.2lf%%)", n, 100*double(n)/total);
}



void makeCutTables(const char* thermalTreeDataGlob="data/makeThermalTree_*.root", const char*  thermalTreeMcGlob="mc/makeThermalTree_mc_*.root"){

  gROOT->ProcessLine("#include \"FFTtools.h\"");
  
  ThermalChain c(thermalTreeDataGlob);
  ThermalChain c2(thermalTreeMcGlob);

  c.SetUseProof(1);
  
  const TCut allPass("true", "1");
  const TCut mcPreselection("mcPreselection", "realTime >= 1418895837");
  std::vector<TString> quickNames;
  std::vector<TCut> waisCuts;
  std::vector<TCut> mcNuCuts;
  std::vector<TCut> dataCutsVPol;
  std::vector<TCut> dataCutsHPol;

  const int n = 10;
  waisCuts.reserve(n);
  mcNuCuts.reserve(n);
  dataCutsVPol.reserve(n);
  dataCutsHPol.reserve(n);
  

  /// 0. pre-selection
  dataCutsHPol.push_back(ThermalTree::analysisSample + TCut("pol==0"));
  dataCutsVPol.push_back(ThermalTree::analysisSample + TCut("pol==1"));
  waisCuts.push_back(ThermalTree::isTaggedAsWaisPulser);
  mcNuCuts.push_back(mcPreselection);
  quickNames.push_back("None");

  
  /// 1. reco efficiency
  dataCutsVPol.push_back(allPass);
  dataCutsHPol.push_back(allPass);  
  waisCuts.push_back(ThermalTree::closeToWais);
  mcNuCuts.push_back(ThermalTree::closeToMC);
  quickNames.push_back("Close to truth");

  /// 2. quality cuts
  dataCutsVPol.push_back(ThermalTree::passAllQualityCuts);
  dataCutsHPol.push_back(ThermalTree::passAllQualityCuts);  
  waisCuts.push_back(ThermalTree::passAllQualityCuts);
  mcNuCuts.push_back(ThermalTree::passAllQualityCuts);
  quickNames.push_back("Quality cuts");  
  
  /// 3. thermal cut
  dataCutsHPol.push_back(ThermalTree::fisherCut);
  dataCutsVPol.push_back(ThermalTree::fisherCut);  
  waisCuts.push_back(ThermalTree::fisherCut);
  mcNuCuts.push_back(ThermalTree::fisherCut);
  quickNames.push_back("Pass thermal cut");

  /// 4. hical cut
  dataCutsHPol.push_back(TCut(!ThermalTree::closeToHiCal));
  dataCutsVPol.push_back(TCut(!ThermalTree::closeToHiCal));  
  waisCuts.push_back(TCut(!ThermalTree::closeToHiCal));
  mcNuCuts.push_back(TCut(!ThermalTree::closeToHiCal));
  quickNames.push_back("Not HiCal");


  bool doVPol = false;
  if(doVPol){
    std::vector<Double_t> dataSeqVPol, dataIfOnlyVPol, dataIfNotVPol;
    getCutTableParams(c.getChain(), true, dataCutsVPol, dataSeqVPol, dataIfOnlyVPol, dataIfNotVPol);
    for(int i=0; i < dataCutsVPol.size(); i++){
      const char* label = i < quickNames.size() ? quickNames[i].Data() : "???";
      std::cout << " | " << label << " | " << int(dataSeqVPol[i]) << " | " << int(dataIfOnlyVPol[i]) << " | " << int(dataIfNotVPol[i]) << " | " << std::endl;    
    }
    return;
  }


  bool doHPol = true;
  if(doHPol){  
    std::vector<Double_t> dataSeqHPol, dataIfOnlyHPol, dataIfNotHPol;
    getCutTableParams(c.getChain(), true, dataCutsHPol, dataSeqHPol, dataIfOnlyHPol, dataIfNotHPol);
    for(int i=0; i < dataCutsHPol.size(); i++){
      const char* label = i < quickNames.size() ? quickNames[i].Data() : "???";
      std::cout << " | " << label << " | " << int(dataSeqHPol[i]) << " | " << int(dataIfOnlyHPol[i]) << " | " << int(dataIfNotHPol[i]) << " | " << std::endl;        
    }
    return;
  }  

  std::vector<Double_t> mcNuSeq, mcNuIfOnly, mcNuIfNot;
  getCutTableParams(c2.getChain(), true, mcNuCuts, mcNuSeq, mcNuIfOnly, mcNuIfNot);
  for(int i=0; i < mcNuSeq.size(); i++){
    const char* label = i < quickNames.size() ? quickNames[i].Data() : "???";
    std::cout << " | " << label << " | " << mcNuSeq[i] << " | " << mcNuIfOnly[i] << " | " << mcNuIfNot[i] << " | " << std::endl;    
  }
  return;
  std::vector<Double_t> waisSeq, waisIfOnly, waisIfNot;
  getCutTableParams(c.getChain(), false, waisCuts, waisSeq, waisIfOnly, waisIfNot);
  std::cout << "cut | is seq | if only | if not |" << std::endl;
  for(int i=0; i < waisSeq.size(); i++){
    const char* label = i < quickNames.size() ? quickNames[i].Data() : "???";
    std::cout << " | " << label << " | " << waisSeq[i] << " | " << waisIfOnly[i] << " | " << waisIfNot[i] << " | " << std::endl;
  }
  
  return;
  
  
  // Double_t w = numPassing(c, "");
  // Double_t n = numPassing(c2, "");
  Int_t nData0 = numPassing(c.getChain(), ThermalTree::analysisSample);
  Int_t nWais0 = numPassing(c.getChain(), ThermalTree::isTaggedAsWaisPulser);
  Double_t nMc0 = numPassing(c2.getChain(), "");
  
  bool detailedQualityCuts = false;
  if(detailedQualityCuts){
    Int_t nData1 = numPassing(c.getChain(), ThermalTree::analysisSample       + ThermalTree::passAllQualityCuts);
    Int_t nWais1 = numPassing(c.getChain(), ThermalTree::isTaggedAsWaisPulser + ThermalTree::passAllQualityCuts);
    Double_t nMc1 = numPassing(c2.getChain(),                           ThermalTree::passAllQualityCuts);
    return;
  }
  
  
  
  Int_t nData1 = numPassing(c.getChain(), ThermalTree::analysisSample       + ThermalTree::passAllQualityCuts);  
  Int_t nWais1 = numPassing(c.getChain(), ThermalTree::isTaggedAsWaisPulser + ThermalTree::passAllQualityCuts);
  Double_t nMc1 = numPassing(c2.getChain(),                           ThermalTree::passAllQualityCuts);




  Int_t nData2 = nData1; //numPassing(c.getChain(), ThermalTree::analysisSample + thermalTree::passAllQualityCuts);
  Int_t nWais2 = numPassing(c.getChain(), ThermalTree::isTaggedAsWaisPulser + ThermalTree::passAllQualityCuts + ThermalTree::closeToWais);
  Double_t nMc2 = numPassing(c2.getChain(),                           ThermalTree::passAllQualityCuts + ThermalTree::closeToMC);

  Int_t nData3 = numPassing(c.getChain(), ThermalTree::analysisSample       + ThermalTree::passAllQualityCuts                            + ThermalTree::fisherCut);
  Int_t nWais3 = numPassing(c.getChain(), ThermalTree::isTaggedAsWaisPulser + ThermalTree::passAllQualityCuts + ThermalTree::closeToWais + ThermalTree::fisherCut);
  Double_t nMc3 = numPassing(c2.getChain(),                           ThermalTree::passAllQualityCuts + ThermalTree::closeToMC   + ThermalTree::fisherCut);

  Int_t nData4 = numPassing(c.getChain(), ThermalTree::analysisSample       + ThermalTree::passAllQualityCuts                            + ThermalTree::fisherCut + !ThermalTree::closeToHiCal);
  Int_t nWais4 = numPassing(c.getChain(), ThermalTree::isTaggedAsWaisPulser + ThermalTree::passAllQualityCuts + ThermalTree::closeToWais + ThermalTree::fisherCut + !ThermalTree::closeToHiCal);
  Double_t nMc4 = numPassing(c2.getChain(),                           ThermalTree::passAllQualityCuts + ThermalTree::closeToMC   + ThermalTree::fisherCut + !ThermalTree::closeToHiCal);


  Int_t nData5 = numPassing(c.getChain(), ThermalTree::analysisSample       + ThermalTree::passAllQualityCuts                            + ThermalTree::fisherCut + !ThermalTree::closeToHiCal + ThermalTree::continentNotIceShelf);
  Int_t nWais5 = numPassing(c.getChain(), ThermalTree::isTaggedAsWaisPulser + ThermalTree::passAllQualityCuts + ThermalTree::closeToWais + ThermalTree::fisherCut + !ThermalTree::closeToHiCal + ThermalTree::continentNotIceShelf);
  Double_t nMc5 = numPassing(c2.getChain(),                           ThermalTree::passAllQualityCuts + ThermalTree::closeToMC   + ThermalTree::fisherCut + !ThermalTree::closeToHiCal + ThermalTree::continentNotIceShelf);
  
  
  // reconstuction efficiency
  std::cout << "Cut\t"                       << "MC (%)"           << "\t" << "Wais (%)"             << "\t" << "Analysis Sample"      << std::endl;
  std::cout << "All\t"                       << result(nMc0, nMc0) << "\t" << result(nWais0, nWais0) << "\t" << result(nData0, nData0) << std::endl;
  std::cout << "QualityCuts\t"               << result(nMc1, nMc0) << "\t" << result(nWais1, nWais0) << "\t" << result(nData1, nData0) << std::endl;
  std::cout << "Reconstruction efficiency\t" << result(nMc2, nMc0) << "\t" << result(nWais2, nWais0) << "\t" << "n/a" << std::endl;
  std::cout << "Thermal Cut\t"               << result(nMc3, nMc0) << "\t" << result(nWais3, nWais0) << "\t" << result(nData3, nData0) << std::endl;
  std::cout << "HiCal Cut\t"                 << result(nMc4, nMc0) << "\t" << result(nWais4, nWais0) << "\t" << result(nData4, nData0) << std::endl;
  std::cout << "Continent & not ice shelf\t" << result(nMc5, nMc0) << "\t" << result(nWais5, nWais0) << "\t" << result(nData5, nData0) << std::endl;    

  
  
}
