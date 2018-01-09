#include <iostream>
#include "SummarySet.h"
#include "DrawStrings.h"
#include "ProgressBar.h"
#include "CutTreeSelector.h"
#include "OutputConvention.h"


using namespace Acclaim;

int main(int argc, char* argv[]){

  if(!(argc == 3 || argc == 2)){
    std::cerr << argv[0] << " 'glob' " << std::endl;
    std::cerr << argv[0] << " 'glob' run " << std::endl;
    return 1;
  }
  const char* outFileName = argv[0];
  const char* glob = argv[1];
  TString runStr = argc == 3 ? argv[2] : "-1";
  Int_t run = atoi(runStr.Data());

  std::vector<const char*> formulas;
  formulas.push_back("run");
  formulas.push_back("eventNumber");
  formulas.push_back("realTime");

  // formulas.push_back("sum.weight()");
  formulas.push_back(SumTree::weight);

  bool mcInput = TString(glob).Contains("_mc_");
  formulas.push_back("mc.phi");
  formulas.push_back("mc.theta");
  formulas.push_back("mc.energy");

  formulas.push_back("peak[][].value");
  formulas.push_back("peak[][].phi");
  formulas.push_back("peak[][].theta");
  formulas.push_back("Iteration$");
  formulas.push_back(SumTree::pol);
  formulas.push_back(SumTree::peakInd);
  // formulas.push_back("sum.highestPeak().value");
  // formulas.push_back("sum.highestPeak().phi");
  // formulas.push_back("sum.highestPeak().theta");
 
  // formulas.push_back("sum.highestPolAsInt()");
  // formulas.push_back("sum.highestPeakInd()");
  // formulas.push_back("sum.highestPeakInd()");
  // formulas.push_back("sum.highestPeakInd()");
  
  formulas.push_back(SumTree::coherent_filtered_fracPowerWindowGradient);
  formulas.push_back(SumTree::deconvolved_filtered_fracPowerWindowGradient);
  // formulas.push_back("sum.highestDeconvolvedFiltered().fracPowerWindowGradient()");
  // formulas.push_back("sum.highestCoherentFiltered().fracPowerWindowGradient()");

  formulas.push_back("coherent_filtered[][].impulsivityMeasure");
  formulas.push_back("deconvolved_filtered[][].impulsivityMeasure");  
  // formulas.push_back("sum.highestDeconvolvedFiltered().impulsivityMeasure");
  // formulas.push_back("sum.highestCoherentFiltered().impulsivityMeasure");
  
  formulas.push_back("coherent_filtered[][].peakHilbert");
  formulas.push_back("deconvolved_filtered[][].peakHilbert");

  formulas.push_back("coherent_filtered[][].peakTime");
  formulas.push_back("deconvolved_filtered[][].peakTime");

  // formulas.push_back("sum.highestDeconvolvedFiltered().peakHilbert");
  // formulas.push_back("sum.highestCoherentFiltered().peakHilbert");

  formulas.push_back("coherent_filtered[][].snr");
  formulas.push_back("deconvolved_filtered[][].snr");

  formulas.push_back("flags.middleOrBottomPower[0]");
  formulas.push_back("flags.middleOrBottomPower[1]");
  formulas.push_back("flags.topPower[0]");
  formulas.push_back("flags.topPower[1]");
  formulas.push_back("flags.maxBottomToTopRatio[0]");
  formulas.push_back("flags.maxBottomToTopRatio[1]");
  formulas.push_back("flags.pulser");
  formulas.push_back("flags.isRF");
  formulas.push_back("flags.isPayloadBlast");
  formulas.push_back("flags.isVarner");
  formulas.push_back("flags.isVarner2");

  std::vector<const TCut *> cuts;
  TCut runCut(TString::Format("run == %d", run).Data());  
  if(run >= 0){
    cuts.push_back(&runCut);
  }
  
  cuts.push_back(&SumTree::highestPeak); // Best peak?
  
  TString hashString = "";
  for(UInt_t i=0; i < formulas.size(); i++){
    hashString += formulas[i];
  }
  for(UInt_t i=0; i < cuts.size(); i++){
    hashString += cuts[i]->GetTitle();
  }
  
  TString proofFileName = TString::Format("/tmp/thermalFile%u.root", hashString.Hash());


  const TString treeName = "thermalTree";
  CutTreeSelector ct(proofFileName, treeName);
  for(UInt_t i=0; i < cuts.size(); i++){
    ct.addCut(cuts.at(i));
  }

  ct.setFormulaStrings(formulas);





  ProgressBar p(1);
  SummarySet ss(glob);
  std::cout << "ss.N() = " << ss.N() << std::endl;
  // ss.SetUseProof(true);
  ss.SetUseProof(false);
  ss.Process(&ct);

  TString theRootPwd = gDirectory->GetPath();
  TFile* thermalFile = TFile::Open(proofFileName);
  TTree* thermalTree = (TTree*) thermalFile->Get(treeName);
  std::cout << treeName << " has " << thermalTree->GetEntries() << std::endl;

  std::vector<const char*> fakeArgv;
  fakeArgv.push_back(outFileName);
  
  if(mcInput){
    fakeArgv.push_back("mc");    
  }
  if(run >= 0){
    fakeArgv.push_back(runStr.Data());  
  }
  OutputConvention oc(fakeArgv.size(), (char**)&fakeArgv[0]);
  TFile* fOut = oc.makeFile();

  fOut->cd();
  TTree* tempTree = (TTree*) thermalTree->CloneTree();
  tempTree->Write();
  delete tempTree;
  
  gDirectory->cd(theRootPwd);

  fOut->Write();
  fOut->Close();
  p++;
  
  return 0;
}
