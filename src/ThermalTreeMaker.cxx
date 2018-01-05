#include "ThermalTreeMaker.h"
#include "TString.h"
#include "CutTreeSelector.h"
#include "SummarySet.h"
#include "OutputConvention.h"


Acclaim::ThermalTreeMaker::ThermalTreeMaker(const char* glob, const char* fileName){

  fSummaryGlob = glob;
  fFileName = fileName;
}







void Acclaim::ThermalTreeMaker::produce(const std::vector<const char*>& formulas, const std::vector<const TCut*> cuts){

  TString hashString = "";
  for(UInt_t i=0; i < formulas.size(); i++){
    hashString += formulas[i];
  }
  for(UInt_t i=0; i < cuts.size(); i++){
    hashString += cuts[i]->GetTitle();
  }

  
  TString proofFileName = TString::Format("/tmp/thermalFile%u.root", hashString.Hash());
  std::cout << proofFileName << std::endl;
  SummarySet ss(fSummaryGlob);
  ss.SetUseProof(true);

  const TString treeName = "thermalTree";
  CutTreeSelector ct(proofFileName, treeName);
  for(UInt_t i=0; i < cuts.size(); i++){
    ct.addCut(cuts.at(i));
  }

  ct.setFormulaStrings(formulas);

  ss.Process(&ct);
  

  TString theRootPwd = gDirectory->GetPath();
  TFile* thermalFile = TFile::Open(proofFileName);
  TTree* thermalTree = (TTree*) thermalFile->Get(treeName);
  std::cout << treeName << " has " << thermalTree->GetEntries() << std::endl;

  const char* fakeArgv[1] = {fFileName.Data()};
  OutputConvention oc(1, (char**)fakeArgv);
  TFile* fOutFile = oc.makeFile();  

  fOutFile->cd();
  TTree* tempTree = (TTree*) thermalTree->CloneTree();
  tempTree->Write();
  delete tempTree;
  TTree* acTree = (TTree*) thermalFile->Get("analysisCutTree");
  fOutFile->cd();
  TTree* tempTree2 = (TTree*) acTree->CloneTree();
  tempTree2->SetName("cutTree");
  tempTree2->Write();
  delete tempTree2;
  
  gDirectory->cd(theRootPwd);

  fOutFile->Write();
  fOutFile->Close();

}
