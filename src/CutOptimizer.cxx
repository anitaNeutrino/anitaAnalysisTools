#include "CutOptimizer.h"
#include "OutputConvention.h"
#include "TObject.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"


Acclaim::CutOptimizer::CutOptimizer(const TString& outFileName, const TString& signalFiles, const TString& backgroundFiles)
    : fOutFileName(outFileName), fSignalName(signalFiles), fBackgroundName(backgroundFiles),  fOutFile(NULL),
      fSignalFile(NULL), fSignalTree(NULL), fBackgroundFile(NULL), fBackgroundTree(NULL) {

  getSignalAndBackgroundTrees();

}

void Acclaim::CutOptimizer::makeOutputFile(){

  int argc = 1;
  
  if(fOutFileName==""){
    fOutFileName = "CutOptimizer";
  }
  
  const char* fNameChar = fOutFileName.Data();
  const char** argv = &fNameChar;

  OutputConvention oc(argc, (char**)argv);
  fOutFile = oc.makeFile();
}


void Acclaim::CutOptimizer::getSignalAndBackgroundTrees(){

  TString treeName = "sumTree";

  fSignalFile = TFile::Open(fSignalName);
  if(!fSignalFile){
    std::cerr << "Error in " <<__PRETTY_FUNCTION__ << ", unable to open signal file " << fSignalName << std::endl;
  }
  
  fSignalTree = (TTree*) fSignalFile->Get("sumTree");
  if(!fSignalTree){
    std::cerr << "Error in " <<__PRETTY_FUNCTION__ << ", unable to find " << treeName << " in " << fSignalName << std::endl;
  }
  


  fBackgroundFile = TFile::Open(fBackgroundName);
  if(!fBackgroundFile){
    std::cerr << "Error in " <<__PRETTY_FUNCTION__ << ", unable to open signal file " << fBackgroundName << std::endl;
  }
  
  fBackgroundTree = (TTree*) fBackgroundFile->Get("sumTree");
  if(!fBackgroundTree){
    std::cerr << "Error in " <<__PRETTY_FUNCTION__ << ", unable to find " << treeName << " in " << fBackgroundName << std::endl;
  }
  
  
}


void Acclaim::CutOptimizer::optimize(){

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,10,0)
 
 makeOutputFile();
  
  TMVA::Factory factory("ThermalCut", fOutFile, "V");
  TMVA::DataLoader dl;

  dl.AddSignalTree(fSignalTree);
  dl.AddBackgroundTree(fBackgroundTree);

  dl.AddVariable("peak[0][0].value");

  factory.BookMethod(&dl, TMVA::Types::EMVA::kFisher, "Fisher");
      
  factory.TrainAllMethods();
  factory.TestAllMethods();
  factory.EvaluateAllMethods();

  fOutFile->Close();
  std::cout << "==> wrote root file " << fOutFile->GetName() << std::endl;
  std::cout << "==> TMVAnalysis is done!" << std::endl;

#else
  std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", TVMA changed a lot in recent ROOT versions. ";
  std::cerr << "This class requires ROOT version 6.10" << std::endl;
#endif

}
