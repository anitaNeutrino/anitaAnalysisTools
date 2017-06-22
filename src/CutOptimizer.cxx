#include "CutOptimizer.h"
#include "OutputConvention.h"
#include "TObject.h"


#if ROOT_VERSION_CODE >= ROOT_VERSION(6,10,0)
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodFisher.h"
#include "TXMLEngine.h"
#endif

bool debug = false;

static void setDebug(bool db){
  debug = db;
}


TH1D* Acclaim::CutOptimizer::FisherResult::makeHist(int nBinsX, const TString& histName, const TString& histTitle, TTree* t, EColor col) const{


  TString command;

  for(int i=0; i < fWeights.size(); i++){
    WeightMap::const_iterator wit = fWeights.find(i);
    if(wit!=fWeights.end()){
      TString w = TString::Format("%lf", wit->second);
      if(i==0){
        command = w;
      }
      else{
        ExpressionMap::const_iterator eit = fExpressions.find(i-1);
        if(eit!=fExpressions.end()){
          command += "+(" + w + "*" + eit->second + ")";
        }
        else{
          std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", couldn't find expression " << i-1 << std::endl;
        }
      }
    }
  }

  if(debug){
    std::cerr << command << std::endl;
  }
  TH1D* h = new TH1D(histName, histTitle, nBinsX, 0, 1);

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
    h->SetCanExtend(TH1::kAllAxes);
#else
    h->SetBit(TH1::kCanRebin);
#endif
  command += ">>" + histName;
  t->Draw(command, "" ,"goff");
  h->SetLineColor(col);
  
  return h;
}

void Acclaim::CutOptimizer::FisherResult::getResultFromXML(const char* filename){

  if(debug){
    std::cerr << "In " << __PRETTY_FUNCTION__ << ",  trying to read back in the xml file: " << filename << std::endl;
  }
  
  // First create engine
  TXMLEngine* xml = new TXMLEngine;
  // Now try to parse xml file
  // Only file with restricted xml syntax are supported
  XMLDocPointer_t xmldoc = xml->ParseFile(filename);
  if (xmldoc==0) {
    delete xml;
    return;
  }
  // take access to main node
  XMLNodePointer_t mainnode = xml->DocGetRootElement(xmldoc);
  // display recursively all nodes and subnodes
  parseNode(xml, mainnode, 1);
  // Release memory before exit
  xml->FreeDoc(xmldoc);
  delete xml;


  if(debug){
    std::cout << "I found expressions..." << std::endl;
    ExpressionMap::const_iterator it;
    for(it = fExpressions.begin(); it != fExpressions.end();  ++it){
      std::cout << it->first << "\t" << it->second << std::endl;
    }
   
    std::cout << "I found weights..." << std::endl;
    WeightMap::const_iterator it2;
    for(it2 = fWeights.begin(); it2 != fWeights.end();  ++it2){
      std::cout << it2->first << "\t" << it2->second << std::endl;
    }
  }
}



void Acclaim::CutOptimizer::FisherResult::parseNode(TXMLEngine* xml, XMLNodePointer_t node, Int_t level){

  // this function display all accessible information about xml node and its children
  // printf("%*c node: %s\n",level,' ', xml->GetNodeName(node));
  // // display namespace

  TString nodeName = xml->GetNodeName(node);

  if(nodeName=="Variable" || nodeName == "Coefficient"){

    TString varIndex;
    TString expression;

    TString index;
    TString weight;

    // XMLNsPointer_t ns = xml->GetNS(node);
    // if (ns!=0){
    //   printf("%*c namespace: %s refer: %s\n",level+2,' ', xml->GetNSName(ns), xml->GetNSReference(ns));
    // }
    
    // display attributes
    XMLAttrPointer_t attr = xml->GetFirstAttr(node);
    while (attr!=0) {
      TString attName(xml->GetAttrName(attr));
      TString attValue(xml->GetAttrValue(attr));      
      // printf("%*c attr: %s value: %s\n",level+2,' ', attName.Data(), attValue.Data());
      attr = xml->GetNextAttr(attr);

      if(attName == "VarIndex"){
        varIndex = attValue;
      }
      else if(attName == "Index"){
        index = attValue;
      }
      else if(attName == "Value"){
        weight = attValue;
      }
      else if(attName == "Expression"){
        expression = attValue;
      }      
    }

    if(varIndex!="" && expression != ""){
      int vi = atoi(varIndex.Data());
      fExpressions[vi] = expression;
    }
    else if (index!="" && weight != ""){
      int vi = atoi(index.Data());
      fWeights[vi] = atof(weight.Data());
    }
    
    // display content (if exists)
    // const char* content = xml->GetNodeContent(node);

    // if (content!=0){
    //   printf("%*c cont: %s\n",level+2,' ', content);
    // }
  }
  
  // parse all child nodes
  XMLNodePointer_t child = xml->GetChild(node);
  while (child!=0) {
    parseNode(xml, child, level+2);
    child = xml->GetNext(child);
  }
}












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

  TString factoryName = "ThermalCut";
  TString option = debug ? "V" : "silent";
  TMVA::Factory factory(factoryName, fOutFile, option);
  TString dlName = fOutFile->GetName();
  dlName.ReplaceAll(".root", "");
  TMVA::DataLoader dl(dlName);

  dl.AddSignalTree(fSignalTree);
  dl.AddBackgroundTree(fBackgroundTree);

  const int nWI = 4;
  const TString waveformInfos[nWI] = {"sum.higherCoherent()", "sum.higherCoherentFiltered()", "sum.higherDeconvolved()", "sum.higherDeconvolvedFiltered()"};  

  dl.AddVariable("sum.higherPeak().value");
  for(int wi=0; wi < nWI; wi++){
    dl.AddVariable(waveformInfos[wi] + ".peakHilbert");
  }

  TString methodTitle = "tc";

  factory.BookMethod(&dl, TMVA::Types::EMVA::kFisher, methodTitle, "");
  factory.TrainAllMethods();
  factory.TestAllMethods();
  factory.EvaluateAllMethods();

  TString weightFileName = dlName + "/weights/" + factoryName + "_" + methodTitle + ".weights.xml";

  FisherResult result(weightFileName.Data());
  fOutFile->cd();
  const int nBins = 1024;
  TH1D* hSignal = result.makeHist(nBins,  "hSignal", "Signal", fSignalTree, kRed);
  TH1D* hBackground = result.makeHist(nBins, "hBackground", "hBackground", fBackgroundTree, kBlue);  
  hSignal->Write();
  hBackground->Write();  
  delete hSignal;
  delete hBackground;
      
  fOutFile->Close();
  
  std::cout << "==> wrote root file " << fOutFile->GetName() << std::endl;
  std::cout << "==> TMVAnalysis is done!" << std::endl;

#else
  std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", TVMA changed a lot in recent ROOT versions. ";
  std::cerr << "This class requires ROOT version 6.10" << std::endl;
#endif

}



