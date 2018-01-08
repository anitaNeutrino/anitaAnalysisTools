#include "CutOptimizer.h"
#include "OutputConvention.h"
#include "SummaryDraw.h"
#include "SummarySet.h"
#include "ProgressBar.h"
#include "AnitaEventSummary.h"
#include "TTreeFormula.h"
#include "TApplication.h"
#include "TXMLEngine.h"
#include "TMath.h"
#include "TBranch.h"
#include "RootTools.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TEfficiency.h"
#include "TDirectory.h"
#include "TProof.h"

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,10,0)
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodFisher.h"
#endif

#include "CutTreeSelector.h"

ClassImp(Acclaim::CutOptimizer::FisherResult)

bool debug = false;


void Acclaim::CutOptimizer::setDebug(bool db){
  debug = db;
}

void Acclaim::CutOptimizer::FisherResult::getExpressions(std::vector<TString>& expressions) const {
  for(ExpressionMap::const_iterator it = fExpressions.begin(); it!= fExpressions.end(); it++){
    expressions.push_back(it->second);
  }
}

double Acclaim::CutOptimizer::FisherResult::getWeight(const char* expression){

  TString expr(expression);
  for(ExpressionMap::const_iterator it = fExpressions.begin(); it!= fExpressions.end(); it++){
    if(it->second == expr){
      return fWeights[it->first];
    }
  }
  return 0;
}


TString Acclaim::CutOptimizer::FisherResult::getFisherFormula() const{

  TString command;
  for(int i=0; i < (int)fWeights.size(); i++){
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
  return command;
}


void Acclaim::CutOptimizer::FisherResult::Print(Option_t* opt) const{
  (void) opt;
  std::cout << "FisherResult:" << std::endl;
  for(int i=0; i < (int)fWeights.size(); i++){
    WeightMap::const_iterator wit = fWeights.find(i);
    if(wit!=fWeights.end()){
      TString w = TString::Format("%lf", wit->second);
      if(i==0){
        std::cout << "Variable " << i << " = constant, weight = " << wit->second << std::endl;
      }
      else{
        ExpressionMap::const_iterator eit = fExpressions.find(i-1);
        if(eit!=fExpressions.end()){
          std::cout << "Variable " << i << " = " << eit->second << ", weight = " << w << std::endl;
        }
        else{
          std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", couldn't find expression " << i-1 << std::endl;
        }
      }
    }
  }
}



TH2D* Acclaim::CutOptimizer::FisherResult::makeTimeHist(int nBinsX, int nBinsY, TTree* t, EColor col, int varInd) const{  

  TString histName;
  TString histTitle;
  TString command;

  if(varInd==0){
    histName = TString::Format("h%sFisher", t->GetName());
    histTitle = TString::Format("%s: Fisher Discriminant", t->GetName());
    command = getFisherFormula();
  }
  else {
    WeightMap::const_iterator wit = fWeights.find(varInd);
    TString w = TString::Format("%lf", wit->second);
    ExpressionMap::const_iterator eit = fExpressions.find(varInd-1);
    if(eit!=fExpressions.end()){
      std::cout << "Variable " << varInd << " = " << eit->second << ", weight = " << w << std::endl;
    }
    else{
      std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", couldn't find expression " << varInd-1 << std::endl;
    }
    histName = TString::Format("h%s_%s", t->GetName(), eit->second.Data());
    histTitle = eit->second + " (weight " + w + ")";
    command = eit->second;
  }
  if(command==""){
    return NULL;
  }

  TH2D* h = new TH2D(histName, histTitle, nBinsX, RootTools::a3StartTime, RootTools::a3EndTime, nBinsY, 0, 1);

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
    h->SetCanExtend(TH1::kAllAxes);
#else
    h->SetBit(TH1::kCanRebin);
#endif
  command += ":realTime>>" + histName;
  t->Draw(command, "" ,"goff");
  std::cerr << h->Integral() << std::endl;
  h->SetLineColor(col);
  h->SetFillColor(col);

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
  SetTitle(getFisherFormula());
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










/** 
 * Constructor
 * 
 * @param summaryFileGlob, a glob expression for getting files containing trees of AnitaEventSummary
 * @param save_trees is a boolian which determines whether or not the intermediate TTrees fed into TMVA.
 */
Acclaim::CutOptimizer::CutOptimizer(const char* signalGlob, const char* backgroundGlob, bool doAllPeaks, bool save_trees)
  : fSignalGlob(signalGlob), fBackgroundGlob(backgroundGlob ? backgroundGlob : ""), fOutFile(NULL),
    fSignalTree(NULL), fBackgroundTree(NULL), fRejectedSignalTree(NULL), fRejectedBackgroundTree(NULL), fDoAllPeaks(doAllPeaks), fSaveTrees(save_trees)
{

}





/** 
 * Destructor
 */
Acclaim::CutOptimizer::~CutOptimizer()
{
  
}






/** 
 * Generate the output file
 * 
 * @param outFileName name to give the file
 * 
 * @return pointer to the newly created TFile
 */
TFile* Acclaim::CutOptimizer::makeOutputFile(const char* outFileName){

  int argc = 1;
  TString ofName = outFileName;
  if(ofName==""){
    ofName = "CutOptimizer";
  }
  
  const char* fNameChar = ofName.Data();
  const char** argv = &fNameChar;

  OutputConvention oc(argc, (char**)argv);
  TFile* f = oc.makeFile();
  return f;
}








/** 
 * Constructor for the FormulaHolder
 * 
 * @param c pointer to the TChain on which the formulae will be evaluated
 */
Acclaim::CutOptimizer::FormulaHolder::FormulaHolder(TChain* c)
    : fChain(c)
{
  if(fChain){
    fChain->SetNotify(this);
  }
}

/** 
 * TChain can only notify one object when the underlying tree changes.
 * The point of this class is to inherit from TObject and overload Notify(),
 * then pass that on to all the TFormulas I'm interested in.
 * 
 * @return always returns true.
 */
Bool_t Acclaim::CutOptimizer::FormulaHolder::Notify(){
  for(UInt_t i=0; i < fForms.size(); i++){
    fForms.at(i)->Notify();
  }
  return true;
}


/** 
 * Destructor: Delete all the contained TFormula and unset the TChain Notification
 */
Acclaim::CutOptimizer::FormulaHolder::~FormulaHolder(){
  for(UInt_t i=0; i < fForms.size(); i++){
    delete fForms.at(i);
  }
  if(fChain){
    fChain->SetNotify(NULL);
  }
  
}





/** 
 * Creates a TFormula object, and put it in the vector
 * 
 * @param formulaStr is the string from which to generate the formula
 * 
 * @return the number of contained TFormula 
 */
size_t Acclaim::CutOptimizer::FormulaHolder::add(const char* formulaStr){
  const int nForms = fForms.size();
  TString formName = TString::Format("form%d", nForms);
  TTreeFormula* form = new TTreeFormula(formName, formulaStr, fChain);

  // apprently this is the signature of an error
  // https://root-forum.cern.ch/t/check-the-status-of-a-ttreeformula/12596/2
  if(form->GetNdim()==0){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", not using bad expression "
              << formulaStr << ", will not included in the output tree" << std::endl;
    delete form;                                                                       
  }
  else{
    fForms.push_back(form);
    fFormStrs.push_back(formulaStr);
  }
  
  return fForms.size();
}






/** 
 * Remove function/object like things from the formula string but make a recognisable branch name
 * 
 * @param formStr is the c-string string containig the fromula
 * 
 * @return the branch name
 */
TString Acclaim::CutOptimizer::branchifyName(const char* formStr){
  TString bName = formStr;


  
  
  bName.ReplaceAll("sum.", ""); // remove the sum.
  bName.ReplaceAll("()", ""); // remove function call bracket
  bName.ReplaceAll("(", "_"); // remove standalone open paren...
  bName.ReplaceAll(")", "_"); // ... and close paren
  bName.ReplaceAll("[]", ""); // remove implicit loop iteration thingy
  bName.ReplaceAll("[", "_"); // or explicit array entry...
  bName.ReplaceAll("]", "_"); // ...
  bName.ReplaceAll(".", "_");  // remove remaining dots
  bName.ReplaceAll("TMath::", "");  // Remove TMath function namespace
  bName.ReplaceAll("FFTtools::", "");  // Remove TMath function namespace
  bName.ReplaceAll(",", "_");  // Remove TMath function namespace    
  bName.ReplaceAll("$", "");  // Remove special identifier
  
  //  bName.ReplaceAll("AnitaPol::", "");  // Remove AnitaPol function namespace

  bName.ReplaceAll("/", "_over_");  // wordify arithmetic/logical operators
  bName.ReplaceAll("+", "_plus_");  // wordify arithmetic/logical operators
  bName.ReplaceAll("-", "_minus_"); // wordify arithmetic/logical operators
  bName.ReplaceAll("*", "_times_"); // wordify arithmetic/logical operators
  bName.ReplaceAll(">=", "_ge_"); // wordify arithmetic/logical operators
  bName.ReplaceAll("<=", "_le_"); // wordify arithmetic/logical operators
  bName.ReplaceAll(">", "_gt_"); // wordify arithmetic/logical operators
  bName.ReplaceAll("<", "_lt_"); // wordify arithmetic/logical operators
  bName.ReplaceAll("%", "_modulo_"); // wordify arithmetic/logical operators
  bName.ReplaceAll(" ", "");  // Remove any remaining spaces  
  bName.ReplaceAll("__", "_");  // Remove double underscore if there is one
  
  bName.Remove(TString::kLeading, '_'); // Remove leading underscore if there is one
  bName.Remove(TString::kTrailing, '_'); // Remove trailing underscore if there is one




  // some special cases...
  
  if(bName=="10_0_times__minus_0_2_times_coherent_filtered_fracPowerWindowEnds_0_minus_coherent_filtered_fracPowerWindowBegins_0__minus_0_1_times_coherent_filtered_fracPowerWindowEnds_1_minus_coherent_filtered_fracPowerWindowBegins_1__plus_0_1_times_coherent_filtered_fracPowerWindowEnds_3_minus_coherent_filtered_fracPowerWindowBegins_3__plus_0_2_times_coherent_filtered_fracPowerWindowEnds_4_minus_coherent_filtered_fracPowerWindowBegins_4"){
    bName = "coherent_filtered_fracPowerWindowGradient";
  }
  else if(bName=="10_0_times__minus_0_2_times_deconvolved_filtered_fracPowerWindowEnds_0_minus_deconvolved_filtered_fracPowerWindowBegins_0__minus_0_1_times_deconvolved_filtered_fracPowerWindowEnds_1_minus_deconvolved_filtered_fracPowerWindowBegins_1__plus_0_1_times_deconvolved_filtered_fracPowerWindowEnds_3_minus_deconvolved_filtered_fracPowerWindowBegins_3__plus_0_2_times_deconvolved_filtered_fracPowerWindowEnds_4_minus_deconvolved_filtered_fracPowerWindowBegins_4"){
    bName = "deconvolved_filtered_fracPowerWindowGradient";
  }
  else if(bName=="_Abs_peak_hwAngle_lt_Abs_peak_hwAngleXPol__times_Abs_peak_hwAngle_plus_Abs_peak_hwAngle_ge_Abs_peak_hwAngleXPol__times_Abs_peak_hwAngleXPol"){
    bName = "minAbsHwAngle";
  }
  else if(bName=="wrap_peak_phi_minus_mc_phi_360_0"){
    bName = "dPhi_mc";
  }
  else if(bName=="_peak_theta_plus_mc_theta"){
    bName = "dTheta_mc";
  }  
  else if(bName=="wrap_peak_phi_minus_wais_phi_360_0"){
    bName = "dPhi_wais";
  }
  else if(bName=="wrap_peak_phi_minus_sun_phi_360_0"){
    bName = "dPhi_sun";
  }
  else if(bName=="floor_Iteration_over_5"){
    bName = "pol";
x  }
  else if(bName=="Iteration_modulo_5"){
    bName = "peakInd";
  }
  
  return bName;
}





/** 
 * Correctly assign the correct type to the branch from the TTreeFormula
 * 
 * @param t is the tree on which to assign the branch
 * @param f is the formula 
 * @param formulaString is the formula string, used to name the branch
 * @param intPtr is the address of the float we will right to if it's an integer
 * @param floatPtr is the address of the float we will right to if it's not an integer
 * 
 * @return 1 if an integer, 0 if a float
 */

Acclaim::CutOptimizer::BranchType Acclaim::CutOptimizer::setBranchFromFormula(TTree* t, const TTreeFormula* f, const char* formulaString, Int_t* intPtr, Float_t* floatPtr){
  TString bName = branchifyName(formulaString);

  BranchType b = f->IsInteger() ? kInt : kFloat;  
  if(b == kInt){
    t->Branch(bName, intPtr);
    // t->Branch(bName, intPtr, bName + "/I");
  }
  else{
    t->Branch(bName, floatPtr);
    // t->Branch(bName, floatPtr, bName + "/F");
  }
  if(debug){
    std::cerr << "Set branch " << bName << " in tree " << t->GetName() << " with status " << t->GetBranchStatus(bName) << std::endl;
  }
  return b;
}



void Acclaim::CutOptimizer::generateSignalAndBackgroundTreesProof(const std::vector<const TCut*>& signalSelection, const std::vector<const TCut*>& backgroundSelection, const std::vector<FormulaString>& formulaStrings){


  std::vector<const char*> fs;
  fs.reserve(formulaStrings.size());
  for(UInt_t i=0; i < formulaStrings.size(); i++){
    fs.push_back(formulaStrings[i].first);
  }


  // now we need to figure out if we have separate globs for signal and background
  // or just one for both...
  std::vector<TString> globs;
  globs.push_back(fSignalGlob);
  if(fBackgroundGlob != ""){
    globs.push_back(fBackgroundGlob);
  }
  else{
    globs.push_back(fSignalGlob); // again...
  }

  const std::vector<const TCut*>* selections[2] = {&signalSelection,
						   &backgroundSelection};
  const char* proofFileNames[2] = {"/tmp/signalCuts.root", "/tmp/backgroundCuts.root"};
  const char* proofTreeNames[2] = {"signalCuts", "backgroundCuts"};
  // for(UInt_t g=0; g < fSpecGlobs.size(); g++){
  //   globs.push_back(fSpecGlobs[g]);
  // }
  // // Then load the master AnitaEventSummary chain using my SummarySet class
  const int nG = globs.size();
  for(int g=0; g < nG; g++){
    SummarySet ss(globs[g]);
    ss.SetUseProof(true);

    CutTreeSelector ct(proofFileNames[g], proofTreeNames[g]);
    for(UInt_t i=0; i < selections[g]->size(); i++){
      ct.addCut(selections[g]->at(i));
    }

    ct.setFormulaStrings(fs);

    ss.Process(&ct);
  }

  TString theRootPwd = gDirectory->GetPath();
  TFile* signalFile = TFile::Open("/tmp/signalCuts.root");
  fSignalTree = (TTree*) signalFile->Get("signalCuts");
  std::cout << "signalTree has " << fSignalTree->GetEntries() << std::endl;

  if(fOutFile){
    fOutFile->cd();
    TTree* tempTree = (TTree*) fSignalTree->CloneTree();
    tempTree->Write();
    delete tempTree;
    TTree* acTree = (TTree*) signalFile->Get("analysisCutTree");
    fOutFile->cd();
    TTree* tempTree2 = (TTree*) acTree->CloneTree();
    tempTree2->SetName("signalAnalysisCutTree");
    tempTree2->Write();
    delete tempTree2;
  }
  
  TFile* backgroundFile = TFile::Open("/tmp/backgroundCuts.root");
  fBackgroundTree = (TTree*) backgroundFile->Get("backgroundCuts");
  std::cout << "backgroundTree has " << fBackgroundTree->GetEntries() << std::endl;

  if(fOutFile){
    fOutFile->cd();
    TTree* tempTree = (TTree*) fBackgroundTree->CloneTree();
    tempTree->Write();
    delete tempTree;
    TTree* acTree = (TTree*) backgroundFile->Get("analysisCutTree");
    fOutFile->cd();
    TTree* tempTree2 = (TTree*) acTree->CloneTree();
    tempTree2->SetName("backgroundAnalysisCutTree");    
    tempTree2->Write();
    delete tempTree2;
  }
  
  gDirectory->cd(theRootPwd);

}




void Acclaim::CutOptimizer::optimize(const std::vector<const TCut*>& signalSelection,
				     const std::vector<const TCut*>& backgroundSelection,
				     const std::vector<FormulaString>& formulaStrings, const char* fileName){

  // Someone with a bit more time can do the backwards compatibility
#if ROOT_VERSION_CODE < ROOT_VERSION(6,10,0)
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", TVMA changed a lot in recent ROOT versions. ";
    std::cerr << "This class requires ROOT version at least 6.10, you only have " << ROOT_VERSION_CODE << std::endl;
#else
  
    fOutFile = makeOutputFile(fileName);

    // generateSignalAndBackgroundTrees(signalSelection, backgroundSelection, formulaStrings);
    generateSignalAndBackgroundTreesProof(signalSelection, backgroundSelection, formulaStrings);
    // return;

    if(!fSignalTree || !fBackgroundTree){
      std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", Non-existent signal or background tree. "
                << "Can't optimize! Aborting!" << std::endl;
      return;
    }
    if(fSignalTree->GetEntries() == 0 || fBackgroundTree->GetEntries() == 0){    
      std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", "
                << "signal tree contains " << fSignalTree->GetEntries() << " events and "
                << "background tree contains " << fBackgroundTree->GetEntries() << " events. "
                << "Can't optimize! Aborting!" << std::endl;
      return;
    }
    
    // fOutFile->cd();
    // fSignalTree->Write();
    // fBackgroundTree->Write();

    TString factoryName = "ThermalCut";
    TString option = debug ? "V" : "silent";
    TMVA::Factory factory(factoryName, fOutFile, option);
    TString dlName = fOutFile->GetName();
    dlName.ReplaceAll(".root", "");
    TMVA::DataLoader dl(dlName);

    dl.AddSignalTree(fSignalTree);
    dl.AddBackgroundTree(fBackgroundTree);

    // should have identical for signal/background
    TObjArray* backBranches = fBackgroundTree->GetListOfBranches();
    for(int i=0; i < backBranches->GetEntries(); i++){
      TBranch* b = (TBranch*) backBranches->At(i);
      TString bName = b->GetName();

      if(bName == "weight"){
        dl.SetWeightExpression("weight");
      }
      else{
        for(std::vector<FormulaString>::const_iterator it = formulaStrings.begin(); it != formulaStrings.end(); ++it){
          // if the name matches and the FormulaString has true, then add
          if(bName == branchifyName(it->first) && it->second){
            dl.AddVariable(bName);
          }
        }
      }    
    }

    TString methodTitle = "tc";

    factory.BookMethod(&dl, TMVA::Types::EMVA::kFisher, methodTitle, "");

    // If an inf/nan ends up in the training sample TMVA will throw an exception.
    // We want to still write the signal/background trees to find out what went wrong,
    // so let's wrap the attempt in a try/catch block
    try{
      factory.TrainAllMethods();
      factory.TestAllMethods();
      factory.EvaluateAllMethods();

      TString weightFileName = dlName + "/weights/" + factoryName + "_" + methodTitle + ".weights.xml";

      FisherResult result(weightFileName.Data());
      fOutFile->cd();

      const int nBinsX = 1024;
      const int nBinsY = 1024;
      const int nT = 2;
      TTree* ts[nT] = {fSignalTree, fBackgroundTree};
    
      for(int t=0; t < nT; t++){
        EColor col = t == 0 ? kRed : kBlue;
        TObjArray* bs = ts[t]->GetListOfBranches();
        int nActive = 0;
        for(int b=0; b < bs->GetEntries(); b++){
          TBranch* br = (TBranch*) bs->At(b);

          int bs = ts[t]->GetBranchStatus(br->GetName());
          if(bs > 0){
            nActive++;
            TH2D* h = result.makeTimeHist(nBinsX, nBinsY, ts[t], col, nActive);
            if(h){
              h->Write();
              delete h;
            }
          }

          ts[t]->SetBranchStatus(br->GetName(), 1);                  
        }
      }

      TH2D* hSignal2 = result.makeTimeHist(nBinsX, nBinsY, fSignalTree, kRed);
      TH2D* hBackground2 = result.makeTimeHist(nBinsX, nBinsY, fBackgroundTree, kBlue);

      TH1D* hSignal = hSignal2->ProjectionY();
      TH1D* hBackground = hBackground2->ProjectionY();

      TH1D* hSigInt = new TH1D("hSigInt", "hSigInt", nBinsY, hSignal->GetXaxis()->GetBinLowEdge(1), hSignal->GetXaxis()->GetBinUpEdge(nBinsY));
      TH1D* hBackInt = new TH1D("hBackInt", "hBackInt", nBinsY, hBackground->GetXaxis()->GetBinLowEdge(1), hBackground->GetXaxis()->GetBinUpEdge(nBinsY));
      hSignal->Scale(1./hSignal->Integral());
      hBackground->Scale(1./hBackground->Integral()); 

      double cumulativeSignal = 1; //hSignal->Integral();
      double cumulativeBackground = 1; //hBackground->Integral();

      for(int bx=1; bx <= nBinsY; bx++){
        hSigInt->SetBinContent(bx,  cumulativeSignal);
        hBackInt->SetBinContent(bx, cumulativeBackground);

        cumulativeSignal -= hSignal->GetBinContent(bx);
        cumulativeBackground -= hBackground->GetBinContent(bx);

      }
      hSignal->Write();
      hSigInt->Write();
      hBackground->Write();
      hBackInt->Write();


      for(unsigned i=0; i < signalSelection.size(); i++){
        // fSignalEffs[kSNR][kIfFirst][i]->SetDirectory(gDirectory);
        // fSignalEffs[kEnergy][kIfFirst][i]->SetDirectory(gDirectory);
        // fSignalEffs[kSNR][kInSequence][i]->SetDirectory(gDirectory);
        // fSignalEffs[kEnergy][kInSequence][i]->SetDirectory(gDirectory);


        // fSignalEffs[kSNR][kIfFirst][i]->Write();
        // fSignalEffs[kEnergy][kIfFirst][i]->Write();
        // fSignalEffs[kSNR][kInSequence][i]->Write();
        // fSignalEffs[kEnergy][kInSequence][i]->Write();


        // delete fSignalEffs[kSNR][kIfFirst][i];
        // delete fSignalEffs[kEnergy][kIfFirst][i];
        // delete fSignalEffs[kSNR][kInSequence][i];
        // delete fSignalEffs[kEnergy][kInSequence][i];
      }

      result.Write();
  
      delete hSignal;
      delete hSigInt;
      delete hBackground;
      delete hBackInt;
  
    
    }
    catch(...){
      std::cerr << "Exception in " << __PRETTY_FUNCTION__ << " can't optimize for some reason!" << std::endl;
    }
    // Activate all the branches!
    // Branches unused by TMVA like run/eventNumber
    // are switched off by the TMVA classes, which is slightly
    // annoying when looking through the output trees.

    const int nT = 2;
    TTree* ts[nT] = {fSignalTree, fBackgroundTree};
    for(int t=0; t < nT; t++){
      TObjArray* bs = ts[t]->GetListOfBranches();
      for(int b=0; b < bs->GetEntries(); b++){
        TBranch* br = (TBranch*) bs->At(b);
        ts[t]->SetBranchStatus(br->GetName(), 1);
      }
    }


  
    fOutFile->Write();
    fOutFile->Close();
  
    std::cout << "==> wrote root file " << fOutFile->GetName() << std::endl;
    std::cout << "==> TMVAnalysis is done!" << std::endl;

#endif

  }
