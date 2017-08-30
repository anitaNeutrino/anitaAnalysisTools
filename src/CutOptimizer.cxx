#include "CutOptimizer.h"
#include "OutputConvention.h"
#include "AnalysisCuts.h"
#include "SummarySet.h"
#include "ProgressBar.h"
#include "AnitaEventSummary.h"
#include "TTreeFormula.h"
#include "TApplication.h"
#include "TXMLEngine.h"
#include "TMath.h"
#include "TBranch.h"
#include "AnalysisPlot.h"
#include "RootTools.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TEfficiency.h"
#include "TDirectory.h"

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,10,0)
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodFisher.h"
#endif

ClassImp(Acclaim::CutOptimizer::FisherResult)

bool debug = false;


void Acclaim::CutOptimizer::setDebug(bool db){
  debug = db;
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
    : fSignalGlob(signalGlob), fBackgroundGlob(backgroundGlob ? backgroundGlob : ""),
      fSignalTree(NULL), fBackgroundTree(NULL), fDoAllPeaks(doAllPeaks), fSaveTrees(save_trees)
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
  bName.ReplaceAll(".", "_");  // remove remaining dots
  bName.ReplaceAll("TMath::", "");  // Remove TMath function namespace  
  //  bName.ReplaceAll("AnitaPol::", "");  // Remove AnitaPol function namespace

  bName.ReplaceAll("/", "_over_");  // wordify arithmetic/logical operators
  bName.ReplaceAll("+", "_plus_");  // wordify arithmetic/logical operators
  bName.ReplaceAll("-", "_minus_"); // wordify arithmetic/logical operators
  bName.ReplaceAll("*", "_times_"); // wordify arithmetic/logical operators
  bName.ReplaceAll(">=", "_ge_"); // wordify arithmetic/logical operators 
  bName.ReplaceAll("<=", "_le_"); // wordify arithmetic/logical operators 
  bName.ReplaceAll(">", "_gt_"); // wordify arithmetic/logical operators 
  bName.ReplaceAll("<", "_lt_"); // wordify arithmetic/logical operators 
  bName.ReplaceAll("__", "_");  // Remove double underscore if there is one
  bName.Remove(TString::kTrailing, '_'); // Remove trailing underscore if there is one
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



/** 
 * This function does the hard work of generating a signal tree and a background tree for the TMVA to work on
 * It creates two new TTrees (signal+background) and generates a branch for each of the formula strings passed.
 * Then it loops through the AnitaEventSummaries, applying the selection cuts.
 * For each event passing the signal/background cuts it adds an entry in the tree by evaluating each of the formulae.
 * 
 * @param signalSelection is a set of my custom cut class defining the signal selection for TMVA training.
 * @param backgroundSelection is a set of my custom cut class defining the background selection for TMVA training.
 * @param formulaStrings are a set of formulae like one would pass to TTree::Drawm which get evaluated for signal/background events
 */
void Acclaim::CutOptimizer::generateSignalAndBackgroundTrees(const std::vector<const Acclaim::AnalysisCut*>& signalSelection, const std::vector<const Acclaim::AnalysisCut*>& backgroundSelection, const std::vector<FormulaString>& formulaStrings){

  // First generate the new signal and background trees
  if(fSignalTree){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", pre-existing signal tree. Deleting and regenerating." << std::endl;
    delete fSignalTree;
    fSignalTree = NULL;
  }
  fSignalTree = new TTree("signalTree",  "signalTree");
  if(!fSaveTrees){
    fSignalTree->SetDirectory(0);
  }

  if(fBackgroundTree){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", pre-existing background tree. Deleting and regenerating." << std::endl;
    delete fBackgroundTree;
    fBackgroundTree = NULL;
  }
  fBackgroundTree = new TTree("backgroundTree",  "backgroundTree");
  if(!fSaveTrees){
    fBackgroundTree->SetDirectory(0);
  }
  
  // now we need to figure out if we have separate globs for signal and background
  // or just one for both...
  std::vector<TString> globs;
  globs.push_back(fSignalGlob);
  if(fBackgroundGlob != ""){
    globs.push_back(fBackgroundGlob);
  }

  // Branches could either be integers or floats, depending on the formula return type
  // We will try to figure it out later from the TTreeFormula, but for now assign
  // arrays of both ints and floats
  const int nFormStr = formulaStrings.size(); // The actual number of formula could be less in the case of a bad expression
  fSignalFloatVals.resize(nFormStr, 0);
  fBackgroundFloatVals.resize(nFormStr, 0);
  fSignalIntVals.resize(nFormStr, 0);
  fBackgroundIntVals.resize(nFormStr, 0);

  // and store the types we do figure out here...
  std::vector<BranchType> assignedSignalBranches(nFormStr, kUnassigned);
  std::vector<BranchType> assignedBackgroundBranches(nFormStr, kUnassigned);

  // book efficiencies
  fSignalEffs[kSNR].reserve(signalSelection.size());
  fSignalEffs[kEnergy].reserve(signalSelection.size());
  bool defSumw2 = TH1::GetDefaultSumw2();
  TH1::SetDefaultSumw2(true); // force tefficiency to have sumw2?
  for(unsigned i=0; i < signalSelection.size(); i++){
    TString cutName = signalSelection[i]->getName();
    TString snrName = "eff_" + cutName + "_vs_SNR";
    fSignalEffs[kSNR].push_back(new TEfficiency(snrName, snrName, 200, 0, 100));

    TString energyName = "eff_" + cutName + "_vs_Energy";
    fSignalEffs[kEnergy].push_back(new TEfficiency(energyName, energyName, 100, 15, 25));
  }
  TH1::SetDefaultSumw2(defSumw2);

  for(UInt_t g=0; g < globs.size(); g++){

    Bool_t doSignal = g == 0 ? true : false;
    Bool_t doBackground = globs.size() == 2 && g == 0 ? false : true;

    // Then load the master AnitaEventSummary chain using my SummarySet class    
    SummarySet ss(globs[g]);
    Long64_t nEntries = ss.N(); // > 10000 ? 10000 : ss.N();

    // Use my custom FolderHolder to handle notification subtleties
    FormulaHolder forms(ss.getChain());
    for(int i=0; i < nFormStr; i++){
      forms.add(formulaStrings.at(i).first);
    }

    // Now we loop over the AnitaEventSummaries
    ProgressBar p(nEntries); // For prettiness  
    for(Long64_t entry=0; entry < nEntries; entry++){
      ss.getEntry(entry);

      AnitaEventSummary* sum = ss.summary();

      // Is this a background event?
      // (Do this first since we probably have more background)
      if(doBackground){
        Bool_t matchBackgroundSelection = true;
        for(UInt_t i=0; i < backgroundSelection.size(); i++){
          int thisCutVal = backgroundSelection.at(i)->apply(sum);
          matchBackgroundSelection = matchBackgroundSelection && thisCutVal!=0;
          if(!matchBackgroundSelection){
            break;
          }
        }

        if(matchBackgroundSelection){
          for(UInt_t i=0; i < forms.N(); i++){

            if(assignedBackgroundBranches.at(i)==kUnassigned){
              assignedBackgroundBranches.at(i) = setBranchFromFormula(fBackgroundTree, forms.at(i), forms.str(i), &fBackgroundIntVals.at(i), &fBackgroundFloatVals.at(i));
            }
            if(assignedBackgroundBranches.at(i)==kInt){
              fBackgroundIntVals.at(i) = forms.at(i)->EvalInstance();
            }
            else{
              fBackgroundFloatVals.at(i) = forms.at(i)->EvalInstance();
            }
          }
          fBackgroundTree->Fill();
        }
      }
      else if(doSignal){
        // If it's not background, it might be signal...
        Bool_t matchSignalSelection = true;
        for(UInt_t i=0; i < signalSelection.size(); i++){
          int thisCutVal = signalSelection.at(i)->apply(sum);

          // std::cerr << signalSelection.at(i)->getName() << "\t" << thisCutVal << std::endl;
          
          fSignalEffs[kSNR][i]->Fill(thisCutVal, sum->mcDeconvolved().snr, sum->weight());
          fSignalEffs[kEnergy][i]->Fill(thisCutVal, TMath::Log10(sum->mc.energy), sum->weight());

          matchSignalSelection = matchSignalSelection && (thisCutVal != 0);
        }

        if(matchSignalSelection){
          for(UInt_t i=0; i < forms.N(); i++){
            if(assignedSignalBranches.at(i)==kUnassigned){
              assignedSignalBranches.at(i) = setBranchFromFormula(fSignalTree, forms.at(i), forms.str(i), &fSignalIntVals.at(i), &fSignalFloatVals.at(i));
            }
            if(assignedSignalBranches.at(i)==kInt){
              fSignalIntVals.at(i) = forms.at(i)->EvalInstance();
            }
            else{
              fSignalFloatVals.at(i) = forms.at(i)->EvalInstance();
            }
          }
          fSignalTree->Fill();
        }
      }
      else{
        std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", not doing signal or background?" << std::endl;
      }
      p.inc(entry, nEntries);
    }
    if(doSignal){
      std::cout << "Created " << fSignalTree->GetName() << " with " << fSignalTree->GetEntries() << " events" << std::endl;
      if(debug){fSignalTree->Print();}
    }
    if(doBackground){
      std::cout << "Created " << fBackgroundTree->GetName() << " with " << fBackgroundTree->GetEntries() << " events" << std::endl;
      if(debug){fBackgroundTree->Print();}
    }
  }
}


void Acclaim::CutOptimizer::optimize(const std::vector<const Acclaim::AnalysisCut*>& signalSelection, const std::vector<const Acclaim::AnalysisCut*>& backgroundSelection, const std::vector<FormulaString>& formulaStrings, const char* fileName){

  // Someone with a bit more time can do the backwards compatibility
#if ROOT_VERSION_CODE < ROOT_VERSION(6,10,0)
  std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", TVMA changed a lot in recent ROOT versions. ";
  std::cerr << "This class requires ROOT version at least 6.10, you only have " << ROOT_VERSION_CODE << std::endl;
#else
  
  TFile* fOutFile = makeOutputFile(fileName);

  generateSignalAndBackgroundTrees(signalSelection, backgroundSelection, formulaStrings);  

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
      fSignalEffs[kSNR][i]->SetDirectory(gDirectory);
      fSignalEffs[kEnergy][i]->SetDirectory(gDirectory);
      
      fSignalEffs[kSNR][i]->Write();
      fSignalEffs[kEnergy][i]->Write();

      delete fSignalEffs[kSNR][i];
      delete fSignalEffs[kEnergy][i];
      
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
