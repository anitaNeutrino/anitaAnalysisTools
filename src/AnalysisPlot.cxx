#include "AnalysisPlot.h"
#include <iostream>
#include "TH2D.h"
#include "TH1D.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TRegexp.h"

ClassImp(Acclaim::AnalysisPlot);
ClassImp(Acclaim::AnalysisProf);


/** 
 * Useful contructor
 * 
 * @param name is the name of the object, all contained histograms will have derived names
 * @param title is the title of the object, all contained histograms will have derived titles
 * @param nx is the number of x-axis bins
 * @param xMin is lower x-axis limit
 * @param xMax is the higher x-axis limit
 * @param ny is the number of y-axis bins (if 0 then contained histograms will be of type TH1D, otherwise TH2D)
 * @param yMin is the lower y-axis limit
 * @param yMax is the higher y-axis limit
 */
Acclaim::AnalysisPlot::AnalysisPlot(const char* name, const char* title,
                                    int nBinsX, double xMin, double xMax,
                                    int nBinsY, double yMin, double yMax)
    : TNamed(name, title),
      fNx(nBinsX), fMinX(xMin), fMaxX(xMax),
      fNy(nBinsY), fMinY(yMin), fMaxY(yMax){

}



/** 
 * Destructor, deletes all histograms.
 */
Acclaim::AnalysisPlot::~AnalysisPlot(){

  for(UInt_t i=0; i < hs.size(); i++){
    delete hs[i];
    hs[i] = NULL;
  }

  for(UInt_t i=0; i < hDummies.size(); i++){
    delete hDummies[i];
    hDummies[i] = NULL;
  }
  
  hs.clear();
  ns.clear();
  indexMultipliers.clear();
  cutNames.clear();
  analysisCuts.clear();
  hDummies.clear();
}




/** 
 * Add a cut to the AnalysisPlot
 * 
 * This function adds a whole bunch of contained histograms, enough so that all combinations of listed cuts can be drawn.
 * The total number of histograms is the product of nRetVals for all added cuts.
 * See the AnalysisCuts namespace for example cuts.
 * 
 * @param AnalysisCut is a pointer to a function that actually implements the cut, returns an integer encoding the result
 * @param nameStr is appended to histogram names (along with the return value)
 * @param titleStr is appended to histogram titles (along with the return value)
 * 
 * @return The number of cuts stored the AnalysisPlot object will proccess on a Fill() command 
 */
size_t Acclaim::AnalysisPlot::addCut(int(*AnalysisCut)(const AnitaEventSummary*), const char* nameStr, const char* titleStr){

  const int nRetVals = AnalysisCut(nullptr);

  // if we haven't allocated any histograms yes
  if(hs.size()==0){
    for(int retVal=0; retVal < nRetVals; retVal++){
      TString name = TString::Format("%s_%s%d", fName.Data(), nameStr, retVal);
      TString title = TString::Format("%s %s=%d", fTitle.Data(), titleStr, retVal);
        
      TH1* h = makeHist(name, title);
      h->SetDirectory(0);

      hs.push_back(h);
    }
  }
  else{

    const UInt_t oldNumHists = hs.size();
    const UInt_t newNumHists = hs.size()*nRetVals;

    while(hs.size() < newNumHists){
      hs.push_back(nullptr);
    }
    
    for(UInt_t i=oldNumHists; i < newNumHists; i++){
      int j = i % oldNumHists;
      hs[i] = (TH1*) hs[j]->Clone();
      hs[i]->SetDirectory(0);
    }

    for(UInt_t i=0; i < newNumHists; i++){
      
      int retVal = i/oldNumHists;

      TString oldName = hs[i]->GetName();
      TString newName = oldName + TString::Format("_%s%d", nameStr, retVal);
      hs[i]->SetName(newName);

      TString oldTitle = hs[i]->GetTitle();
      TString newTitle = oldTitle + TString::Format(" %s=%d", titleStr, retVal);
      hs[i]->SetTitle(newTitle);
    }
  }
  ns.push_back(nRetVals);
  indexMultipliers.push_back(ns.size() == 1 ? 1 : ns.at(ns.size()-1));  
  cutNames.push_back(nameStr);
  analysisCuts.push_back(AnalysisCut);

  return analysisCuts.size();
}




/** 
 * Execute all the cut functions and get the corresponding histogram index
 * 
 * @param sum is the AnitaEventSummary to process
 * 
 * @return the histogram index
 */
int Acclaim::AnalysisPlot::getIndexFromCuts(const AnitaEventSummary* sum){
  int index = 0;
  for(UInt_t i=0; i < analysisCuts.size(); i++){

    int retVal = analysisCuts.at(i)(sum);

    if(retVal < 0 || retVal >= ns.at(i)){
      std::cerr << "Error in " << __PRETTY_FUNCTION__
                << ", unexpected return value from AnalysisCut function "
                << cutNames.at(i) << ". Setting AnalysisCut return to 0."
                << std::endl;
      retVal = 0;
    }

    index += indexMultipliers.at(i)*retVal;
  }
  return index;
}


/** 
 * Apply all the cut functions added with addCut(), and fill the appropriate histogram
 * 
 * @param sum the AnitaEventSummary associated with the event
 * @param x is the value to fill along the x-axis
 * @param y is an event weight in the TH1D case, or is the value to fill along the y-axis in the TH2D case
 * @param z does nothing in the TH1D case, or is an event weight in the TH2D case
 * 
 * @return the index of the bin filled, as returned by TH1D::Fill or TH2D::Fill
 */
int Acclaim::AnalysisPlot::Fill(const AnitaEventSummary* sum, double x, double y, double z){

  int index = getIndexFromCuts(sum);
  
  if(fNy==0){
    TH1D* h1 = dynamic_cast<TH1D*>(hs[index]);
    return h1->Fill(x, y);
  }
  else{
    TH2D* h2 = dynamic_cast<TH2D*>(hs[index]);
    return h2->Fill(x, y, z);
  }
}




/** 
 * Default draw function, just draws all the histograms
 * 
 * @param opt is the draw option passed to the histogram
 */
void Acclaim::AnalysisPlot::Draw(Option_t* opt){
  // will assume you just want to plot everything
  TString selection = "*";
  Draw(opt, selection);
}



/** 
 * Print out some useful information to the terminal
 * 
 * @param opt is an optional string, if it contains "v" prints all the histogram names
 */
void Acclaim::AnalysisPlot::Print(Option_t* opt) const{
  TNamed::Print();
  std::cout << fName << " Cuts:" << std::endl;
  for(UInt_t i=0; i < cutNames.size(); i++){
    std::cout << cutNames.at(i) << "\t" << ns.at(i) << std::endl;
  }
  TString s(opt);
  if(s.Contains("v")){
    std::cout << fName << " Histograms:" << std::endl;
    for(UInt_t i=0; i < hs.size(); i++){
      std::cout << hs.at(i)->GetName() << std::endl;
    }
  }
}




/** 
 * Draw the information contained inside the histograms
 * 
 * @param opt is the ROOT draw option
 * @param selection is a regexp to select the names of the histogram
 */
void Acclaim::AnalysisPlot::Draw(Option_t* opt, const TString& selection){

  TString name = TString::Format("%s_Draw%lu", fName.Data(), hDummies.size());
  TString title = selection == "*" ? fTitle : selection;
  TH1* h = makeHist(name, selection);
  hDummies.push_back(h);
  
  TString boundedSelection = "*" + selection + "*";
  TRegexp reggie(boundedSelection.Data(), true);
  
  for(UInt_t i=0; i < hs.size(); i++){
    TString histName = hs[i]->GetName();
    Int_t index = histName.Index(reggie);
    // std::cout << selection << "\t" << hs[i]->GetName() << "\t" << i << "\t" << index << std::endl;
    if(index!=-1){
      // std::cout << "here!" << std::endl;
      h->Add(hs[i]);
    }      
  }  
  h->Draw(opt);
}




/** 
 * Utility function to construct another contained histogram 
 * 
 * @param name 
 * @param title 
 * 
 * @return the newly constructed histogram cast as a pointer to a TH1
 */
TH1* Acclaim::AnalysisPlot::makeHist(const char* name, const char* title) const{
  TH1* h = NULL;
  if(fNy){
    h = new TH2D(name, title, fNx, fMinX, fMaxX, fNy, fMinY, fMaxY);
  }
  else{
    h = new TH1D(name, title, fNx, fMinX, fMaxX);
  }
  return h;
}






/** 
 * Useful contructor for the derived AnalysisProf class
 * 
 * @param name is the name of the object, all contained profiles will have derived names
 * @param title is the title of the object, all contained profiles will have derived titles
 * @param nx is the number of x-axis bins
 * @param xMin is lower x-axis limit
 * @param xMax is the higher x-axis limit
 * @param ny is the number of y-axis bins (if 0 then contained profiles will be of type TProfile, otherwise TProfile2D)
 * @param yMin is the lower y-axis limit
 * @param yMax is the higher y-axis limit
 */
Acclaim::AnalysisProf::AnalysisProf(const char* name, const char* title,
                                    int nBinsX, double xMin, double xMax,
                                    int nBinsY, double yMin, double yMax)
    : AnalysisPlot(name, title, nBinsX, xMin, xMax, nBinsY, yMin, yMax){
  
}




/** 
 * Utility function to construct another contained profile 
 * 
 * @param name 
 * @param title 
 * 
 * @return the newly constructed profile cast as a pointer to a TH1
 */
TH1* Acclaim::AnalysisProf::makeHist(const char* name, const char* title) const{
  TH1* h = NULL;
  if(fNy){
    h = new TProfile2D(name, title, fNx, fMinX, fMaxX, fNy, fMinY, fMaxY);
  }
  else{
    h = new TProfile(name, title, fNx, fMinX, fMaxX);
  }
  return h;
}




/** 
 * Apply all the cut functions added with addCut(), and fill the appropriate profile
 * 
 * @param sum the AnitaEventSummary associated with the event
 * @param x is the value to fill along the x-axis
 * @param y is the value to profile in the TProfile case, the y-axis value in the TProfile2D case
 * @param z does nothing in the TProfile case, or is the value to profile in the TProfile2D case
 * 
 * @return the index of the bin filled, as returned by TProfile::Fill or TProfile2D::Fill
 */
int Acclaim::AnalysisProf::Fill(const AnitaEventSummary* sum, double x, double y, double z){

  int index = getIndexFromCuts(sum);
  
  if(fNy==0){
    TProfile* h1 = dynamic_cast<TProfile*>(hs[index]);
    return h1->Fill(x, y);
  }
  else{
    TProfile2D* h2 = dynamic_cast<TProfile2D*>(hs[index]);
    return h2->Fill(x, y, z);
  }
}
















