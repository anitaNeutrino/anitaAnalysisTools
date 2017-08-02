#include "AnalysisPlot.h"
#include <iostream>
#include "TH2D.h"
#include "TH1D.h"
#include "TRegexp.h"

ClassImp(Acclaim::AnalysisPlot);


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
                                    int nx, double xMin, double xMax,
                                    int ny, double yMin, double yMax)
    : TNamed(name, title),
      fNx(nx), fMinX(xMin), fMaxX(xMax),
      fNy(ny), fMinY(yMin), fMaxY(yMax){

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
  cutFuncs.clear();
  hDummies.clear();
}




/** 
 * Add a cut to the AnalysisPlot
 * 
 * This function adds a whole bunch of contained histograms, enough so that all combinations of listed cuts can be drawn.
 * The total number of histograms is the product of nRetVals for all added cuts.
 * 
 * @param nRetVals is the number of possible return values, should be small i.e. 2,3,4
 * @param nameStr is appended to histogram names (along with the return value)
 * @param titleStr is appended to histogram titles (along with the return value)
 * @param AnalysisCut is a pointer to a function taking a AnitaEventSummary, and returning an int.
 * 
 * @return The number of cuts to be applied
 */
size_t Acclaim::AnalysisPlot::addCut(int nRetVals, const char* nameStr, const char* titleStr, int(*AnalysisCut)(AnitaEventSummary*)){

  // if we haven't allocated any histograms yes
  if(hs.size()==0){
    for(int retVal=0; retVal < nRetVals; retVal++){
      TString name = TString::Format("%s_%s%d", fName.Data(), nameStr, retVal);
      TString title = TString::Format("%s %s=%d", fTitle.Data(), titleStr, retVal);
        
      TH1* h = makeHist(name, title);
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
  cutFuncs.push_back(AnalysisCut);

  return cutFuncs.size();
}


int Acclaim::AnalysisPlot::Fill(AnitaEventSummary* sum, double xVal, double yVal, double zVal){

  int index = 0;
  for(UInt_t i=0; i < cutFuncs.size(); i++){

    int retVal = cutFuncs.at(i)(sum);

    if(retVal < 0 || retVal > ns.at(i)){
      std::cerr << "Error in " << __PRETTY_FUNCTION__
                << ", unexpected return value from AnalysisCut function "
                << cutNames.at(i) << ". Setting AnalysisCut return to 0."
                << std::endl;
      retVal = 0;
    }

    index += indexMultipliers.at(i)*retVal;
  }

  if(fNy==0){
    TH1D* h1 = dynamic_cast<TH1D*>(hs[index]);
    return h1->Fill(xVal, yVal);
  }
  else{
    TH2D* h2 = dynamic_cast<TH2D*>(hs[index]);
    return h2->Fill(xVal, yVal, zVal);
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
 */
void Acclaim::AnalysisPlot::Print(Option_t* opt) const{
  TNamed::Print();
  std::cout << "Cuts:" << std::endl;
  for(UInt_t i=0; i < cutNames.size(); i++){
    std::cout << cutNames.at(i) << "\t" << ns.at(i) << std::endl;
  }
  TString s(opt);
  if(s.Contains("v")){
    std::cout << "Histograms:" << std::endl;
    for(UInt_t i=0; i < hs.size(); i++){
      std::cout << hs.at(i)->GetName() << std::endl;
    }
  }
}



void Acclaim::AnalysisPlot::Draw(Option_t* opt, const TString& selection){

  TString name = TString::Format("%s_Draw%lu", fName.Data(), hDummies.size());
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
TH1* Acclaim::AnalysisPlot::makeHist(const char* name, const char* title){
  TH1* h = NULL;
  if(fNy){
    h = new TH2D(name, title, fNx, fMinX, fMaxX, fNy, fMinY, fMaxY);
  }
  else{
    h = new TH1D(name, title, fNx, fMinX, fMaxX);
  }
  return h;
}

