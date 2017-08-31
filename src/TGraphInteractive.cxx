#include "TGraphInteractive.h"
#include "TVirtualPad.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "RootTools.h"

Acclaim::TGraphInteractive::TGraphInteractive(const TGraph* gr)
    : TGraphAligned(gr->GetN(), gr->GetX(), gr->GetY()),
      fParent(NULL)
{
  SetLineColor(gr->GetLineColor());
  SetLineStyle(gr->GetLineStyle());
  SetFillColor(0);
  SetFillStyle(gr->GetFillStyle());
  SetName(gr->GetName());
  
  GetXaxis()->SetTitle(gr->GetXaxis()->GetTitle());
  GetYaxis()->SetTitle(gr->GetYaxis()->GetTitle());
  TString title = gr->GetTitle();
  if(title.Contains(";")){
    Ssiz_t p = title.First(';');
    title = title(0, p);
  }
  SetTitle(title);


  GetXaxis()->SetRange(gr->GetXaxis()->GetFirst(), gr->GetXaxis()->GetLast());
  
}



Acclaim::TGraphInteractive::TGraphInteractive(const TGraphInteractive* gr)
    : TGraphInteractive((const TGraph*)gr)
{
  for(UInt_t i=0; i < gr->fChildren.size(); i++){
    if(gr->fChildren[i]!=NULL){
      add(gr->fChildren[i]);
    }
  }
}



Acclaim::TGraphInteractive::~TGraphInteractive(){

  // Remove from parent's reference
  if(fParent){
    fParent->removeReference(this);
  }

  // delete all children
  for(unsigned i=0; i < fChildren.size(); i++){
    if(fChildren[i]!=NULL){
      delete fChildren[i];
      // I think the children should then have zero'd the parent pointer
      // in their own destructor?
    }
  }
  
}




size_t Acclaim::TGraphInteractive::add(const TGraph* gr){

  // create a copy, that this object owns
  TGraphInteractive* grChild = new TGraphInteractive(gr);
 // I am your father
  grChild->fParent = this;
  fChildren.push_back(grChild);  

  return fChildren.size();
}


void Acclaim::TGraphInteractive::removeReference(TGraphInteractive* gr){
  UInt_t id = gr->GetUniqueID();
  for(unsigned i=0; i < fChildren.size(); i++){
    if(fChildren[i]->GetUniqueID() == id){
      fChildren[i] = NULL;
      break;
    }
  }
}





Acclaim::TGraphInteractive* Acclaim::TGraphInteractive::getParent(){
  return fParent;
}



/** 
 * Go up the parent tree until we find a TGraphInteractive with fParent = NULL
 * 
 * @return the top level TGraphInteractive
 */
Acclaim::TGraphInteractive* Acclaim::TGraphInteractive::findOriginator(){
  
  TGraphInteractive* current = this;
  TGraphInteractive* next = this->getParent();

  while(next!=NULL){
    current = next;
    next = current->getParent();
  }
  return current;  
}


void Acclaim::TGraphInteractive::Draw(Option_t* opt){

  TGraphAligned::Draw(opt);
  
  TString drawOpt = opt;
  drawOpt.ReplaceAll("same", "");
  drawOpt.ReplaceAll("a", "");

  drawOpt += "same";

  // gPad->Update();
  int firstBin = GetXaxis()->GetFirst();
  int lastBin = GetXaxis()->GetLast();
  double lowerLimit = GetXaxis()->GetBinLowEdge(firstBin);
  double upperLimit = GetXaxis()->GetBinUpEdge(lastBin);

  double yMax=0, yMin=0, xMax=0, xMin=0;
  RootTools::getMaxMinWithinLimits(this, yMax, xMax, yMin, xMin, lowerLimit, upperLimit);

  for(unsigned i=0; i < fChildren.size(); i++){
    if(fChildren[i]){
      fChildren[i]->Draw(drawOpt);

      double childMax, childMin;
      RootTools::getMaxMinWithinLimits(fChildren[i], childMax, xMax, childMin, xMin, lowerLimit, upperLimit);
      // std::cerr << lowerLimit << "\t" << upperLimit << "\t" << firstBin << "\t" << lastBin << "\t" << childMax << "\t" << childMin << std::endl;

      if(childMax > yMax){
        yMax = childMax;
      }
      if(childMin < yMin){
        yMin = childMin;
      }
    }
  }

  double padding = 0.1*(yMax - yMin);
  SetMaximum(yMax + padding);
  SetMinimum(yMin - padding);
}


void Acclaim::TGraphInteractive::ExecuteEvent(int event, int x, int y){
  
  if(event == kButton1Double){

    TGraphInteractive* originator = findOriginator();
    TString drawOpt = originator->GetDrawOption();

    // force options for drawing a new axis
    drawOpt.ReplaceAll("same", "");
    if(!drawOpt.Contains("a")){
      drawOpt += "a"; 
    }

    TGraphInteractive* grNew = new TGraphInteractive(originator);
    grNew->SetBit(kCanDelete); // now the canvas will own it
    
    TCanvas* c1 = new TCanvas();
    grNew->Draw(drawOpt);

    TLegend* l1 = c1->BuildLegend();
    l1->SetBit(kCanDelete);    
  }  

  // pass on up the chain
  TGraphAligned::ExecuteEvent(event, x, y);  
}
