#include "TGraphInteractive.h"
#include "TVirtualPad.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "RootTools.h"



Acclaim::GuiParent::~GuiParent(){
  // delete all children
  for(unsigned i=0; i < fChildren.size(); i++){
    if(fChildren[i]!=NULL){
      delete fChildren[i];
      // I think the children should then have zero'd the parent pointer
      // in their own destructor?
    }
  }
}


size_t Acclaim::GuiParent::addGuiChild(const TGraph& gr, Option_t* drawOpt){
  // create a copy, that this object owns
  TGraphInteractive* grChild = new TGraphInteractive(&gr, drawOpt);
  return addGuiChild(grChild);
}


size_t Acclaim::GuiParent::addGuiChild(TGraphInteractive* grPtr){
  grPtr->fParent = this;
  fChildren.push_back(grPtr);  

  return fChildren.size();
}


size_t Acclaim::GuiParent::copyChildren(const GuiParent* that){
  for(UInt_t i=0; i < that->fChildren.size(); i++){
    this->addGuiChild(*that->fChildren[i], that->fChildren[i]->fDrawOpt); // use copy add function
  }
  return fChildren.size();
}



void Acclaim::GuiParent::removeReference(TGraphInteractive* gr){
  UInt_t hash = gr->Hash();
  for(unsigned i=0; i < fChildren.size(); i++){
    if(fChildren[i]->Hash() == hash){
      fChildren[i] = NULL;
      break;
    }
  }
}


void Acclaim::GuiParent::DrawGroup(Option_t* opt){
  Draw(opt);
  for(UInt_t i=0; i < fChildren.size(); i++){
    fChildren[i]->Draw(); // draw with saved draw option
  }
}


void Acclaim::GuiParent::ExecuteEvent(int event, int x, int y){
  (void) x;
  (void) y;
  if(event == kButton1Double){

    TCanvas* c1 = new TCanvas();
    DrawGroup("colz");

    TLegend* l1 = c1->BuildLegend();
    l1->SetBit(kCanDelete);    
  }  
}























Acclaim::TGraphInteractive::TGraphInteractive(const TGraph* gr, Option_t* drawOpt)
    : TGraphAligned(gr->GetN(), gr->GetX(), gr->GetY()),
      fParent(NULL), fDrawOpt(drawOpt)
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




Acclaim::TGraphInteractive::~TGraphInteractive(){

  // Remove from parent's reference
  if(fParent){
    fParent->removeReference(this);
  }
}




/** 
 * Go up the parent tree until we find a GuiParent where getParent() returns NULL
 * 
 * @return the top level Graphical object
 */
Acclaim::GuiParent* Acclaim::TGraphInteractive::findOriginator() const{
  
  const GuiParent* current = static_cast<const GuiParent*>(this);
  const GuiParent* next = this->getParent();

  while(next!=NULL){
    current = next;
    next = current->getParent();
  }
  return const_cast<GuiParent*>(current);
}




void Acclaim::TGraphInteractive::Draw(Option_t* opt){

  TString passedOpt(opt);
  if(passedOpt!=""){
    fDrawOpt = opt;
  }

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


