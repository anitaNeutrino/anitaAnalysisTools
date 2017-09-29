#include "TGraphInteractive.h"
#include "TVirtualPad.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "RootTools.h"



Acclaim::GuiParent::~GuiParent(){
  // delete all children
  deleteChildren();
}


void Acclaim::GuiParent::deleteChildren(){
  for(unsigned i=0; i < fChildren.size(); i++){
    if(fChildren[i]){
      delete fChildren[i];
      fChildren[i] = NULL;
    }
  }
  fChildren.resize(0);
}



size_t Acclaim::GuiParent::addGuiChild(const TGraph& gr, Option_t* drawOpt){
  // create a copy, that this object owns
  TGraphInteractive* grChild = new TGraphInteractive(&gr, drawOpt);
  return addGuiChild(grChild);
}




size_t Acclaim::GuiParent::addGuiChild(TGraphInteractive* grPtr){
  if(grPtr){
    grPtr->fParent = this;
    fChildren.push_back(grPtr);
  }
  return fChildren.size();
}




size_t Acclaim::GuiParent::copyChildren(const GuiParent* that){
  for(UInt_t i=0; i < that->fChildren.size(); i++){
    if(that->fChildren[i]){
      this->addGuiChild(*that->fChildren[i], that->fChildren[i]->fDrawOpt); // use copy add function
    }
  }
  return fChildren.size();
}




void Acclaim::GuiParent::removeReference(TGraphInteractive* gr){
  UInt_t hash = gr->Hash();
  for(unsigned i=0; i < fChildren.size(); i++){
    if(fChildren[i] && fChildren[i]->Hash() == hash){
      fChildren[i] = NULL;
      break;
    }
  }
}



void Acclaim::GuiParent::DrawGroup(Option_t* opt){
  this->Draw(opt);
  for(UInt_t i=0; i < fChildren.size(); i++){
    if(fChildren[i]){
      TString thisDrawOpt = fChildren[i]->fDrawOpt + "same";
      fChildren[i]->Draw(thisDrawOpt);
    }
  }
}


/** 
 * Get a child TGraphInteractive with name matching name
 * Will return the first graph matching name.
 * The graph pointed to is still owned by guiParent, do not delete!
 * 
 * @param name is the name to match
 * 
 * @return pointer to the first graph in fChildren called name
 */
const Acclaim::TGraphInteractive* Acclaim::GuiParent::findChild(const char* name){

  const TGraphInteractive* grChild = NULL;

  for(UInt_t i=0; i < fChildren.size(); i++){
    if(fChildren[i] && strcmp(fChildren[i]->GetName(), name)==0){
      return grChild = fChildren[i];
      break;
    }
  }
  return grChild;
}
























Acclaim::TGraphInteractive::TGraphInteractive(const TGraph* gr, Option_t* drawOpt)
    : TGraphAligned(gr->GetN(), gr->GetX(), gr->GetY()),
      fParent(NULL), fDrawOpt(drawOpt)
{
  SetLineColor(gr->GetLineColor());
  SetLineStyle(gr->GetLineStyle());
  SetLineWidth(gr->GetLineWidth());
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


Acclaim::TGraphInteractive::TGraphInteractive(int n, const double* x, const double* y, Option_t* drawOpt)
    : TGraphAligned(n, x, y), fParent(NULL), fDrawOpt(drawOpt)
{
}


// Acclaim::TGraphInteractive::TGraphInteractive(const TGraphInteractive* gr)
//     : TGraphInteractive((const TGraph*)gr, gr->fDrawOpt)
// {
//   copyChildren(gr);
// }




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




void Acclaim::TGraphInteractive::DrawGroup(Option_t* opt){

  TString passedOpt(opt);
  if(passedOpt!=""){
    fDrawOpt = opt;
  }

  Draw(opt);
  
  TString drawOpt = opt;
  // drawOpt.ReplaceAll("same", "");
  // drawOpt.ReplaceAll("a", "");

  // drawOpt += "same";

  // gPad->Update();

  if(fChildren.size() > 0){
    int firstBin = GetXaxis()->GetFirst();
    int lastBin = GetXaxis()->GetLast();
    double lowerLimit = GetXaxis()->GetBinLowEdge(firstBin);
    double upperLimit = GetXaxis()->GetBinUpEdge(lastBin);

    double yMax=0, yMin=0, xMax=0, xMin=0;
    RootTools::getMaxMinWithinLimits(this, yMax, xMax, yMin, xMin, lowerLimit, upperLimit);

    for(unsigned i=0; i < fChildren.size(); i++){
      if(fChildren[i]){
        // fChildren[i]->Draw(drawOpt);

        TString drawOpt = fChildren[i]->GetDrawOpt();
        drawOpt.ReplaceAll("same", "");
        drawOpt.ReplaceAll("a", "");
        drawOpt += "same";

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
}


void Acclaim::TGraphInteractive::ExecuteEvent(int event, int px, int py){
  (void) px;
  (void) py;
  if(event==kButton1Double){
    GuiParent* p = findOriginator();

    new TCanvas();
    TGraphInteractive* grP = dynamic_cast<TGraphInteractive*>(p);
    if(grP){
      TGraphInteractive* copy = new TGraphInteractive(grP, grP->fDrawOpt);
      copy->copyChildren(grP);
      copy->SetBit(kCanDelete);      
      copy->DrawGroup();
    }
    else{
      p->DrawGroup();
    }
    
  }
  TGraphAligned::ExecuteEvent(event, px, py);
}
