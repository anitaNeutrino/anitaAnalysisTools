#ifndef TGRAPH_INTERACTIVE
#define TGRAPH_INTERACTIVE

#include "TGraphAligned.h"
#include "TCanvas.h"

namespace Acclaim{


/** 
 * @class A minimalistic extension to TGraphAligned for some GUI bells and whistles
 * 
 * I'm thinking of using this primarily with the AnalysisReco class
 * TODO: 
 *      Try to incorporate some of the features of the now redundant TGraphFB class,
 *      inside the now redundant FourierBuffer class.
 * 
 * There is a parent/child model to bundle graphs together
 * Each TGraphInteractive can have only one parent.
 * However, they can have any number of children.
 * (There's a crisis of family values in the TGraphInteractive community)
 * 
 */

class TGraphInteractive;

/** 
 * Skeletal class to inherit from if you want to have a bunch of TGraphs which follow you around a ROOT GUI
 * Draw the whole set (self + children) with DrawGroup(Option_t*)
 */
class GuiParent {
 public: 

  // Must be defined so derived class can DrawGroup (just make it reference base class draw)
  // e.g. virtual void Draw(Option_t* opt){TObject::Draw(opt);}
  virtual void Draw(Option_t* opt) = 0;
  
  // Must be overloaded by children (TGraphInteractive*) to return pointer to parent
  virtual GuiParent* getParent() const {return NULL;}



  void DrawGroup(Option_t* opt=""); // Draw self and loop over/draw TGraphInteractive children
  size_t addGuiChild(TGraphInteractive* grPtr); // Does not copy, takes ownership of heap object
  size_t addGuiChild(const TGraph& grRef,  Option_t* drawOpt); // Copies and owns the copy
  size_t copyChildren(const GuiParent* that); // copy all children from that, and add to this  

  const TGraphInteractive* findChild(const char* name);

  GuiParent(){;}
  virtual ~GuiParent(); // Delete children on destruction in here
  void deleteChildren();
  
 private:
  
  friend class TGraphInteractive;
  void removeReference(TGraphInteractive* grChild);  
  std::vector<TGraphInteractive*> fChildren;
};



class TGraphInteractive : public TGraphAligned, public GuiParent {
 public:
  TGraphInteractive() {;}
  TGraphInteractive(int n, const double* x, const double* y, Option_t* drawOpt = "");
  TGraphInteractive(const TGraph* gr, Option_t* drawOpt);
  virtual ~TGraphInteractive();

  // Satisfy pure virtual overload of GuiParent
  virtual void Draw(Option_t* opt = ""){
    TGraphAligned::Draw(opt);
  }

  virtual void ExecuteEvent(int event, int px, int py);

  // Walk up parent tree until a NULL is returned
  GuiParent* findOriginator() const;


  virtual void DrawGroup(Option_t* opt = "");

 public:
  friend class GuiParent;
  virtual GuiParent* getParent() const {return fParent;}
 private:
  GuiParent* fParent;
  TString fDrawOpt; // keep it somewhere
  
};

}
#endif
