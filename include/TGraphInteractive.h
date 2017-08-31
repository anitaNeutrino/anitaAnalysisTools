#ifndef TGRAPH_INTERACTIVE
#define TGRAPH_INTERACTIVE

#include "TGraphAligned.h"

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

class TGraphInteractive : public TGraphAligned {

 public:
  TGraphInteractive() {;}
  TGraphInteractive(const TGraph* gr);
  TGraphInteractive(const TGraphInteractive* gr);
  virtual ~TGraphInteractive();
  void ExecuteEvent(int event, int px, int py);
  
  size_t add(const TGraph*);
  TGraphInteractive* getParent();
  TGraphInteractive* findOriginator();
  void Draw(Option_t* opt);
  
 protected:
  std::vector<TGraphInteractive*> fChildren;
  TGraphInteractive* fParent;

  void removeReference(TGraphInteractive* gr);
  
};

}
#endif
