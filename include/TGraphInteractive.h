#ifndef TGRAPH_INTERACTIVE
#define TGRAPH_INTERACTIVE

#include "TGraphAligned.h"
#include "TCanvas.h"


namespace Acclaim{

  class TGraphInteractive;
  
  /**
   * @class GuiParent
   * @brief Inherit from this to draw interactive TGraphs on top of you
   * 
   * Skeletal class to inherit from if you want to have a bunch of TGraphs which follow you around a ROOT GUI
   * Draw the whole set (self + children) with DrawGroup(Option_t*)
   */
  
  class GuiParent
  {
  public: 


    /** 
     * Default constructor
     */
    GuiParent(){;}

    virtual ~GuiParent();

    void deleteChildren();
    
    /** 
     * Must be defined so derived class can DrawGroup (just make it reference base class draw)
     * 
     * For example, the following would suffice: 
     * virtual void Draw(Option_t* opt){TGraph::Draw(opt);}
     * 
     * @param opt is the draw option to use
     */
    virtual void Draw(Option_t* opt) = 0;
    virtual GuiParent* getParent() const {return NULL;} /// *Must* be overloaded by children (TGraphInteractive*) to return pointer to parent

    virtual void DrawGroup(Option_t* opt=""); 
    size_t addGuiChild(TGraphInteractive* grPtr);
    size_t addGuiChild(const TGraph& grRef,  Option_t* drawOpt);
    size_t copyChildren(const GuiParent* that);
    const TGraphInteractive* findChild(const char* name);
  
  private:
  
    friend class TGraphInteractive;
    void removeReference(TGraphInteractive* grChild);  
    std::vector<TGraphInteractive*> fChildren;
  };




  /**
   * @class TGraphInteractive
   * @brief A minimalistic extension to TGraphAligned for some GUI bells and whistles
   * 
   * I'm thinking of using this primarily with the AnalysisReco class
   * TODO: 
   *      Try to incorporate some of the features of the now redundant TGraphFB class,
   *      inside the now redundant FourierBuffer class.
   * 
   * There is a parent/child model to bundle graphs together.
   * Each TGraphInteractive can have only one parent.
   * However, they can have any number of children.
   * (There's a crisis of family values in TGraphInteractive land)
   */

  class TGraphInteractive : public TGraphAligned, public GuiParent {
  public:
    TGraphInteractive() {;} /// Default constructor
    TGraphInteractive(int n, const double* x, const double* y, Option_t* drawOpt = "");
    TGraphInteractive(const TGraph* gr, Option_t* drawOpt);
    virtual ~TGraphInteractive();

    virtual void ExecuteEvent(int event, int px, int py);

    GuiParent* findOriginator() const;
    virtual void DrawGroup(Option_t* opt = "");

    
    /** 
     * Call TGraphAligned draw Satisfy pure virtual overload of GuiParent
     * 
     * @param opt is the ROOT drawing option to use
     */
    virtual void Draw(Option_t* opt = ""){ 
      TGraphAligned::Draw(opt);
    }

    /** 
     * Set the default draw option
     * This value is used if Draw() is called with default arg
     * 
     * @param drawOpt is the default draw option to set
     */
    void SetDrawOpt(Option_t* drawOpt){
      fDrawOpt = drawOpt;
    }

    /** 
     * Retrieve default draw option
     * 
     * @return the default draw option(fDrawOpt)
     */
    Option_t* GetDrawOpt(){
      return fDrawOpt;
    }

  public:
    friend class GuiParent;
    virtual GuiParent* getParent() const {return fParent;}	/// Returns pointer to parent
  private:
    GuiParent* fParent;						/// Pointer to parent
    TString fDrawOpt;						/// Internal copy of ROOT draw option, since ROOT stores it strangely 
  
  };

}
#endif
