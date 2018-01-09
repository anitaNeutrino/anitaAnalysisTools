/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             Read in trees of ANITA event summaries and try to separate signal and background...
***********************************************************************************************************/

#ifndef CUT_OPTIMIZER_H
#define CUT_OPTIMIZER_H


#include "TString.h"
#include "TChain.h"
#include "TXMLEngine.h"


class TH2D;
class TEfficiency;

namespace Acclaim{

  class SummarySet;

  /** 
   * @class CutOptimizer
   * @brief A class to parse AnitaEventSummary trees, extract quantities of interest,
   *        feed them into a TMVA framework and separate them
   */
  class CutOptimizer {

  public:

    // typedef std::pair<const char*, bool> FormulaString;

    static void setDebug(bool db);
    static TString branchifyName(const char* formStr);
  
    CutOptimizer(const char* signalGlob, const char* backgroundGlob = NULL);
    virtual ~CutOptimizer();
    void optimize(const std::vector<const TCut*>& signalSelection,
		  const std::vector<const TCut*>& backgroundSelection,
		  const std::vector<TString>& formulaStrings,
		  const char* outFileName = "");


    // /** 
    //  * Add to a list of extra tree for which you wish to calculate the MVA score
    //  * 
    //  * @param treeName is the name of the tree to create
    //  * @param glob is the wildcard filenames for the ROOT files containing AnitaEventSummary trees (named sunTree)
    //  * @param spectatorSelection is a vector of TCut variables used to choose the spectators 
    //  * 
    //  * @return the number of spectator trees to be prepared
    //  */
    // size_t addSpectatorTree(const char* treeName, const char* glob, const std::vector<const TCut*>& spectatorSelection){
    //   fSpecTreeNames.push_back(treeName);
    //   fSpecGlobs.push_back(glob);
    //   fSpecSelections.push_back(spectatorSelection);
    //   return fSpecSelections.size();
    // }

  protected:

    enum BranchType{
      kUnassigned,
      kInt,
      kFloat
    };

    TFile* makeOutputFile(const char* outFileName);
    // void generateSignalAndBackgroundTreesProof(const std::vector<const TCut*>& signalSelection,
    // 					       const std::vector<const TCut*>& backgroundSelection,
    // 					       const std::vector<FormulaString>& treeVars);
    // BranchType setBranchFromFormula(TTree* t, const TTreeFormula* f, const char* formulaString, Int_t* intPtr, Float_t* floatPtr);

    TString fSignalGlob;
    TString fBackgroundGlob;
    TString fOutFileName;
    TFile* fOutFile;
    TChain* fSignalTree;
    TChain* fBackgroundTree;
    // TTree* fRejectedSignalTree;
    // TTree* fRejectedBackgroundTree;
    // Bool_t fDoAllPeaks;
    // Bool_t fSaveTrees;
    // std::vector<Float_t> fSignalFloatVals;
    // std::vector<Float_t> fBackgroundFloatVals;
    // std::vector<Int_t> fSignalIntVals;
    // std::vector<Int_t> fBackgroundIntVals;

    // std::vector<TString> fSpecTreeNames;
    // std::vector<TString> fSpecGlobs;
    // std::vector<std::vector<const TCut*> > fSpecSelections;
    // std::vector<TTree*> fSpecTrees;

    static const int numEffVars = 2;
    enum{
      kSNR,
      kEnergy
    };
    static const int numCutOrders = 2;
    enum {
      kInSequence,
      kIfFirst
      // kIfLast
    };

    std::vector<TEfficiency*>fSignalEffs[numEffVars][numCutOrders];









    // /**
    //  * @class FormulaHolder
    //  * @brief Contains the TTreeFormula and trick TChain into notifying all the formulas
    //  */

    // class FormulaHolder : public TObject {
    // public:
    //   FormulaHolder(TChain* c);
    //   virtual ~FormulaHolder();
    //   virtual Bool_t Notify();
    //   virtual size_t add(const char* formulaString);
    //   TTreeFormula* at(UInt_t i) {return fForms.at(i);}
    //   const char* str(UInt_t i) const {return fFormStrs.at(i);}
    //   size_t N(){return fForms.size();}
    // protected:
    //   TChain* fChain;
    //   std::vector<TTreeFormula*> fForms;
    //   std::vector<const char*> fFormStrs;
    // };







    

  public:

    /**
     * @class FisherResult
     * @brief Get the results of the Fisher Discriminant into a more useful form
     */

    class FisherResult : public TNamed {

      typedef std::map<int, double> WeightMap;
      typedef std::map<int, TString> ExpressionMap;

    public:
      /** 
       * Default constructor
       * 
       * @param fileName is the fileName containing the xml
       */
      FisherResult(const char* fileName = "") : TNamed("FisherResult", fileName){
	getResultFromXML(fileName);
      }
      TString getFisherFormula() const;
      TH2D* makeTimeHist(int nBinsX, int nBinsY, TTree* t, EColor col, int varInd=0, const char* extraName = "") const;
      virtual void Print(Option_t* opt = "") const;
      void getExpressions(std::vector<TString>& expressions) const;
      double getWeight(const char* expression);

    protected:
      void getResultFromXML(const char* filename);
      void parseNode(TXMLEngine* xml, XMLNodePointer_t node, Int_t level);

      ExpressionMap fExpressions; /// List of variable names
      WeightMap fWeights; /// List of variable weights

      ClassDef(FisherResult, 1);
    };
    

    
    
  
  };

}
#endif //CUT_OPTIMIZER
