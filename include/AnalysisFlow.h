/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             Applies the analysis in sequence
***********************************************************************************************************/


#ifndef ANALYSIS_FLOW_H
#define ANALYSIS_FLOW_H

#include "BlindDataset.h"
#include "AnalysisReco.h"
#include "TFile.h"
#include "FilterStrategy.h"


namespace Acclaim
{
  class AnalysisFlow {

  public:

    enum selection{
      kAll             = 0,
      kWaisPulser      = 1,    
      kDecimated       = 2,
      kQuietTime       = 3
    };



    AnalysisFlow(const char* outFileBaseName, int run, selection selection, FilterStrategy* filterStrat=NULL, BlindDataset::strategy blindStrat=BlindDataset::kDefault, int theDivision=0, int theNumDivisions=1);
    ~AnalysisFlow();  

    void doAnalysis();  
    void prepareDataSet();
    void prepareOutputFiles();
    Bool_t shouldIDoThisEvent(RawAnitaHeader* header, UsefulAdu5Pat* usefulPat);

  private:
    int fRun;
    selection fSelection;
    BlindDataset::strategy fBlindStrat;
    int fDivision;
    int fNumDivisions;


    BlindDataset* fData;
    Long64_t fFirstEntry;
    Long64_t fLastEntry;


    AnalysisReco* fReco;  
    FilterStrategy* fFilterStrat;
  
    TString fOutFileBaseName;
    TFile* fOutFile;
    TTree* fSumTree;

    
  };

}


#endif
