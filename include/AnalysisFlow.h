/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             Applies the analysis in sequence
***********************************************************************************************************/


#ifndef ANALYSIS_FLOW_H
#define ANALYSIS_FLOW_H

#include "AnitaDataset.h"
#include "AnalysisReco.h"
#include "AnalysisSettings.h"
#include "TFile.h"
#include "FilterStrategy.h"


namespace Acclaim
{

/**
 * @class AnalysisFlow
 * @brief High level interface to join up the various bits of the analysis
 *
 * AnalysisFlow takes care of event selection, loading input data, saving output data, and contains a template function to loop over said data while performing the analysis reconstruction
 */
  class AnalysisFlow {

  public:

    enum selection{
      kAll             = 0,
      kWaisPulser      = 1,
      kDecimated       = 2,
      kQuietTime       = 3
    };



    AnalysisFlow(const char* outFileBaseName, int run, selection selection, FilterStrategy* filterStrat=NULL, AnitaDataset::BlindingStrategy blindStrat=AnitaDataset::kDefault, int theDivision=0, int theNumDivisions=1);
    ~AnalysisFlow();

    void doAnalysis();  
    void prepareDataSet();
    void prepareOutputFiles();
    Bool_t shouldIDoThisEvent(RawAnitaHeader* header, UsefulAdu5Pat* usefulPat);
    
  protected:
    int fRun;
    selection fSelection;
    AnitaDataset::BlindingStrategy fBlindStrat;
    int fDivision;
    int fNumDivisions;

    AnitaDataset* fData;
    Long64_t fFirstEntry;
    Long64_t fLastEntry;

    AnalysisReco* fReco;
    FilterStrategy* fFilterStrat;
    AnalysisSettings* fSettings;
  
    TString fOutFileBaseName;
    TFile* fOutFile;
    TTree* fSumTree;

    ANALYSIS_SETTING(Int_t, DoAll);
    ANALYSIS_SETTING(Double_t, NoiseTimeScaleSeconds);
    ANALYSIS_SETTING(Int_t, NoiseEvenWaveforms);
  };

}


#endif
