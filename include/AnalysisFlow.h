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

class AnitaEventSummary;
class NoiseMonitor;

namespace Acclaim
{

/**
 * @class AnalysisFlow
 * @brief High level interface to join up the various bits of the analysis
 *
 * AnalysisFlow takes care of event selection, loading input data, saving output data, and contains a template function to loop over said data while performing the analysis reconstruction
 */
class AnalysisFlow : public TObject{

  public:

    enum selection{
      kAll                = 0,
      kWaisPulser         = 1,
      kDecimated          = 2,
      kQuietTime          = 3,
      kWaisPulserAndNonRF = 4 // I need min bias events for SNR estimates with the noise monitor class
    };

    AnalysisFlow(const char* outFileBaseName, int run, selection selection, FilterStrategy* filterStrat=NULL, AnitaDataset::BlindingStrategy blindStrat=AnitaDataset::kDefault, int theDivision=0, int theNumDivisions=1);
    ~AnalysisFlow();

    void doAnalysis();
    AnitaEventSummary* doEntry(Long64_t entry);
    AnitaEventSummary* doEvent(UInt_t eventNumber);
    AnalysisReco* getReco(){return fReco;}

 protected:

    void prepareDataSet();
    void prepareOutputFiles();
    Bool_t shouldIDoThisEvent(RawAnitaHeader* header, UsefulAdu5Pat* usefulPat);
    Bool_t isPulserWAIS(RawAnitaHeader* header, UsefulAdu5Pat* usefulPat);
    Bool_t isPulserLDB(RawAnitaHeader* header, UsefulAdu5Pat* usefulPat);    
    void setPulserFlags(RawAnitaHeader* header, UsefulAdu5Pat* usefulPat, AnitaEventSummary* sum);  

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
    AnitaEventSummary* fEventSummary;
    NoiseMonitor* fNoiseMonitor;
    UInt_t fLastEventConsidered;
  
    TString fOutFileBaseName;
    TFile* fOutFile;
    TTree* fSumTree;

    void prepareEverything();

    ANALYSIS_SETTING(Int_t, DoAll);
    ANALYSIS_SETTING(Double_t, NoiseTimeScaleSeconds);
    ANALYSIS_SETTING(Int_t, NoiseEvenWaveforms);
    ANALYSIS_SETTING(Int_t, OutFileCompressionLevel);

    ClassDef(AnalysisFlow, 0);
  };

}


#endif
