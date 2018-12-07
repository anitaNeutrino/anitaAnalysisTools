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

  class CmdLineArgs;


  /**
   * @class AnalysisFlow
   * @brief High level interface to join up the various bits of the analysis
   *
   * AnalysisFlow takes care of event selection, loading input data, saving output data, and contains a template function to loop over said data while performing the analysis reconstruction
   */
  class AnalysisFlow : public TObject{

  public:

    /**
     * @enum selection
     * @brief A human readable method of selecting the set of events to process, this can be expanded as required, @see #shouldIDoThisEvent(RawAnitaHeader* header, UsefulAdu5Pat* usefulPat) where the logic lives
     */
    enum selection{
      kAll,
      kWaisPulser,
      kDecimated
    };


    /** 
     * Translate the AnalysisFlow::selection enum to a descriptive string
     * 
     * @param sel is the AnalysisFlow event selection
     * 
     * @return c-style string, pointer to const char* array
     */
    static const char* selectionAsString(selection sel){
      switch(sel){
      case kWaisPulser: return "wais";
      case kDecimated: return "decimated";
      case kAll:
      default:
	return "all";
      }
    }


    /** 
     * @brief Preferred constructor
     *
     * Also searches for the environment variable SGE_TASK_ID, which indicates this code is running on the hoffman2 cluster.
     * If this is the case, then some of the passed variables (#fRun,  #fDivision, #fNumDivisions) are overwritten with information endoded in SGE_TASK_ID.
     * This allows for easy cluster scripting.
     * 
     * @param args contains the command line arguments
     * @param filterStrat is the filter stragtegy to use
     */    
    AnalysisFlow(CmdLineArgs* args, FilterStrategy* filterStrat=NULL);


    /** 
     * @brief Destructor
     */    
    ~AnalysisFlow();


    /** 
     * Does the main analysis loop for all specified events
     * 
     * @param startAtThisEntry optional, (default = -1), if you wish force the analysis loop to start at a specific tree entry
     */
    void doAnalysis(Long64_t startAtThisEntry=-1);


    /** 
     * Does my analysis on a single event in AnitaDataset member #fData trees, referenced by entry 
     * 
     * @param entry is the entry to process
     * 
     * @return the generated AnitaEventSummary, it is the caller's responsibility to delete this.
     */    
    AnitaEventSummary* doEntry(Long64_t entry);

    /** 
     * Does my analysis on a single event in AnitaDataset member #fData trees, referenced by RawAnitaHeader::eventNumber
     * 
     * @param eventNumber is the event to process, must be in the run and event selection!
     * 
     * @return the generated AnitaEventSummary, it is the caller's responsibility to delete this.
     */    
    AnitaEventSummary* doEvent(UInt_t eventNumber);

    /** 
     * Get a pointer to the AnalysisReco member, #fReco, which does the hard work of event reconstruction
     * @return pointer to #fData
     */
    AnalysisReco* getReco(){return fReco;} 

    /** 
     * Get a pointer to the FilteredAnitaEvent, #fEv (this is mostly for MagicDisplay)
     * @return A poiner to the #fEv
     */
    FilteredAnitaEvent* getEvent() {return fEv;}

    /** 
     * Get a pointer to the contained AnitaDataset #fData member
     * @return pointer to #fData
     */
    const AnitaDataset* getData(){return fData;}

    /** 
     * Get the first entry to be fetched in the #fData
     * @return the value of #fFirstEntry
     */
    Long64_t firstEntry(){return fFirstEntry;}

    /** 
     * Get the final entry to be processed in the #fData
     * 
     * This is not the last proccessed entry, but the final entry to be processed in #fData.
     * @return the value of #fLastEntry
     */
    Long64_t lastEntry(){return fLastEntry;}

  protected:
    
    selection fSelection;		///< Which event selection is applied? @see AnalysisFlow::selection
    int fRun;                           ///< Which run in the data set should be is being processed?
    int fNumDivisions;			///< For parallel processing the run is divided into this many portions
    int fDivision;			///< Which portion of the run to process (runs from 0 to fNumDivisions)

    AnitaDataset* fData;		///< Rootified data handler
    Long64_t fFirstEntry;		///< The first entry to process (derived from fDivision and fNumDivisions)
    Long64_t fLastEntry;		///< The final entry to process (derived from fDivision and fNumDivisions)

    AnalysisReco* fReco;		///< The reconstruction class, fills an AnitaEventSummary
    FilterStrategy* fFilterStrat;	///< Which filter strategy is applied to the events
    AnalysisSettings* fSettings;	///< Contains configurable numbers parsed from an Acclaim analysis settings file.
    AnitaEventSummary* fEventSummary;	///< The most recently filled AnitaEventSummary, updated after each event processed
    FilteredAnitaEvent* fEv;		///< The most recently produced FilteredAnitaEvent, updated after each event processed
    NoiseMonitor* fNoiseMonitor;	///< Measures the noise
    UInt_t fLastEventConsidered;	///< Tracks the event numbers processed by the class

    TString fOutFileBaseName;		///< The meat of the output file name
    TFile* fOutFile;			///< the output file, will contain TTree of AnitaEventSummary
    TTree* fSumTree;			///< The produced TTree of AnitaEventSummary

    /** 
     * @brief Create the data set if not already done
     *
     * Creates an instance of the AnitaDataset class and finds the first/last entries to process using the division/numDivision member variables.
     */
    void prepareDataSet();


    /** 
     * @brief Coax the OutputConvention class into creating appropriately named output files
     *
     * The OutputConvention class was written some time ago to convert the cpp default main arguments (argc/argv) into an output file with a helpful name.
     * This function dances around the original implementation intention to generate some fake argc/argv variables and pass them to an OutputConvention object.
     */    
    void prepareOutputFiles();

    /** 
     * @brief Applies high level event selection
     * 
     * Look at the (useful) GPS and header information and implement the event selection logic for the AnalysisFlow::selection enum.
     * 
     * @param header is the RawAnitaHeader for this event
     * @param usefulPat is the UsefulAdu5Pat for the event
     * 
     * @return true is event satisfies selection criteria, false otherwise
     */    
    Bool_t shouldIDoThisEvent(RawAnitaHeader* header, UsefulAdu5Pat* usefulPat);



    /** 
     * Applies my WAIS pulser selection, which works for ANITA-3.
     * 
     * @param header is the event header
     * @param usefulPat is a usefulAdu5Pat object
     * 
     * @return true if the event matches the timing criteria
     */    
    Bool_t isPulserWAIS(RawAnitaHeader* header, UsefulAdu5Pat* usefulPat);



    /** 
     * Applies Linda's LDB pulser selection from the thesis(!) version of my analysis. Should still work, but I've not tested this in a while.
     * 
     * @param header is the event header
     * @param usefulPat is a usefulAdu5Pat object
     * 
     * @return true if the event matches the timing criteria
     */    
    Bool_t isPulserLDB(RawAnitaHeader* header, UsefulAdu5Pat* usefulPat);


    /**
     * Set the pulser flags in the passed AnitaEventSummary
     *
     * @param header is the event header
     * @param usefulPat is the ANITA gps data
     * @param sum is the AnitaEventSummary in which to set the flag
     */    
    void setPulserFlags(RawAnitaHeader* header, UsefulAdu5Pat* usefulPat, AnitaEventSummary* sum);


    /** 
     * @brief Initialize all input and output objects.
     */
    void prepareEverything(const char* preferredSettingsFileName = NULL);

    /**
     * @brief Look for an environmental variable called SGE_TASK_ID for running on the hoffman2 cluster.
     * 
     * If it is found this function *overrides* the fRun, fDivision, fNumDivision members.
     * A task array on the hoffman2 cluster gives each process a different SGE_TASK_ID so this is useful
     * to make the analysis run in parallel.
     * 
     * @return true if the variable is present.
     */
    Bool_t checkForSgeTaskId();


    /** 
     * @fn GetDebug
     * @brief Get the value of #fDebug
     */
    /** 
     * @fn SetDebug
     * @brief Set the value of #fDebug
     */
    /** 
     * @var fDebug
     * @brief Controls the printing of debug info, feel free to wrap noisy stderr messages with if(fDebug),
     * @see AnalysisSettings
     */        
    ANALYSIS_SETTING(Int_t, Debug);


    /** 
     * @fn GetDoAll
     * @brief Get the value of #fDoAll
     */
    /** 
     * @fn SetDoAll
     * @brief Set the value of #fDoAll
     */
    /** 
     * @var fDoAll
     * @brief Process all events (that won't crash the analysis software) regardless of quality 
     * @see AnalysisSettings
     */
    ANALYSIS_SETTING(Int_t, DoAll);


    /** 
     * @fn GetEnumBlindStrat
     * @brief Get the value of #fBlindStrat as AnalysisFlow::selection enum
     * @see AnalysisSettings
     */
    /** 
     * @fn GetBlindStrat
     * @brief Get the value of #fBlindStrat as an integer
     * @see AnalysisSettings
     */    
    /** 
     * @fn SetBlindStrat(Int_t val)
     * @brief Set the value of #fBlindStrat from an integer value
     * @see AnalysisSettings
     */
    /** 
     * @fn SetBlindStrat(AnitaDataset::BlindingStrategy val)
     * @brief Set the value of #fBlindStrat from an enum value
     * @see AnalysisSettings
     */    
    /** 
     * @var fBlindStrat
     * @brief The blinding strategy with which to initialize #fData
     * @see AnalysisSettings
     */    
    ENUM_ANALYSIS_SETTING(AnitaDataset::BlindingStrategy, BlindStrat);




    /** 
     * @fn GetUseNoiseMonitor
     * @brief Get the value of #fUseNoiseMonitor
     * @see AnalysisSettings
     */
    /** 
     * @fn SetUseNoiseMonitor
     * @brief Set the value of #fUseNoiseMonitor
     * @see AnalysisSettings
     */
    /** 
     * @var fUseNoiseMonitor
     * @brief Whether or not to use the noise monitor
     * @see AnalysisSettings
     */
    ANALYSIS_SETTING(Int_t, UseNoiseMonitor);

    
    /** 
     * @fn GetNoiseTimeScaleSeconds
     * @brief Get the value of #fNoiseTimeScaleSeconds
     * @see AnalysisSettings
     */
    /** 
     * @fn SetNoiseTimeScaleSeconds
     * @brief Set the value of #fNoiseTimeScaleSeconds
     * @see AnalysisSettings
     */
    /** 
     * @var fNoiseTimeScaleSeconds
     * @brief The N value for SNR measurements, comes from waveform RMS in nearby min-bias triggers, this variable controls how long to average over in seconds
     * @see AnalysisSettings
     */
    ANALYSIS_SETTING(Double_t, NoiseTimeScaleSeconds);


    /** 
     * @fn GetNoiseEvenWaveforms
     * @brief Get the value of #fNoiseEvenWaveforms
     * @see AnalysisSettings
     */
    /** 
     * @fn SetNoiseEvenWaveforms
     * @brief Set the value of #fNoiseEvenWaveforms
     * @see AnalysisSettings
     */
    /** 
     * @var fNoiseEvenWaveforms
     * @brief Derived the noise value for RMS from the even or the uneven waveforms? Which choice is sensible will depend on the filter strategy, #fFilterStrat.
     * @see AnalysisSettings
     */    
    ANALYSIS_SETTING(Int_t, NoiseEvenWaveforms);


    /** 
     * @fn GetOutFileCompressionLevel
     * @brief Get the value of #fOutFileCompressionLevel
     * @see AnalysisSettings
     */
    /** 
     * @fn SetOutFileCompressionLevel
     * @brief Set the value of #fOutFileCompressionLevel
     * @see AnalysisSettings
     */
    /** 
     * @var fOutFileCompressionLevel
     * @brief What compression level for the ROOT output files? There's no reason not to set this to the maximum possible value
     * @see AnalysisSettings
     */    
    ANALYSIS_SETTING(Int_t, OutFileCompressionLevel);

    /** 
     * @fn GetOutFileCompressionAlgo
     * @brief Get the value of #fOutFileCompressionAlgo
     * @see AnalysisSettings
     */
    /** 
     * @fn SetOutFileCompressionAlgo
     * @brief Set the value of #fOutFileCompressionAlgo
     * @see AnalysisSettings
     */
    /** 
     * @var fOutFileCompressionAlgo
     * @brief Which compression algorithm for the ROOT output files?
     * @see AnalysisSettings
     */    
    ANALYSIS_SETTING(Int_t, OutFileCompressionAlgo);


    /**
     * @brief Required for ROOT awareness of getter/setter functions.
     * @see AnalysisSettings
     */
    ClassDef(AnalysisFlow, 0);
  };

}


#endif
