#include "AnalysisFlow.h"
#include "OutputConvention.h"
#include "ProgressBar.h"
#include "FilterStrategy.h"
#include "QualityCut.h"
#include "AcclaimFilters.h"
#include "NoiseMonitor.h"

ClassImp(Acclaim::AnalysisFlow);

/** 
 * @brief Constructor
 *
 * Also searches for the environment variable SGE_TASK_ID, which indicates this code is running on the hoffman2 cluster.
 * If this is the case, then some of the passed variables (run,  division, numDivisions) are overwritten with information endoded in SGE_TASK_ID.
 * This allows for easy cluster scripting.
 * 
 * @param run is the run to analyse
 * @param selection is the event selection to apply
 * @param blindStrat is the blinding strategy to use
 * @param division selects which sub division of the run to do. The run is divied into numDivisions divisions, so division goes from 0 -> numDivisions-1, default is 0.
 * @param numDivisions divides the run into pieces. 
 */

Acclaim::AnalysisFlow::AnalysisFlow(const char* outFileBaseName, int run, Acclaim::AnalysisFlow::selection selection, FilterStrategy* filterStrat, AnitaDataset::BlindingStrategy blindStrat, int division, int numDivisions){

  fOutFileBaseName = outFileBaseName ? TString::Format("%s", outFileBaseName) : "";
  fSelection = selection;
  fFilterStrat = filterStrat;
  fBlindStrat = blindStrat;


  // on the Hoffman2 cluster, this is the task ID
  // we're going to use it to try and figure out how to set the division if it wasn't explicitly set
  // assuming this wasn't explicitly set...
  const char* sgeTaskId = "SGE_TASK_ID";
  const char* sgeTaskIdEnv = getenv(sgeTaskId);
  if(sgeTaskIdEnv && numDivisions==1){    
    // we have runs 130-439
    // I'm going to set the job array like:
    // 130-439 for 1 job per run = 309 jobs
    // 1300-4390 for 10 jobs per run = 3090 jobs
    // 13000-43900 for 100 jobs per run = 30900 jobs (that's probably excessive)
    Int_t jobArrayIndex = atoi(sgeTaskIdEnv);

    if(jobArrayIndex < 1e3){
      fNumDivisions = 1;
      fRun = jobArrayIndex;
      fDivision = 0;
    }
    else if(jobArrayIndex < 1e4){
      fNumDivisions = 10;
      fRun = jobArrayIndex/fNumDivisions;
      fDivision = jobArrayIndex%10;
    }
    else if(jobArrayIndex < 1e5){
      fNumDivisions = 100;
      fRun = jobArrayIndex/fNumDivisions;      
      fDivision = jobArrayIndex%100;
    }
    else{
      std::cerr << "Error in " << __FILE__ << " couldn't figure out the run/divisions from " << sgeTaskId << ". Giving up." << std::endl;
      exit(1);
    }

    std::cout << "Info in " << __PRETTY_FUNCTION__ << ". "
              << "Found " << sgeTaskId << " " << sgeTaskIdEnv
	      << ", so the jobArrayIndex = " << jobArrayIndex
	      << ", set fRun = " << fRun
	      << ", fNumDivisions = " << fNumDivisions
	      << ", fDivision = " << fDivision << "." << std::endl;
  }
  else{
    fDivision = division;
    fNumDivisions = numDivisions;
    fRun = run;
  }


  fSumTree = NULL;
  fData = NULL;
  fReco = NULL;
  fOutFile = NULL;
  fSettings = NULL;
  fEventSummary = NULL;

  fFirstEntry=0;
  fLastEntry=0;
  fLastEventConsidered = 0;

  prepareEverything();
}





/** 
 * @brief Destructor
 * 
 */
Acclaim::AnalysisFlow::~AnalysisFlow(){

  if(fData){
    delete fData;
    fData = NULL;
  }

  if(fSettings){
    if(fOutFile){
      fSettings->write(fOutFile);
    }
    delete fSettings;
  }

  if(fReco){
    delete fReco;
    fReco = NULL;
  }

  if(fFilterStrat){
    delete fFilterStrat;
    fFilterStrat = NULL;
  }

  if(fOutFile){
    fOutFile->Write();
    fOutFile->Close();
    fOutFile = NULL;
  }
}




/** 
 * @brief Create the data set if not already done
 *
 * Creates an instance of the AnitaDataset class and finds the first/last entries to process using the division/numDivision member variables.
 */
void Acclaim::AnalysisFlow::prepareDataSet(){

  if(fData==NULL){
    bool doDecimated = fSelection == kDecimated ? true : false;
    fData = new AnitaDataset(fRun, doDecimated, WaveCalType::kDefault, AnitaDataset::ANITA_ROOT_DATA, fBlindStrat);

    Long64_t numEntriesInRun = fData->N();

    Long64_t numEventsPerDivision = numEntriesInRun/fNumDivisions;

    fFirstEntry = numEventsPerDivision*fDivision;
    fLastEntry = fFirstEntry + numEventsPerDivision;

    if(fDivision==fNumDivisions-1){
      Long64_t numRemainder = numEntriesInRun%fNumDivisions;
      fLastEntry += numRemainder;
    }
  }
}



/** 
 * @brief Coax the OutputConvention class into creating appropriately named output files
 *
 * The OutputConvention class was written some time ago to convert the cpp default main arguments (argc/argv) into an output file with a helpful name.
 * This function goes around the houses to generate some fake argc/argv variables and passes them to an OutputConvention object. 
 */
void Acclaim::AnalysisFlow::prepareOutputFiles(){


  if(fOutFile==NULL && fOutFileBaseName!=""){

    // The output convention was originally designed to take argc and argv from the start of main
    // so here I reconstruct what they might have been.
    // Previously I used the argv to specify the run and subdivision in the file name, so I do that again here.
    
    std::vector<char*> fakeArgv;
    fakeArgv.push_back((char*) fOutFileBaseName.Data());
    TString runStr = TString::Format("%d", fRun);
    fakeArgv.push_back((char*) runStr.Data());
    
    TString extra;
    if(fNumDivisions > 1){
      
      Int_t numDigitsTotal = floor(TMath::Log10(fNumDivisions) + 1);
      Int_t numDigitsThis = fDivision == 0 ? 1 : floor(TMath::Log10(fDivision) + 1);

      std::cout << fNumDivisions << "\t" << numDigitsTotal << std::endl;
      std::cout << fDivision << "\t" << numDigitsThis << std::endl;

      for(int i=0; i < numDigitsTotal - numDigitsThis; i++){
	std::cout << i << "\t" << numDigitsTotal - numDigitsThis << std::endl;
	
	extra += TString::Format("0");
      }
      extra += TString::Format("%d", fDivision);

      std::cout << extra << std::endl;
      
      fakeArgv.push_back((char*) extra.Data());
    }
    Int_t fakeArgc = (Int_t) fakeArgv.size();

    // std::cout << fakeArgc << ", " << &fakeArgv[0] << std::endl;
    // for(int i=0; i < fakeArgc; i++){
    //   std::cout << "\t:" << i << ", " << fakeArgv[i] << std::endl;    
    // }

    OutputConvention oc(fakeArgc, &fakeArgv[0]);

    fOutFile = oc.makeFile();
    fOutFile->SetCompressionLevel(fOutFileCompressionLevel);

    // std::cout << fOutFile << "\t" << fOutFile->GetName() << std::endl;
  }  
  
}



/** 
 * @brief Applies high level event selection
 * 
 * @param header is the RawAnitaHeader for this event
 * @param usefulPat is the UsefulAdu5Pat for the event
 * 
 * @return true is event satisfies selection criteria, false otherwise
 */
Bool_t Acclaim::AnalysisFlow::shouldIDoThisEvent(RawAnitaHeader* header, UsefulAdu5Pat* usefulPat){

  Bool_t doEvent = false;

  switch (fSelection){

    case kWaisPulser:
      doEvent = isPulserWAIS(header, usefulPat);
      break;

    case kDecimated:
      doEvent = true;
      break;

    case kAll:
      doEvent = true;
      break;

    case kWaisPulserAndNonRF:
      doEvent = isPulserWAIS(header, usefulPat) || !(header->getTriggerBitRF());
      break;
      
    default:
      std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unknown event selection." << std::endl;
      doEvent = false;
      break;
  }

  return doEvent;
}



/** 
 * Applies my WAIS pulser selection
 * 
 * @param header is the event header
 * @param usefulPat is a usefulAdu5Pat object
 * 
 * @return true if the event matches the timing criteria
 */
Bool_t Acclaim::AnalysisFlow::isPulserWAIS(RawAnitaHeader* header, UsefulAdu5Pat* usefulPat){

  const UInt_t triggerTimeNsExpected = usefulPat->getWaisDivideTriggerTimeNs();
  const int maxSeparationMeters = 1e6; // 1000 km
  const double c_ns = 0.3; // speed of light in nano-seconds
  Bool_t isWais = false;
  if(triggerTimeNsExpected < maxSeparationMeters/c_ns){
  
    const Double_t maxDeltaTriggerTimeNs = 1200;
    UInt_t triggerTimeNs = header->triggerTimeNs;
    Int_t deltaTriggerTimeNs = Int_t(triggerTimeNs) - Int_t(triggerTimeNsExpected);

    isWais = TMath::Abs(deltaTriggerTimeNs) < maxDeltaTriggerTimeNs;
    // std::cerr << __PRETTY_FUNCTION__ << "\t" << isWais << "\t" << triggerTimeNs << "\t" << triggerTimeNsExpected << std::endl;
  }
  return isWais;
}



/** 
 * Applies Linda's LDB pulser selection, I've not actually tested this in a while
 * 
 * @param header is the event header
 * @param usefulPat is a usefulAdu5Pat object
 * 
 * @return true if the event matches the timing criteria
 */
Bool_t Acclaim::AnalysisFlow::isPulserLDB(RawAnitaHeader* header, UsefulAdu5Pat* usefulPat){

  const UInt_t triggerTimeNsExpected = usefulPat->getLDBTriggerTimeNs();
  const int maxSeparationMeters = 1e6; // 1000 km
  const double c_ns = 0.3; // speed of light in nano-seconds
  Bool_t isLDB = false;

  if(triggerTimeNsExpected < maxSeparationMeters/c_ns){
  
    const double cutTimeNs = 1200;
    Int_t delay=0;
  
    if (fRun<150){ // LDB VPOL 
      delay =  25000000; // V-POL pulse at 25 ms runs 145-149
    }
    else if (fRun<154){ // LDB HPOL
      delay =  50000000; // H-POL pulse at 50 ms
    }
    else if (fRun<172){ // LDB VPOL
      delay =  50000000; // V-POL pulse at 50 ms runs 154-171    
    }

    const Int_t delayGenerator = 199996850; // delay generator was not exactly 200 ms
    const Int_t constdelay = 500;

    Int_t deltaTriggerTimeNs    = Int_t(header->triggerTimeNs) - Int_t(triggerTimeNsExpected);
    deltaTriggerTimeNs    = deltaTriggerTimeNs%(delayGenerator) - delay - constdelay;

    isLDB = TMath::Abs(deltaTriggerTimeNs) < cutTimeNs;
  }
  return isLDB;
}





/**
 * Set the pulser flags in the AnitaEventSummary
 *
 * @param header is the event header
 * @param usefulPat is the ANITA gps data
 * @param sum is the AnitaEventSummary in which to set the flag
 */
void Acclaim::AnalysisFlow::setPulserFlags(RawAnitaHeader* header, UsefulAdu5Pat* usefulPat, AnitaEventSummary* sum){

  if(isPulserWAIS(header, usefulPat)){
    sum->flags.pulser = AnitaEventSummary::EventFlags::WAIS;
  }
  else if(isPulserLDB(header, usefulPat)){
    sum->flags.pulser = AnitaEventSummary::EventFlags::LDB;    
  }
  // std::cerr << __PRETTY_FUNCTION__ << ": I just set the flags to be " << sum->flags.pulser << std::endl;
}



/** 
 * Set up all I/O
 * 
 */
void Acclaim::AnalysisFlow::prepareEverything(){

  if(fRun >= 257 && fRun <= 263){
    if(AnitaVersion::get()==3){
      std::cerr << "No ANITA-3 data for run " << fRun << " so won't do analysis" << std::endl;
      return;
    }
  }
  
  if(!fData){
    prepareDataSet();
  }

  if(!fSettings){
    fSettings = new AnalysisSettings();
    fSettings->apply(dynamic_cast<TObject*>(this));
  }

  if(!fReco){
    fReco = new AnalysisReco();
    fSettings->apply(fReco);
  }

  if(!fOutFile){
    prepareOutputFiles();
  }

  fEventSummary = NULL;
  if(fOutFile && !fSumTree){
    fSumTree = new TTree("sumTree", "Tree of AnitaEventSummaries");
    fSumTree->Branch("sum", &fEventSummary);    
  }

  if(fFilterStrat){
    // this doesn't necessarily mean the output will be saved
    // it will only be saved if operations were added with the optional enable_output bool = true
    fFilterStrat->attachFile(fOutFile);
  }
  else{
    // empty strategy does nothing
    fFilterStrat = new FilterStrategy();
  }  
  
  fNoiseMonitor = new NoiseMonitor(fNoiseTimeScaleSeconds, NoiseMonitor::kUneven, fOutFile);
  fLastEventConsidered = 0;
}



/** 
 * Does my analysis on a single eventNumber in the run
 * 
 * @param eventNumber is the event to process, must be in the run and event selection!
 * 
 * @return the generated AnitaEventSummary, it is the caller's responsibility to delete this.
 */
AnitaEventSummary* Acclaim::AnalysisFlow::doEvent(UInt_t eventNumber){
  int entry = fData->getEvent(eventNumber);
  return entry > 0 ? doEntry(entry) : NULL;
}



/** 
 * Does my analysis on a single entry in the run
 * 
 * @param entry is the entry to process
 * 
 * @return the generated AnitaEventSummary, it is the caller's responsibility to delete this.
 */
AnitaEventSummary* Acclaim::AnalysisFlow::doEntry(Long64_t entry){

  fData->getEntry(entry);
  RawAnitaHeader* header = fData->header();
  UsefulAnitaEvent* usefulEvent = fData->useful();

  // make FourierBuffer filters behave as if we were considering sequential events
  // this is useful for the decimated data set...
  if(fLastEventConsidered + 1 != header->eventNumber){
    Filters::makeFourierBuffersLoadHistoryOnNextEvent(fFilterStrat);
  }

  Adu5Pat* pat = fData->gps();
  UsefulAdu5Pat usefulPat(pat);

  Bool_t needToReconstruct = shouldIDoThisEvent(header, &usefulPat);

  fEventSummary = NULL;
  
  if(needToReconstruct){
    fEventSummary = new AnitaEventSummary(header, &usefulPat);
    Bool_t isGoodEvent = QualityCut::applyAll(usefulEvent, fEventSummary);

    if(isGoodEvent || fDoAll){

      setPulserFlags(header, &usefulPat, fEventSummary);
        
      // since we now have rolling averages make sure the filter strategy is sees every event before deciding whether or not to reconstruct
      FilteredAnitaEvent filteredEvent(usefulEvent, fFilterStrat, pat, header, false);

      fNoiseMonitor->update(&filteredEvent);
        
      fReco->process(&filteredEvent, &usefulPat, fEventSummary, fNoiseMonitor);
    }

    if(fSumTree){
      fSumTree->Fill();
    }
  }
    
  fLastEventConsidered = header->eventNumber;
  return fEventSummary;

}



/** 
 * Does the main analysis loop for all specified events
 */
void Acclaim::AnalysisFlow::doAnalysis(){

  fLastEventConsidered = 0;

  const Long64_t numEntries = fLastEntry-fFirstEntry;
  ProgressBar p(numEntries);
  
  for(Long64_t entry = fFirstEntry; entry < fLastEntry; entry++){
    AnitaEventSummary* sum = doEntry(entry);
    if(sum){
      delete sum;
    }
    p.inc(entry, numEntries);
  }
}
