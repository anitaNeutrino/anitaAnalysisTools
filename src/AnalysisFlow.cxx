#include "AnalysisFlow.h"
#include "OutputConvention.h"
#include "ProgressBar.h"
#include "FilterStrategy.h"
#include "QualityCut.h"
#include "AcclaimFilters.h"

/** 
 * Constructor, sets up some of the options for the analysis 
* 
 * @param run is the run to analyse
 * @param selection is the event selection to apply
 * @param blindStrat is the blinding strategy to use
 * @param division selects which subset of the run to do, goes from 0 -> numDivisions -1, default is 0.
 * @param numDivisions is the number of bits we're splitting the event into
 */
Acclaim::AnalysisFlow::AnalysisFlow(const char* outFileBaseName, int run, Acclaim::AnalysisFlow::selection selection, FilterStrategy* filterStrat, BlindDataset::strategy blindStrat, int division, int numDivisions){

  fOutFileBaseName = TString::Format("%s", outFileBaseName);
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

    std::cout << "Found " << sgeTaskId << " " << sgeTaskIdEnv
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
}





/** 
 * Destructor
 * 
 */
Acclaim::AnalysisFlow::~AnalysisFlow(){

  if(fData){
    delete fData;
    fData = NULL;
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
 * Create the data set.
 * 
 */
void Acclaim::AnalysisFlow::prepareDataSet(){

  if(fData==NULL){
    bool doDecimated = fSelection == kDecimated ? true : false;
    fData = new BlindDataset(fBlindStrat, fRun, doDecimated);

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
 * Coax the OutputConvention class into making an appropriately named output files
 */
void Acclaim::AnalysisFlow::prepareOutputFiles(){


  if(fOutFile==NULL){

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

    // std::cout << fOutFile << "\t" << fOutFile->GetName() << std::endl;
  }  
  
}



/** 
 * Applies event selection
 * 
 * @param header is the RawAnitaHeader for this event
 * @param usefulPat is the UsefulAdu5Pat for the event
 * 
 * @return true is event satisfies selection criteria, false otherwise
 */
Bool_t Acclaim::AnalysisFlow::shouldIDoThisEvent(RawAnitaHeader* header, UsefulAdu5Pat* usefulPat){

  Bool_t doEvent = false;

  // for the WAIS selection
  const Double_t maxDeltaTriggerTimeNs = 1200;  
  UInt_t triggerTimeNsExpected = usefulPat->getWaisDivideTriggerTimeNs();
  UInt_t triggerTimeNs = header->triggerTimeNs;
  Int_t deltaTriggerTimeNs = Int_t(triggerTimeNs) - Int_t(triggerTimeNsExpected);

  switch (fSelection){

  case kWaisPulser:
    doEvent = TMath::Abs(deltaTriggerTimeNs) < maxDeltaTriggerTimeNs ? true : false;
    break;

  case kDecimated:
    doEvent = true;
    break;

  case kAll:
    doEvent = true;
    break;
    
  default:
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unknown event selection." << std::endl;
    doEvent = false;
    break;
  }

  return doEvent;
}







/** 
 * Does the main analysis loop
 */
void Acclaim::AnalysisFlow::doAnalysis(){
    
  if(!fData){
    prepareDataSet();
  }

  if(!fReco){
    fReco = new AnalysisReco();
  }

  if(!fOutFile){
    prepareOutputFiles();
  }

  if(!fSumTree){
    fSumTree = new TTree("sumTree", "Tree of AnitaEventSummaries");
  }

  AnitaEventSummary* eventSummary = NULL;
  fSumTree->Branch("sum", &eventSummary);

  const Long64_t numEntries = fLastEntry-fFirstEntry;
  ProgressBar p(numEntries);

  if(fFilterStrat){
    // this doesn't necessarily mean the output will be saved
    // it will only be saved if operations were added with the optional enable_output bool = true
    fFilterStrat->attachFile(fOutFile);
  }
  else{
    // empty strategy does nothing
    fFilterStrat = new FilterStrategy();
  }

  UInt_t lastEventConsidered = 0;
  for(Long64_t entry = fFirstEntry; entry < fLastEntry; entry++){

    fData->getEntry(entry);
    RawAnitaHeader* header = fData->header();
    UsefulAnitaEvent* usefulEvent = fData->useful();

    // make FourierBuffer filters behave as if we were considering sequential events
    // this is useful for the decimated data set...
    if(lastEventConsidered + 1 != header->eventNumber){
      Filters::makeFourierBuffersLoadHistoryOnNextEvent(fFilterStrat);
    }


    Adu5Pat* pat = fData->gps();
    UsefulAdu5Pat usefulPat(pat);

    Bool_t needToReconstruct = shouldIDoThisEvent(header, &usefulPat);

    if(needToReconstruct){
      eventSummary = new AnitaEventSummary(header, &usefulPat);
      Bool_t isGoodEvent = QualityCut::applyAll(usefulEvent, eventSummary);

      if(!isGoodEvent){
	// since we now have rolling averages make sure the filter strategy is sees every event before deciding whether or not to reconstruct
	FilteredAnitaEvent filteredEvent(usefulEvent, fFilterStrat, pat, header, false);
	// fReco->reconstructEvent(&filteredEvent, usefulPat, eventSummary);
	fReco->process(&filteredEvent, &usefulPat, eventSummary);
      }

      fSumTree->Fill();
      delete eventSummary;
      eventSummary = NULL;
    }
    
    lastEventConsidered = header->eventNumber;
    p.inc(entry, numEntries);
  }
}


