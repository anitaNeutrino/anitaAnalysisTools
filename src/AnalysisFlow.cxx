#include "AnalysisFlow.h"
#include "OutputConvention.h"
#include "ProgressBar.h"
#include "FilterStrategy.h"

/** 
 * Constructor, sets up some of the options for the analysis
 * 
 * @param run is the run to analyse
 * @param selection is the event selection to apply
 * @param blindStrat is the blinding strategy to use
 * @param division selects which subset of the run to do, goes from 0 -> numDivisions -1, default is 0.
 * @param numDivisions is the number of bits we're splitting the event into
 */
AnalysisFlow::AnalysisFlow(const char* outFileBaseName, int run, AnalysisFlow::selection selection, BlindDataset::strategy blindStrat, int division, int numDivisions){

  fOutFileBaseName = TString::Format("%s", outFileBaseName);
  fSelection = selection;
  fBlindStrat = blindStrat;
  fDivision = division;
  fNumDivisions = numDivisions;
  fRun = run;

  fData = NULL;
  fCrossCorr = NULL;
  fOutFile = NULL;
}



/** 
 * Destructor
 * 
 */
AnalysisFlow::~AnalysisFlow(){

  if(fData){
    delete fData;
  }

  if(fCrossCorr){
    delete fCrossCorr;
    fCrossCorr = NULL;
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
void AnalysisFlow::prepareDataSet(){

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
 * Coax the OutputConvention class into making an appropriately named output file/
 */
void AnalysisFlow::prepareOutputFiles(){


  if(fOutFile==NULL){

    // The output convention was originally designed to take argc and argv from the start of main
    // so here I reconstruct what they might have been.
    // Previously I used the argv to specify the run and subdivision in the file name, so I do that again here.
    
    std::vector<char*> fakeArgv;
    fakeArgv.push_back((char*) fOutFileBaseName.Data());
    TString runStr = TString::Format("%d", fRun);
    fakeArgv.push_back((char*) runStr.Data());


    if(fNumDivisions > 1){
      
      Int_t numDigitsTotal = TMath::Log10(fNumDivisions);
      Int_t numDigitsThis = fDivision > 0 ? TMath::Log10(fDivision) : 1;

      TString extra;
      for(int i=0; i < numDigitsTotal - numDigitsThis; i++){
	extra += TString::Format("0");
      }
      extra += TString::Format("%d", fDivision);

      fakeArgv.push_back((char*) extra.Data());
    }
    Int_t fakeArgc = (Int_t) fakeArgv.size();

    OutputConvention oc(fakeArgc, &fakeArgv[0]);

    fOutFile = oc.makeFile();  
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
Bool_t AnalysisFlow::shouldIDoThisEvent(RawAnitaHeader* header, UsefulAdu5Pat* usefulPat){

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
void AnalysisFlow::doAnalysis(){
    
  if(!fData){
    prepareDataSet();
  }

  if(!fCrossCorr){
    fCrossCorr = new CrossCorrelator();
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

  FilterStrategy filterStrat(fOutFile);

  for(Long64_t entry = fFirstEntry; entry < fLastEntry; entry++){

    fData->getEntry(entry);
    RawAnitaHeader* header = fData->header();
    Adu5Pat* pat = fData->gps();
    UsefulAdu5Pat usefulPat(pat);

    Bool_t selectedEvent = shouldIDoThisEvent(header, &usefulPat);

    if(selectedEvent){
    
      UsefulAnitaEvent* usefulEvent = fData->useful();

      
      FilteredAnitaEvent filteredEvent(usefulEvent, &filterStrat, pat, header, false);

      eventSummary = new AnitaEventSummary(header, &usefulPat);
      fCrossCorr->reconstructEvent(&filteredEvent, usefulPat, eventSummary);

      fSumTree->Fill();
      delete eventSummary;
      eventSummary = NULL;
    }
    
    p.inc(entry, numEntries);
  }
}  
