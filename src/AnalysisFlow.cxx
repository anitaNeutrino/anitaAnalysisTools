#include "AnalysisFlow.h"
#include "OutputConvention.h"
#include "ProgressBar.h"
#include "FilterStrategy.h"
#include "QualityCut.h"
#include "AcclaimFilters.h"
#include "NoiseMonitor.h"
#include "AcclaimCmdLineArgs.h"
#include "UsefulAnitaEvent.h"
#include "Compression.h"

ClassImp(Acclaim::AnalysisFlow);


Acclaim::AnalysisFlow::AnalysisFlow(Acclaim::CmdLineArgs* args, FilterStrategy* filterStrat){

  fOutFileBaseName = args->output_filename;
  
  fSelection = (AnalysisFlow::selection) args->event_selection;
  fFilterStrat = filterStrat;
  
  if(!checkForSgeTaskId()){
    fDivision = args->division;
    fNumDivisions = args->numdivisions;
    fRun = args->run;
  }



  fDebug = 0;
  fOutFileCompressionLevel = 9;
  fOutFileCompressionAlgo = 2;
  fBlindStrat = AnitaDataset::kNoBlinding;

  prepareEverything(args->settings_filename);
}






Bool_t Acclaim::AnalysisFlow::checkForSgeTaskId(){
  // on the Hoffman2 cluster, this is the task ID
  // we're going to use it to try and figure out the run, division and numdivision
  // assuming this wasn't explicitly set...
  const char* sgeTaskId = "SGE_TASK_ID";
  const char* sgeTaskIdEnv = getenv(sgeTaskId);

  if(sgeTaskIdEnv){    
    // we have runs 130-439
    // I'm going to set the job array like:
    // 130-439 for 1 job per run = 309 jobs
    // 1300-4399 for 10 jobs per run = 3099 jobs
    // 13000-43999 for 100 jobs per run = 30999 jobs (that's probably excessive)
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
    return true;
  }
  else{
    return false;
  }
}



Acclaim::AnalysisFlow::~AnalysisFlow(){

  if(fEv){
    delete fEv;
    fEv = NULL;
  }

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

    // std::cerr << "Info in " << __PRETTY_FUNCTION__ << "Current blinding strategy: " << fData->getDescription(fData->getStrategy()) << std::endl;

  }
}


void Acclaim::AnalysisFlow::prepareOutputFiles(){


  if(fOutFile==NULL && fOutFileBaseName!=""){

    // The output convention was originally designed to take argc and argv from the start of main
    // so here I reconstruct what they might have been.
    // Previously I used the argv to specify the run and subdivision in the file name, so I do that again here.
    
    std::vector<char*> fakeArgv;
    fakeArgv.push_back((char*) fOutFileBaseName.Data());
    TString runStr = TString::Format("%d", fRun);
    fakeArgv.push_back((char*) runStr.Data());
    
    TString extraDigits;
    if(fNumDivisions > 1){
      
      Int_t numDigitsTotal = floor(TMath::Log10(fNumDivisions) + 1);
      Int_t numDigitsThis = fDivision == 0 ? 1 : floor(TMath::Log10(fDivision) + 1);

      // std::cout << fNumDivisions << "\t" << numDigitsTotal << std::endl;
      // std::cout << fDivision << "\t" << numDigitsThis << std::endl;

      for(int i=0; i < numDigitsTotal - numDigitsThis; i++){
	// std::cout << i << "\t" << numDigitsTotal - numDigitsThis << std::endl;
	
	extraDigits += TString::Format("0");
      }
      extraDigits += TString::Format("%d", fDivision);

      // std::cout << extraDigits << std::endl;
      
      fakeArgv.push_back((char*) extraDigits.Data());
    }
    Int_t fakeArgc = (Int_t) fakeArgv.size();

    // std::cout << fakeArgc << ", " << &fakeArgv[0] << std::endl;
    // for(int i=0; i < fakeArgc; i++){
    //   std::cout << "\t:" << i << ", " << fakeArgv[i] << std::endl;    
    // }

    OutputConvention oc(fakeArgc, &fakeArgv[0]);

    fOutFile = oc.makeFile();
    if(fOutFileCompressionLevel >= 0){
      fOutFile->SetCompressionLevel(fOutFileCompressionLevel);
      // kUseGlobalCompressionSetting,
      //                           kZLIB,
      //                           kLZMA,
      //                           kOldCompressionAlgo
      fOutFile->SetCompressionAlgorithm((ROOT::ECompressionAlgorithm) fOutFileCompressionAlgo);
    }

    // std::cout << fOutFile << "\t" << fOutFile->GetName() << std::endl;
  }  
  
}



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

    // case kWaisPulserAndNonRF:
    //   doEvent = isPulserWAIS(header, usefulPat) || !(header->getTriggerBitRF());
    //   break;
      
    default:
      std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unknown event selection." << std::endl;
      doEvent = false;
      break;
  }

  return doEvent;
}



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

    Int_t deltaTriggerTimeNs = Int_t(header->triggerTimeNs) - Int_t(triggerTimeNsExpected);
    deltaTriggerTimeNs = deltaTriggerTimeNs%(delayGenerator) - delay - constdelay;

    isLDB = TMath::Abs(deltaTriggerTimeNs) < cutTimeNs;
  }
  return isLDB;
}





void Acclaim::AnalysisFlow::setPulserFlags(RawAnitaHeader* header, UsefulAdu5Pat* usefulPat, AnitaEventSummary* sum){

  if(isPulserWAIS(header, usefulPat)){
    sum->flags.pulser = AnitaEventSummary::EventFlags::WAIS;
  }
  else if(isPulserLDB(header, usefulPat)){
    sum->flags.pulser = AnitaEventSummary::EventFlags::LDB;    
  }
  // std::cerr << __PRETTY_FUNCTION__ << ": I just set the flags to be " << sum->flags.pulser << std::endl;
}



void Acclaim::AnalysisFlow::prepareEverything(const char* preferredSettingsFileName){

  if(fRun >= 257 && fRun <= 263){
    if(AnitaVersion::get()==3){
      std::cerr << "No ANITA-3 data for run " << fRun << " so won't do analysis" << std::endl;
      return;
    }
  }
  
  if(!fSettings){
    fSettings = new AnalysisSettings(preferredSettingsFileName);
    fSettings->apply(dynamic_cast<TObject*>(this));
  }

  if(!fReco){
    fReco = new AnalysisReco();
    fSettings->apply(fReco);
  }

  if(!fData){ // do this after reading blinding strategy this->fBlindStrat from config file
    prepareDataSet();
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

  if(fUseNoiseMonitor){
    fNoiseMonitor = new NoiseMonitor(fFilterStrat);
  }

  fLastEventConsidered = 0;
}



AnitaEventSummary* Acclaim::AnalysisFlow::doEvent(UInt_t eventNumber){
  int entry = fData->getEvent(eventNumber);
  return entry > 0 ? doEntry(entry) : NULL;
}



AnitaEventSummary* Acclaim::AnalysisFlow::doEntry(Long64_t entry){

  fData->getEntry(entry);

  RawAnitaHeader* header = fData->header();
  UsefulAnitaEvent* usefulEvent = fData->useful();

  if(fDebug){
    std::cerr << "Debug in " << __PRETTY_FUNCTION__ << " doing entry " << entry << " for selection " << fSelection
              << ", header->eventNumber = " << header->eventNumber << ", usefulEvent->eventNumber = " << usefulEvent->eventNumber
              << std::endl;
  }

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
    fEventSummary = new AnitaEventSummary(header, &usefulPat, fData->truth());
    Bool_t isGoodEvent = QualityCut::applyAll(usefulEvent, fEventSummary);

    // varner2 contains a check that all channels have lots of points, too few can crash interpolation
    if(isGoodEvent || (fDoAll && !fEventSummary->flags.isVarner2)){

      setPulserFlags(header, &usefulPat, fEventSummary);
        
      // since we now have rolling averages make sure the filter strategy is sees every event before deciding whether or not to reconstruct
      if(fEv){
        delete fEv;
      }
      fEv = new FilteredAnitaEvent(usefulEvent, fFilterStrat, pat, header, false);
      fReco->process(fEv, fEventSummary, fNoiseMonitor);
    }

    if(fSumTree){
      fSumTree->Fill();
    }
  }
    
  fLastEventConsidered = header->eventNumber;
  return fEventSummary;

}



void Acclaim::AnalysisFlow::doAnalysis(Long64_t startAtThisEntry){

  fLastEventConsidered = 0;

  std::cout << "Info in " << __PRETTY_FUNCTION__ << ", doing run " << fRun
	    << " from entry " << fFirstEntry << " to " << fLastEntry << std::endl;

  if(startAtThisEntry >= fFirstEntry && startAtThisEntry < fLastEntry){
    std::cout << "Info in " << __PRETTY_FUNCTION__ << ", starting at entry " << startAtThisEntry
	      << " instead of entry " << fFirstEntry << std::endl;
    fFirstEntry = startAtThisEntry;
  }
  else if(startAtThisEntry!=-1){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << " got requested start entry "
	      << startAtThisEntry << ", which lies outside range " << fFirstEntry
	      << " to " << fLastEntry << ". Ignoring request!" << std::endl;
  }

  const Long64_t numEntries = fLastEntry-fFirstEntry;
  ProgressBar p(numEntries);
  for(Long64_t entryInd = 0; entryInd < numEntries; entryInd++){
    Long64_t entry = fFirstEntry + entryInd;
    AnitaEventSummary* sum = doEntry(entry);

    if(sum){
      delete sum;
    }
    p.inc(entryInd, numEntries);
  }
}
