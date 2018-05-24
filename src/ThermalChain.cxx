#include "ThermalChain.h"
#include "SummarySet.h"
#include "TChain.h"
#include "TEntryList.h"
#include "ProgressBar.h"
#include "TROOT.h"
#include "AnitaDataset.h"


Acclaim::ThermalChain::ThermalChain(const char* glob, const char* treeName){

  fChain = new TChain(treeName);
  fChain->Add(glob);

  TString hiCalGlob(glob);
  hiCalGlob.ReplaceAll("makeThermalTree", "makeHiCalTree");

  fFriendChain1 = new TChain("hiCalTree");
  fFriendChain1->Add(hiCalGlob);

  TString surfaceGlob(glob);
  surfaceGlob.ReplaceAll("makeThermalTree", "makeSurfaceTree");

  fFriendChain2 = new TChain("surfaceTree");
  fFriendChain2->Add(surfaceGlob);

  // std::cout << fFriendChain2->GetEntries() << std::endl;
  fChain->AddFriend(fFriendChain1);
  fChain->AddFriend(fFriendChain2);

  gROOT->ProcessLine("#include \"FFTtools.h\""); // hack to get various delta phi wrap Draw things to work in stand alone executables

  fCut = "";  
  fEntryListDirty = true;
  fEntryList = NULL;
  fUseProof = false;
  fAnitaVersion = 0;
  setBranches();
}

Acclaim::ThermalChain::~ThermalChain(){

  if(fEntryList){
    fChain->SetEntryList(0);
    delete fEntryList;
    fEntryList = NULL;
  }

  // Pretty sure ROOT takes care of this...
  fFriendChain1 = NULL;
  fFriendChain2 = NULL;
  // if(fFriendChain){
  //   delete fFriendChain;
  //   fFriendChain = NULL;
  // }  
  
  if(fChain){
    delete fChain;
    fChain = NULL;
  }
}


void Acclaim::ThermalChain::setBranches(){

  // these ones need a type conversion before public consumption...
  fChain->SetBranchAddress("pol", &polFloat);
  fChain->SetBranchAddress("peakInd", &peakIndFloat);
  fChain->SetBranchAddress("eventNumber", &eventNumberInt);
  fChain->SetBranchAddress("realTime", &realTimeInt);
  
  fChain->SetBranchAddress("run", &run);
  fChain->SetBranchAddress("peak_phi", &peak_phi);
  fChain->SetBranchAddress("peak_theta", &peak_theta);
  fChain->SetBranchAddress("anitaLocation_longitude", &anita_longitude);
  fChain->SetBranchAddress("anitaLocation_latitude", &anita_latitude);
  fChain->SetBranchAddress("anitaLocation_altitude", &anita_altitude);
  fChain->SetBranchAddress("anitaLocation_heading", &anita_heading);
  fChain->SetBranchAddress("coherent_filtered_snr", &coherent_filtered_snr);
  fChain->SetBranchAddress("weight", &weight);
  fChain->SetBranchAddress("mc_energy", &mc_energy);

  fChain->SetBranchAddress("peak_value", &peak_value);
  fChain->SetBranchAddress("coherent_filtered_peakHilbert", &coherent_filtered_peakHilbert);
  fChain->SetBranchAddress("deconvolved_filtered_peakHilbert", &deconvolved_filtered_peakHilbert);
  fChain->SetBranchAddress("coherent_filtered_impulsivityMeasure", &coherent_filtered_impulsivityMeasure);  
  fChain->SetBranchAddress("deconvolved_filtered_impulsivityMeasure", &deconvolved_filtered_impulsivityMeasure);  
  fChain->SetBranchAddress("coherent_filtered_fracPowerWindowGradient", &coherent_filtered_fracPowerWindowGradient);
  fChain->SetBranchAddress("deconvolved_filtered_fracPowerWindowGradient", &deconvolved_filtered_fracPowerWindowGradient);


  fFriendChain1->SetBranchAddress("duringHiCal", &duringHiCal);
  fFriendChain1->SetBranchAddress("hiCalPhi", &hiCalPhi);
  fFriendChain1->SetBranchAddress("hiCalTheta", &hiCalTheta);
  fFriendChain1->SetBranchAddress("eventNumber2", &eventNumber2);
  fFriendChain1->SetBranchAddress("run2", &run2);


  fFriendChain2->SetBranchAddress("longitude", &longitude);
  fFriendChain2->SetBranchAddress("latitude", &latitude);
  fFriendChain2->SetBranchAddress("altitude", &altitude);
  fFriendChain2->SetBranchAddress("thetaAdjustmentRequired", &thetaAdjustmentRequired);
  fFriendChain2->SetBranchAddress("onContinent", &onContinent);
  fFriendChain2->SetBranchAddress("onIceShelf", &onIceShelf);
  fFriendChain2->SetBranchAddress("iceThickness", &iceThickness);
  fFriendChain2->SetBranchAddress("eventNumber3", &eventNumber3);
  fFriendChain2->SetBranchAddress("run3", &run3);
  
}





/** 
 * Overwrites the current cut with a new cut
 * 
 * @param cut the new cut to apply to the chain of thermalTrees
 */
void Acclaim::ThermalChain::setCut(const TCut& cut){
  fCut = "";
  addCut(cut);
}

/** 
 * Overwrites the current cut with a new cut
 * 
 * @param cut the new cut to apply to the chain of thermalTrees
 */
void Acclaim::ThermalChain::setCut(const char* cut){
  fCut = "";
  addCut(cut);
}


/** 
 * Adds to the current cut with a new cut
 * (invalidates the cached entrylist)
 * 
 * @param cut the new cut to add to the cuts already applied to the chain of thermalTrees
 */
void Acclaim::ThermalChain::addCut(const char* cut){
  TCut cut2(cut);
  addCut(cut2);
}


/** 
 * Adds to the current cut with a new cut
 * (invalidates the cached entrylist)
 * 
 * @param cut the new cut to add to the cuts already applied to the chain of thermalTrees
 */
void Acclaim::ThermalChain::addCut(const TCut& cut){
  TCut oldCut;
  fCut += cut;
  if(fCut != oldCut){
    fEntryListDirty = true;
  }
}


void Acclaim::ThermalChain::SetUseProof(bool useProof){
  if(useProof){
    SummarySet::startProof();
  }
  fChain->SetProof(useProof);
  fUseProof = useProof;
}


/** 
 * Work horse function to update the entrylist 
 * 
 */
void Acclaim::ThermalChain::makeSelection() const {
  if(fEntryListDirty){

    if(fEntryList){
      fChain->SetEntryList(0);
      delete fEntryList;
      fEntryList = NULL;
    }

    fChain->SetProof(false);    

    std::cout << "Info in " << __PRETTY_FUNCTION__ << ", updating fEntryList..." << std::endl;
    ProgressBar p(1);
    fChain->Draw(">>fEntryList", fCut, "entrylist");
    fEntryList = dynamic_cast<TEntryList*>(gROOT->FindObject("fEntryList"));
    if(!fEntryList){
      std::cerr << "Error! couldn't find fEntryList!" << std::endl;
    }
    else{
      fEntryListDirty = false;
      fChain->SetEntryList(fEntryList);      
    }

    // Turn proof back on, if it was on...
    fChain->SetProof(fUseProof);
    
    p++;
  }
}


/** 
 * How many tree entries match the current selection?
 * @return The number of entries in the chain matching the current selection
 */
Long64_t Acclaim::ThermalChain::N() const {
  makeSelection();
  return fEntryList->GetN();
}

Long64_t Acclaim::ThermalChain::getEntry(Long64_t entry){
  
  makeSelection();
  Int_t treeIndex = -1;
  Long64_t treeEntry = fEntryList->GetEntryAndTree(entry,treeIndex);
  Long64_t chainEntry = treeEntry+fChain->GetTreeOffset()[treeIndex];

  Long64_t retVal = fChain->GetEntry(chainEntry);  
  
  doTypeConversions();

  return retVal;
}

void Acclaim::ThermalChain::doTypeConversions(){
  pol = (AnitaPol::AnitaPol_t) polFloat;
  peakInd = (Int_t) peakIndFloat;
  eventNumber = (UInt_t) eventNumberInt;
  realTime = (UInt_t) realTimeInt;

  if(eventNumber2 != eventNumber || run2 != run){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", mismatch in friendChain1!" << std::endl;
  }
  if(eventNumber3 != eventNumber || run3 != run){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", mismatch in friendChain2!" << std::endl;
  }
}




Long64_t Acclaim::ThermalChain::getEvent(UInt_t eventNumber){
  // Rolled my own lookup-by-eventNumber since TChainIndex doesn't seem up to the task
  // maybe because the eventNumber marker isn't as well behaved as we would like?

  // First set the correct AnitaVersion in case I compile with the default to ANITA-4  
  if(!fAnitaVersion){ 
    fChain->GetEntry(0);
    doTypeConversions();
    AnitaVersion::setVersionFromUnixTime(realTime);
    fAnitaVersion = AnitaVersion::get();
  }

  // Look up the run, and find the appropriate file
  Int_t run = AnitaDataset::getRunContainingEventNumber(eventNumber);
  TString desiredFileName = TString::Format("makeThermalTree_%d_", run);
  TIter next(fChain->GetListOfFiles());
  Int_t desiredTreeNumber = -1;
  while(TObject* f = next()){
    desiredTreeNumber++;    
    TString fileName = f->GetTitle();
    if(fileName.Contains(desiredFileName)){
      break;
    }
  }


  // Check boundaries
  if(desiredTreeNumber > -1 && desiredTreeNumber < fChain->GetListOfFiles()->GetEntries()){

    // Get the first entry in file corresponding to run
    Long64_t entry = fChain->GetTreeOffset()[desiredTreeNumber];
    fChain->GetEntry(entry);
    doTypeConversions();

    // If event number was perfectly monotonic with no gaps this would work first time
    // But as it is, need multiple attempts
    const int maxTries = 50; // (Allow a lot of attempts)
    for(int numTries = 0; numTries < maxTries; numTries++){

      // Update entry difference in eventNumber from last event and desired event
      Int_t deltaEntries = eventNumber - this->eventNumber;
      entry += deltaEntries;

      // Get most recent guessed entry
      Long64_t nb = fChain->GetEntry(entry); 
      doTypeConversions();
      
      // If it's a match then we're done, if not then loop again
      if(this->eventNumber == eventNumber){
	// std::cout << this->eventNumber << "\t" << eventNumber << std::endl;
	return nb;
	break;
      }

    }
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", unable to find eventNumber "
	      << eventNumber << " after " << maxTries << " attempts, giving up!" << std::endl;
    return -1;
  }
  std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", unable to find file corresponding to run " << run << std::endl;
  return -1;
}


Adu5Pat Acclaim::ThermalChain::pat(){
  Adu5Pat pat;
  pat.latitude = anita_latitude;
  pat.longitude = anita_longitude;
  pat.altitude = anita_altitude;
  pat.heading = anita_heading;
  pat.realTime = realTime;
  return pat;
}



Double_t Acclaim::ThermalChain::fisherDiscriminant(){
  //////////////////////////////////////////
  // THIS NEEDS TO BE MAINTAINED BY HAND! //
  //////////////////////////////////////////
  
  // from DrawStrings.h
  // const TString fisherDiscriminant = "0.898497+(1.929594*coherent_filtered_fracPowerWindowGradient/deconvolved_filtered_fracPowerWindowGradient)+(-0.195909*deconvolved_filtered_fracPowerWindowGradient)+(5.943355*coherent_filtered_impulsivityMeasure)+(0.826114*deconvolved_filtered_impulsivityMeasure)+(0.021763*coherent_filtered_peakHilbert)+(-0.012670*deconvolved_filtered_peakHilbert)+(-0.394201*peak_value)";

  return 0.898497+(1.929594*coherent_filtered_fracPowerWindowGradient/deconvolved_filtered_fracPowerWindowGradient)+(-0.195909*deconvolved_filtered_fracPowerWindowGradient)+(5.943355*coherent_filtered_impulsivityMeasure)+(0.826114*deconvolved_filtered_impulsivityMeasure)+(0.021763*coherent_filtered_peakHilbert)+(-0.012670*deconvolved_filtered_peakHilbert)+(-0.394201*peak_value);
}
