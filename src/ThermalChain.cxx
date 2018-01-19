#include "ThermalChain.h"
#include "SummarySet.h"
#include "TChain.h"
#include "TEntryList.h"
#include "ProgressBar.h"
#include "TROOT.h"

Acclaim::ThermalChain::ThermalChain(const char* glob, const char* treeName){

  fChain = new TChain(treeName);
  fChain->Add(glob);

  TString hiCalGlob(glob);
  hiCalGlob.ReplaceAll("makeThermalTree", "makeHiCalTree");
  fFriendChain = new TChain("hiCalTree");
  fFriendChain->Add(hiCalGlob);
  std::cout << fFriendChain->GetEntries() << std::endl;
  fChain->AddFriend(fFriendChain);

  
  gROOT->ProcessLine("#include \"FFTtools.h\""); // hack to get various delta phi wrap Draw things to work in stand alone executables

  fCut = "";  
  fEntryListDirty = true;
  fEntryList = NULL;
  fUseProof = false;
  setBranches();
}

Acclaim::ThermalChain::~ThermalChain(){

  if(fEntryList){
    fChain->SetEntryList(0);
    delete fEntryList;
    fEntryList = NULL;
  }

  // Pretty sure ROOT takes care of this...
  fFriendChain = NULL;
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



  fFriendChain->SetBranchAddress("duringHiCal", &duringHiCal);
  fFriendChain->SetBranchAddress("hiCalPhi", &hiCalPhi);
  fFriendChain->SetBranchAddress("hiCalTheta", &hiCalTheta);
  fFriendChain->SetBranchAddress("eventNumber2", &eventNumber2);
  fFriendChain->SetBranchAddress("run2", &run2);
  
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
  pol = (AnitaPol::AnitaPol_t) polFloat;
  peakInd = (Int_t) peakIndFloat;
  eventNumber = (UInt_t) eventNumberInt;
  realTime = (UInt_t) realTimeInt;

  if(eventNumber2 != eventNumber || run2 != run){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", mismatch in chains!" << std::endl;    
  }
  
  return retVal;
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

