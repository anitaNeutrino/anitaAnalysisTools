#include "ThermalChain.h"
#include "SummarySet.h"
#include "TChain.h"
#include "TEntryList.h"
#include "ProgressBar.h"
#include "TROOT.h"

Acclaim::ThermalChain::ThermalChain(const char* glob, const char* treeName){

  fChain = new TChain(treeName);
  fChain->Add(glob);

  fCut = "";  
  fEntryListDirty = true;
  fEntryList = NULL;
}

Acclaim::ThermalChain::~ThermalChain(){

  if(fEntryList){
    fChain->SetEntryList(0);
    delete fEntryList;
    fEntryList = NULL;
  }
  
  if(fChain){
    delete fChain;
    fChain = NULL;
  }
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
  if(fCut != cut){
    fEntryListDirty = true;
  }
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
  return fChain->GetEntry(chainEntry);
}
