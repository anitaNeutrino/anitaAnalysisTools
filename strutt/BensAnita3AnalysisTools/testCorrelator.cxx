/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             Program to test the functionality of the CrossCorrelator class, nothing more complicated than that.
*************************************************************************************************************** */


#include "CrossCorrelator.h"

#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "RawAnitaHeader.h"

void testImageGPUStyle();
void testNewCombinatorics();
void testImageFullStyle();
void writeCorrelationGraphs(CrossCorrelator* cc);

int main(){

  //  testImageGPUStyle();
  //  testNewCombinatorics();
  testImageFullStyle();
  return 0;

}


void testNewCombinatorics(){
  /*
    This function tests the combinatorics for cross-correlating channels within +- 2 phi-sectors.
   */

  char eventFileName[1024];

  Int_t run = 11672;
  sprintf(eventFileName, "../antarctica2014/PrioritizerdCalib/localData/run%d/eventFile%d.root", run, run);
  TFile* eventFile = TFile::Open(eventFileName);
  TTree* eventTree = (TTree*) eventFile->Get("eventTree");

  char rawHeaderFileName[1024];
  sprintf(rawHeaderFileName, "../antarctica2014/PrioritizerdCalib/localData/run%d/headFile%d.root", run, run);
  TFile* rawHeaderFile = TFile::Open(rawHeaderFileName);
  TTree* headTree = (TTree*) rawHeaderFile->Get("headTree");

  RawAnitaEvent* event = NULL;
  eventTree->SetBranchAddress("event", &event);
  RawAnitaHeader* header = NULL;
  headTree->SetBranchAddress("header", &header);

  Long64_t numEntries = 200; //eventTree->GetEntries();
  Int_t entry = 107;

  CrossCorrelator* cc = new CrossCorrelator();

  headTree->GetEntry(entry);
  eventTree->GetEntry(entry);

  UsefulAnitaEvent* realEvent(new UsefulAnitaEvent(event, WaveCalType::kVTBenS, header));
  cc->correlateEvent(realEvent);
  Int_t nc = -1;
  for(Int_t ant1=0; ant1<NUM_SEAVEYS; ant1++){
    for(UInt_t ant2Ind; ant2Ind < cc->ant2s[ant1].size(); ant2Ind++){
      Int_t ant2 = cc->ant2s[ant1].at(ant2Ind);
      //    for(int& ant2 : cc->ant2s[ant1]){
      std::cout << ant1 << " " << ant2 << " "
		<< cc->comboIndices[ant1][ant2] << "  " 
		<< cc->comboIndices[ant2][ant1] << std::endl;
      if(cc->comboIndices[ant1][ant2]>nc){
	nc = cc->comboIndices[ant1][ant2];
      }
    }
  }
  std::cout << nc << " combos" << std::endl; 

  TFile* outFile = new TFile("/tmp/testNewCombinatorics.root","recreate");
  TH2D* hImage = cc->makeImage(AnitaPol::kVertical);
  outFile->Write();
  outFile->Close();

  delete cc;
  
}

void testImageGPUStyle(){
  /*
    This function tests the cross-correlation functions which mimic what the GPU is doing on some level.
   */

  char eventFileName[1024];

  Int_t run = 11672;
  sprintf(eventFileName, "../antarctica2014/PrioritizerdCalib/localData/run%d/eventFile%d.root", run, run);
  TFile* eventFile = TFile::Open(eventFileName);
  TTree* eventTree = (TTree*) eventFile->Get("eventTree");

  char rawHeaderFileName[1024];
  sprintf(rawHeaderFileName, "../antarctica2014/PrioritizerdCalib/localData/run%d/headFile%d.root", run, run);
  TFile* rawHeaderFile = TFile::Open(rawHeaderFileName);
  TTree* headTree = (TTree*) rawHeaderFile->Get("headTree");

  RawAnitaEvent* event = NULL;
  eventTree->SetBranchAddress("event", &event);
  RawAnitaHeader* header = NULL;
  headTree->SetBranchAddress("header", &header);

  TFile* outFile = new TFile("/tmp/testImageGPUStyle.root","recreate");
  TTree* corrTree = new TTree("corrTree","corrTree");
  Double_t imagePeak;
  corrTree->Branch("imagePeak", &imagePeak);

  //  Long64_t numEntries = eventTree->GetEntries();
  Long64_t numEntries = 200; //eventTree->GetEntries();
  numEntries = 108;

  CrossCorrelator* cc = new CrossCorrelator();

  for(Long64_t entry=107; entry<numEntries; entry++){

    imagePeak = -1;

    headTree->GetEntry(entry);
    eventTree->GetEntry(entry);

    UsefulAnitaEvent* realEvent(new UsefulAnitaEvent(event, WaveCalType::kVTBenS, header));
    cc->correlateEventGPU(realEvent);
    TH2D* hImage = cc->makeImageGPU(AnitaPol::kVertical);
    hImage->SetName("hImageFrom3PhiSectorPulsing");
    hImage->SetTitle("Image From 3 Phi-Sector Pulsing");

    delete realEvent;
  }

  outFile->Write();
  outFile->Close();

  delete cc;

}


void testImageFullStyle(){
  /*
    This function tests the full +-2 phi-sector reconstruction.
   */

  char eventFileName[1024];

  Int_t run = 11672;
  sprintf(eventFileName, "../antarctica2014/PrioritizerdCalib/localData/run%d/eventFile%d.root", run, run);
  TFile* eventFile = TFile::Open(eventFileName);
  TTree* eventTree = (TTree*) eventFile->Get("eventTree");

  char rawHeaderFileName[1024];
  sprintf(rawHeaderFileName, "../antarctica2014/PrioritizerdCalib/localData/run%d/headFile%d.root", run, run);
  TFile* rawHeaderFile = TFile::Open(rawHeaderFileName);
  TTree* headTree = (TTree*) rawHeaderFile->Get("headTree");

  RawAnitaEvent* event = NULL;
  eventTree->SetBranchAddress("event", &event);
  RawAnitaHeader* header = NULL;
  headTree->SetBranchAddress("header", &header);

  TFile* outFile = new TFile("/tmp/testImageFullStyle.root","recreate");
  TTree* corrTree = new TTree("corrTree","corrTree");
  Double_t imagePeak;
  corrTree->Branch("imagePeak", &imagePeak);

  //  Long64_t numEntries = eventTree->GetEntries();
  Long64_t numEntries = 200; //eventTree->GetEntries();
  numEntries = 108;

  CrossCorrelator* cc = new CrossCorrelator();

  for(Long64_t entry=107; entry<numEntries; entry++){

    imagePeak = -1;

    headTree->GetEntry(entry);
    eventTree->GetEntry(entry);

    UsefulAnitaEvent* realEvent(new UsefulAnitaEvent(event, WaveCalType::kVTBenS, header));
    cc->correlateEvent(realEvent);

    writeCorrelationGraphs(cc);

    TH2D* hImage = cc->makeImage(AnitaPol::kVertical);
    hImage->SetName("hImageFrom3PhiSectorPulsing");
    hImage->SetTitle("Image From 3 Phi-Sector Pulsing, correlating 5 phi-sectors");

    delete realEvent;
  }

  outFile->Write();
  outFile->Close();

  delete cc;

}

void writeCorrelationGraphs(CrossCorrelator* cc){
  /* for debugging */

  for(Int_t ant1=0; ant1<NUM_SEAVEYS; ant1++){
    char name[1024];
    sprintf(name, "grsInterp_%d", ant1);
    // cc->grsInterp[AnitaPol::kVertical][ant1]->SetName(name);
    // cc->grsInterp[AnitaPol::kVertical][ant1]->SetTitle(name);
    // cc->grsInterp[AnitaPol::kVertical][ant1]->Write();
    for(Int_t ant2=0; ant2<NUM_SEAVEYS; ant2++){
      TGraph* grCorr = cc->getCrossCorrelationGraph(AnitaPol::kVertical, ant1, ant2);
      if(grCorr){
	grCorr->Write();
      }
    }
  }
}
