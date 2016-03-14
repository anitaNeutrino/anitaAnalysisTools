/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             Program to test the functionality of the CrossCorrelator class, nothing more complicated than that. Writes some output to /tmp.

*************************************************************************************************************** */


#include "CrossCorrelator.h"
#include "RootTools.h"
#include "FFTtools.h"

#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include <RawAnitaHeader.h>
#include <UsefulAnitaEvent.h>
#include <CalibratedAnitaEvent.h>

//void testImageGPUStyle();
void testNewCombinatorics();
void testImageFullStyle();
void writeCorrelationGraphs(CrossCorrelator* cc);
void testCoherentlySummedWaveform();
// void testFileWriting();

int main(){

  //  testImageGPUStyle();
  //  testNewCombinatorics();
  // testImageFullStyle();
  testCoherentlySummedWaveform();  
  //  testFileWriting();

  return 0;

}

// void testFileWriting(){

//   CrossCorrelator* c = new CrossCorrelator();
//   c->writeDeltaTsFile();

//   c->readDeltaTsFile();
//   for(int combo=0; combo<NUM_COMBOS; combo++){
//     std::cout << (int)c->deltaTs[combo][0][0] << ", ";
//   }
//   std::cout << std::endl;
// }

void testCoherentlySummedWaveform(){
  /*
    This function tests the full +-2 phi-sector reconstruction.
   */

  char eventFileName[1024];

  Int_t run = 352;
  sprintf(eventFileName, "~/UCL/ANITA/calibratedFlight1415/run%d/calEventFile%d.root", run, run);

  TFile* eventFile = TFile::Open(eventFileName);
  TTree* eventTree = (TTree*) eventFile->Get("eventTree");

  char rawHeaderFileName[1024];
  sprintf(rawHeaderFileName, "~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);

  TFile* rawHeaderFile = TFile::Open(rawHeaderFileName);
  TTree* headTree = (TTree*) rawHeaderFile->Get("headTree");

  CalibratedAnitaEvent* event = NULL;  
  eventTree->SetBranchAddress("event", &event);

  RawAnitaHeader* header = NULL;
  headTree->SetBranchAddress("header", &header);

  TFile* outFile = new TFile("/tmp/testCoherentlySummedWaveform.root","recreate");
  TTree* corrTree = new TTree("corrTree","corrTree");
  Double_t imagePeak;
  Double_t imagePeakZoom;  
  corrTree->Branch("imagePeak", &imagePeak);

  CrossCorrelator* cc = new CrossCorrelator();

  headTree->BuildIndex("eventNumber");
  eventTree->BuildIndex("eventNumber");  
  
  imagePeak = -1;
  imagePeakZoom = -1;    
  Double_t peakPhiDeg = -1;
  Double_t peakThetaDeg = -1;
  Double_t peakPhiDegZoom = -1;
  Double_t peakThetaDegZoom = -1;

  const Long64_t eventNumber = 60832108;
  headTree->GetEntryWithIndex(eventNumber);
  eventTree->GetEntryWithIndex(eventNumber);
  std::cout << header->eventNumber << "\t" << event->eventNumber << std::endl;  

  
  // UsefulAnitaEvent* realEvent(new UsefulAnitaEvent(event, WaveCalType::kDefault,  header));
  UsefulAnitaEvent* realEvent = new UsefulAnitaEvent(event);
  cc->correlateEvent(realEvent);

  // writeCorrelationGraphs(cc);

  std::vector<TH2D*> hImages;
  std::vector<TH2D*> hZoomedImages;    

  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;
  
  UInt_t l3TrigPattern = 0;
  if(pol==AnitaPol::kHorizontal){
    l3TrigPattern = header->l3TrigPatternH;
  }
  else if(pol==AnitaPol::kVertical){
    l3TrigPattern = header->l3TrigPattern;
  }
  else{
    std::cerr << "??????????????????????" << std::endl;
  }
  std::cout << l3TrigPattern << std::endl;
      
  hImages.push_back(cc->makeTriggeredImage(pol, imagePeak, peakPhiDeg, peakThetaDeg, l3TrigPattern));
  hZoomedImages.push_back(cc->makeZoomedImage(pol, imagePeakZoom, peakPhiDegZoom, peakThetaDegZoom,
					      l3TrigPattern, peakPhiDeg, peakThetaDeg));

  Double_t snr = 0;
  TGraph* grCoherent = cc->makeCoherentlySummedWaveform(pol, peakPhiDegZoom, peakThetaDegZoom,
  							l3TrigPattern, snr);

  if(grCoherent){
    grCoherent->Write();
    TGraph* grHilbert = FFTtools::getHilbertEnvelope(grCoherent);


    grHilbert->SetName("grHilbert");
    grHilbert->Write();

    delete grCoherent;
    grCoherent = NULL;
    delete grHilbert;
    grHilbert = NULL;
  }

  
  for(Int_t phiSector=0; phiSector<NUM_PHI; phiSector++){
    Int_t doPhiSector = RootTools::getBit(phiSector, l3TrigPattern);
    if(doPhiSector > 0){
      for(int ring=0; ring<NUM_RING; ring++){
	int ant = phiSector + NUM_PHI*ring;
	TGraph* gr = (TGraph*) cc->grsResampled[pol][ant]->Clone();
	for(Int_t samp=0; samp<gr->GetN(); samp++){
	  gr->GetY()[samp]*=cc->interpRMS[pol][ant];
	}
	gr->SetName(TString::Format("grInterp_%d", ant));
	gr->Write();
      }
    }
  }
      
    
  delete realEvent;

  outFile->Write();
  outFile->Close();

  delete cc;

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

  Int_t entry = 107;

  CrossCorrelator* cc = new CrossCorrelator();

  headTree->GetEntry(entry);
  eventTree->GetEntry(entry);

  UsefulAnitaEvent* realEvent(new UsefulAnitaEvent(event, WaveCalType::kDefault,  header));
  cc->correlateEvent(realEvent);

  TFile* outFile = new TFile("/tmp/testNewCombinatorics.root","recreate");
  Double_t imagePeak = -1;
  Double_t peakPhiDeg = -1;
  Double_t peakThetaDeg = -1;
  
  // TH2D* hImage = cc->makeGlobalImage(AnitaPol::kVertical, imagePeak, peakPhiDeg, peakThetaDeg);
  TH2D* hImage = cc->makeTriggeredImage(AnitaPol::kVertical, imagePeak, peakPhiDeg, peakThetaDeg, header->l3TrigPatternH);  
  hImage->Write();
  delete hImage;

  outFile->Write();
  outFile->Close();

  delete cc;
  
}

void testImageFullStyle(){
  /*
    This function tests the full +-2 phi-sector reconstruction.
   */

  char eventFileName[1024];

  Int_t run = 352;
  sprintf(eventFileName, "~/UCL/ANITA/calibratedFlight1415/run%d/calEventFile%d.root", run, run);

  TFile* eventFile = TFile::Open(eventFileName);
  TTree* eventTree = (TTree*) eventFile->Get("eventTree");

  char rawHeaderFileName[1024];
  sprintf(rawHeaderFileName, "~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);

  TFile* rawHeaderFile = TFile::Open(rawHeaderFileName);
  TTree* headTree = (TTree*) rawHeaderFile->Get("headTree");

  CalibratedAnitaEvent* event = NULL;  
  eventTree->SetBranchAddress("event", &event);

  RawAnitaHeader* header = NULL;
  headTree->SetBranchAddress("header", &header);

  TFile* outFile = new TFile("/tmp/testImageFullStyle.root","recreate");
  TTree* corrTree = new TTree("corrTree","corrTree");
  Double_t imagePeak;
  corrTree->Branch("imagePeak", &imagePeak);

  //  Long64_t numEntries = eventTree->GetEntries();
  Long64_t numEntries = 200; //eventTree->GetEntries();
  numEntries = 408;

  CrossCorrelator* cc = new CrossCorrelator();

  for(Long64_t entry=407; entry<numEntries; entry++){

    imagePeak = -1;
    Double_t peakPhiDeg = -1;
    Double_t peakThetaDeg = -1;

    headTree->GetEntry(entry);
    eventTree->GetEntry(entry);

    // UsefulAnitaEvent* realEvent(new UsefulAnitaEvent(event, WaveCalType::kDefault,  header));
    UsefulAnitaEvent* realEvent = new UsefulAnitaEvent(event);
    cc->correlateEvent(realEvent);

    writeCorrelationGraphs(cc);

    std::vector<TH2D*> hImages;
    for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
      hImages.push_back(cc->makeGlobalImage(AnitaPol::AnitaPol_t(pol), imagePeak, peakPhiDeg, peakThetaDeg));
    }

    // hImages.push_back(cc->makeTriggeredImage(AnitaPol::kHorizontal, imagePeak, peakPhiDeg, peakThetaDeg,
    // 					     header->l3TrigPatternH));  
    // hImages.push_back(cc->makeTriggeredImage(AnitaPol::kVertical, imagePeak, peakPhiDeg, peakThetaDeg,
    // 					     header->l3TrigPattern));  

    
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
