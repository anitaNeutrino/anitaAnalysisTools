/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             Program to test the functionality of the CrossCorrelator class, nothing more complicated than that. Writes some output to /tmp.

*************************************************************************************************************** */


#include <CrossCorrelator.h>
#include <RootTools.h>

#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

#include <RawAnitaHeader.h>
#include <UsefulAnitaEvent.h>
#include <CalibratedAnitaEvent.h>

//void testImageGPUStyle();
void testNewCombinatorics();
void testImageFullStyle();
void writeCorrelationGraphs(CrossCorrelator* cc);
void hackyNormalizationTest();
// void testFileWriting();

int main(){

  //  testImageGPUStyle();
  //  testNewCombinatorics();
  testImageFullStyle();
  // hackyNormalizationTest();
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
  Int_t nc = -1;
  for(Int_t ant1=0; ant1<NUM_SEAVEYS; ant1++){
    for(UInt_t ant2Ind=0; ant2Ind < cc->ant2s[ant1].size(); ant2Ind++){
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
  hImage->Write();
  delete hImage;

  outFile->Write();
  outFile->Close();

  delete cc;
  
}

// void testImageGPUStyle(){
//   /*
//     This function tests the cross-correlation functions which mimic what the GPU is doing on some level.
//    */

//   char eventFileName[1024];

//   Int_t run = 11672;
//   sprintf(eventFileName, "../antarctica2014/PrioritizerdCalib/localData/run%d/eventFile%d.root", run, run);
//   TFile* eventFile = TFile::Open(eventFileName);
//   TTree* eventTree = (TTree*) eventFile->Get("eventTree");

//   char rawHeaderFileName[1024];
//   sprintf(rawHeaderFileName, "../antarctica2014/PrioritizerdCalib/localData/run%d/headFile%d.root", run, run);
//   TFile* rawHeaderFile = TFile::Open(rawHeaderFileName);
//   TTree* headTree = (TTree*) rawHeaderFile->Get("headTree");

//   RawAnitaEvent* event = NULL;
//   eventTree->SetBranchAddress("event", &event);
//   RawAnitaHeader* header = NULL;
//   headTree->SetBranchAddress("header", &header);

//   TFile* outFile = new TFile("/tmp/testImageGPUStyle.root","recreate");
//   TTree* corrTree = new TTree("corrTree","corrTree");
//   Double_t imagePeak;
//   corrTree->Branch("imagePeak", &imagePeak);

//   //  Long64_t numEntries = eventTree->GetEntries();
//   Long64_t numEntries = 108; //200; //eventTree->GetEntries();

//   CrossCorrelator* cc = new CrossCorrelator();

//   for(Long64_t entry=107; entry<numEntries; entry++){

//     imagePeak = -1;

//     headTree->GetEntry(entry);
//     eventTree->GetEntry(entry);

//     UsefulAnitaEvent* realEvent(new UsefulAnitaEvent(event, WaveCalType::kDefault,  header));
//     cc->correlateEventGPU(realEvent);
//     TH2D* hImage = cc->makeImageGPU(AnitaPol::kVertical);
//     hImage->SetName("hImageFrom3PhiSectorPulsing");
//     hImage->SetTitle("Image From 3 Phi-Sector Pulsing");


//     delete realEvent;
//   }

//   outFile->Write();
//   outFile->Close();

//   delete cc;

// }


void testImageFullStyle(){
  /*
    This function tests the full +-2 phi-sector reconstruction.
   */

  char eventFileName[1024];

  // Int_t run = 11672;
  // sprintf(eventFileName, "~/UCL/ANITA/antarctica2014/PrioritizerdCalib/localData/run%d/eventFile%d.root", run, run);
  Int_t run = 352;
  sprintf(eventFileName, "~/UCL/ANITA/calibratedFlight1415/run%d/calEventFile%d.root", run, run);

  TFile* eventFile = TFile::Open(eventFileName);
  TTree* eventTree = (TTree*) eventFile->Get("eventTree");

  // char rawHeaderFileName[1024];
  // sprintf(rawHeaderFileName, "~/UCL/ANITA/antarctica2014/PrioritizerdCalib/localData/run%d/headFile%d.root", run, run);
  char rawHeaderFileName[1024];
  sprintf(rawHeaderFileName, "~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);

  TFile* rawHeaderFile = TFile::Open(rawHeaderFileName);
  TTree* headTree = (TTree*) rawHeaderFile->Get("headTree");

  // RawAnitaEvent* event = NULL;
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
  numEntries = 108;

  CrossCorrelator* cc = new CrossCorrelator(2);

  for(Long64_t entry=107; entry<numEntries; entry++){

    imagePeak = -1;

    headTree->GetEntry(entry);
    eventTree->GetEntry(entry);

    // UsefulAnitaEvent* realEvent(new UsefulAnitaEvent(event, WaveCalType::kDefault,  header));
    UsefulAnitaEvent* realEvent = new UsefulAnitaEvent(event);
    cc->correlateEvent(realEvent);

    writeCorrelationGraphs(cc);

    std::vector<TH2D*> hImages;
    for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
      hImages.push_back(cc->makeImage(AnitaPol::AnitaPol_t(pol)));
    }
    
    delete realEvent;
  }

  outFile->Write();
  outFile->Close();

  delete cc;

}



void hackyNormalizationTest(){
  /*
    Check FFT normalization
   */

  char eventFileName[1024];

  Int_t run = 11672;
  sprintf(eventFileName, "~/UCL/ANITA/antarctica2014/PrioritizerdCalib/localData/run%d/eventFile%d.root", run, run);
  TFile* eventFile = TFile::Open(eventFileName);
  TTree* eventTree = (TTree*) eventFile->Get("eventTree");

  char rawHeaderFileName[1024];
  sprintf(rawHeaderFileName, "~/UCL/ANITA/antarctica2014/PrioritizerdCalib/localData/run%d/headFile%d.root", run, run);
  TFile* rawHeaderFile = TFile::Open(rawHeaderFileName);
  TTree* headTree = (TTree*) rawHeaderFile->Get("headTree");

  RawAnitaEvent* event = NULL;
  eventTree->SetBranchAddress("event", &event);
  RawAnitaHeader* header = NULL;
  headTree->SetBranchAddress("header", &header);

  CrossCorrelator* cc = new CrossCorrelator();

  Long64_t numEntries = 108;
  for(Long64_t entry=107; entry<numEntries; entry++){

    headTree->GetEntry(entry);
    eventTree->GetEntry(entry);

    UsefulAnitaEvent* realEvent(new UsefulAnitaEvent(event, WaveCalType::kDefault, header));
    TGraph* gr1 = realEvent->getGraph(15, AnitaPol::kVertical);
    TGraph* gr2 = realEvent->getGraph(15, AnitaPol::kVertical);

    TGraph* gr1Int = cc->interpolateWithStartTime(gr1, gr1->GetX()[0]);
    TGraph* gr2Int = cc->interpolateWithStartTime(gr2, gr2->GetX()[0]);
    RootTools::normalize(gr1Int);
    RootTools::normalize(gr2Int);

    // cc->correlateEvent(realEvent);
    Double_t* corrs = cc->crossCorrelateFourier(gr1Int, gr2Int);
    std::cout << "Normalization factor is " << corrs[0] << std::endl;
    
    delete corrs;

    delete realEvent;
  }

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
