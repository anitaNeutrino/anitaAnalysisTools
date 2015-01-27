#include "CrossCorrelator.h"

#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "RawAnitaHeader.h"

int main(){

  //  const int startRun = 11672;
  //  const int startRun = 126;
  //  const int startRun = 137;
  const int startRun = 158;
  //  const int endRun = 165; // when anita went out of LOS...
  //  const int endRun = 162; // when anita went out of LOS...
  const int endRun = 159; // when anita went out of LOS...
  //  const int endRun = 11673;   

  TChain* headChain = new TChain("headTree");
  TChain* eventChain = new TChain("eventTree");
  for(int run=startRun; run<endRun; run++){
    char eventFileName[1024];
    sprintf(eventFileName, "root/run%d/eventFile%d.root", run, run);
    eventChain->Add(eventFileName);

    char rawHeaderFileName[1024];
    sprintf(rawHeaderFileName, "root/run%d/eventHeadFile%d.root", run, run);
    headChain->Add(rawHeaderFileName);
  }

  RawAnitaEvent* event = NULL;
  eventChain->SetBranchAddress("event", &event);
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);

  const int nCuts = 1;
  double triggerTimeNsLow[nCuts] = {77.54e6}; //75.3e6, 25.2e6, 50.3e6};
  double triggerTimeNsHigh[nCuts] = {77.8e6}; //77.8e6, 26.8e6, 51.8e6};

  const int numSteps=10;

  TFile* outFile = new TFile("phaseCentrePlots.root","recreate");
  TTree* seaveyTree = new TTree("seaveyTree","seaveyTree");
  Double_t gpuPeakTheta;
  Double_t gpuPeakPhi;
  Double_t gpuPeak;
  Double_t imagePeakThetaV[numSteps];
  Double_t imagePeakPhiV[numSteps];
  Double_t imagePeakV[numSteps];
  Double_t imagePeakThetaH[numSteps];
  Double_t imagePeakPhiH[numSteps];
  Double_t imagePeakH[numSteps];

  Long64_t eventNumberHeader;
  seaveyTree->Branch("imagePeakThetaV", &imagePeakThetaV, "imagePeakThetaV[10]/D");
  seaveyTree->Branch("imagePeakPhiV", &imagePeakPhiV, "imagePeakPhiV[10]/D");
  seaveyTree->Branch("imagePeakV", &imagePeakV, "imagePeakV[10]/D");
  seaveyTree->Branch("imagePeakThetaH", &imagePeakThetaH, "imagePeakThetaH[10]/D");
  seaveyTree->Branch("imagePeakPhiH", &imagePeakPhiH, "imagePeakPhiH[10]/D");
  seaveyTree->Branch("imagePeakH", &imagePeakH, "imagePeakH[10]/D");
  seaveyTree->Branch("gpuPeak", &gpuPeak);
  seaveyTree->Branch("gpuPeakTheta", &gpuPeakTheta);
  seaveyTree->Branch("gpuPeakPhi", &gpuPeakPhi);
  seaveyTree->Branch("eventNumberHeader", &eventNumberHeader);

  Long64_t numEntries = headChain->GetEntries();

  CrossCorrelator* cc = new CrossCorrelator();
  cc->dlPerStep = 0.1; //metres
  cc->angleDeg = -10;


  eventChain->BuildIndex("event->eventNumber");
  int maxEvents = 10;

  int numEvents = 0;
  for(Long64_t entry=0; entry<numEntries; entry++){
    headChain->GetEntry(entry);
    eventChain->GetEntryWithIndex(header->eventNumber);

    if(event){

      int tt1e8 = header->triggerTimeNs%(int(1e8));
      int goodTriggerTimeNsFlag = 0;
      for(int cut=0; cut<nCuts; cut++){
  	if(tt1e8 > triggerTimeNsLow[cut] && tt1e8 < triggerTimeNsHigh[cut]){
  	  goodTriggerTimeNsFlag = 1;
  	}
      }

      if(goodTriggerTimeNsFlag){
  	if(event->eventNumber != header->eventNumber){
  	  std::cerr << "Couldn't find matching RawAnitaEvent for RawAnitaHeader "<< header->eventNumber << ". Skipping event. " << std::endl;
  	  continue;
  	}

  	UsefulAnitaEvent* realEvent = new UsefulAnitaEvent(event, WaveCalType::kVTBenS, header);
	cc->correlateEvent(realEvent);

	/* for debugging */
	//writeCorrelationGraphs(cc);

	for(int step=0; step<numSteps; step++){
	  cc->feedOffsetStep = step;
	  TH2D* hImageV = cc->makeImage(AnitaPol::kVertical, imagePeakV[step], 
					imagePeakThetaV[step], imagePeakPhiV[step]);
	  TH2D* hImageH = cc->makeImage(AnitaPol::kHorizontal, imagePeakH[step],
					imagePeakThetaH[step], imagePeakPhiH[step]);

	  char name[1024];
	  sprintf(name, "h%uV_step%d", event->eventNumber, step);
	  hImageV->SetName(name);
	  sprintf(name, "h%uH_step%d", event->eventNumber, step);
	  hImageH->SetName(name);

	  hImageV->Write();
	  hImageH->Write();

	  delete hImageV;
	  delete hImageH;

	} // step

	gpuPeak = header->getImagePeak();
	gpuPeakTheta = -1*header->getPeakThetaDeg();
	gpuPeakPhi = header->getPeakPhiDeg();
	eventNumberHeader = header->eventNumber;


  	delete realEvent;

  	seaveyTree->Fill();

	numEvents++;
	std::cout << "numEvents = " << numEvents << std::endl;
      }
    }
    else{
      std::cerr << "Can't find entry " << entry << " in eventTree" << std::endl;
    }

    if(numEvents==maxEvents){
      break;
    }

  }
  
  outFile->Write();
  outFile->Close();
  
  delete cc;
}
