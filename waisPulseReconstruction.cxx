#include "CrossCorrelator.h"
#include "FancyTTreeInterpolator.h"

#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "RawAnitaHeader.h"
#include "UsefulAdu5Pat.h"

Double_t getRMS(TGraph* gr);

int main(){

  //  const Int_t startRun = 11672;
  //  const Int_t startRun = 126;
  //  const Int_t startRun = 137;
  //  const Int_t startRun = 158;
  const Int_t startRun = 331;
  //  const Int_t endRun = 165; // when anita went out of LOS...
  //  const Int_t endRun = 162; // when anita went out of LOS...
  //  const Int_t endRun = 159; // when anita went out of LOS...
  const Int_t endRun = 355; // when anita went out of LOS...
  //  const Int_t endRun = 11673;   

  TChain* headChain = new TChain("headTree");
  TChain* eventChain = new TChain("eventTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  for(Int_t run=startRun; run<endRun; run++){
    char eventFileName[1024];
    sprintf(eventFileName, "root/run%d/eventFile%d.root", run, run);
    eventChain->Add(eventFileName);

    char rawHeaderFileName[1024];
    sprintf(rawHeaderFileName, "root/run%d/eventHeadFile%d.root", run, run);
    headChain->Add(rawHeaderFileName);


    char gpsFileName[1024];
    sprintf(gpsFileName, "root/run%d/gpsFile%d.root", run, run);
    gpsChain->Add(gpsFileName);
  }

  RawAnitaEvent* event = NULL;
  eventChain->SetBranchAddress("event", &event);
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);

  FancyTTreeInterpolator gpsInterp(gpsChain, "pat->realTime");
  gpsInterp.add("pat->heading", "pat->heading>-500", 360.0);
  gpsInterp.add("pat->altitude");
  gpsInterp.add("pat->latitude");
  gpsInterp.add("pat->longitude");

  CrossCorrelator* cc = new CrossCorrelator();
  cc->dlPerStep = 0.1; //metres
  cc->angleDeg = -10;
  const Int_t numSteps= cc->numSteps;
  const Double_t minRMS = 1e-6;

  TFile* outFile = new TFile("waisPulseReconstructionPlots.root","recreate");
  gpsInterp.get("pat->heading")->SetName("grHead");
  gpsInterp.get("pat->altitude")->SetName("grAlt");
  gpsInterp.get("pat->latitude")->SetName("grLat");
  gpsInterp.get("pat->longitude")->SetName("grLong");
  
  std::shared_ptr<TGraph> gr = gpsInterp.get("pat->heading");
  
  for(int i=1; i<gr->GetN(); i++){
    double dx = gr->GetX()[i] - gr->GetX()[i-1];
    if(dx < 0 ){
      std::cout << i << " " << gr->GetX()[i] - gr->GetX()[i-1] << " " 
		<< gr->GetX()[i] << " " <<  gr->GetX()[i-1] << std::endl;
    }
  }
  gpsInterp.get("pat->heading")->Write();
  gpsInterp.get("pat->altitude")->Write();
  gpsInterp.get("pat->latitude")->Write();
  gpsInterp.get("pat->longitude")->Write();

  TTree* waisTree = new TTree("waisTree","waisTree");
  Double_t gpuPeakTheta;
  Double_t gpuPeakPhi;
  Double_t gpuPeak;
  Double_t thetaExpected = -1;
  Double_t phiExpected = -1;
  // Double_t imagePeakThetaV; //[numSteps];
  // Double_t imagePeakPhiV; //[numSteps];
  // Double_t imagePeakV; //[numSteps];
  Double_t imagePeakThetaH; //[numSteps];
  Double_t imagePeakPhiH; //[numSteps];
  Double_t imagePeakH; //[numSteps];
  TH2D* hImageH;// = cc->makeImageGPU(AnitaPol::kHorizontal);

  // Int_t droppedPacketsV[16][3] = {{0}};
  Int_t droppedPacketsH[16][3] = {{0}};

  Long64_t eventNumberHeader;
  // waisTree->Branch("imagePeakThetaV", &imagePeakThetaV); //imagePeakThetaV[10]/D");
  // waisTree->Branch("imagePeakPhiV", &imagePeakPhiV); //imagePeakPhiV[10]/D");
  // waisTree->Branch("imagePeakV", &imagePeakV); //imagePeakV[10]/D");

  // waisTree->Branch("hImageH", "TH2D", &hImageH, 32000, 0); //imagePeakThetaH[10]/D");
  waisTree->Branch("imagePeakThetaH", &imagePeakThetaH); //imagePeakThetaH[10]/D");
  waisTree->Branch("imagePeakPhiH", &imagePeakPhiH); //imagePeakPhiH[10]/D");
  waisTree->Branch("imagePeakH", &imagePeakH); //imagePeakH[10]/D");
  waisTree->Branch("gpuPeak", &gpuPeak);
  waisTree->Branch("gpuPeakTheta", &gpuPeakTheta);
  waisTree->Branch("gpuPeakPhi", &gpuPeakPhi);
  waisTree->Branch("thetaExpected", &thetaExpected);
  waisTree->Branch("phiExpected", &phiExpected);
  // waisTree->Branch("droppedPacketsV", &droppedPacketsV, "droppedPacketsV[16][3]/I");
  waisTree->Branch("droppedPacketsH", &droppedPacketsH, "droppedPacketsH[16][3]/I");  
  waisTree->Branch("eventNumberHeader", &eventNumberHeader);

  Long64_t numEntries = headChain->GetEntries();

  Double_t sourceLat = - (79 + (27.93728/60));
  Double_t sourceLon = -(112 + (6.74974/60));
  Double_t sourceAlt = 1813.42;

  Double_t cutTimeNs = 1200; //60e3;

  eventChain->BuildIndex("event->eventNumber");
  gpsChain->BuildIndex("pat->realTime");
  Int_t maxEvents = 10;

  Int_t numEvents = 0;

  std::cout.precision(10);
  for(Long64_t entry=0; entry<numEntries; entry++){
    headChain->GetEntry(entry);
    int eventIndex = eventChain->GetEntryWithIndex(header->eventNumber);
    // int gpsIndex = gpsChain->GetEntryWithIndex(header->realTime);

    // if(eventIndex > 0 && gpsIndex > 0){
    
      // std::cout << header->realTime << " " <<  pat->realTime << std::endl;
      
    // std::cout << gpsInterp.fXmin << " " << header->realTime << std::endl;
    // std::cout << gpsInterp.fXmin - header->realTime << std::endl;
    
    if(eventIndex > 0 && header->realTime >= gpsInterp.fXmin && header->realTime <= gpsInterp.fXmax){

      Int_t triggerTimeNs = header->triggerTimeNs;
      Int_t tt1e8 = header->triggerTimeNs%(Int_t(1e8));
      Int_t goodTriggerTimeNsFlag = 0;
      // for(Int_t cut=0; cut<nCuts; cut++){
      // 	if(tt1e8 > triggerTimeNsLow[cut] && tt1e8 < triggerTimeNsHigh[cut]){
      // 	  goodTriggerTimeNsFlag = 1;
      // 	}
      // }

      Adu5Pat pat2;
      pat2.heading = gpsInterp.interp("pat->heading", header->realTime);
      pat2.latitude = gpsInterp.interp("pat->latitude", header->realTime);
      pat2.longitude = gpsInterp.interp("pat->longitude", header->realTime);
      pat2.altitude = gpsInterp.interp("pat->altitude", header->realTime);

      
      //      UsefulAdu5Pat usefulPat(pat);
      UsefulAdu5Pat usefulPat(&pat2);
      UInt_t triggerTimeNsExpected = usefulPat.getTriggerTimeNsFromSource(sourceLat, sourceLon, sourceAlt);
      usefulPat.getThetaAndPhiWaveAnita3(sourceLon, sourceLat, sourceAlt, thetaExpected, phiExpected);
      thetaExpected *= TMath::RadToDeg();
      phiExpected *= TMath::RadToDeg();
      if(TMath::Abs(Double_t(triggerTimeNsExpected) - Double_t(triggerTimeNs)) < cutTimeNs){
	goodTriggerTimeNsFlag = 1;
      }

      if(goodTriggerTimeNsFlag){
      	if(event->eventNumber != header->eventNumber){
      	  std::cerr << "Couldn't find matching RawAnitaEvent for RawAnitaHeader "<< header->eventNumber << ". Skipping event. " << std::endl;
      	  continue;
      	}

  	UsefulAnitaEvent* realEvent = new UsefulAnitaEvent(event, WaveCalType::kVTBenS, header);

	for(int phi=0; phi<16; phi++){
	  for(int ring=0; ring<3; ring++){
	    // droppedPacketsV[phi][ring] = 0;
	    droppedPacketsH[phi][ring] = 0;
	    int ant = phi + 16*ring;
	    // TGraph* gV = realEvent->getGraph(ant, AnitaPol::kVertical);
	    // if(getRMS(gV) > 0) droppedPacketsV[phi][ring] = 1;
	    TGraph* gH = realEvent->getGraph(ant, AnitaPol::kHorizontal);
	    if(getRMS(gH) < minRMS) droppedPacketsH[phi][ring] = 1;
	    // delete gV;
	    delete gH;
	  }
	}


	cc->correlateEventGPU(realEvent);

	/* for debugging */
	//writeCorrelationGraphs(cc);

	// TH2D* hImageV = cc->makeImageGPU(AnitaPol::kVertical);
	// , imagePeakV[step], 
	//     imagePeakThetaV[step], imagePeakPhiV[step]);
	hImageH = cc->makeImageGPU(AnitaPol::kHorizontal); //, UInt_t(header->l3TrigPattern));
	// , imagePeakH[step],
	//     imagePeakThetaH[step], imagePeakPhiH[step]);

	imagePeakH = cc->findImagePeak(hImageH, imagePeakThetaH, imagePeakPhiH);
	// imagePeakV = cc->findImagePeak(hImageV, imagePeakThetaV, imagePeakPhiV);

	char name[1024];
	// sprintf(name, "h%uV", event->eventNumber);
	// hImageV->SetName(name);
	sprintf(name, "h%uH", event->eventNumber);
	hImageH->SetName(name);

	// hImageV->Write();
	hImageH->Write();

	// delete hImageV;
	delete hImageH;

	gpuPeak = header->getImagePeak();
	gpuPeakTheta = header->getPeakThetaDeg();
	gpuPeakPhi = header->getPeakPhiDeg();
	eventNumberHeader = header->eventNumber;


  	delete realEvent;

  	waisTree->Fill();

	numEvents++;
	std::cout << "numEvents = " << numEvents << std::endl;
      }
    }

    if(numEvents==maxEvents){
      break;
    }

  }
  
  outFile->Write();
  outFile->Close();
  
  delete cc;
}


Double_t getRMS(TGraph* gr){
  Double_t square = 0;
  Double_t mean = 0;
  for(int i=0; i<gr->GetN(); i++){
    square += gr->GetY()[i]*gr->GetY()[i];
    mean += gr->GetY()[i];
  }
  mean /= gr->GetN();
  return TMath::Sqrt(square/gr->GetN() - mean*mean);
}
