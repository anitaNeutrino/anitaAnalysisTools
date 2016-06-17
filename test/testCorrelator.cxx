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
void testNormalization();
void testHackyFilter();
void testOffAxisDelay();
// void testFileWriting();

int main(){

  //  testImageGPUStyle();
  //  testNewCombinatorics();
  // testImageFullStyle();
  // testCoherentlySummedWaveform();
  // testNormalization();
  // testHackyFilter();
  testOffAxisDelay();
  //  testFileWriting();

  return 0;

}

void testOffAxisDelay(){

  CrossCorrelator* cc = new CrossCorrelator();
  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;


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
  
  
  headTree->BuildIndex("eventNumber");
  eventTree->BuildIndex("eventNumber");
  
  const Long64_t eventNumber = 60832108;
  headTree->GetEntryWithIndex(eventNumber);
  eventTree->GetEntryWithIndex(eventNumber);
  std::cout << header->eventNumber << "\t" << event->eventNumber << std::endl;

  UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(event);  
  cc->kUseOffAxisDelay = 1;
  cc->reconstructEvent(usefulEvent, 2, 2);


  TFile* outFile = new TFile("/tmp/testOffAxisDelay.root","recreate");    

  Double_t peakValue, imagePeak, peakPhiDeg, peakThetaDeg;
  TH2D* hCoarse = cc->getMap(pol, peakValue, peakPhiDeg, peakThetaDeg);
  // hCoarse->Write();
  

  for(int i=0; i < 2; i++){
    std::cout << i << "\t" << cc->coarseMapPeakValues[pol][i] << "\t"
	      << cc->coarseMapPeakPhiDegs[pol][i] << "\t" 
	      << cc->coarseMapPeakThetaDegs[pol][i] << std::endl;
  }
  
  cc->kUseOffAxisDelay = 1;  
  TH2D* hWithOffAxisDelay = cc->makeZoomedImage(pol, imagePeak, peakPhiDeg, peakThetaDeg,
						   cc->coarseMapPeakPhiDegs[pol][0],
						   cc->coarseMapPeakThetaDegs[pol][0]);
  hWithOffAxisDelay->SetName("hWithOffAxisDelay");
  hWithOffAxisDelay->SetTitle("With Off-Axis Delay");
  
  cc->kUseOffAxisDelay = 0; 
  TH2D* hWithOutOffAxisDelay = cc->makeZoomedImage(pol, imagePeak, peakPhiDeg, peakThetaDeg,
						   cc->coarseMapPeakPhiDegs[pol][0],
						   cc->coarseMapPeakThetaDegs[pol][0]);

  hWithOutOffAxisDelay->SetName("hWithOutOffAxisDelay");
  hWithOutOffAxisDelay->SetTitle("Without Off-Axis Delay");  



  Int_t ant1 = 1;

  for(int antInd=0; antInd < 2; antInd++){
    Int_t ant2 = 2 + antInd;
    const double dPhi = RootTools::getDeltaAngleDeg(cc->phiArrayDeg[pol].at(ant2),
						    cc->phiArrayDeg[pol].at(ant1));
    const double mp = cc->phiArrayDeg[pol].at(ant1) + 0.5*dPhi;

    std::cout << cc->phiArrayDeg[pol].at(ant1) << "\t" << mp
	      << "\t" << cc->phiArrayDeg[pol].at(ant2) << "\t" << dPhi << std::endl;
    
    TGraph* grOffAxisDelay = new TGraph();
    for(Double_t phiFromMidPoint = -90; phiFromMidPoint <= 90; phiFromMidPoint += 1){
      Double_t delay = cc->relativeOffAxisDelay(pol, ant1, ant2, phiFromMidPoint+mp);
      grOffAxisDelay->SetPoint(grOffAxisDelay->GetN(), phiFromMidPoint, delay);
    }

    TString name = TString::Format("grOffAxisDelay_%d_%d_%d", pol, ant1, ant2);
    grOffAxisDelay->SetName(name);
    TString title = TString::Format("Antennas %d and %d", ant1, ant2);
    grOffAxisDelay->SetTitle(title);
    grOffAxisDelay->Write();

    delete grOffAxisDelay;
  }



  
  
  delete cc;
  
  outFile->Write();
  outFile->Close();  
}




void testNormalization(){
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

  TFile* outFile = new TFile("/tmp/testNormalization.root","recreate");

  CrossCorrelator* cc = new CrossCorrelator();

  headTree->BuildIndex("eventNumber");
  eventTree->BuildIndex("eventNumber");
  
  const Long64_t eventNumber = 60832108;
  headTree->GetEntryWithIndex(eventNumber);
  eventTree->GetEntryWithIndex(eventNumber);
  std::cout << header->eventNumber << "\t" << event->eventNumber << std::endl;
  
  // UsefulAnitaEvent* realEvent(new UsefulAnitaEvent(event, WaveCalType::kDefault,  header));
  UsefulAnitaEvent* realEvent = new UsefulAnitaEvent(event);
  cc->correlateEvent(realEvent);

  if(cc->grsResampled[0][1]){
    delete cc->grsResampled[0][1];
    cc->grsResampled[0][1] = (TGraph*) cc->grsResampled[0][0]->Clone("grCopy");
  }

  cc->doFFTs(AnitaPol::kHorizontal);
  cc->doAllCrossCorrelationsThreaded(AnitaPol::kHorizontal);
  cc->doUpsampledCrossCorrelationsThreaded(AnitaPol::kHorizontal, 0);

  TGraph* grCrossCorr = cc->getCrossCorrelationGraph(AnitaPol::kHorizontal, 0, 1);
  grCrossCorr->Write();

  TGraph* grCrossCorr2 = cc->getUpsampledCrossCorrelationGraph(AnitaPol::kHorizontal, 0, 1);
  grCrossCorr2->Write();
  
  Double_t mean=0, rms=1;
  for(int pol=0; pol < AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      RootTools::getMeanAndRms(cc->grsResampled[pol][ant], mean, rms);
      // std::cout << pol << "\t" << ant << "\t" << mean << "\t" << rms << std::endl;
      std::cout << pol << "\t" << ant << "\t" << mean << "\t" << rms << "\t";
      std::cout << pow(cc->interpRMS[pol][ant], 2) << "\t" << pow(cc->interpRMS2[pol][ant], 2) << "\t"
		<< pow(cc->interpRMS[pol][ant],2)/pow(cc->interpRMS2[pol][ant],2) << "\t"
		<< std::endl;
    }
  }

  Double_t* volts = FancyFFTs::doInvFFT(cc->numSamples, cc->ffts[0][0], true);
  Double_t* voltsUpsampled = FancyFFTs::doInvFFT(cc->numSamplesUpsampled, cc->fftsPadded[0][0], true);

  std::vector<Double_t> times(cc->numSamples, 0);
  for(int samp=0; samp < cc->numSamples; samp++){
    times.at(samp) = cc->nominalSamplingDeltaT*samp;
  }
  
  std::vector<Double_t> timesUpsampled(cc->numSamplesUpsampled, 0);
  for(int samp=0; samp < cc->numSamplesUpsampled; samp++){
    timesUpsampled.at(samp) = cc->correlationDeltaT*samp;    
  }

  TGraph* grVolts = new TGraph(cc->numSamples, &times[0], volts);  
  TGraph* grVoltsUpsampled = new TGraph(cc->numSamplesUpsampled, &timesUpsampled[0], voltsUpsampled);

  {
    Double_t mean, rms;
    RootTools::getMeanAndRms(grVolts, mean, rms);
    std::cout << "grVolts\t" << mean << "\t" << rms << std::endl;
    RootTools::getMeanAndRms(grVoltsUpsampled, mean, rms);    
    std::cout << "grVoltsUpsampled\t" << mean << "\t" << rms << std::endl;    
  }

  grVolts->SetName("grVolts");
  grVoltsUpsampled->SetName("grVoltsUpsampled");

  grVolts->Write();
  grVoltsUpsampled->Write();

  delete [] volts;
  delete [] voltsUpsampled;  
  
  const int numFreqs = FancyFFTs::getNumFreqs(cc->numSamples);
  const int numFreqsPadded = FancyFFTs::getNumFreqs(cc->numSamplesUpsampled);

  for(int pol=0; pol < AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      Double_t sumSquaredPadded = 0;
      for(int freqInd=0; freqInd < numFreqsPadded; freqInd++){
	if(freqInd == numFreqsPadded - 1){
	  sumSquaredPadded += std::norm(cc->fftsPadded[pol][ant][freqInd])/2;
	}
	// if(freqInd == numFreqs - 1){
	//   sumSquaredPadded += std::norm(cc->fftsPadded[pol][ant][freqInd])/2;	  
	// }	
	else{
	  sumSquaredPadded += std::norm(cc->fftsPadded[pol][ant][freqInd]);
	}
      }
      std::cout << pol << "\t" << ant << "\t" << sumSquaredPadded << "\t" << cc->numSamplesUpsampled << "\t" << sumSquaredPadded/cc->numSamplesUpsampled << std::endl;

      Double_t sumSquared = 0;
      for(int freqInd=0; freqInd < numFreqs; freqInd++){
	if(freqInd == numFreqs - 1){
	  sumSquared += std::norm(cc->ffts[pol][ant][freqInd])/2;
	}
	else{
	  sumSquared += std::norm(cc->ffts[pol][ant][freqInd]);
	}
      }
      std::cout << pol << "\t" << ant << "\t" << sumSquared << "\t" << cc->numSamples << "\t" << sumSquared/cc->numSamples << std::endl;
    }
  }

  
  outFile->Write();
  outFile->Close();  
}




void testHackyFilter(){
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

  TFile* outFile = new TFile("/tmp/testHackyFilter.root","recreate");

  CrossCorrelator* cc = new CrossCorrelator();
  cc->kDoSimpleSatelliteFiltering = 1;

  headTree->BuildIndex("eventNumber");
  eventTree->BuildIndex("eventNumber");
  
  const Long64_t eventNumber = 60832108;
  headTree->GetEntryWithIndex(eventNumber);
  eventTree->GetEntryWithIndex(eventNumber);
  std::cout << header->eventNumber << "\t" << event->eventNumber << std::endl;
  
  // UsefulAnitaEvent* realEvent(new UsefulAnitaEvent(event, WaveCalType::kDefault,  header));
  UsefulAnitaEvent* realEvent = new UsefulAnitaEvent(event);
  cc->correlateEvent(realEvent);

  if(cc->grsResampled[0][1]){
    delete cc->grsResampled[0][1];
    cc->grsResampled[0][1] = (TGraph*) cc->grsResampled[0][0]->Clone("grCopy");
  }

  Double_t mean = 0;
  Double_t rms = 0;
  Int_t n = cc->grsResampled[0][1]->GetN();
  for(int samp=0; samp < n; samp++){
    Double_t y = cc->grsResampled[0][1]->GetY()[samp];
    mean += y;
    rms += y*y;
  }

  mean /= n;
  rms  = rms / n - mean*mean;
  std::cout << "mean, rms, n: " << mean << "\t" << rms << "\t" << n << std::endl;

  cc->doFFTs(AnitaPol::kHorizontal);
  cc->doAllCrossCorrelationsThreaded(AnitaPol::kHorizontal);
  cc->doUpsampledCrossCorrelationsThreaded(AnitaPol::kHorizontal, 0);

  for(Int_t polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(Int_t ant=0; ant < NUM_SEAVEYS; ant++){
      cc->simple260MHzSatelliteNotch(pol, ant);
    }
  }
  const int numFreqs = FancyFFTs::getNumFreqs(cc->numSamples);
  std::vector<Double_t> powSpec(numFreqs, 0);
  std::vector<Double_t> freqsMHz(numFreqs, 0);  
  
  for(int freqInd=0; freqInd < numFreqs; freqInd++){
    freqsMHz.at(freqInd) = freqInd*1e3/(cc->nominalSamplingDeltaT*cc->numSamples);
    powSpec.at(freqInd) = std::norm(cc->ffts[0][0][freqInd]);

    if(powSpec.at(freqInd) >= 1e-10){
      powSpec.at(freqInd) = 10*TMath::Log10(powSpec.at(freqInd));
    }
    else{
      powSpec.at(freqInd) = -10;
    }
    std::cout << powSpec.at(freqInd) << "\t"  << cc->ffts[0][0][freqInd] << std::endl;
  }

  TGraph* grPowSpec = new TGraph(numFreqs, &freqsMHz[0], &powSpec[0]);
  grPowSpec->SetName("grPowSpecTest");
  grPowSpec->Write();
  
  outFile->Write();
  outFile->Close();  
}






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
