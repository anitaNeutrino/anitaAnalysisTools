#include "DataQualityMonitor.h"

DataQualityMonitor::DataQualityMonitor(){

}

DataQualityMonitor::DataQualityMonitor(TChain* c){

  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      
      maxAbsSecondDeriv[polInd][ant] = 0;
      maxVolts[polInd][ant] = 0;
      numPoints[polInd][ant] = 0;
    }
  }
  eventNumber = 0;
  
  setBranches(c);
}


void DataQualityMonitor::setBranches(TChain* c){
  dataQualityChain = c; 
  dataQualityChain->SetBranchAddress("maxAbsSecondDeriv", &maxAbsSecondDeriv[0][0]);
  dataQualityChain->SetBranchAddress("maxVolts", &maxVolts[0][0]);
  dataQualityChain->SetBranchAddress("numPoints", &numPoints[0][0]);
  dataQualityChain->SetBranchAddress("eventNumber", &eventNumber);  


  maxVoltsBlastThresh = 400;
  saturationVolts = 1500;  
}

DataQualityMonitor::~DataQualityMonitor(){
  
}



Int_t DataQualityMonitor::processEntry(Long64_t entry, UInt_t eventNumberCheck){

  // std::cerr << dataQualityChain->GetEntry(entry) << std::endl;
  dataQualityChain->GetEntry(entry);
  
  if(eventNumberCheck > 0){
    if(eventNumberCheck != eventNumber){
      std::cerr << "Warning in " << __PRETTY_FUNCTION__ << " for entry " << entry << ", eventNumber mismatch! eventNumber = "
		<< eventNumber << ", but eventNumberCheck = " << eventNumberCheck << std::endl;
      return -1;
    }
  }
  
  
  numChannelsAboveSurfSaturation = 0;
  numAboveVoltsBlastThresh = 0;
  numPhiAboveMaxVoltsBlastThresh = 0;  
  for(int phi=0; phi<NUM_PHI; phi++){
    phiAboveMaxVoltsThresh[phi] = 0;
  }

  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    for(int phi=0; phi<NUM_PHI; phi++){
      phiAboveMaxVoltsThresh[phi] = 0;
    }
      
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      if(maxVolts[polInd][ant] > maxVoltsBlastThresh){
	numAboveVoltsBlastThresh++;
	int phi = ant%NUM_PHI;
	phiAboveMaxVoltsThresh[phi] = 1;
      }

      if(maxVolts[polInd][ant] > saturationVolts){
	numChannelsAboveSurfSaturation++;
      }
    }
  }

  for(int phi=0; phi<NUM_PHI; phi++){
    if(phiAboveMaxVoltsThresh[phi] > 0){
      numPhiAboveMaxVoltsBlastThresh++;
    }
  }
    
  // hNumChannelsAboveMaxVolts->Fill(numAboveVoltsBlastThresh);
  // hNumPhiAboveMaxVolts->Fill(numPhiAboveMaxVoltsBlastThresh);    
  // hAnyChannelsAboveSaturation->Fill(numChannelsAboveSurfSaturation);
  // p++;


  Int_t dataQualityNum = 0;
  if(numAboveVoltsBlastThresh >= 15 || numPhiAboveMaxVoltsBlastThresh >= 9){
    std::cerr << std::endl;
    std::cerr << "blast ? " << eventNumber << "\t" << numAboveVoltsBlastThresh << "\t" << numPhiAboveMaxVoltsBlastThresh << std::endl;
    dataQualityNum = 1;
  }
  else if(numChannelsAboveSurfSaturation >= 3){
    std::cerr << std::endl;
    std::cerr << "three or more saturating channels ? " << eventNumber << "\t" << numChannelsAboveSurfSaturation << std::endl;
    dataQualityNum = 2;    
  }

  return dataQualityNum;
}
