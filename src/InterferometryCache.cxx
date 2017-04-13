#include "InterferometryCache.h"
#include "CrossCorrelator.h"
#include "AnalysisReco.h"

Acclaim::InterferometryCache::InterferometryCache(){
  initialized = false;
  numCombos = 0;
  nCoarseBinsPhi = 0;
  nCoarseBinsTheta = 0;
}


Acclaim::InterferometryCache::InterferometryCache(CrossCorrelator* cc, const AnalysisReco* reco){
  initialized = false;
  numCombos = 0;
  nCoarseBinsPhi = 0;
  nCoarseBinsTheta = 0;
  populateCache(cc, reco);
  populateFineCache(cc, reco);  
}


void Acclaim::InterferometryCache::init(CrossCorrelator* cc, const AnalysisReco* reco, bool forceCacheRecalculation){
  if(!initialized || forceCacheRecalculation){
    populateCache(cc, reco);
    populateFineCache(cc, reco);
  }
  initialized = true;
}


void Acclaim::InterferometryCache::populateCache(CrossCorrelator* cc, const AnalysisReco* reco){


  kUseOffAxisDelay = reco->kUseOffAxisDelay;
  // std::cout << "in cache kUseOffAxisDelay " << kUseOffAxisDelay << std::endl;
  numCombos = cc->numCombos;
  
  const std::vector<Double_t> coarsePhiBinEdges = InterferometricMap::getCoarseBinEdgesPhi();
  nCoarseBinsPhi = coarsePhiBinEdges.size()-1;

  const std::vector<Double_t> coarseThetaBinEdges = InterferometricMap::getCoarseBinEdgesTheta();
  nCoarseBinsTheta = coarseThetaBinEdges.size()-1;
  

  std::vector<Double_t> phiWaveLookup(nCoarseBinsPhi);
  for(Int_t phiIndex=0; phiIndex < nCoarseBinsPhi; phiIndex++){
    Double_t phiDeg = coarsePhiBinEdges.at(phiIndex);
    Double_t phiWave = TMath::DegToRad()*phiDeg;
    phiWaveLookup[phiIndex] = phiWave;
  }

  std::vector<Double_t> thetaWaves(nCoarseBinsTheta);

  for(Int_t thetaIndex=0; thetaIndex < nCoarseBinsTheta; thetaIndex++){
    Double_t thetaWaveDeg = coarseThetaBinEdges.at(thetaIndex); 
    Double_t thetaWave = thetaWaveDeg*TMath::DegToRad();
    thetaWaves[thetaIndex] = thetaWave;
  }

  
  int numDeltaTs = AnitaPol::kNotAPol*numCombos*nCoarseBinsTheta*nCoarseBinsPhi;
  deltaTs.resize(numDeltaTs); // reserve memory and set the correct size so we use .at() in the loop
  
  for(Int_t polInd=0; polInd<AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;

    for(Int_t combo=0; combo<numCombos; combo++){
      Int_t ant1 = cc->comboToAnt1s.at(combo);
      Int_t ant2 = cc->comboToAnt2s.at(combo);

      for(Int_t phiBin = 0; phiBin < nCoarseBinsPhi; phiBin++){
	Double_t phiWave = phiWaveLookup[phiBin];

	for(Int_t thetaBin = 0; thetaBin < nCoarseBinsTheta; thetaBin++){
	  Double_t thetaWave = thetaWaves[thetaBin];
	  deltaTs.at(coarseIndex(pol, combo, phiBin, thetaBin)) = reco->getDeltaTExpected(pol, ant1, ant2, phiWave, thetaWave);
	}
      }
    }
  }
}






void Acclaim::InterferometryCache::populateFineCache(CrossCorrelator* cc, const AnalysisReco* reco){

  // std::cerr << "here" << std::endl;
  
  // Double_t minThetaDegZoom = MIN_THETA - THETA_RANGE_ZOOM/2;
  // Double_t minPhiDegZoom = InterferometricMap::getBin0PhiDeg() - PHI_RANGE_ZOOM/2;

  

  // Double_t minThetaDegZoom = MIN_THETA - THETA_RANGE_ZOOM/2;
  // Double_t minPhiDegZoom = InterferometricMap::getBin0PhiDeg() - PHI_RANGE_ZOOM/2;

  std::vector<double> fineBinsTheta = InterferometricMap::getFineBinEdgesTheta();
  std::vector<double> fineBinsPhi = InterferometricMap::getFineBinEdgesPhi();  

  Double_t minThetaDegZoom = fineBinsTheta.at(0);
  Double_t minPhiDegZoom = fineBinsPhi.at(0);
  
  nFineBinsTheta = fineBinsTheta.size() - 1;
  nFineBinsPhi = fineBinsPhi.size() - 1;
  

  // reserve cache and populate with 0s
  zoomedPhiWaveLookup.resize(nFineBinsPhi);
  for(unsigned phiIndex=0; phiIndex < zoomedPhiWaveLookup.size(); phiIndex++){
    Double_t phiDeg = fineBinsPhi.at(phiIndex);
    // Double_t phiDeg = minPhiDegZoom + phiIndex*ZOOM_BIN_SIZE_PHI;
    // zoomedPhiWaveLookup[phiIndex] = phiDeg*TMath::DegToRad();
    zoomedPhiWaveLookup.at(phiIndex) = phiDeg*TMath::DegToRad();
  }
  
  zoomedThetaWaves.resize(nFineBinsTheta);
  zoomedTanThetaWaves.resize(nFineBinsTheta);
  zoomedCosThetaWaves.resize(nFineBinsTheta);
  dtFactors.resize(nFineBinsTheta);
  

  for(unsigned thetaIndex=0; thetaIndex < zoomedThetaWaves.size(); thetaIndex++){
    // Double_t thetaWaveDeg = minThetaDegZoom + thetaIndex*ZOOM_BIN_SIZE_THETA;
    Double_t thetaWaveDeg = fineBinsTheta.at(thetaIndex);
    Double_t thetaWave = -1*thetaWaveDeg*TMath::DegToRad();
    zoomedThetaWaves[thetaIndex] = thetaWave;
    zoomedTanThetaWaves[thetaIndex] = tan(thetaWave);
    zoomedCosThetaWaves[thetaIndex] = cos(thetaWave);
    dtFactors[thetaIndex] = zoomedCosThetaWaves[thetaIndex]/(SPEED_OF_LIGHT_NS*cc->correlationDeltaT);
  }


  // std::cerr << "here" << std::endl;
  
  partBAsZoom.resize(AnitaPol::kNotAPol*numCombos*nFineBinsTheta);
  for(Int_t pol=0; pol<AnitaPol::kNotAPol; pol++){
    for(Int_t combo=0; combo < NUM_COMBOS; combo++){
      Int_t ant1 = cc->comboToAnt1s.at(combo);
      Int_t ant2 = cc->comboToAnt2s.at(combo);
      for(Int_t thetaIndex=0; thetaIndex < nFineBinsTheta; thetaIndex++){
	// partBAsZoom[pol][combo][thetaIndex] = zoomedTanThetaWaves[thetaIndex]*(zArray[pol].at(ant2)-zArray[pol].at(ant1));
	// partBAsZoom[pol][combo][thetaIndex] = zoomedTanThetaWaves[thetaIndex]*(reco->zArray[pol].at(ant2)-reco->zArray[pol].at(ant1));		
	partBAsZoom[partBAsIndex(pol, combo, thetaIndex)] = zoomedTanThetaWaves[thetaIndex]*(reco->zArray[pol].at(ant2)-reco->zArray[pol].at(ant1));
	// partBAsZoom[pol][combo][thetaIndex] = zoomedTanThetaWaves[thetaIndex]*(reco->zArray[pol].at(ant2)-reco->zArray[pol].at(ant1));
      }
    }
  }

  // std::cerr << "here" << std::endl;  

  zoomedCosPartLookup.resize(AnitaPol::kNotAPol*NUM_SEAVEYS*nFineBinsPhi);

         offAxisDelays.resize(AnitaPol::kNotAPol*numCombos*nFineBinsPhi);
  offAxisDelaysDivided.resize(AnitaPol::kNotAPol*numCombos*nFineBinsPhi);  
           part21sZoom.resize(AnitaPol::kNotAPol*numCombos*nFineBinsPhi);  


  for(Int_t pol=0; pol<AnitaPol::kNotAPol; pol++){
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      for(Int_t phiIndex=0; phiIndex < nFineBinsPhi; phiIndex++){
	// Double_t phiDeg = minPhiDegZoom + phiIndex*ZOOM_BIN_SIZE_PHI;
	Double_t phiDeg = fineBinsPhi.at(phiIndex);
	Double_t phiWave = phiDeg*TMath::DegToRad();
	// zoomedCosPartLookup[pol][ant][phiIndex] = rArray[pol].at(ant)*cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol].at(ant));
	// std::cout << pol << "\t" << ant << "\t" << phiIndex << "\t" << zoomedCosPartIndex(pol, ant, phiIndex) << std::endl;
	zoomedCosPartLookup.at(zoomedCosPartIndex(pol, ant, phiIndex)) = reco->rArray[pol].at(ant)*cos(phiWave-TMath::DegToRad()*reco->phiArrayDeg[pol].at(ant));
      }
    }

    for(Int_t combo=0; combo < NUM_COMBOS; combo++){
      Int_t ant1 = cc->comboToAnt1s.at(combo);
      Int_t ant2 = cc->comboToAnt2s.at(combo);

      for(Int_t phiIndex=0; phiIndex < nFineBinsPhi; phiIndex++){
	// Double_t phiWave = TMath::DegToRad()*phiIndex*ZOOM_BIN_SIZE_PHI;
	Double_t phiDeg = fineBinsPhi.at(phiIndex);
	// Double_t offAxisDelay = getOffAxisDelay((AnitaPol::AnitaPol_t)pol, ant1, ant2, phiWave);
	// offAxisDelays[pol][combo][phiIndex] = offAxisDelay;
	Double_t offAxisDelay = reco->relativeOffAxisDelay((AnitaPol::AnitaPol_t)pol, ant2, ant1, phiDeg);

	int p21 = part21sIndex(pol, combo, phiIndex);
	offAxisDelays[p21] = offAxisDelay;
	offAxisDelaysDivided[p21] = offAxisDelay/cc->correlationDeltaT;
	part21sZoom[p21] = zoomedCosPartLookup.at((zoomedCosPartIndex(pol, ant2, phiIndex))) - zoomedCosPartLookup.at((zoomedCosPartIndex(pol, ant1, phiIndex)));
      }
    }
  }

  // std::cout << "done" << std::endl;
}
  
 
