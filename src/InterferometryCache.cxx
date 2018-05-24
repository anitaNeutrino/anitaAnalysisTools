#include "InterferometryCache.h"
#include "CrossCorrelator.h"
#include "AnalysisReco.h"

Acclaim::InterferometryCache::InterferometryCache(){
  fInitialized = false;
  fNumCombos = 0;
  fNCoarseBinsPhi = 0;
  fNCoarseBinsTheta = 0;
}

Acclaim::InterferometryCache::InterferometryCache(CrossCorrelator* cc, const AnalysisReco* reco){
  fInitialized = false;
  fNumCombos = 0;
  fNCoarseBinsPhi = 0;
  fNCoarseBinsTheta = 0;
  populateCache(cc, reco);
  populateFineCache(cc, reco);  
}

void Acclaim::InterferometryCache::init(CrossCorrelator* cc, const AnalysisReco* reco, bool forceCacheRecalculation){
  if(!fInitialized || forceCacheRecalculation){
    populateCache(cc, reco);
    populateFineCache(cc, reco);
  }
  fInitialized = true;
}

void Acclaim::InterferometryCache::populateCache(CrossCorrelator* cc, const AnalysisReco* reco){

  fUseOffAxisDelay = reco->GetUseOffAxisDelay();
  // std::cout << "in cache fUseOffAxisDelay " << fUseOffAxisDelay << std::endl;
  fNumCombos = cc->numCombos;
  
  const std::vector<Double_t> coarsePhiBinEdges = InterferometricMap::getCoarseBinEdgesPhi();
  fNCoarseBinsPhi = coarsePhiBinEdges.size()-1;

  const std::vector<Double_t> coarseThetaBinEdges = InterferometricMap::getCoarseBinEdgesTheta();
  fNCoarseBinsTheta = coarseThetaBinEdges.size()-1;
  

  std::vector<Double_t> phiWaveLookup(fNCoarseBinsPhi);
  for(Int_t phiIndex=0; phiIndex < fNCoarseBinsPhi; phiIndex++){
    Double_t phiDeg = coarsePhiBinEdges.at(phiIndex);
    Double_t phiWave = TMath::DegToRad()*phiDeg;
    phiWaveLookup[phiIndex] = phiWave;
  }

  std::vector<Double_t> thetaWaves(fNCoarseBinsTheta);

  for(Int_t thetaIndex=0; thetaIndex < fNCoarseBinsTheta; thetaIndex++){
    Double_t thetaWaveDeg = coarseThetaBinEdges.at(thetaIndex); 
    Double_t thetaWave = thetaWaveDeg*TMath::DegToRad();
    thetaWaves[thetaIndex] = thetaWave;
  }

  int numDeltaTs = AnitaPol::kNotAPol*fNumCombos*fNCoarseBinsTheta*fNCoarseBinsPhi;
  fDeltaTs.resize(numDeltaTs); // reserve memory and set the correct size so we use .at() in the loop
  
  for(Int_t polInd=0; polInd<AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;

    for(Int_t combo=0; combo<fNumCombos; combo++){
      Int_t ant1 = cc->comboToAnt1s.at(combo);
      Int_t ant2 = cc->comboToAnt2s.at(combo);

      for(Int_t phiBin = 0; phiBin < fNCoarseBinsPhi; phiBin++){
	Double_t phiWave = phiWaveLookup[phiBin];

	for(Int_t thetaBin = 0; thetaBin < fNCoarseBinsTheta; thetaBin++){
	  Double_t thetaWave = thetaWaves[thetaBin];
	  fDeltaTs.at(coarseIndex(pol, combo, phiBin, thetaBin)) = reco->getDeltaTExpected(pol, ant1, ant2, phiWave, thetaWave);
	}
      }
    }
  }
}

void Acclaim::InterferometryCache::populateFineCache(CrossCorrelator* cc, const AnalysisReco* reco){

  std::vector<double> fineBinsTheta = InterferometricMap::getFineBinEdgesTheta();
  std::vector<double> fineBinsPhi = InterferometricMap::getFineBinEdgesPhi();  

  fNFineBinsTheta = fineBinsTheta.size() - 1;
  fNFineBinsPhi = fineBinsPhi.size() - 1;

  // reserve cache and populate with 0s
  fZoomedPhiWaveLookup.resize(fNFineBinsPhi);
  for(unsigned phiIndex=0; phiIndex < fZoomedPhiWaveLookup.size(); phiIndex++){
    Double_t phiDeg = fineBinsPhi.at(phiIndex);
    // Double_t phiDeg = minPhiDegZoom + phiIndex*ZOOM_BIN_SIZE_PHI;
    // zoomedPhiWaveLookup[phiIndex] = phiDeg*TMath::DegToRad();
    fZoomedPhiWaveLookup.at(phiIndex) = phiDeg*TMath::DegToRad();
  }
  
  fZoomedThetaWaves.resize(fNFineBinsTheta);
  fZoomedTanThetaWaves.resize(fNFineBinsTheta);
  fZoomedCosThetaWaves.resize(fNFineBinsTheta);
  fDtFactors.resize(fNFineBinsTheta);
  

  for(unsigned thetaIndex=0; thetaIndex < fZoomedThetaWaves.size(); thetaIndex++){
    // Double_t thetaWaveDeg = minThetaDegZoom + thetaIndex*ZOOM_BIN_SIZE_THETA;
    Double_t thetaWaveDeg = fineBinsTheta.at(thetaIndex);
    Double_t thetaWave = -1*thetaWaveDeg*TMath::DegToRad();
    fZoomedThetaWaves[thetaIndex] = thetaWave;
    fZoomedTanThetaWaves[thetaIndex] = tan(thetaWave);
    fZoomedCosThetaWaves[thetaIndex] = cos(thetaWave);
    fDtFactors[thetaIndex] = fZoomedCosThetaWaves[thetaIndex]/(SPEED_OF_LIGHT_NS*cc->correlationDeltaT);
  }


  // std::cerr << "here" << std::endl;
  
  fPartBAsZoom.resize(AnitaPol::kNotAPol*fNumCombos*fNFineBinsTheta);
  for(Int_t pol=0; pol<AnitaPol::kNotAPol; pol++){
    for(Int_t combo=0; combo < NUM_COMBOS; combo++){
      Int_t ant1 = cc->comboToAnt1s.at(combo);
      Int_t ant2 = cc->comboToAnt2s.at(combo);
      for(Int_t thetaIndex=0; thetaIndex < fNFineBinsTheta; thetaIndex++){
	fPartBAsZoom[partBAsIndex(pol, combo, thetaIndex)] = fZoomedTanThetaWaves[thetaIndex]*(reco->fZArray[pol].at(ant2)-reco->fZArray[pol].at(ant1));
      }
    }
  }

  fZoomedCosPartLookup.resize(AnitaPol::kNotAPol*NUM_SEAVEYS*fNFineBinsPhi);

  fOffAxisDelays.resize(AnitaPol::kNotAPol*fNumCombos*fNFineBinsPhi);
  fOffAxisDelaysDivided.resize(AnitaPol::kNotAPol*fNumCombos*fNFineBinsPhi);  
  fPart21sZoom.resize(AnitaPol::kNotAPol*fNumCombos*fNFineBinsPhi);  


  for(Int_t pol=0; pol<AnitaPol::kNotAPol; pol++){
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      for(Int_t phiIndex=0; phiIndex < fNFineBinsPhi; phiIndex++){
	Double_t phiDeg = fineBinsPhi.at(phiIndex);
	Double_t phiWave = phiDeg*TMath::DegToRad();
	fZoomedCosPartLookup.at(zoomedCosPartIndex(pol, ant, phiIndex)) = reco->fRArray[pol].at(ant)*cos(phiWave-TMath::DegToRad()*reco->fPhiArrayDeg[pol].at(ant));
      }
    }

    for(Int_t combo=0; combo < NUM_COMBOS; combo++){
      Int_t ant1 = cc->comboToAnt1s.at(combo);
      Int_t ant2 = cc->comboToAnt2s.at(combo);

      for(Int_t phiIndex=0; phiIndex < fNFineBinsPhi; phiIndex++){
	Double_t phiDeg = fineBinsPhi.at(phiIndex);
	Double_t offAxisDelay = reco->relativeOffAxisDelay((AnitaPol::AnitaPol_t)pol, ant2, ant1, phiDeg);

	int p21 = part21sIndex(pol, combo, phiIndex);
	fOffAxisDelays[p21] = offAxisDelay;
	fOffAxisDelaysDivided[p21] = offAxisDelay/cc->correlationDeltaT;
	fPart21sZoom[p21] = fZoomedCosPartLookup.at((zoomedCosPartIndex(pol, ant2, phiIndex))) - fZoomedCosPartLookup.at((zoomedCosPartIndex(pol, ant1, phiIndex)));
      }
    }
  }
}
  
 
