#include "InterferometryCache.h"
#include "AnitaConventions.h"
#include "CrossCorrelator.h"
#include "InterferometricMapMaker.h"

InterferometryCache::InterferometryCache(){
  numCombos = 0;
  nCoarseBinsPhi = 0;
  nCoarseBinsTheta = 0;
}


InterferometryCache::InterferometryCache(CrossCorrelator* cc, InterferometricMapMaker* mm){
  numCombos = 0;
  nCoarseBinsPhi = 0;
  nCoarseBinsTheta = 0;
  populateCache(cc, mm);
}



void InterferometryCache::populateCache(CrossCorrelator* cc, InterferometricMapMaker* mm){

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
	  deltaTs.at(coarseIndex(pol, combo, phiBin, thetaBin)) = mm->getDeltaTExpected(pol, ant1, ant2, phiWave, thetaWave);
	}
      }
    }
  }
}
  
 
