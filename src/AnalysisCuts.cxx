#include "AnalysisCuts.h"



AnalysisCuts::Status_t AnalysisCuts::applyBottomToTopRingPeakToPeakRatioCut(AnitaPol::AnitaPol_t pol, Double_t* peakToPeak, Double_t& maxRatio){

  Status_t status = kPass;
  maxRatio = 0;
  for(int phi=0; phi < NUM_PHI; phi++){
    if(pol==AnitaPol::kVertical && phi==7){
      continue;
    }
    if(pol==AnitaPol::kHorizontal && phi==4){
      continue;
    }

    Double_t ratio = peakToPeak[phi+2*NUM_PHI]/peakToPeak[phi];
    if(ratio > maxRatio){
      maxRatio = ratio;
    }
  }
  if(maxRatio > ratioCutHigh || maxRatio < ratioCutLow){
    // std::cerr << eventNumberDQ << "\t" << maxRatio << std::endl;
    status = kFail;
  }

  return status;
}


AnalysisCuts::Status_t AnalysisCuts::applyBottomToTopRingPeakToPeakRatioCut(Double_t maxPeakToPeakRatio){

  Status_t status = kPass;
  if(maxPeakToPeakRatio > ratioCutHigh || maxPeakToPeakRatio < ratioCutLow){
    // std::cerr << eventNumberDQ << "\t" << maxRatio << std::endl;
    status = kFail;
  }

  return status;
}





AnalysisCuts::Status_t AnalysisCuts::L3TriggerDirectionCut(AnitaPol::AnitaPol_t pol, RawAnitaHeader* header, Double_t recoPhiDeg, Int_t& deltaPhiSect){

  Status_t status = kPass;

  // CUT FLOW
  // Step 2: cut phi-sector angle triggers
  const double aftForeOffset = 45; // in cc it's
  const double bin0PhiDeg = -aftForeOffset + DEGREES_IN_CIRCLE;
  double angleThroughPhiSectors = recoPhiDeg - bin0PhiDeg;
  angleThroughPhiSectors += angleThroughPhiSectors < 0 ? DEGREES_IN_CIRCLE : 0;
  angleThroughPhiSectors -= angleThroughPhiSectors >= DEGREES_IN_CIRCLE ? DEGREES_IN_CIRCLE : 0;
  if(angleThroughPhiSectors < 0 || angleThroughPhiSectors >= DEGREES_IN_CIRCLE){
    std::cerr << "you moron " << angleThroughPhiSectors << "\t" << recoPhiDeg << std::endl;
  }
  // const double phiRelativeToPhiSector0 = RootTools::getDeltaAngleDeg(recoPhiDeg, -aftForeOffset);
  Int_t phiSectorOfPeak = Int_t(angleThroughPhiSectors/PHI_RANGE);
  if(phiSectorOfPeak < 0 || phiSectorOfPeak >= NUM_PHI){
    std::cerr << "you idiot again  " << phiSectorOfPeak << "\t" << recoPhiDeg << std::endl;
  }

  int isMinBias = header->getTriggerBitSoftExt() | header->getTriggerBitG12() | header->getTriggerBitADU5();

  // UShort_t phiTrigMask = pol==AnitaPol::kVertical ? header->phiTrigMask : header->phiTrigMaskH;
  bool peakIsInMaskedPhiSector = false;

  bool wasAnL3Trigger = false;
  // Int_t deltaPhiSect = NUM_PHI/2;
  for(int phi=0; phi<NUM_PHI; phi++){

    // check against masked phi-sectors
    UInt_t phiMaskBit = header->isInL1Mask(phi, pol) || header->isInPhiMask(phi, pol);
    if(phiMaskBit > 0 && phiSectorOfPeak==phi){
      peakIsInMaskedPhiSector = true;
    }

    // check against L3 trigger
    UInt_t phiBit = RootTools::getBit(phi, header->getL3TrigPattern(pol));
    if(phiBit > 0){
      Int_t dPhiSect = phiSectorOfPeak - phi;
      if(dPhiSect < -NUM_PHI/2){
	dPhiSect += NUM_PHI;
      }
      else if(dPhiSect >= NUM_PHI/2){
	dPhiSect -= NUM_PHI;
      }
      if(TMath::Abs(dPhiSect) < deltaPhiSect){
	deltaPhiSect = TMath::Abs(dPhiSect);
	wasAnL3Trigger = true;
      }
    }
  }
  if(wasAnL3Trigger == true && deltaPhiSect >= NUM_PHI/2){
    std::cerr << "You bloody fool of a took" << wasAnL3Trigger << "\t" << deltaPhiSect << std::endl;
  }
  // if(deltaPhiSect >= NUM_PHI/2){
  // 	std::cerr << header->run << "\t" << header->eventNumber << std::endl;
  // }
  if(peakIsInMaskedPhiSector==true){
    status = kFail;
  }
  else if(isMinBias==0 && TMath::Abs(deltaPhiSect) > maxAbsDeltaPhiSect){
    status = kFail;
  }

  return status;
}



AnalysisCuts::Status_t AnalysisCuts::applySunPointingCut(Double_t deltaSolarPhiDeg){
  Status_t status = kPass;
  if(deltaSolarPhiDeg < -DEGREES_IN_CIRCLE/2 || deltaSolarPhiDeg >= DEGREES_IN_CIRCLE/2){
    std::cerr << "Warning! in " << __PRETTY_FUNCTION__
	      << " You've not normalized deltaSolarPhiDeg to lie in the range "
	      << -DEGREES_IN_CIRCLE/2 << " to " <<  DEGREES_IN_CIRCLE/2 << " degrees." << std::endl;
  }

  if(TMath::Abs(deltaSolarPhiDeg) < deltaSolarPhiDegCut){
    status = kFail;
  }
  return status;
}



AnalysisCuts::Status_t AnalysisCuts::applyThermalBackgroundCut(Double_t imagePeak, Double_t hilbertPeak, Double_t& fisher){

  Status_t status = kPass;
  fisher = fisherWeights[0] + imagePeak*fisherWeights[1] + hilbertPeak*fisherWeights[2];

  if(fisher < fisherCutVal){
    status = kFail;
  }

  return status;
}



AnalysisCuts::Status_t AnalysisCuts::applyImagePeakRatioCut(Double_t p1, Double_t p2, Double_t& peakRatio){

  Status_t status = kPass;
  peakRatio = p2/p1;
  if(peakRatio > maxImagePeakRatio){
    status = kFail;
  }

  return status;
}



AnalysisCuts::Status_t AnalysisCuts::applyThetaAngleCut(Double_t thetaDeg){
  // +ve theta is up in my convention

  Status_t status = kPass;
  if(thetaDeg > maxThetaDeg || thetaDeg < minThetaDeg){
    status = kFail;
  }

  return status;
}





AnalysisCuts::Status_t AnalysisCuts::applySurfSaturationCut(Double_t maxVolts[][NUM_SEAVEYS], Double_t minVolts[][NUM_SEAVEYS], Double_t& maxMaxVolts, Double_t& minMinVolts, Double_t& absSumMaxMin){
  Status_t status = kPass;

  maxMaxVolts = 0;
  minMinVolts = 0;
  for(int pol=0; pol <AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      if(maxVolts[pol][ant] > maxMaxVolts){
	maxMaxVolts = maxVolts[pol][ant];
      }

      if(minVolts[pol][ant] < minMinVolts){
	minMinVolts = minVolts[pol][ant];
      }
    }
  }
  absSumMaxMin = TMath::Abs(maxMaxVolts + minMinVolts);

  if(maxMaxVolts > maxVoltsLimit || minMinVolts < minVoltsLimit || absSumMaxMin > absMaxMinSumLimit){
    status = kFail;
  }

  return status;

}

AnalysisCuts::Status_t AnalysisCuts::applySurfSaturationCutBetter(Double_t theMaxVolts, Double_t theMinVolts, Double_t absSumMaxMin){
  Status_t status = kPass;

  if(theMaxVolts > maxVoltsLimit || theMinVolts < minVoltsLimit || absSumMaxMin > absMaxMinSumLimit){
    status = kFail;
  }

  return status;

}
