#include "InterferometricMapMaker.h"
#include "InterferometricMap.h"
#include "InterferometryCache.h"


InterferometricMapMaker::InterferometricMapMaker(){
  initializeInternals();
}

InterferometricMapMaker::~InterferometricMapMaker(){
  if(cc){
    delete cc;
  }
}

void InterferometricMapMaker::process(const FilteredAnitaEvent * ev, UsefulAdu5Pat* pat ,AnitaEventSummary * summary) const{
  std::cout << "Blank for now... " << ev << "\t" << pat << "\t" << summary << std::endl;
}



void InterferometricMapMaker::initializeInternals(){

  cc = new CrossCorrelator();

  
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  for(Int_t pol=0; pol < AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant<NUM_SEAVEYS; ant++){
      rArray[pol].push_back(geom->getAntR(ant, AnitaPol::AnitaPol_t(pol)));
      zArray[pol].push_back(geom->getAntZ(ant, AnitaPol::AnitaPol_t(pol)));
      phiArrayDeg[pol].push_back(geom->getAntPhiPositionRelToAftFore(ant, AnitaPol::AnitaPol_t(pol))*TMath::RadToDeg());
    }
  }

  coarseMaps[AnitaPol::kHorizontal] = NULL;//new InterferometricMap("h0H", "h0H", InterferometricMap::getBin0PhiDeg());
  coarseMaps[AnitaPol::kVertical] = NULL; //new InterferometricMap("h0V", "h0V", InterferometricMap::getBin0PhiDeg());
  
  fillDeltaTLookup();

  for(Int_t pol=0; pol < AnitaPol::kNotAPol; pol++){
    eventNumber[pol] = 0;    
    for(Int_t peakInd=0; peakInd < MAX_NUM_PEAKS; peakInd++){
      coarseMapPeakValues[pol][peakInd] = -9999;
      coarseMapPeakPhiDegs[pol][peakInd] = -9999;
      coarseMapPeakThetaDegs[pol][peakInd] = -9999;

      fineMapPeakValues[pol][peakInd] = -9999;
      fineMapPeakPhiDegs[pol][peakInd] = -9999;
      fineMapPeakThetaDegs[pol][peakInd] = -9999;
    }
  }

  mapModeNames[kGlobal] = "Global";
  mapModeNames[kTriggered] = "Triggered";
  zoomModeNames[kZoomedOut] = "";
  zoomModeNames[kZoomedIn] = "Zoom";


  maxDPhiDeg = 0;
  kUseOffAxisDelay = 1;
  coherentDeltaPhi = 0;


}





AnitaPol::AnitaPol_t InterferometricMapMaker::reconstructEventPeakPol(FilteredAnitaEvent* usefulEvent, Int_t numFinePeaks ,Int_t numCoarsePeaks){

  for(Int_t polInd = AnitaPol::kHorizontal; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t)polInd;

    cc->correlateEvent(usefulEvent, pol);    

    reconstruct(pol);
    coarseMaps[pol]->findPeakValues(numCoarsePeaks, coarseMapPeakValues[pol],
				    coarseMapPeakPhiDegs[pol], coarseMapPeakThetaDegs[pol]);

  }
  AnitaPol::AnitaPol_t peakPol = AnitaPol::kVertical;
  if(coarseMapPeakValues[AnitaPol::kHorizontal][0] > coarseMapPeakValues[AnitaPol::kVertical][0]){
    peakPol = AnitaPol::kHorizontal;
  }


  if(numFinePeaks > 0 && numFinePeaks <= numCoarsePeaks){
    for(Int_t peakInd=numFinePeaks-1; peakInd >= 0; peakInd--){

      reconstructZoom(peakPol, fineMapPeakValues[peakPol][peakInd],
		      fineMapPeakPhiDegs[peakPol][peakInd], fineMapPeakThetaDegs[peakPol][peakInd],
		      coarseMapPeakPhiDegs[peakPol][peakInd], coarseMapPeakThetaDegs[peakPol][peakInd],
		      peakInd);
    }
  }
  return peakPol;
}




void InterferometricMapMaker::reconstructEvent(FilteredAnitaEvent* usefulEvent, UsefulAdu5Pat& usefulPat, AnitaEventSummary* eventSummary){

  // reconstructEvent(usefulEvent, MAX_NUM_PEAKS, MAX_NUM_PEAKS);
  // const int thisNumPeaks = MAX_NUM_PEAKS;
  const int thisNumPeaks = 1;
  reconstructEvent(usefulEvent, thisNumPeaks, thisNumPeaks);

  for(Int_t polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;

    for(Int_t peakInd=0; peakInd < thisNumPeaks; peakInd++){

      Double_t coarsePeak, coarsePhi, coarseTheta;
      getCoarsePeakInfo(pol, peakInd, coarsePeak, coarsePhi, coarseTheta);

      getFinePeakInfo(pol, peakInd,
		      eventSummary->peak[pol][peakInd].value,
		      eventSummary->peak[pol][peakInd].phi,
		      eventSummary->peak[pol][peakInd].theta);

      // std::cout << eventSummary->eventNumber << "\t" << eventSummary->peak[pol][peakInd].value << "\t" << eventSummary->peak[pol][peakInd].phi << "\t" << eventSummary->peak[pol][peakInd].theta << std::endl;

      eventSummary->peak[pol][peakInd].dphi_rough = eventSummary->peak[pol][peakInd].phi - coarsePhi;
      eventSummary->peak[pol][peakInd].dtheta_rough = eventSummary->peak[pol][peakInd].theta - coarseTheta;


      Double_t minY = 0;
      TGraph* grZ0 = makeUpsampledCoherentlySummedWaveform(pol,
							   eventSummary->peak[pol][peakInd].phi,
							   eventSummary->peak[pol][peakInd].theta,
							   coherentDeltaPhi,
							   eventSummary->peak[pol][peakInd].snr);

      TGraph* grZ0Hilbert = FFTtools::getHilbertEnvelope(grZ0);

      RootTools::getMaxMin(grZ0Hilbert, eventSummary->coherent[pol][peakInd].peakHilbert, minY);

      delete grZ0;
      delete grZ0Hilbert;

      Double_t phiWave = eventSummary->peak[pol][peakInd].phi*TMath::DegToRad();
      Double_t thetaWave = -1*eventSummary->peak[pol][peakInd].theta*TMath::DegToRad();
      Double_t sourceLat, sourceLon, sourceAlt;
      int success = usefulPat.getSourceLonAndLatAtAlt(phiWave, thetaWave,
						      sourceLon, sourceLat, sourceAlt);

      if(success==1){
	eventSummary->peak[pol][peakInd].latitude = sourceLat;
	eventSummary->peak[pol][peakInd].longitude = sourceLon;
	eventSummary->peak[pol][peakInd].altitude = sourceAlt;
	eventSummary->peak[pol][peakInd].distanceToSource = SPEED_OF_LIGHT_NS*usefulPat.getTriggerTimeNsFromSource(sourceLat, sourceLon, sourceAlt);
      }
      else{
	eventSummary->peak[pol][peakInd].latitude = -9999;
	eventSummary->peak[pol][peakInd].longitude = -9999;
	eventSummary->peak[pol][peakInd].altitude = -9999;
	eventSummary->peak[pol][peakInd].distanceToSource = -9999;
      }
    }
  }
}


void InterferometricMapMaker::reconstructEvent(FilteredAnitaEvent* usefulEvent, Int_t numFinePeaks ,Int_t numCoarsePeaks){
  
  for(Int_t polInd = AnitaPol::kHorizontal; polInd < AnitaPol::kNotAPol; polInd++){
    eventNumber[polInd] = usefulEvent->eventNumber;
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t)polInd;

    cc->correlateEvent(usefulEvent, pol);

    reconstruct(pol);
    coarseMaps[pol]->findPeakValues(numCoarsePeaks, coarseMapPeakValues[pol],
				    coarseMapPeakPhiDegs[pol], coarseMapPeakThetaDegs[pol]);

    if(numFinePeaks > 0 && numFinePeaks <= numCoarsePeaks){
      for(Int_t peakInd=numFinePeaks-1; peakInd >= 0; peakInd--){

	// std::cerr << peakInd << "\t" << numFinePeaks << "\t" << numCoarsePeaks << "\t"
	// 	  << coarseMapPeakPhiDegs[pol][peakInd] << "\t" << coarseMapPeakThetaDegs[pol][peakInd]
	// 	  << std::endl;


	// reconstructZoom(pol, fineMapPeakValues[pol][peakInd],
	// 		fineMapPeakPhiDegs[pol][peakInd], fineMapPeakThetaDegs[pol][peakInd],
	// 		0, 0,
	// 		peakInd);

	reconstructZoom(pol, fineMapPeakValues[pol][peakInd],
			fineMapPeakPhiDegs[pol][peakInd], fineMapPeakThetaDegs[pol][peakInd],
			coarseMapPeakPhiDegs[pol][peakInd], coarseMapPeakThetaDegs[pol][peakInd],
			peakInd);

	// std::cerr << fineMapPeakValues[pol][peakInd] << "\t" << fineMapPeakPhiDegs[pol][peakInd] << "\t"
	// 	  << fineMapPeakThetaDegs[pol][peakInd] << std::endl << std::endl;
      }
    }
  }
}



void InterferometricMapMaker::getCoarsePeakInfo(AnitaPol::AnitaPol_t pol, Int_t peakIndex,
						Double_t& value, Double_t& phiDeg, Double_t& thetaDeg){

  if(peakIndex < MAX_NUM_PEAKS){
    value = coarseMapPeakValues[pol][peakIndex];
    phiDeg = coarseMapPeakPhiDegs[pol][peakIndex];
    thetaDeg = coarseMapPeakThetaDegs[pol][peakIndex];
  }
  else{
    std::cerr << "Warning in "<< __PRETTY_FUNCTION__ << " in " << __FILE__ << "."
 	      << "Requested peak info with index too large. peakIndex = "
	      << peakIndex << ", MAX_NUM_PEAKS = " << MAX_NUM_PEAKS << "." << std::endl;
    value = -999;
    phiDeg = -999;
    thetaDeg = -999;
  }
}



void InterferometricMapMaker::getFinePeakInfo(AnitaPol::AnitaPol_t pol, Int_t peakIndex,
					      Double_t& value, Double_t& phiDeg, Double_t& thetaDeg){

  if(peakIndex < MAX_NUM_PEAKS){
    value = fineMapPeakValues[pol][peakIndex];
    phiDeg = fineMapPeakPhiDegs[pol][peakIndex];
    thetaDeg = fineMapPeakThetaDegs[pol][peakIndex];
  }
  else{
    std::cerr << "Warning in "<< __PRETTY_FUNCTION__ << " in " << __FILE__ << "."
 	      << "Requested peak info with index too large. peakIndex = "
	      << peakIndex << ", MAX_NUM_PEAKS = " << MAX_NUM_PEAKS << "." << std::endl;
    value = -999;
    phiDeg = -999;
    thetaDeg = -999;
  }
}




Double_t InterferometricMapMaker::singleAntennaOffAxisDelay(Double_t deltaPhiDeg) {


  // These are the numbers from Linda's fits...
  //  FCN=2014.5 FROM HESSE     STATUS=NOT POSDEF     40 CALLS         407 TOTAL
  //  EDM=5.13972e-17    STRATEGY= 1      ERR MATRIX NOT POS-DEF
  //  EXT PARAMETER                APPROXIMATE        STEP         FIRST
  //  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
  //  1  p0           0.00000e+00     fixed
  //  2  p1           0.00000e+00     fixed
  //  3  p2          -1.68751e-05   1.90656e-07   4.07901e-10   1.39356e-03
  //  4  p3           0.00000e+00     fixed
  //  5  p4           2.77815e-08   9.38089e-11   7.26358e-14   2.34774e+01
  //  6  p5           0.00000e+00     fixed
  //  7  p6          -8.29351e-12   1.78286e-14   7.64605e-18  -1.72486e+05
  //  8  p7           0.00000e+00     fixed
  //  9  p8           1.15064e-15   1.78092e-18   6.93019e-22  -1.31237e+09
  //  10  p9           0.00000e+00     fixed
  //  11  p10         -7.71170e-20   1.63489e-22   6.05470e-26   4.32831e+13
  //  12  p11          0.00000e+00     fixed
  //  13  p12          1.99661e-24   9.79818e-27   1.84698e-29  -6.15528e+16

  // ... in a const array
  const Int_t numPowers = 13;
  const Double_t params[numPowers] = {0.00000e+00, 0.00000e+00, -1.68751e-05, 0.00000e+00,
				      2.77815e-08,  0.00000e+00, -8.29351e-12,  0.00000e+00,
				      1.15064e-15,  0.00000e+00, -7.71170e-20,  0.00000e+00,
				      1.99661e-24};

  // Sum up the powers in off boresight angle.
  Double_t offBoresightDelay = 0;
  for(int power=0; power < numPowers; power++){
    offBoresightDelay += params[power]*TMath::Power(deltaPhiDeg, power);
  }

  return offBoresightDelay;
}

Double_t InterferometricMapMaker::relativeOffAxisDelay(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
					       Double_t phiDeg) {

  Double_t deltaPhiDeg1 = RootTools::getDeltaAngleDeg(phiArrayDeg[pol].at(ant1), phiDeg);
  Double_t deltaPhiDeg2 = RootTools::getDeltaAngleDeg(phiArrayDeg[pol].at(ant2), phiDeg);

  // std::cout << phiDeg << "\t" << deltaPhiDeg1 << "\t" << deltaPhiDeg2 << std::endl;

  // Double_t leftPhi=x[0]+11.25;
  // Double_t rightPhi=x[0]-11.25;

  //leftPhi*=-1;
  //  rightPhi*=-1;
  // Double_t leftDelay=singleDelayFuncMod(&leftPhi,par);
  // Double_t rightDelay=singleDelayFuncMod(&rightPhi,par);

  Double_t delay1 = singleAntennaOffAxisDelay(deltaPhiDeg1);
  Double_t delay2 = singleAntennaOffAxisDelay(deltaPhiDeg2);

  // return (leftDelay-rightDelay);
  return (delay1-delay2);
}


Double_t InterferometricMapMaker::getDeltaTExpected(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t phiWave, Double_t thetaWave){

  // Double_t tanThetaW = tan(thetaWave);
  Double_t tanThetaW = tan(-1*thetaWave);
  Double_t part1 = zArray[pol].at(ant1)*tanThetaW - rArray[pol].at(ant1)*cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol].at(ant1));
  Double_t part2 = zArray[pol].at(ant2)*tanThetaW - rArray[pol].at(ant2)*cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol].at(ant2));
  Double_t tdiff = 1e9*((cos(thetaWave) * (part2 - part1))/SPEED_OF_LIGHT); // Returns time in ns

  return tdiff;
}



void InterferometricMapMaker::insertPhotogrammetryGeometry(){
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->usePhotogrammetryNumbers(1);
  for(Int_t pol=0; pol < AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant<NUM_SEAVEYS; ant++){
      rArray[pol].at(ant) = geom->getAntR(ant, AnitaPol::AnitaPol_t(pol));
      zArray[pol].at(ant) = geom->getAntZ(ant, AnitaPol::AnitaPol_t(pol));
      phiArrayDeg[pol].at(ant) = geom->getAntPhiPositionRelToAftFore(ant, AnitaPol::AnitaPol_t(pol))*TMath::RadToDeg();
    }
  }

  AnitaEventCalibrator* cal = AnitaEventCalibrator::Instance();
  for(int surf=0; surf < NUM_SURF; surf++){
    for(int chan=0; chan < NUM_CHAN; chan++){
      cal->relativePhaseCenterToAmpaDelays[surf][chan] = 0;
    }
  }

  fillDeltaTLookup();
  geom->usePhotogrammetryNumbers(0);

}







void InterferometricMapMaker::fillDeltaTLookup(){

  dtCache.populateCache(cc, this);
  dtCache.populateFineCache(cc, this);  
  std::cout << "done populate fine cache" << std::endl;

  
  minThetaDegZoom = MIN_THETA - THETA_RANGE_ZOOM/2;
  minPhiDegZoom = InterferometricMap::getBin0PhiDeg() - PHI_RANGE_ZOOM/2;

}




InterferometricMap* InterferometricMapMaker::getMap(AnitaPol::AnitaPol_t pol){

  return coarseMaps[pol];
}





TH2D* InterferometricMapMaker::getZoomMap(AnitaPol::AnitaPol_t pol, Int_t peakInd){

  TString name = "hZoom";
  name += pol == AnitaPol::kVertical ? "ImageV" : "ImageH";
  name += TString::Format("%u", eventNumber[pol]);

  TString title = TString::Format("Event %u ", eventNumber[pol]);
  title += (pol == AnitaPol::kVertical ? "VPOL" : "HPOL");
  title += " Zoomed In";
  title += " Map";

  // Here I'm hacking my socks off to make the map elevation = -1*theta,
  // where theta is the internal class representation of the vertical angle
  TH2D* hImage = new TH2D(name, title,
			  NUM_BINS_PHI_ZOOM, zoomPhiMin[pol], zoomPhiMin[pol] + PHI_RANGE_ZOOM,
			  NUM_BINS_THETA_ZOOM, zoomThetaMin[pol], zoomThetaMin[pol] + THETA_RANGE_ZOOM);

  hImage->GetXaxis()->SetTitle("Azimuth (Degrees)");
  hImage->GetYaxis()->SetTitle("Elevation (Degrees)");

  for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA_ZOOM; thetaBin++){
    for(Int_t phiBin = 0; phiBin < NUM_BINS_PHI_ZOOM; phiBin++){
      hImage->SetBinContent(phiBin+1, thetaBin+1, fineMap[pol][peakInd][thetaBin][phiBin]);
    }
  }

  return hImage;
}







TH2D* InterferometricMapMaker::makeZoomedImage(AnitaPol::AnitaPol_t pol,
				       Double_t& imagePeak, Double_t& peakPhiDeg, Double_t& peakThetaDeg,
				       Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg){

  reconstructZoom(pol, imagePeak, peakPhiDeg, peakThetaDeg,
		  zoomCenterPhiDeg, zoomCenterThetaDeg);
  return getZoomMap(pol);

}

TH2D* InterferometricMapMaker::makeZoomedImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
				       Double_t& peakThetaDeg, UShort_t l3TrigPattern,
				       Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg){
  reconstructZoom(pol, imagePeak, peakPhiDeg, peakThetaDeg,
		  zoomCenterPhiDeg, zoomCenterThetaDeg);
  return getZoomMap(pol);

}


TH2D* InterferometricMapMaker::makeZoomedImage(AnitaPol::AnitaPol_t pol, UShort_t l3TrigPattern,
				       Double_t zoomCenterPhiDeg,Double_t zoomCenterThetaDeg){

  Double_t imagePeak, peakPhiDeg, peakThetaDeg;
  // return makeZoomedImage(pol, imagePeak, peakPhiDeg, peakThetaDeg,
  // 			 l3TrigPattern, zoomCenterPhiDeg, zoomCenterThetaDeg);
  reconstructZoom(pol, imagePeak, peakPhiDeg, peakThetaDeg,
		  zoomCenterPhiDeg, zoomCenterThetaDeg);
  return getZoomMap(pol);
}

void InterferometricMapMaker::reconstruct(AnitaPol::AnitaPol_t pol){

  makerOwnsMap[pol] = false;
  if(!makerOwnsMap[pol] || coarseMaps[pol]==NULL)
  {
    // I don't own the map or there isn't one, so I'll make a new one
    TString name = "h";
    name += pol == AnitaPol::kVertical ? "ImageV" : "ImageH";
    name += TString::Format("%u", eventNumber[pol]);

    TString title = TString::Format("Event %u ", eventNumber[pol]);
    title += (pol == AnitaPol::kVertical ? "VPOL" : "HPOL");
    title += " Map";

    coarseMaps[pol] = new InterferometricMap(name, title);
  }  
  else{// if(makerOwnsMap[pol]){  {
    // I own the map so I can overwrite it  
    for(int phiBin=1; phiBin<=coarseMaps[pol]->GetNbinsPhi(); phiBin++){
      for(int thetaBin=1; thetaBin<=coarseMaps[pol]->GetNbinsTheta(); thetaBin++){
	coarseMaps[pol]->SetBinContent(phiBin, thetaBin, 0);	
      }
    }
  }

  coarseMaps[pol]->Fill(pol, cc, &dtCache);
}


void InterferometricMapMaker::reconstructZoom(AnitaPol::AnitaPol_t pol, Double_t& imagePeak,
					      Double_t& peakPhiDeg, Double_t& peakThetaDeg,
					      Double_t zoomCenterPhiDeg,
					      Double_t zoomCenterThetaDeg,
					      Int_t peakIndex){

  // Some kind of sanity check here due to the unterminating while loop inside RootTools::getDeltaAngleDeg
  if(zoomCenterPhiDeg < -500 || zoomCenterThetaDeg < -500 ||
     zoomCenterPhiDeg >= 500 || zoomCenterThetaDeg >= 500){

    std::cerr << "Warning in "<< __PRETTY_FUNCTION__ << " in " << __FILE__
	      << ". zoomCenterPhiDeg = " << zoomCenterPhiDeg
	      << " and zoomCenterThetaDeg = " << zoomCenterThetaDeg
	      << " ...these values look suspicious so I'm skipping this reconstruction."
	      << " eventNumber = " << eventNumber[pol] << std::endl;
    imagePeak = -9999;
    peakPhiDeg = -9999;
    peakThetaDeg = -9999;
    return;
  }

  
  Double_t deltaPhiDegPhi0 = RootTools::getDeltaAngleDeg(zoomCenterPhiDeg, InterferometricMap::getBin0PhiDeg());
  deltaPhiDegPhi0 = deltaPhiDegPhi0 < 0 ? deltaPhiDegPhi0 + DEGREES_IN_CIRCLE : deltaPhiDegPhi0;

  Int_t phiSector = floor(deltaPhiDegPhi0/PHI_RANGE);
  cc->doUpsampledCrossCorrelations(pol, phiSector);


  zoomCenterPhiDeg = (TMath::Nint(zoomCenterPhiDeg/ZOOM_BIN_SIZE_PHI))*ZOOM_BIN_SIZE_PHI;
  zoomCenterThetaDeg = (TMath::Nint(zoomCenterThetaDeg/ZOOM_BIN_SIZE_THETA))*ZOOM_BIN_SIZE_THETA;

  zoomPhiMin[pol] = zoomCenterPhiDeg - PHI_RANGE_ZOOM/2;
  zoomThetaMin[pol] = zoomCenterThetaDeg - THETA_RANGE_ZOOM/2;

  std::vector<Int_t>* combosToUse = &cc->combosToUseGlobal[phiSector];

  Int_t phiZoomBase = TMath::Nint((zoomPhiMin[pol] - minPhiDegZoom)/ZOOM_BIN_SIZE_PHI);
  Int_t thetaZoomBase = TMath::Nint((zoomThetaMin[pol] - minThetaDegZoom)/ZOOM_BIN_SIZE_THETA);

  for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA_ZOOM; thetaBin++){
    for(Int_t phiBin = 0; phiBin < NUM_BINS_PHI_ZOOM; phiBin++){
      fineMap[pol][peakIndex][thetaBin][phiBin]=0;
    }
  }

  const Int_t offset = cc->numSamplesUpsampled/2;
  for(UInt_t comboInd=0; comboInd<combosToUse->size(); comboInd++){
    Int_t combo = combosToUse->at(comboInd);
    if(cc->kOnlyThisCombo >= 0 && combo!=cc->kOnlyThisCombo){
      continue;
    }
    int ant1 = cc->comboToAnt1s[combo];
    int ant2 = cc->comboToAnt2s[combo];
    for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA_ZOOM; thetaBin++){
      Int_t zoomThetaInd = thetaZoomBase + thetaBin;
      // Double_t zoomThetaWave = zoomedThetaWaves[zoomThetaInd];
      // Double_t partBA = partBAsZoom[pol][combo][zoomThetaInd];
      Double_t partBA = dtCache.partBAsZoom[dtCache.partBAsIndex(pol, combo, zoomThetaInd)]; //)[pol][combo][zoomThetaInd];      
      // Double_t dtFactor = dtFactors[zoomThetaInd];
      Double_t dtFactor = dtCache.dtFactors[zoomThetaInd];      

      for(Int_t phiBin = 0; phiBin < NUM_BINS_PHI_ZOOM; phiBin++){
	Int_t zoomPhiInd = phiZoomBase + phiBin;
	// Double_t zoomPhiWave = zoomedPhiWaveLookup[zoomPhiInd];

	int p21 = dtCache.part21sIndex(pol, combo, zoomPhiInd);
	Double_t offsetLowDouble = dtFactor*(partBA - dtCache.part21sZoom[p21]);//[pol][combo][zoomPhiInd]);		

	// Double_t offsetLowDouble = dtFactor*(partBA - part21sZoom[pol][combo][zoomPhiInd]);
	// Double_t offsetLowDouble = dtFactor*(partBA - part21sZoom[pol][combo][zoomPhiInd]);	
	// offsetLowDouble += kUseOffAxisDelay > 0 ? offAxisDelaysDivided[pol][combo][zoomPhiInd] : 0;
	offsetLowDouble += kUseOffAxisDelay > 0 ? dtCache.offAxisDelaysDivided[p21] : 0;
	
	offsetLowDouble += cc->startTimes[pol][ant1]/cc->correlationDeltaT;
	offsetLowDouble -= cc->startTimes[pol][ant2]/cc->correlationDeltaT;
	

	// hack for floor()
	Int_t offsetLow = (int) offsetLowDouble - (offsetLowDouble < (int) offsetLowDouble);
	// Double_t deltaT = getDeltaTExpected(pol, ant1, ant2, zoomPhiWave, zoomThetaWave);

	// Int_t offsetLow = floor(deltaT/cc->correlationDeltaT);
	// offsetLow += offset;

	// deltaT -= offsetLow*cc->correlationDeltaT;

	Double_t deltaT = (offsetLowDouble - offsetLow);
	offsetLow += offset;
	Double_t c1 = cc->crossCorrelationsUpsampled[pol][combo][offsetLow];
	Double_t c2 = cc->crossCorrelationsUpsampled[pol][combo][offsetLow+1];
	Double_t cInterp = deltaT*(c2 - c1) + c1;

	// std::cout << pol << "\t" << peakIndex << "\t"
	// 	  << thetaBin << "\t" << phiBin << "\t"
	// 	  << zoomThetaWave << "\t" << zoomPhiWave << "\t"
	// 	  << deltaT << "\t" << c1 << std::endl;
	// fineMap[pol][peakIndex][thetaBin][phiBin] += c1;
	fineMap[pol][peakIndex][thetaBin][phiBin] += cInterp;
      }
    }
  }
  
  // const Int_t offset = cc->numSamplesUpsampled/2;
  // for(UInt_t comboInd=0; comboInd<combosToUse->size(); comboInd++){
  //   Int_t combo = combosToUse->at(comboInd);
  //   if(cc->kOnlyThisCombo >= 0 && combo!=cc->kOnlyThisCombo){
  //     continue;
  //   }
  //   for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA_ZOOM; thetaBin++){
  //     Int_t zoomThetaInd = thetaZoomBase + thetaBin;
  //     Double_t partBA = partBAsZoom[pol][combo][zoomThetaInd];
  //     // Double_t dtFactor = zoomedCosThetaWaves[zoomThetaInd]/SPEED_OF_LIGHT_NS;
  //     Double_t dtFactor = dtFactors[zoomThetaInd];
  //     for(Int_t phiBin = 0; phiBin < NUM_BINS_PHI_ZOOM; phiBin++){
  // 	Int_t zoomPhiInd = phiZoomBase + phiBin;
  // 	Double_t offsetLowDouble = dtFactor*(partBA - part21sZoom[pol][combo][zoomPhiInd]);
  // 	offsetLowDouble += kUseOffAxisDelay > 0 ? offAxisDelaysDivided[pol][combo][zoomPhiInd] : 0;
  // 	// hack for floor()
  // 	Int_t offsetLow = (int) offsetLowDouble - (offsetLowDouble < (int) offsetLowDouble);

  // 	Double_t deltaT = (offsetLowDouble - offsetLow);
  // 	offsetLow += offset;
  // 	Double_t c1 = cc->crossCorrelationsUpsampled[pol][combo][offsetLow];
  // 	Double_t c2 = cc->crossCorrelationsUpsampled[pol][combo][offsetLow+1];
  // 	Double_t cInterp = deltaT*(c2 - c1) + c1;

  // 	fineMap[pol][peakIndex][thetaBin][phiBin] += cInterp;
  //     }
  //   }
  // }

  Double_t normFactor = cc->kOnlyThisCombo < 0 && combosToUse->size() > 0 ? combosToUse->size() : 1;
  // absorb the removed inverse FFT normalization
  normFactor*=(cc->numSamples*cc->numSamples);

  // set peak finding variables
  imagePeak = -DBL_MAX;
  Int_t peakPhiBin = -1;
  Int_t peakThetaBin = -1;
  
  for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA_ZOOM; thetaBin++){
    for(Int_t phiBin = 0; phiBin < NUM_BINS_PHI_ZOOM; phiBin++){
      fineMap[pol][peakIndex][thetaBin][phiBin] /= normFactor;
      // std::cout << pol << "\t" << thetaBin << "\t" << phiBin << "\t" << fineMap[pol][peakIndex][thetaBin][phiBin]  << std::endl;

      if(fineMap[pol][peakIndex][thetaBin][phiBin] > imagePeak){	
	imagePeak = fineMap[pol][peakIndex][thetaBin][phiBin];
	peakPhiBin = phiBin;
	peakThetaBin = thetaBin;
      }
    }
  }

  peakPhiDeg = zoomPhiMin[pol] + peakPhiBin*ZOOM_BIN_SIZE_PHI;
  peakThetaDeg = zoomThetaMin[pol] + peakThetaBin*ZOOM_BIN_SIZE_THETA;
}




Int_t InterferometricMapMaker::directlyInsertGeometry(TString pathToLindasFile, AnitaPol::AnitaPol_t pol){

  // Since I am simulataneously testing many of Linda's geometries on lots of different files
  // I need the help of a machine to check I'm testing the geometry I think I'm testing.

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->usePhotogrammetryNumbers(0); // i.e. definitely use the numbers I am inserting.
  AnitaEventCalibrator* cal = AnitaEventCalibrator::Instance();


  std::ifstream lindasNums(pathToLindasFile.Data());
  if(lindasNums.is_open()==0){
    return 1; // This is an error
  }

  Int_t ant;
  Double_t dr, dPhiRad, dz, dt;

  while(lindasNums >> ant >> dr >> dz >> dPhiRad >> dt){

    Int_t surf, chan, ant2;
    geom->getSurfChanAntFromRingPhiPol(AnitaRing::AnitaRing_t (ant/NUM_PHI), ant%NUM_PHI, pol,
				       surf, chan, ant2);

    Double_t newR = geom->rPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] + dr;
    Double_t newPhi = geom->azPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] + dPhiRad;
    Double_t newZ = geom->zPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] + dz;
    Double_t newT = dt;

    if(newPhi >= TMath::TwoPi()){
      newPhi -= TMath::TwoPi();
    }
    else if(newPhi < 0){
      newPhi += TMath::TwoPi();
    }

    geom->rPhaseCentreFromVerticalHorn[ant][pol] = newR;
    geom->azPhaseCentreFromVerticalHorn[ant][pol] = newPhi;
    geom->zPhaseCentreFromVerticalHorn[ant][pol] = newZ;
    cal->relativePhaseCenterToAmpaDelays[surf][chan] = newT;

    // std::cout << ant << "\t" << dr << "\t" << dz << "\t" << dPhiRad << "\t" << dt << std::endl;

  }
  if(ant != NUM_SEAVEYS - 1){
    // then you didn't make it to the end of the file!
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << " in " << __FILE__ << std::endl;
    std::cerr << "It looks like you didn't read to the end of Linda's geometry file!" << std::endl;
    return 1;
  }

  return 0;
}

Int_t InterferometricMapMaker::getPhiSectorOfAntennaClosestToPhiDeg(AnitaPol::AnitaPol_t pol, Double_t phiDeg){
  Int_t phiSectorOfPeak = -1;
  Double_t bestDeltaPhiOfPeakToAnt = DEGREES_IN_CIRCLE;
  for(int ant=0; ant < NUM_SEAVEYS; ant++){
    Double_t phiOfAnt = phiArrayDeg[pol].at(ant);
    Double_t deltaPhiOfPeakToAnt = TMath::Abs(RootTools::getDeltaAngleDeg(phiOfAnt, phiDeg));
    if(deltaPhiOfPeakToAnt < bestDeltaPhiOfPeakToAnt){
      bestDeltaPhiOfPeakToAnt = deltaPhiOfPeakToAnt;
      phiSectorOfPeak = (ant % NUM_PHI);
    }
  }
  return phiSectorOfPeak;
}




TGraph* InterferometricMapMaker::makeUpsampledCoherentlySummedWaveform(AnitaPol::AnitaPol_t pol, Double_t phiDeg,
							       Double_t thetaDeg, Int_t maxDeltaPhiSect,
							       Double_t& snr){
  return makeCoherentWorker(pol, phiDeg, thetaDeg, maxDeltaPhiSect, snr, cc->numSamplesUpsampled);
}


TGraph* InterferometricMapMaker::makeCoherentlySummedWaveform(AnitaPol::AnitaPol_t pol, Double_t phiDeg,
						      Double_t thetaDeg, Int_t maxDeltaPhiSect,
						      Double_t& snr){
  return makeCoherentWorker(pol, phiDeg, thetaDeg, maxDeltaPhiSect, snr, cc->numSamples);
}


TGraph* InterferometricMapMaker::makeCoherentWorker(AnitaPol::AnitaPol_t pol, Double_t phiDeg,
					    Double_t thetaDeg, Int_t maxDeltaPhiSect,
					    Double_t& snr,
					    Int_t nSamp){
  return NULL;
}
