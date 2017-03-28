#include "InterferometricMapMaker.h"


InterferometricMapMaker::InterferometricMapMaker(){
  initializeInternals();
}

InterferometricMapMaker::~InterferometricMapMaker(){
  
}

void InterferometricMapMaker::process(const FilteredAnitaEvent * ev, UsefulAdu5Pat* pat ,AnitaEventSummary * summary) const{
  std::cout << "Blanks for now... " << ev << "\t" << pat << "\t" << summary << std::endl;
}



void InterferometricMapMaker::initializeInternals(){

  cc = new CrossCorrelator();    
  
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  for(Int_t pol=0; pol < NUM_POL; pol++){
    for(int ant=0; ant<NUM_SEAVEYS; ant++){
      rArray[pol].push_back(geom->getAntR(ant, AnitaPol::AnitaPol_t(pol)));
      zArray[pol].push_back(geom->getAntZ(ant, AnitaPol::AnitaPol_t(pol)));
      phiArrayDeg[pol].push_back(geom->getAntPhiPositionRelToAftFore(ant, AnitaPol::AnitaPol_t(pol))*TMath::RadToDeg());
    }
  }


  aftForeOffset = geom->aftForeOffsetAngleVertical*TMath::RadToDeg(); //phiArrayDeg[0].at(0);

  fillDeltaTLookup();

  for(Int_t pol=0; pol < NUM_POL; pol++){
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




Double_t InterferometricMapMaker::getBin0PhiDeg(){

  Double_t phi0 = -aftForeOffset;
  if(phi0 < -DEGREES_IN_CIRCLE/2){
    phi0+=DEGREES_IN_CIRCLE;
  }
  else if(phi0 >= DEGREES_IN_CIRCLE/2){
    phi0-=DEGREES_IN_CIRCLE;
  }
  return phi0 - PHI_RANGE/2;
}






void InterferometricMapMaker::findPeakValues(AnitaPol::AnitaPol_t pol, Int_t numPeaks, Double_t* peakValues,
				     Double_t* phiDegs, Double_t* thetaDegs){

  // In this function I want to find numPeak peaks and set an exclusion zone around each peak
  // so that the next peak isn't just a neighbour of a true peak.

  if(numPeaks > MAX_NUM_PEAKS){
    // You can have numPeaks less than or requal MAX_NUM_PEAKS, but not greater than.
    std::cerr << "Warning in "<< __PRETTY_FUNCTION__ << " in " << __FILE__
	      << ". numPeaks = " << numPeaks << ", InterferometricMapMaker compiled with MAX_NUM_PEAKS  = "
	      << MAX_NUM_PEAKS << ", setting numPeaks = " << MAX_NUM_PEAKS << std::endl;
    numPeaks = MAX_NUM_PEAKS;
  }

  // Set not crazy, but still debug visible values for peak values/location
  // Believe -DBL_MAX was causing me some while loop issues in RootTools::getDeltaAngleDeg(...)
  // which has an unrestricted while loop inside.
  for(Int_t peakInd=0; peakInd < numPeaks; peakInd++){
    peakValues[peakInd] = -999;
    phiDegs[peakInd] = -999;
    thetaDegs[peakInd] = -999;
  }

  // As we start, all regions are allowed.
  Int_t allowedBins[NUM_BINS_PHI*NUM_PHI][NUM_BINS_THETA];
  for(Int_t phiBin=0; phiBin<NUM_BINS_PHI*NUM_PHI; phiBin++){
    for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA; thetaBin++){
      allowedBins[phiBin][thetaBin] = 1;
    }
  }

  for(Int_t peakInd=0; peakInd < numPeaks; peakInd++){
    Int_t gotHere = 0;
    for(Int_t phiBin=0; phiBin<NUM_BINS_PHI*NUM_PHI; phiBin++){
      for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA; thetaBin++){

	if(allowedBins[phiBin][thetaBin] > 0){
	  // then we are in an allowed region

	  Double_t mapValue = coarseMap[pol][phiBin][thetaBin];
	  if(mapValue > peakValues[peakInd]){
	    // higher than the current highest map value

	    gotHere = 1;
	    // Time for something a little tricky to read...
	    // I want to stop bins that are on the edge of allowed regions registering as peaks
	    // if the bins just inside the disallowed regions are higher valued.
	    // To ensure a local maxima so I am going to ensure that all neighbouring
	    // bins are less than the current value.
	    // Checking the 3x3 grid around the current bin (wrapping in phi).

	    // Wrap phi edges
	    Int_t lastPhiBin = phiBin != 0 ? phiBin - 1 : (NUM_BINS_PHI*NUM_PHI) - 1;
	    Int_t nextPhiBin = phiBin != (NUM_BINS_PHI*NUM_PHI) - 1 ? phiBin + 1 : 0;


	    // Is current bin higher than phi neighbours in same theta row?
	    if(coarseMap[pol][lastPhiBin][thetaBin] < mapValue &&
	       coarseMap[pol][nextPhiBin][thetaBin] < mapValue){

	      gotHere = 2;
	      // So far so good, now check theta neigbours
	      // This doesn't wrap, so I will allow edge cases to pass as maxima
	      Int_t nextThetaBin = thetaBin != (NUM_BINS_THETA) - 1 ? thetaBin+1 : thetaBin;
	      Int_t lastThetaBin = thetaBin != 0 ? thetaBin-1 : thetaBin;

	      // Check theta bins below
	      gotHere = 3;
	      if(thetaBin == 0 || (coarseMap[pol][lastPhiBin][lastThetaBin] < mapValue &&
				   coarseMap[pol][phiBin][lastThetaBin] < mapValue &&
				   coarseMap[pol][nextPhiBin][lastThetaBin] < mapValue)){

		// Check theta bins above
		gotHere = 4;
		if(thetaBin == NUM_BINS_THETA-1 || (coarseMap[pol][lastPhiBin][nextThetaBin] < mapValue &&
						    coarseMap[pol][phiBin][nextThetaBin] < mapValue &&
						    coarseMap[pol][nextPhiBin][nextThetaBin] < mapValue)){

		  // Okay okay okay... you're a local maxima, you can be my next peak.
		  peakValues[peakInd] = coarseMap[pol][phiBin][thetaBin];
		  phiDegs[peakInd] = phiWaveLookup[phiBin]*TMath::RadToDeg();
		  thetaDegs[peakInd] = thetaWaves[thetaBin]*TMath::RadToDeg();
		}
	      }
	    }
	  }
	}
      } // thetaBin
    } // phiBin

    // Still inside the peakInd loop

    // Here I set the exclusion zone around the peak bin.
    if(peakValues[peakInd] >= -999){ // Checks that a peak was found, probably unnecessary

      for(Int_t phiBin=0; phiBin<NUM_BINS_PHI*NUM_PHI; phiBin++){

	Double_t phiDeg = phiWaveLookup[phiBin]*TMath::RadToDeg();
	Double_t absDeltaPhi = TMath::Abs(RootTools::getDeltaAngleDeg(phiDegs[peakInd], phiDeg));
	if(absDeltaPhi < PEAK_PHI_DEG_RANGE){

	  for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA; thetaBin++){

	    Double_t thetaDeg = thetaWaves[thetaBin]*TMath::RadToDeg();
	    Double_t absDeltaTheta = TMath::Abs(RootTools::getDeltaAngleDeg(thetaDegs[peakInd], thetaDeg));
	    if(absDeltaTheta < PEAK_THETA_DEG_RANGE){

	      // Now this region in phi/theta is disallowed on the next loop through...
	      allowedBins[phiBin][thetaBin] = 0;
	    }
	  }
	}
      }
    }
    if(peakValues[peakInd] < 0){
      std::cerr << "Peak " << peakInd << " = " << peakValues[peakInd] << "\t" << gotHere << std::endl;
      for(int pi=0; pi <= peakInd; pi++){
    	std::cerr << peakValues[pi] << "\t" << phiDegs[pi] << "\t" << thetaDegs[pi] << std::endl;
      }
    }
  }
}




template <class NiceAnitaEvent> // needs eventNumber member and getGraph(int ant, AnitaPol::AnitaPol_t pol)
AnitaPol::AnitaPol_t InterferometricMapMaker::reconstructEventPeakPol(NiceAnitaEvent* usefulEvent, Int_t numFinePeaks ,Int_t numCoarsePeaks){

  for(Int_t polInd = AnitaPol::kHorizontal; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t)polInd;

    cc->correlateEvent(usefulEvent, pol);

    reconstruct(pol, coarseMapPeakValues[pol][0], coarseMapPeakPhiDegs[pol][0], coarseMapPeakThetaDegs[pol][0]);


    findPeakValues(pol, numCoarsePeaks, coarseMapPeakValues[pol],
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




template <class NiceAnitaEvent> // needs eventNumber member and getGraph(int ant, AnitaPol::AnitaPol_t pol)
void InterferometricMapMaker::reconstructEvent(NiceAnitaEvent* usefulEvent, UsefulAdu5Pat& usefulPat, AnitaEventSummary* eventSummary){

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



template <class NiceAnitaEvent> // needs eventNumber member and getGraph(int ant, AnitaPol::AnitaPol_t pol)
void InterferometricMapMaker::reconstructEvent(NiceAnitaEvent* usefulEvent, Int_t numFinePeaks ,Int_t numCoarsePeaks){

  for(Int_t polInd = AnitaPol::kHorizontal; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t)polInd;

    cc->correlateEvent(usefulEvent, pol);

    reconstruct(pol, coarseMapPeakValues[pol][0], coarseMapPeakPhiDegs[pol][0], coarseMapPeakThetaDegs[pol][0]);    
    findPeakValues(pol, numCoarsePeaks, coarseMapPeakValues[pol],
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
  for(Int_t pol=0; pol < NUM_POL; pol++){
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

  Double_t phi0 = getBin0PhiDeg();
  const Double_t phiBinSize = Double_t(PHI_RANGE)/NUM_BINS_PHI;
  for(Int_t phiIndex=0; phiIndex < NUM_BINS_PHI*NUM_PHI; phiIndex++){
    Double_t phiDeg = phi0 + phiIndex*phiBinSize;
    Double_t phiWave = TMath::DegToRad()*phiDeg;
    phiWaveLookup[phiIndex] = phiWave;
  }

  // const Double_t thetaBinSize = (Double_t(THETA_RANGE)/NUM_BINS_THETA);
  const Double_t thetaBinSize = (Double_t(MAX_THETA - MIN_THETA)/NUM_BINS_THETA);
  for(Int_t thetaIndex=0; thetaIndex < NUM_BINS_THETA; thetaIndex++){
    Double_t thetaWaveDeg = MIN_THETA + thetaIndex*thetaBinSize;
    // Double_t thetaWaveDeg = (thetaIndex-NUM_BINS_THETA/2)*thetaBinSize;
    Double_t thetaWave = thetaWaveDeg*TMath::DegToRad();
    thetaWaves[thetaIndex] = thetaWave;
  }

  int numCombos = cc->numCombos;

  for(Int_t polInd=0; polInd<NUM_POL; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(Int_t combo=0; combo<numCombos; combo++){
      Int_t ant1 = cc->comboToAnt1s.at(combo);
      Int_t ant2 = cc->comboToAnt2s.at(combo);

      for(Int_t phiSector = 0; phiSector<NUM_PHI; phiSector++){
	for(Int_t phiInd = 0; phiInd < NUM_BINS_PHI; phiInd++){
	  Int_t phiBin = phiSector*NUM_BINS_PHI + phiInd;
	  Double_t phiWave = phiWaveLookup[phiBin];

	  for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA; thetaBin++){
	    Double_t thetaWave = thetaWaves[thetaBin];
	    Double_t deltaT = getDeltaTExpected(pol, ant1, ant2, phiWave, thetaWave);
  	    // Int_t offset = TMath::Nint(deltaT/nominalSamplingDeltaT);
  	    // deltaTs[pol][phiBin][thetaBin][combo] = deltaT; //offset;
  	    // deltaTs[pol][combo][thetaBin][phiBin] = deltaT; //offset;
	    Int_t offsetLow = floor(deltaT/cc->nominalSamplingDeltaT);
	    offsetLows[pol][combo][phiBin][thetaBin] = offsetLow;
	    // offsetLows[pol][phiBin][combo][thetaBin] = offsetLow;
	    Double_t dt1 = offsetLow*cc->nominalSamplingDeltaT;
	    interpPreFactors[pol][combo][phiBin][thetaBin] = (deltaT - dt1)/cc->nominalSamplingDeltaT;
	    // interpPreFactors[pol][phiBin][combo][thetaBin] = (deltaT - dt1)/nominalSamplingDeltaT;

	    // Here we account for the fact that we are now time ordering the correlations
	    offsetLows[pol][combo][phiBin][thetaBin]+=cc->numSamples/2;
	    // offsetLows[pol][phiBin][combo][thetaBin]+=numSamples/2;
	  }
  	}
      }
    }
  }

  // minThetaDegZoom = -78.5;
  // minPhiDegZoom = -59.25;
  // minThetaDegZoom = -THETA_RANGE/2 - THETA_RANGE_ZOOM/2;
  minThetaDegZoom = MIN_THETA - THETA_RANGE_ZOOM/2;
  minPhiDegZoom = getBin0PhiDeg() - PHI_RANGE_ZOOM/2;

  for(Int_t thetaIndex=0; thetaIndex < NUM_BINS_THETA_ZOOM_TOTAL; thetaIndex++){
    Double_t thetaWaveDeg = minThetaDegZoom + thetaIndex*ZOOM_BIN_SIZE_THETA;
    Double_t thetaWave = -1*thetaWaveDeg*TMath::DegToRad();
    zoomedThetaWaves[thetaIndex] = thetaWave;
    zoomedTanThetaWaves[thetaIndex] = tan(thetaWave);
    zoomedCosThetaWaves[thetaIndex] = cos(thetaWave);
    dtFactors[thetaIndex] = zoomedCosThetaWaves[thetaIndex]/(SPEED_OF_LIGHT_NS*cc->correlationDeltaT);
  }

  for(Int_t pol=0; pol<AnitaPol::kNotAPol; pol++){
    for(Int_t combo=0; combo < NUM_COMBOS; combo++){
      Int_t ant1 = cc->comboToAnt1s.at(combo);
      Int_t ant2 = cc->comboToAnt2s.at(combo);
      for(Int_t thetaIndex=0; thetaIndex < NUM_BINS_THETA_ZOOM_TOTAL; thetaIndex++){
	partBAsZoom[pol][combo][thetaIndex] = zoomedTanThetaWaves[thetaIndex]*(zArray[pol].at(ant2)-zArray[pol].at(ant1));
      }
    }
  }

  for(Int_t pol=0; pol<AnitaPol::kNotAPol; pol++){
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      for(Int_t phiIndex=0; phiIndex < NUM_BINS_PHI_ZOOM_TOTAL; phiIndex++){
	Double_t phiDeg = minPhiDegZoom + phiIndex*ZOOM_BIN_SIZE_PHI;
	Double_t phiWave = phiDeg*TMath::DegToRad();
	zoomedCosPartLookup[pol][ant][phiIndex] = rArray[pol].at(ant)*cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol].at(ant));
      }
    }

    for(Int_t combo=0; combo < NUM_COMBOS; combo++){
      Int_t ant1 = cc->comboToAnt1s.at(combo);
      Int_t ant2 = cc->comboToAnt2s.at(combo);

      for(Int_t phiIndex=0; phiIndex < NUM_BINS_PHI_ZOOM_TOTAL; phiIndex++){
	// Double_t phiWave = TMath::DegToRad()*phiIndex*ZOOM_BIN_SIZE_PHI;
	Double_t phiDeg = minPhiDegZoom + phiIndex*ZOOM_BIN_SIZE_PHI;
	zoomedPhiWaveLookup[phiIndex] = phiIndex*ZOOM_BIN_SIZE_PHI;

	// Double_t offAxisDelay = getOffAxisDelay((AnitaPol::AnitaPol_t)pol, ant1, ant2, phiWave);
	// offAxisDelays[pol][combo][phiIndex] = offAxisDelay;
	Double_t offAxisDelay = relativeOffAxisDelay((AnitaPol::AnitaPol_t)pol, ant2, ant1, phiDeg);
	offAxisDelays[pol][combo][phiIndex] = offAxisDelay;
	offAxisDelaysDivided[pol][combo][phiIndex] = offAxisDelay/cc->correlationDeltaT;

	part21sZoom[pol][combo][phiIndex] = zoomedCosPartLookup[pol][ant2][phiIndex] - zoomedCosPartLookup[pol][ant1][phiIndex];

      }
    }
  }
}




TH2D* InterferometricMapMaker::getMap(AnitaPol::AnitaPol_t pol, Double_t& peakValue,
			      Double_t& peakPhiDeg, Double_t& peakThetaDeg,
			      UShort_t l3TrigPattern){

  TString name = "h";
  name += pol == AnitaPol::kVertical ? "ImageV" : "ImageH";
  name += TString::Format("%u", eventNumber[pol]);

  TString title = TString::Format("Event %u ", eventNumber[pol]);
  title += (pol == AnitaPol::kVertical ? "VPOL" : "HPOL");
  title += " Map";

  Double_t phiMin = getBin0PhiDeg();
  Double_t phiMax = phiMin + DEGREES_IN_CIRCLE;
  Double_t thetaMin = MIN_THETA;
  Double_t thetaMax = MAX_THETA;

  TH2D* hImage = new TH2D(name, title,
			  NUM_BINS_PHI*NUM_PHI, phiMin, phiMax,
			  NUM_BINS_THETA, thetaMin, thetaMax);
  hImage->GetXaxis()->SetTitle("Azimuth (Degrees)");
  hImage->GetYaxis()->SetTitle("Elevation (Degrees)");

  peakValue = -2;
  peakPhiDeg = -9999;
  peakThetaDeg = -9999;

  for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA; thetaBin++){
    for(Int_t phiSector=0; phiSector<NUM_PHI; phiSector++){
      Int_t doPhiSector = RootTools::getBit(phiSector, l3TrigPattern);
      if(doPhiSector){
	for(Int_t phiBin = phiSector*NUM_BINS_PHI; phiBin < NUM_BINS_PHI*(phiSector+1); phiBin++){
	  hImage->SetBinContent(phiBin+1, thetaBin+1, coarseMap[pol][phiBin][thetaBin]);
	  if(coarseMap[pol][phiBin][thetaBin] > peakValue){
	    peakValue = coarseMap[pol][phiBin][thetaBin];
	    peakPhiDeg = hImage->GetXaxis()->GetBinLowEdge(phiBin+1);
	    peakThetaDeg = hImage->GetYaxis()->GetBinLowEdge(thetaBin+1);
	  }
	}
      }
    }
  }

  return hImage;
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







TH2D* InterferometricMapMaker::makeGlobalImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
				       Double_t& peakThetaDeg){

  return getMap(pol, imagePeak, peakPhiDeg, peakThetaDeg);
}



TH2D* InterferometricMapMaker::makeGlobalImage(AnitaPol::AnitaPol_t pol){
  Double_t imagePeak, peakPhiDeg, peakThetaDeg;
  return makeGlobalImage(pol, imagePeak, peakPhiDeg, peakThetaDeg);
}

TH2D* InterferometricMapMaker::makeTriggeredImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak,
					  Double_t& peakPhiDeg, Double_t& peakThetaDeg,
					  UShort_t l3TrigPattern){
  return getMap(pol, imagePeak, peakPhiDeg, peakThetaDeg, l3TrigPattern);
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

void InterferometricMapMaker::reconstruct(AnitaPol::AnitaPol_t pol, Double_t& imagePeak,
				  Double_t& peakPhiDeg, Double_t& peakThetaDeg){

  imagePeak = -DBL_MAX;
  peakPhiDeg = -9999;
  peakThetaDeg = -9999;
  
  Int_t peakPhiBin = -1;
  Int_t peakThetaBin = -1;

  // zero internal map
  for(Int_t phiSector = 0; phiSector < NUM_PHI; phiSector++){
    Int_t startPhiBin = phiSector*NUM_BINS_PHI;
    Int_t endPhiBin = (phiSector+1)*NUM_BINS_PHI;
    // std::cout << startPhiBin << "\t" << endPhiBin << std::endl;
    for(Int_t phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){
      for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA; thetaBin++){
	coarseMap[pol][phiBin][thetaBin] = 0;
      }
    }
  }

  std::vector<Int_t>* combosToUse = NULL;
  for(Int_t phiSector = 0; phiSector < NUM_PHI; phiSector++){
    combosToUse = &cc->combosToUseGlobal[phiSector];

    Int_t startPhiBin = phiSector*NUM_BINS_PHI;
    Int_t endPhiBin = (phiSector+1)*NUM_BINS_PHI;
    for(UInt_t comboInd=0; comboInd<combosToUse->size(); comboInd++){
      Int_t combo = combosToUse->at(comboInd);
      if(cc->kOnlyThisCombo >= 0 && combo!=cc->kOnlyThisCombo){
	continue;
      }
      for(Int_t phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){

	for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA; thetaBin++){
	  Int_t offsetLow = offsetLows[pol][combo][phiBin][thetaBin];
	  Double_t c1 = cc->crossCorrelations[pol][combo][offsetLow];
	  Double_t c2 = cc->crossCorrelations[pol][combo][offsetLow+1];
	  Double_t cInterp = interpPreFactors[pol][combo][phiBin][thetaBin]*(c2 - c1) + c1;

	  coarseMap[pol][phiBin][thetaBin] += cInterp;
	}
      }
    }
  }

  for(Int_t phiSector = 0; phiSector < NUM_PHI; phiSector++){
    combosToUse = &cc->combosToUseGlobal[phiSector];

    Double_t normFactor = cc->kOnlyThisCombo < 0 && combosToUse->size() > 0 ? combosToUse->size() : 1;
    // absorb the removed inverse FFT normalization
    normFactor*=(cc->numSamples*cc->numSamples);

    Int_t startPhiBin = phiSector*NUM_BINS_PHI;
    Int_t endPhiBin = (phiSector+1)*NUM_BINS_PHI;
    for(Int_t phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){
      for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA; thetaBin++){

	coarseMap[pol][phiBin][thetaBin]/=normFactor;
	if(coarseMap[pol][phiBin][thetaBin] > imagePeak){
	  imagePeak = coarseMap[pol][phiBin][thetaBin];
	  peakPhiBin = phiBin;
	  peakThetaBin = thetaBin;
	}
      }
    }
  }

  peakPhiDeg = phiWaveLookup[peakPhiBin]*TMath::RadToDeg();
  peakThetaDeg = thetaWaves[peakThetaBin]*TMath::RadToDeg();

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

  Double_t deltaPhiDegPhi0 = RootTools::getDeltaAngleDeg(zoomCenterPhiDeg, getBin0PhiDeg());
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
    for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA_ZOOM; thetaBin++){
      Int_t zoomThetaInd = thetaZoomBase + thetaBin;
      Double_t partBA = partBAsZoom[pol][combo][zoomThetaInd];
      // Double_t dtFactor = zoomedCosThetaWaves[zoomThetaInd]/SPEED_OF_LIGHT_NS;
      Double_t dtFactor = dtFactors[zoomThetaInd];
      for(Int_t phiBin = 0; phiBin < NUM_BINS_PHI_ZOOM; phiBin++){
	Int_t zoomPhiInd = phiZoomBase + phiBin;
	Double_t offsetLowDouble = dtFactor*(partBA - part21sZoom[pol][combo][zoomPhiInd]);
	offsetLowDouble += kUseOffAxisDelay > 0 ? offAxisDelaysDivided[pol][combo][zoomPhiInd] : 0;
	// hack for floor()
	Int_t offsetLow = (int) offsetLowDouble - (offsetLowDouble < (int) offsetLowDouble);

	Double_t deltaT = (offsetLowDouble - offsetLow);
	offsetLow += offset;
	Double_t c1 = cc->crossCorrelationsUpsampled[pol][combo][offsetLow];
	Double_t c2 = cc->crossCorrelationsUpsampled[pol][combo][offsetLow+1];
	Double_t cInterp = deltaT*(c2 - c1) + c1;

	fineMap[pol][peakIndex][thetaBin][phiBin] += cInterp;
      }
    }
  }

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

  // I can get an out of range error (I think) in the
  if(phiDeg < -500 || thetaDeg < -500 || phiDeg >= 500 || thetaDeg >= 500){
    phiDeg = 0;
    thetaDeg = 0;

    // std::cerr << "Requesting coherently summed waveform with implausible incoming angle." << std::endl;
    // return NULL;
  }

  Double_t theDeltaT = nSamp == cc->numSamples ? cc->nominalSamplingDeltaT : cc->correlationDeltaT;

  Double_t phiRad = phiDeg*TMath::DegToRad();
  Double_t thetaRad = thetaDeg*TMath::DegToRad();
  Int_t numAnts = 0;

  Int_t centerPhiSector = getPhiSectorOfAntennaClosestToPhiDeg(pol, phiDeg);

  const Int_t firstAnt = centerPhiSector;


  std::pair<Int_t, Int_t> key(nSamp, 0);
  // vArray is actually internal memory managed by FancyFFTs... don't delete this!!!
  Double_t* vArray = FancyFFTs::getRealArray(key);
  // std::complex<double>* cArray = FancyFFTs::getComplexArray(key);
  // Double_t* rmsArray = NULL;
  // if(nSamp==numSamples){
  //   FancyFFTs::doInvFFT(nSamp, ffts[pol][firstAnt], false);
  //   rmsArray = interpRMS[pol];
  // }
  // else{
  //   // Here is the fix, do something hacky manually copying the waveform fft into the longer internal array
  //   const int numFreqs = FancyFFTs::getNumFreqs(numSamples);
  //   for(int freqInd=0; freqInd < numFreqs; freqInd++){
  //     cArray[freqInd] = ffts[pol][firstAnt][freqInd];
  //   }
  //   const int moreFreqs = FancyFFTs::getNumFreqs(numSamplesUpsampled);
  //   for(int freqInd=numFreqs; freqInd < moreFreqs; freqInd++){
  //     cArray[freqInd] = 0;
  //   }
  //   FancyFFTs::doInvFFT(nSamp, cArray, false);
  //   rmsArray = interpRMS2[pol];
  //   //rmsArray = interpRMS[pol];
  // }

  // double lengthFftNorm = double(nSamp)/numSamples;

  std::vector<Double_t> tArray(nSamp, 0);
  Double_t t0 = cc->grsResampled[pol][firstAnt]->GetX()[0];
  for(Int_t samp=0; samp<nSamp; samp++){
    double time = t0 + samp*theDeltaT;
    tArray.at(samp) = time;
    // vArray[samp] *= interpRMS[pol][firstAnt]*lengthFftNorm; // Undo the normalization.
    vArray[samp] = cc->grsResampled[pol][firstAnt]->Eval(time);
  }

  // sum of rms of first few ns of each waveform
  Double_t rms = 0;
  Int_t numSampRms =0;

  // Factor of two here to drop the zero padding at the back of the waveform
  // which was used during correlations.
  TGraph* grCoherent = new TGraph(nSamp/2, &tArray[0], &vArray[0]);

  for(Int_t deltaPhiSect=-maxDeltaPhiSect; deltaPhiSect<=maxDeltaPhiSect; deltaPhiSect++){
  // for(int i=0; i < 0; i++){
    // int deltaPhiSect = 0;
    Int_t phiSector = deltaPhiSect + centerPhiSector;
    phiSector = phiSector < 0 ? phiSector + NUM_PHI : phiSector;
    phiSector = phiSector >= NUM_PHI ? phiSector - NUM_PHI : phiSector;

    for(Int_t ring=0; ring<NUM_RING; ring++){
    // for(Int_t ring=0; ring<1; ring++){
      Int_t ant= phiSector + ring*NUM_PHI;


      if(firstAnt!=ant){ // Don't do the first antenna twice

	Double_t dt = getDeltaTExpected(pol, firstAnt, ant, phiRad, thetaRad);

	Double_t t0 = cc->grsResampled[pol][ant]->GetX()[0];

	for(Int_t samp=0; samp<grCoherent->GetN(); samp++){
	  double time = t0 + samp*theDeltaT + dt;

	  // grCoherent->GetY()[samp] += grsResampled[pol][ant]->Eval(time - deltaT); //interpRMS[pol][ant];
	  // grCoherent->GetY()[samp] += grsResampled[pol][ant]->Eval(time + deltaT); //interpRMS[pol][ant];
	  grCoherent->GetY()[samp] += cc->grsResampled[pol][ant]->Eval(time + dt); //interpRMS[pol][ant];

	}
      }

      // Here we look at the RMS of the first few ns of the uninterpolated waveforms
      const Double_t timeToEvalRms = 10; // ns
      // for(Int_t samp3=0; samp3 < cc->grs[pol][ant]->GetN(); samp3++){
      // start time of the interpolated waveform... but this could be front padded with zero...
      Double_t t0 = cc->grsResampled[pol][ant]->GetX()[0];
      // this is when we've got the good stuff..., we really want to start counting from here
      Double_t t0Good = cc->grs[pol][ant]->GetX()[0];
      for(Int_t samp3=0; samp3 < nSamp; samp3++){
	Double_t t = t0 + samp3*theDeltaT;
	if(t >= t0Good && t < t0Good + timeToEvalRms){
	  // rms += cc->grs[pol][ant]->GetY()[samp3]*cc->grs[pol][ant]->GetY()[samp3];
	  double V = cc->grsResampled[pol][ant]->GetY()[samp3];
	  rms += V*V;
	  numSampRms++;

	  // std::cout << ant << "\t" << t << "\t" << samp3 << "\t" << V << "\t" << numSampRms << std::endl;
	}
      }
      numAnts++;
    }
  }

  // Normalize
  if(numAnts > 0){
    TString name = nSamp == cc->numSamples ? "grCoherent" : "grInterpCoherent";
    TString title;
    for(Int_t samp=0; samp<grCoherent->GetN(); samp++){
      grCoherent->GetY()[samp]/=numAnts;
    }

    if(pol==AnitaPol::kHorizontal){
      // name += TString::Format("H_%u", eventNumber[pol]);
      name += TString::Format("H_%u", eventNumber[pol]);
      title = "HPOL ";
    }
    else{
      // name += TString::Format("V_%u", eventNumber[pol]);
      name += TString::Format("V_%u", eventNumber[pol]);
      title = "VPOL ";
    }

    title += TString::Format("Coherently Summed Waveform for arrival direction elevation %4.2lf (Deg) and azimuth %4.2lf (Deg); Time (ns); Voltage (mV)", thetaDeg, phiDeg);

    grCoherent->SetName(name);
    grCoherent->SetTitle(title);

    rms/=numSampRms;
    rms = TMath::Sqrt(rms);

    Double_t maxY = 0;
    Double_t maxX = 0;
    Double_t minY = 0;
    Double_t minX = 0;
    RootTools::getLocalMaxToMin(grCoherent, maxY, maxX, minY, minX);
    snr = (maxY - minY)/(2*rms);
    // std::cout << std::endl << snr << "\t" << maxY - maxX << "\t" << rms << "\t" << numSampRms << std::endl;
  }
  else{
    snr = -9999;
  }
  return grCoherent;
}







template void InterferometricMapMaker::reconstructEvent<UsefulAnitaEvent>(UsefulAnitaEvent* usefulEvent, Int_t numFinePeaks=MAX_NUM_PEAKS, Int_t numCoarsePeaks=MAX_NUM_PEAKS);
template void InterferometricMapMaker::reconstructEvent<FilteredAnitaEvent>(FilteredAnitaEvent* usefulEvent, Int_t numFinePeaks=MAX_NUM_PEAKS, Int_t numCoarsePeaks=MAX_NUM_PEAKS);

template AnitaPol::AnitaPol_t InterferometricMapMaker::reconstructEventPeakPol<UsefulAnitaEvent>(UsefulAnitaEvent* usefulEvent, Int_t numFinePeaks=MAX_NUM_PEAKS, Int_t numCoarsePeaks=MAX_NUM_PEAKS);
template AnitaPol::AnitaPol_t InterferometricMapMaker::reconstructEventPeakPol<FilteredAnitaEvent>(FilteredAnitaEvent* usefulEvent, Int_t numFinePeaks=MAX_NUM_PEAKS, Int_t numCoarsePeaks=MAX_NUM_PEAKS);


template void InterferometricMapMaker::reconstructEvent<UsefulAnitaEvent>(UsefulAnitaEvent* usefulEvent, UsefulAdu5Pat& usefulPat, AnitaEventSummary* eventSummary);
template void InterferometricMapMaker::reconstructEvent<FilteredAnitaEvent>(FilteredAnitaEvent* usefulEvent, UsefulAdu5Pat& usefulPat, AnitaEventSummary* eventSummary);

