#include "CrossCorrelator.h"




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Constructor
 */
CrossCorrelator::CrossCorrelator(){
  initializeVariables();
}






//---------------------------------------------------------------------------------------------------------
/**
 * @brief Destructor
 */
CrossCorrelator::~CrossCorrelator(){
  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    deleteAllWaveforms((AnitaPol::AnitaPol_t)pol);
  }

  if(corrInterp!=NULL){
    delete corrInterp;
    corrInterp = NULL;
  }
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Workhorse function to set internal variables.
 */
void CrossCorrelator::initializeVariables(){

  kDeltaPhiSect = 2;
  multiplyTopRingByMinusOne = 0;

  corrInterp = NULL;
  coherentDeltaPhi = 0;

  numThreads = 1;
  while((Int_t)threadImagePeak.size() < numThreads){
    threadImagePeak.push_back(0);
    threadPeakPhiDeg.push_back(0);
    threadPeakThetaDeg.push_back(0);
    threadImagePeakZoom.push_back(0);
    threadPeakPhiDegZoom.push_back(0);
    threadPeakThetaDegZoom.push_back(0);
    threadPeakPhiBinZoom.push_back(0);
    threadPeakThetaBinZoom.push_back(0);
  }

  // Initialize with NULL otherwise very bad things will happen with gcc
  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      grs[pol][ant] = NULL;
      grsResampled[pol][ant] = NULL;
      interpRMS[pol][ant] = 0;
    }
  }

  maxDPhiDeg = 0;
  kOnlyThisCombo = -1;
  kUseOffAxisDelay = 1;
  numSamples = PAD_FACTOR*NUM_SAMPLES; // Factor of two for padding. Turns circular xcor into linear xcor.
  numSamplesUpsampled = numSamples*UPSAMPLE_FACTOR; // For upsampling

  nominalSamplingDeltaT = NOMINAL_SAMPLING_DELTAT;
  correlationDeltaT = nominalSamplingDeltaT/UPSAMPLE_FACTOR;

  do5PhiSectorCombinatorics();

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

  threadPol = AnitaPol::kHorizontal;
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Prints some summary information about the CrossCorrelator to stdout.
 */
void CrossCorrelator::printInfo(){
  std::cerr << __PRETTY_FUNCTION__ << std::endl;
  std::cerr << "\tupsample factor = " << UPSAMPLE_FACTOR << std::endl;
  // std::cerr << "\tBin size theta (deg) = " << Double_t(THETA_RANGE)/NUM_BINS_THETA << std::endl;
  std::cerr << "\tBin size theta (deg) = " << Double_t(MAX_THETA - MIN_THETA)/NUM_BINS_THETA << std::endl;
  std::cerr << "\tBin size phi (deg) = " << Double_t(PHI_RANGE)/NUM_BINS_PHI << std::endl;
  std::cerr << "\tdeltaTs array size = "
	    << sizeof(Double_t)*NUM_POL*numCombos*NUM_PHI*NUM_BINS_PHI*NUM_BINS_THETA
	    << " bytes" << std::endl;
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Single function to get the angle of the first bin of the interferometric histogram.
 *
 * @returns position of antenna 0 in ADU5Pat coordinates, offset by half a phi-sector.
 */
Double_t CrossCorrelator::getBin0PhiDeg(){

  Double_t phi0 = -aftForeOffset;
  if(phi0 < -DEGREES_IN_CIRCLE/2){
    phi0+=DEGREES_IN_CIRCLE;
  }
  else if(phi0 >= DEGREES_IN_CIRCLE/2){
    phi0-=DEGREES_IN_CIRCLE;
  }
  return phi0 - PHI_RANGE/2;
}






//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates an interpolated TGraph with zero mean.
 *
 * @param grIn is a TGraph containing the waveform to interpolate / pad
 * @param startTime is the time to begin the interpolation
 * @param dt is the step size to interpolate with
 * @param nSamp is the number of samples to pad to.
 *
 * @return a new TGraph containing the zero meaned, interpolated TGraph.
 */
TGraph* CrossCorrelator::interpolateWithStartTimeAndZeroMean(TGraph* grIn, Double_t startTime,
							     Double_t dt, Int_t nSamp){

  // I want to interpolate the waveform such that the mean of the post interpolation waveform is zero...

  std::vector<Double_t> newTimes = std::vector<Double_t>(nSamp, 0); // for the interp volts
  std::vector<Double_t> newVolts = std::vector<Double_t>(nSamp, 0); // for the interp times
  std::vector<Int_t> isPadding = std::vector<Int_t>(nSamp, 1); // whether or not the value is padding
  Double_t thisStartTime = grIn->GetX()[0];
  Double_t lastTime = grIn->GetX()[grIn->GetN()-1];

  // // Quantizes the start and end times so data points lie at integer multiples of nominal sampling
  // // perhaps not really necessary if all waveforms have the same start time...
  // startTime = dt*TMath::Nint(startTime/dt + 0.5);
  // lastTime = dt*TMath::Nint(lastTime/dt - 0.5);

   //ROOT interpolator object constructor takes std::vector objects
  std::vector<Double_t> tVec(grIn->GetX(), grIn->GetX() + grIn->GetN());
  std::vector<Double_t> vVec(grIn->GetY(), grIn->GetY() + grIn->GetN());

  // This is ROOT's interpolator object

  const int minGlsInterpPoints = 5;
  if(tVec.size() < minGlsInterpPoints || vVec.size() < minGlsInterpPoints){
    std::cerr << "Warning! Very short waveform! Make sure you have a data quality cut on number of samples" << std::endl;
    // not really random, but a hack to make sure GSL doesn't crash
    // this is why you should do data quality checks BEFORE reconstruction...
    const int randomPadLength = 20;
    for(int i=0; i < randomPadLength; i++){
      vVec.push_back(0);
      tVec.push_back(nominalSamplingDeltaT + (tVec.at(tVec.size()-1)));
    }
  }

  ROOT::Math::Interpolator chanInterp(tVec,vVec,ROOT::Math::Interpolation::kAKIMA);

  // Put new data into arrays
  Double_t sum = 0;
  Double_t time = startTime;
  Int_t meanCount = 0;
  for(Int_t samp = 0; samp < nSamp; samp++){
    newTimes.at(samp) = time;
    if(time >= thisStartTime && time <= lastTime){
      Double_t V = chanInterp.Eval(time);
      newVolts.at(samp) = V;
      sum += V;
      isPadding.at(samp) = 0;
      meanCount++;
    }
    else{
      newVolts.at(samp) = 0;
    }
    time += dt;
  }

  // so now we have the new arrays, time to normalize these bad boys
  Double_t mean = sum/meanCount; // the mean of the non-padded values
  for(Int_t samp=0; samp < nSamp; samp++){
    if(isPadding.at(samp)==0){
      newVolts.at(samp) -= mean;
    }
  }

  return new TGraph(nSamp, &newTimes[0], &newVolts[0]);

}



//---------------------------------------------------------------------------------------------------------
/**
 * @brief Makes evenly re-sampled, normalized waveform graphs from a UsefulAnitaEvent or FilteredAnitaEvent.
 *
 * @param usefulEvent points to the UsefulAnitaEvent of interest.
 * @param pol is the polarization of interest.
 */
// Allow functions to take UsefulAnitaEvent and/or FilteredAnitaEvent.
template <class NiceAnitaEvent> // needs eventNumber member and getGraph(int ant, AnitaPol::AnitaPol_t pol)
void CrossCorrelator::getNormalizedInterpolatedTGraphs(NiceAnitaEvent* usefulEvent,
						       AnitaPol::AnitaPol_t pol){

  // now that MagicDisplay allows filtering on the fly, remove all lastEventNumber crap
  // still want this somewhere, so put it here for now...
  eventNumber[pol] = usefulEvent->eventNumber;

  // Delete any old waveforms (at start rather than end to leave in memory to be examined if need be)
  deleteAllWaveforms(pol);


  std::vector<Double_t> earliestStart(NUM_POL, 100); // ns
  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){

    grs[pol][ant] = usefulEvent->getGraph(ant, (AnitaPol::AnitaPol_t)pol);

    if(multiplyTopRingByMinusOne > 0 && ant < NUM_PHI){ // top ring
      for(int samp=0; samp < grs[pol][ant]->GetN(); samp++){
	grs[pol][ant]->GetY()[samp] *= -1;
      }
    }

    // Find the start time of all waveforms
    if(grs[pol][ant]->GetX()[0]<earliestStart.at(pol)){
      earliestStart.at(pol) = grs[pol][ant]->GetX()[0];
    }
  }

  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
    grsResampled[pol][ant] = interpolateWithStartTimeAndZeroMean(grs[pol][ant], earliestStart.at(pol),
								 nominalSamplingDeltaT, numSamples);

    Double_t sumOfVSquared = 0;
    for(int samp=0; samp < grsResampled[pol][ant]->GetN(); samp++){
      Double_t V = grsResampled[pol][ant]->GetY()[samp];
      sumOfVSquared += V*V;
    }
    interpRMS[pol][ant] = TMath::Sqrt(sumOfVSquared/numSamples);
    interpRMS2[pol][ant] = TMath::Sqrt(sumOfVSquared/numSamplesUpsampled);
    // for(int samp=0; samp < grsResampled[pol][ant]->GetN(); samp++){
    // 	grsResampled[pol][ant]->GetY()[samp]/=interpRMS[pol][ant];
    // }
  }

  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){

    TString name, title;
    TString name2, title2;
    if(pol==AnitaPol::kHorizontal){
      name = TString::Format("gr%dH_%u", ant, usefulEvent->eventNumber);
      name2 = TString::Format("grInterp%dH_%u", ant, usefulEvent->eventNumber);
      title = TString::Format("Antenna %d HPOL eventNumber %u", ant, usefulEvent->eventNumber);
      title2 = TString::Format("Antenna %d HPOL eventNumber %u (evenly resampled)",
			       ant, usefulEvent->eventNumber);
    }
    else if(pol==AnitaPol::kVertical){
      name = TString::Format("gr%dV_%u", ant, usefulEvent->eventNumber);
      name2 = TString::Format("grInterp%dV_%u", ant, usefulEvent->eventNumber);
      title = TString::Format("Antenna %d VPOL eventNumber %u", ant, usefulEvent->eventNumber);
      title2 = TString::Format("Antenna %d VPOL eventNumber %u (evenly resampled)",
			       ant, usefulEvent->eventNumber);
    }
    title += "; Time (ns); Voltage (mV)";
    grs[pol][ant]->SetName(name);
    grsResampled[pol][ant]->SetName(name2);
    grs[pol][ant]->SetTitle(title);
    grsResampled[pol][ant]->SetTitle(title2);
  }
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Takes FFTs of the normalized, evenly resampled waveforms and puts them in memory
 * @param pol is which polarization to process
 *
 * Now also creates the zero padded FFTs used for upsampled cross-correlations.
 */
void CrossCorrelator::doFFTs(AnitaPol::AnitaPol_t pol){

  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
    FancyFFTs::doFFT(numSamples, grsResampled[pol][ant]->GetY(), ffts[pol][ant]);
  }

}















//---------------------------------------------------------------------------------------------------------
/**
 * @brief Goes through the coarseMap and finds the top N values.
 *
 * @param pol is the polarization of interest.
 * @param numPeaks is the number of peaks you wish to find, and the assumed size arrays pointed to by the following 3 variables.
 * @param peakValues is a pointer to the first element of an array of length numPeaks and will be filled with the values of the coarse map at the peaks.
 * @param phiDegs is a pointer to the first element of an array of length numPeaks and will be filled with the values of phi (degrees relative to ADU5-aft-fore) at the peaks.
 * @param thetaDegs is a pointer to the first element of an array of length numPeaks and will be filled with the values of theta (degrees relative to ADU5-aft-fore) at the peaks.
 *
 * Assumes the event has been correlated and reconstructed.
 * These peaks must be separated by the PEAK_PHI_DEG_RANGE in phi (Degrees) and PEAK_THETA_DEG_RANGE in theta (Degrees)
 */
void CrossCorrelator::findPeakValues(AnitaPol::AnitaPol_t pol, Int_t numPeaks, Double_t* peakValues,
				     Double_t* phiDegs, Double_t* thetaDegs){

  // In this function I want to find numPeak peaks and set an exclusion zone around each peak
  // so that the next peak isn't just a neighbour of a true peak.

  if(numPeaks > MAX_NUM_PEAKS){
    // You can have numPeaks less than or requal MAX_NUM_PEAKS, but not greater than.
    std::cerr << "Warning in "<< __PRETTY_FUNCTION__ << " in " << __FILE__
	      << ". numPeaks = " << numPeaks << ", CrossCorrelator compiled with MAX_NUM_PEAKS  = "
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








//---------------------------------------------------------------------------------------------------------
/**
 * @brief Reconstruct event but only do interpolated peak finding on the polarization with the highest interpolated peak finding
 *
 * @param usefulEvent is the event to process
 * @param numFinePeaks is the number of fine peaks to reconstruct (should be >= numCoarsePeaks but < MAX_NUM_PEAKS)
 * @param numCoarsePeaks is the number of coarse peaks to reconstruct (should < MAX_NUM_PEAKS)
 *
 * Wraps the key reconstruction algorithms and puts the results in internal memory.
 */
template <class NiceAnitaEvent> // needs eventNumber member and getGraph(int ant, AnitaPol::AnitaPol_t pol)
AnitaPol::AnitaPol_t CrossCorrelator::reconstructEventPeakPol(NiceAnitaEvent* usefulEvent, Int_t numFinePeaks ,Int_t numCoarsePeaks){

  for(Int_t polInd = AnitaPol::kHorizontal; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t)polInd;

    // now calls reconstruct inside correlate event
    correlateEvent(usefulEvent, pol);

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

      // std::cerr << fineMapPeakValues[pol][peakInd] << "\t" << fineMapPeakPhiDegs[pol][peakInd] << "\t"
      // 	  << fineMapPeakThetaDegs[pol][peakInd] << std::endl << std::endl;
    }
  }
  return peakPol;
}











//---------------------------------------------------------------------------------------------------------
/**
 * @brief Reconstruct event
 *
 * @param usefulEvent is the event to process
 * @param numFinePeaks is the number of fine peaks to reconstruct (should be >= numCoarsePeaks but < MAX_NUM_PEAKS)
 * @param numCoarsePeaks is the number of coarse peaks to reconstruct (should < MAX_NUM_PEAKS)
 *
 * Wraps the key reconstruction algorithms and puts the results in internal memory.
 */
template <class NiceAnitaEvent> // needs eventNumber member and getGraph(int ant, AnitaPol::AnitaPol_t pol)
void CrossCorrelator::reconstructEvent(NiceAnitaEvent* usefulEvent, UsefulAdu5Pat& usefulPat, AnitaEventSummary* eventSummary){

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

























//---------------------------------------------------------------------------------------------------------
/**
 * @brief Reconstruct event
 *
 * @param usefulEvent is the event to process
 * @param numFinePeaks is the number of fine peaks to reconstruct (should be >= numCoarsePeaks but < MAX_NUM_PEAKS)
 * @param numCoarsePeaks is the number of coarse peaks to reconstruct (should < MAX_NUM_PEAKS)
 *
 * Wraps the key reconstruction algorithms and puts the results in internal memory.
 */
template <class NiceAnitaEvent> // needs eventNumber member and getGraph(int ant, AnitaPol::AnitaPol_t pol)
void CrossCorrelator::reconstructEvent(NiceAnitaEvent* usefulEvent, Int_t numFinePeaks ,Int_t numCoarsePeaks){

  for(Int_t polInd = AnitaPol::kHorizontal; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t)polInd;

    // now calls reconstruct inside correlate event
    correlateEvent(usefulEvent, pol);

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







//---------------------------------------------------------------------------------------------------------
/**
 * @brief Gets the results from the coarse reconstruction
 *
 * @param pol is the polarization of interest
 * @param peakIndex runs from 0 to MAX_NUM_PEAKS and indexes the peak (0 is the largest).
 * @param value is the bin content of the image peak.
 * @param phiDeg is the phi coordinate in degrees.
 * @param thetaDeg is the theta coordinate in degrees.
 */

void CrossCorrelator::getCoarsePeakInfo(AnitaPol::AnitaPol_t pol, Int_t peakIndex,
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





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Gets the results from the fine reconstruction
 *
 * @param pol is the polarization of interest
 * @param peakIndex runs from 0 to MAX_NUM_PEAKS and indexes the peak (0 is the largest).
 * @param value is the bin content of the image peak.
 * @param phiDeg is the phi coordinate in degrees.
 * @param thetaDeg is the theta coordinate in degrees.
 */

void CrossCorrelator::getFinePeakInfo(AnitaPol::AnitaPol_t pol, Int_t peakIndex,
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





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Correlate the event
 *
 * First step in interferometry is to pass the pointer to the UsefulAnitaEvent in here.
 * This populates the internal arrays of normalized, interpolated waveforms and set of cross correlations.
 * These data are then used as look up tables for generating interferometic images.
 *
 * @param usefulEvent is the event to process.
 */
template <class NiceAnitaEvent> // needs eventNumber member and getGraph(int ant, AnitaPol::AnitaPol_t pol)
void CrossCorrelator::correlateEvent(NiceAnitaEvent* usefulEvent){

  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    correlateEvent(usefulEvent, (AnitaPol::AnitaPol_t)pol);
  }
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Correlate the event
 *
 * First step in interferometry is to pass the pointer to the UsefulAnitaEvent/FilteredAnitaEvent in here.
 * This populates the internal arrays of normalized, interpolated waveforms and set of cross correlations.
 * These data are then used as look up tables for generating interferometic images.
 *
 * @param usefulEvent is the event to process.
 * @param pol tells CrossCorrelator to only do this polarization.
 */
template <class NiceAnitaEvent> // needs eventNumber member and getGraph(int ant, AnitaPol::AnitaPol_t pol)
void CrossCorrelator::correlateEvent(NiceAnitaEvent* usefulEvent, AnitaPol::AnitaPol_t pol){

  // Read TGraphs from events into memory (also deletes old TGraphs)
  getNormalizedInterpolatedTGraphs(usefulEvent, pol);

  // Now cross correlate those already FFT'd waveforms

  // Generate set of ffts for cross correlation (each waveform only needs to be done once)
  doFFTs(pol);


  doAllCrossCorrelationsThreaded(pol);

  // reconstruct
  reconstruct(pol, coarseMapPeakValues[pol][0], coarseMapPeakPhiDegs[pol][0], coarseMapPeakThetaDegs[pol][0]);

  // fprintf(stderr, "%lf\t%lf\t%lf\n", coarseMapPeakValues[pol][0], coarseMapPeakPhiDegs[pol][0], coarseMapPeakThetaDegs[pol][0]);

  // Safety check to make sure we don't do any hard work twice.
}



//---------------------------------------------------------------------------------------------------------
/**
 * @brief Wrapper function which launches the threads for generating all the cross correlations
 *
 * @param pol tells CrossCorrelator to only do this polarization.
 */
void CrossCorrelator::doAllCrossCorrelationsThreaded(AnitaPol::AnitaPol_t pol){

  // Set variable for use in threads
  threadPol = pol;

  if(numThreads > 1){
    std::vector<TThread*> corrThreads(numThreads, NULL);
    for(long threadInd=0; threadInd<numThreads; threadInd++){
      TString name = TString::Format("threadCorr%ld", threadInd);
      CrossCorrelator::threadArgs threadArgVals;
      threadArgVals.threadInd = threadInd;
      threadArgVals.ptr = this;

      corrThreads.at(threadInd) = new TThread(name.Data(),
					      CrossCorrelator::doSomeCrossCorrelationsThreaded,
					      (void*)&threadArgVals);
    }

    for(long threadInd=0; threadInd<numThreads; threadInd++){
      corrThreads.at(threadInd)->Run();
    }

    for(long threadInd=0; threadInd<numThreads; threadInd++){
      corrThreads.at(threadInd)->Join();
    }

    for(long threadInd=0; threadInd<numThreads; threadInd++){
      delete corrThreads.at(threadInd);
    }
    corrThreads.clear();
  }
  else{
    CrossCorrelator::threadArgs threadArgVals;
    threadArgVals.threadInd = 0;
    threadArgVals.ptr = this;
    CrossCorrelator::doSomeCrossCorrelationsThreaded((void*) &threadArgVals);
  }
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Static member function which generates the finely binned set of cross correlations from the FFTs held in memory.
 *
 * @param voidPtrArgs contains a pointer to a CrossCorrelator::threadArgs struct
 */
void* CrossCorrelator::doSomeUpsampledCrossCorrelationsThreaded(void* voidPtrArgs){

  // Disgusting hacks to get ROOT threading to compile inside a class.
  CrossCorrelator::threadArgs* args = (CrossCorrelator::threadArgs*) voidPtrArgs;
  Long_t threadInd = args->threadInd;
  CrossCorrelator* ptr = args->ptr;
  AnitaPol::AnitaPol_t pol = ptr->threadPol;
  Int_t phiSector = ptr->threadPhiSector;
  Int_t numThreads = ptr->numThreads;
  const std::vector<Int_t>* combosToUse = &ptr->combosToUseGlobal[phiSector];
  Int_t numCombosAllThreads = combosToUse->size();
  Int_t numCorrPerThread = numCombosAllThreads/numThreads;
  Int_t numRemainder = numCombosAllThreads%numThreads;

  Int_t startComboInd = threadInd*numCorrPerThread;

  // Double_t stash[NUM_SAMPLES*UPSAMPLE_FACTOR*PAD_FACTOR];
  Double_t* ccInternalArray = FancyFFTs::getRealArray(std::pair<int, int>(ptr->numSamplesUpsampled, threadInd));

  for(int comboInd=startComboInd; comboInd<startComboInd+numCorrPerThread; comboInd++){
    Int_t combo = combosToUse->at(comboInd);

    Int_t ant1 = ptr->comboToAnt1s.at(combo);
    Int_t ant2 = ptr->comboToAnt2s.at(combo);

    FancyFFTs::crossCorrelatePadded(ptr->numSamples,
				    UPSAMPLE_FACTOR,
				    ptr->ffts[pol][ant2],
				    ptr->ffts[pol][ant1],
				    ptr->crossCorrelationsUpsampled[pol][combo],
				    false,
				    threadInd, false);

    const int offset = ptr->numSamplesUpsampled/2;

    // copies first half of original array (times >= 0) into second half of internal storage
    for(Int_t samp=0; samp < ptr->numSamplesUpsampled/2; samp++){
      ptr->crossCorrelationsUpsampled[pol][combo][samp+offset] = ccInternalArray[samp];
    }
    // copies second half of original array (times < 0) into first half of internal storage
    for(Int_t samp=ptr->numSamplesUpsampled/2; samp < ptr->numSamplesUpsampled; samp++){
      ptr->crossCorrelationsUpsampled[pol][combo][samp-offset] = ccInternalArray[samp];
    }
  }

  if(threadInd < numRemainder){
    Int_t numDoneInAllThreads = numThreads*numCorrPerThread;
    Int_t comboInd = numDoneInAllThreads + threadInd;
    Int_t combo = combosToUse->at(comboInd);
    Int_t ant1 = ptr->comboToAnt1s.at(combo);
    Int_t ant2 = ptr->comboToAnt2s.at(combo);

    FancyFFTs::crossCorrelatePadded(ptr->numSamples,
				    UPSAMPLE_FACTOR,
				    ptr->ffts[pol][ant2],
				    ptr->ffts[pol][ant1],
				    ptr->crossCorrelationsUpsampled[pol][combo],
				    false,
				    threadInd, false);
    const int offset = ptr->numSamplesUpsampled/2;

    // copies first half of original array (times >= 0) into second half of internal storage
    for(Int_t samp=0; samp < ptr->numSamplesUpsampled/2; samp++){
      ptr->crossCorrelationsUpsampled[pol][combo][samp+offset] = ccInternalArray[samp];
    }
    // copies second half of original array (times < 0) into first half of internal storage
    for(Int_t samp=ptr->numSamplesUpsampled/2; samp < ptr->numSamplesUpsampled; samp++){
      ptr->crossCorrelationsUpsampled[pol][combo][samp-offset] = ccInternalArray[samp];
    }
  }
  return 0;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Static member function which generates the coarsely binned set of cross correlations from the FFTs held in memory.
 *
 * @param voidPtrArgs contains a pointer to a CrossCorrelator::threadArgs struct
 */
void* CrossCorrelator::doSomeCrossCorrelationsThreaded(void* voidPtrArgs){

  // Disgusting hacks to get ROOT threading to compile inside a class.
  CrossCorrelator::threadArgs* args = (CrossCorrelator::threadArgs*) voidPtrArgs;
  Long_t threadInd = args->threadInd;
  CrossCorrelator* ptr = args->ptr;
  // TThread::Lock();
  // std::cerr << __PRETTY_FUNCTION__ << "\t" << threadInd << std::endl;
  // TThread::UnLock();
  AnitaPol::AnitaPol_t pol = ptr->threadPol;
  Int_t numThreads = ptr->numThreads;
  Int_t numCorrPerThread = NUM_COMBOS/numThreads;
  Int_t startCombo = threadInd*numCorrPerThread;

  // Double_t stash[NUM_SAMPLES*PAD_FACTOR];

  Double_t* ccInternalArray = FancyFFTs::getRealArray(std::pair<int, int>(ptr->numSamples, threadInd));

  for(int combo=startCombo; combo<startCombo+numCorrPerThread; combo++){
    Int_t ant1 = ptr->comboToAnt1s.at(combo);
    Int_t ant2 = ptr->comboToAnt2s.at(combo);
    FancyFFTs::crossCorrelatePadded(ptr->numSamples,
				    1,
				    ptr->ffts[pol][ant2],
				    ptr->ffts[pol][ant1],
				    ptr->crossCorrelations[pol][combo],
				    false,
				    threadInd,
				    false);
    const int offset = ptr->numSamples/2;

    // copies first half of original array (times >= 0) into second half of internal storage
    for(Int_t samp=0; samp < ptr->numSamples/2; samp++){
      ptr->crossCorrelations[pol][combo][samp+offset] = ccInternalArray[samp];
      // ptr->crossCorrelations[pol][combo][samp+offset] = stash[samp];
    }
    // copies second half of original array (times < 0) into first half of internal storage
    for(Int_t samp=ptr->numSamples/2; samp < ptr->numSamples; samp++){
      ptr->crossCorrelations[pol][combo][samp-offset] = ccInternalArray[samp];
      // ptr->crossCorrelations[pol][combo][samp-offset] = stash[samp];
    }
  }
  return 0;
}







//---------------------------------------------------------------------------------------------------------
/**
 * @brief Function which launches threads
 *
 * @param pol tells CrossCorrelator to only do this polarization.
 * @param phiSector is used to figure out which finely binned cross correlations are required.
 */
void CrossCorrelator::doUpsampledCrossCorrelationsThreaded(AnitaPol::AnitaPol_t pol, Int_t phiSector){

  threadPhiSector = phiSector;
  threadPol = pol;

  if(numThreads > 1){
    std::vector<TThread*> upsampledCorrThreads(numThreads, NULL);
    for(long threadInd=0; threadInd<numThreads; threadInd++){
      TString name = TString::Format("threadUpsampledCorr%ld", threadInd);
      CrossCorrelator::threadArgs threadArgVals;
      threadArgVals.threadInd = threadInd;
      threadArgVals.ptr = this;
      upsampledCorrThreads.at(threadInd) = new TThread(name.Data(),
						       CrossCorrelator::doSomeUpsampledCrossCorrelationsThreaded,
						       (void*)&threadArgVals);
    }

    for(long threadInd=0; threadInd<numThreads; threadInd++){
      upsampledCorrThreads.at(threadInd)->Run();
    }

    for(long threadInd=0; threadInd<numThreads; threadInd++){
      upsampledCorrThreads.at(threadInd)->Join();
    }

    for(long threadInd=0; threadInd<numThreads; threadInd++){
      delete upsampledCorrThreads.at(threadInd);
    }
    upsampledCorrThreads.clear();
  }
  else{
    CrossCorrelator::threadArgs threadArgVals;
    threadArgVals.threadInd = 0;
    threadArgVals.ptr = this;
    CrossCorrelator::doSomeUpsampledCrossCorrelationsThreaded((void*)&threadArgVals);
  }

}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Get the maximum correlation time and value from a polarization and antenna combo index
 *
 * @param pol is the polarization of the antenna combo
 * @param combo is the index of the combination (see comboIndices[NUM_SEAVEYS][NUM_SEAVEYS])
 * @param time is reference variable, updated with the maximum correlation time.
 * @param value is a reference variable, updated with the correlation coefficient at time.
 *
 * You probably want to call getMaxCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t& time, Double_t& value) instead.
 */
void CrossCorrelator::getMaxCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t combo,
						 Double_t& time, Double_t& value){

  Int_t maxInd = RootTools::getIndexOfMaximum(numSamples, crossCorrelations[pol][combo]);
  Int_t offset = maxInd - numSamples/2;
  value = crossCorrelations[pol][combo][maxInd];
  time = offset*nominalSamplingDeltaT;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Get the maximum correlation time and value from a polarization and pair of antennas
 *
 * @param pol is the polarization of the antenna combo
 * @param ant1 is the first antenna
 * @param ant2 is the second antenna
 * @param time is reference variable, updated with the maximum correlation time.
 * @param value is a reference variable, updated with the correlation coefficient at time.
 */
void CrossCorrelator::getMaxCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
						 Double_t& time, Double_t& value){
  Int_t combo = comboIndices[ant1][ant2];
  getMaxCorrelationTimeValue(pol, combo, time, value);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Get the maximum upsampled correlation time and value from a polarization and antenna combo index
 *
 * @param pol is the polarization of the antenna combo
 * @param combo is the index of the combination (see comboIndices[NUM_SEAVEYS][NUM_SEAVEYS])
 * @param time is reference variable, updated with the maximum correlation time.
 * @param value is a reference variable, updated with the correlation coefficient at time.
 */
void CrossCorrelator::getMaxUpsampledCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t combo,
							  Double_t& time, Double_t& value){

  Int_t maxInd = RootTools::getIndexOfMaximum(numSamplesUpsampled, crossCorrelationsUpsampled[pol][combo]);
  Int_t offset = maxInd - numSamplesUpsampled/2;
  value = crossCorrelationsUpsampled[pol][combo][maxInd];
  time = offset*correlationDeltaT;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Get the maximum upsampled correlation time and value from a polarization and pair of antennas
 *
 * @param pol is the polarization of the antenna combo
 * @param ant1 is the first antenna
 * @param ant2 is the second antenna
 * @param time is reference variable, updated with the maximum correlation time.
 * @param value is a reference variable, updated with the correlation coefficient at time.
 */
void CrossCorrelator::getMaxUpsampledCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
							  Double_t& time, Double_t& value){
  Int_t combo = comboIndices[ant1][ant2];
  getMaxUpsampledCorrelationTimeValue(pol, combo, time, value);
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Get the off-axis delay for an off boresight angle.
 *
 * @param deltaPhiDeg is the angle (Degrees) of the plane wave relative to an antenna boresight.
 * @return the delay relative to boresight in nano-seconds.
 *
 */
Double_t CrossCorrelator::singleAntennaOffAxisDelay(Double_t deltaPhiDeg) {


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





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Get the relative off-axis delay between an antenna pair.
 *
 * @param pol is the polarization of the plane wave.
 * @param ant1 is the first antenna.
 * @param ant2 is the second antenna.
 * @param phiDeg is the angle of the planewave in Degrees relative to ADU5-aft-fore.
 * @return the difference between the off-axis delays of the two antennas.
 */
Double_t CrossCorrelator::relativeOffAxisDelay(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
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





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Function to calculate the time taken by a plane wave to cross the payload
 *
 * @param pol is the polarization of the antennas
 * @param ant1 is the first antenna
 * @param ant2 is the second antenna
 * @param phiWave is the phi direction of the incoming wave in radians (relative to the ADU5)
 * @param thetaWave is the theta direction of the incoming wave in radians (relative to the ADU5)
 * @returns time at ant2 - time at ant1 in ns.
 */
Double_t CrossCorrelator::getDeltaTExpected(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t phiWave, Double_t thetaWave){

  // Double_t tanThetaW = tan(thetaWave);
  Double_t tanThetaW = tan(-1*thetaWave);
  Double_t part1 = zArray[pol].at(ant1)*tanThetaW - rArray[pol].at(ant1)*cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol].at(ant1));
  Double_t part2 = zArray[pol].at(ant2)*tanThetaW - rArray[pol].at(ant2)*cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol].at(ant2));
  Double_t tdiff = 1e9*((cos(thetaWave) * (part2 - part1))/SPEED_OF_LIGHT); // Returns time in ns

  return tdiff;
}






//---------------------------------------------------------------------------------------------------------
/**
 * @brief Inserts the photogrammetry geometry into the CrossCorrelator internals *AND* the AnitaEventCalibrator internals
 */
void CrossCorrelator::insertPhotogrammetryGeometry(){
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









//---------------------------------------------------------------------------------------------------------
/**
 * @brief Function to index all possible antenna pairs for use in reconstuction
 *
 * In order to avoid using separate regions of memory many of the interals rely on knowing NUM_COMBOS.
 * This function counts up numCombos and expects it to equal NUM_COMBOS.
 * This class will probably break if they are not equal (especially if numCombos > NUM_COMBOS).
 * Prints a warning if they are not equal.
 */
void CrossCorrelator::do5PhiSectorCombinatorics(){
  // For checking later...
  for(Int_t ant1=0; ant1 < NUM_SEAVEYS; ant1++){
    for(Int_t ant2=0; ant2 < NUM_SEAVEYS; ant2++){
      comboIndices[ant1][ant2] = -1;
    }
  }

  numCombos=0;
  for(Int_t ant1=0; ant1<NUM_SEAVEYS; ant1++){
    Double_t phiSect1 = ant1%NUM_PHI;
    for(Int_t ant2=ant1+1; ant2<NUM_SEAVEYS; ant2++){
      Double_t phiSect2 = ant2%NUM_PHI;

      if(TMath::Abs(phiSect1 - phiSect2) <= DELTA_PHI_SECT || TMath::Abs(phiSect1 - phiSect2) >= (NUM_PHI-DELTA_PHI_SECT)){
	comboIndices[ant1][ant2] = numCombos;
	comboIndices[ant2][ant1] = numCombos;
	comboToAnt1s.push_back(ant1);
	comboToAnt2s.push_back(ant2);
	numCombos++;
      }
    }
  }
  if(numCombos != NUM_COMBOS){
    std::cerr << "Warning! in = " << __PRETTY_FUNCTION__ << std::endl;
    std::cerr << "\tnumCombos = " << numCombos << "." << std::endl;
    std::cerr << "\tNUM_COMBOS = " << NUM_COMBOS << "." << std::endl;
    std::cerr << "\tExpecting NUM_COMBOS == numCombos." << std::endl;
    std::cerr << "\tSeriously bad things are probably about to happen." << std::endl;
    std::cerr << "\tCheck the combinatorics!" << std::endl;
  }

  fillCombosToUse();
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Fills in an array of cached deltaTs between antenna pairs as a function of arrival direction
 */
void CrossCorrelator::fillDeltaTLookup(){

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

  for(Int_t polInd=0; polInd<NUM_POL; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(Int_t combo=0; combo<numCombos; combo++){
      Int_t ant1 = comboToAnt1s.at(combo);
      Int_t ant2 = comboToAnt2s.at(combo);

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
	    Int_t offsetLow = floor(deltaT/nominalSamplingDeltaT);
	    offsetLows[pol][combo][phiBin][thetaBin] = offsetLow;
	    // offsetLows[pol][phiBin][combo][thetaBin] = offsetLow;
	    Double_t dt1 = offsetLow*nominalSamplingDeltaT;
	    interpPreFactors[pol][combo][phiBin][thetaBin] = (deltaT - dt1)/nominalSamplingDeltaT;
	    // interpPreFactors[pol][phiBin][combo][thetaBin] = (deltaT - dt1)/nominalSamplingDeltaT;

	    // Here we account for the fact that we are now time ordering the correlations
	    offsetLows[pol][combo][phiBin][thetaBin]+=numSamples/2;
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
    dtFactors[thetaIndex] = zoomedCosThetaWaves[thetaIndex]/(SPEED_OF_LIGHT_NS*correlationDeltaT);
  }

  for(Int_t pol=0; pol<AnitaPol::kNotAPol; pol++){
    for(Int_t combo=0; combo < NUM_COMBOS; combo++){
      Int_t ant1 = comboToAnt1s.at(combo);
      Int_t ant2 = comboToAnt2s.at(combo);
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
      Int_t ant1 = comboToAnt1s.at(combo);
      Int_t ant2 = comboToAnt2s.at(combo);

      for(Int_t phiIndex=0; phiIndex < NUM_BINS_PHI_ZOOM_TOTAL; phiIndex++){
	// Double_t phiWave = TMath::DegToRad()*phiIndex*ZOOM_BIN_SIZE_PHI;
	Double_t phiDeg = minPhiDegZoom + phiIndex*ZOOM_BIN_SIZE_PHI;
	zoomedPhiWaveLookup[phiIndex] = phiIndex*ZOOM_BIN_SIZE_PHI;

	// Double_t offAxisDelay = getOffAxisDelay((AnitaPol::AnitaPol_t)pol, ant1, ant2, phiWave);
	// offAxisDelays[pol][combo][phiIndex] = offAxisDelay;
	Double_t offAxisDelay = relativeOffAxisDelay((AnitaPol::AnitaPol_t)pol, ant2, ant1, phiDeg);
	offAxisDelays[pol][combo][phiIndex] = offAxisDelay;
	offAxisDelaysDivided[pol][combo][phiIndex] = offAxisDelay/correlationDeltaT;

	part21sZoom[pol][combo][phiIndex] = zoomedCosPartLookup[pol][ant2][phiIndex] - zoomedCosPartLookup[pol][ant1][phiIndex];

      }
    }
  }
}

//---------------------------------------------------------------------------------------------------------
/**
 * @brief This function encodes whether a pair of antennas should be used in a particular phi sector.
 *
 * @param ant1 is the first antenna
 * @param ant2 is the second antenna
 * @param phiSector is the phi-sector being considered
 * @param deltaPhiSect is the range of phi-sectors to allow in the reconstruction. Note that deltaPhiSect < 0 implies that one antenna must be in phiSector, but allows a range out to abs(deltaPhiSect).
 *
 * This function gets used when comparing different reconstruction strategies by member variable kDeltaPhiSect.
 */
Bool_t CrossCorrelator::useCombo(Int_t ant1, Int_t ant2, Int_t phiSector, Int_t deltaPhiSect){

  // I want to be able to choose whether or not require one of the antennas to be the phi-sector
  // of interest or just to have both in range of deltaPhiSect.
  // I'm going to use the sign of deltaPhiSect to do this.
  // If deltaPhiSect < 0 this implies one of the antenna pairs must be in an l3Triggered phi-sector

  Int_t absDeltaPhiSect = TMath::Abs(deltaPhiSect);

  // Require that the difference in phi-sectors be <= absDeltaPhiSect
  Int_t phiSectorOfAnt1 = ant1%NUM_PHI;
  Bool_t ant1InRange = TMath::Abs(phiSector - (phiSectorOfAnt1))<=absDeltaPhiSect;

  // Takes account of wrapping around payload (e.g. antennas in phi-sectors 1 and 16 are neighbours)
  ant1InRange = ant1InRange || TMath::Abs(phiSector - (phiSectorOfAnt1))>=(NUM_PHI-absDeltaPhiSect);

  Int_t phiSectorOfAnt2 = ant2%NUM_PHI;
  Bool_t ant2InRange = TMath::Abs(phiSector - (phiSectorOfAnt2))<=absDeltaPhiSect;
  ant2InRange = ant2InRange || TMath::Abs(phiSector - (phiSectorOfAnt2))>=(NUM_PHI-absDeltaPhiSect);

  // See rather rant above.
  Bool_t extraCondition = true;
  if(deltaPhiSect < 0){
    extraCondition = (phiSectorOfAnt1==phiSector || phiSectorOfAnt2==phiSector);
  }

  return (ant1InRange && ant2InRange && extraCondition);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates vectors of antenna combo indices and puts them in the combosToUseGlobal map
 */
void CrossCorrelator::fillCombosToUse(){

  for(Int_t phiSector = 0; phiSector<NUM_PHI; phiSector++){
    if(combosToUseGlobal[phiSector].size() == 0){
      for(Int_t combo=0; combo<numCombos; combo++){
	Int_t ant1 = comboToAnt1s.at(combo);
	Int_t ant2 = comboToAnt2s.at(combo);
	if(useCombo(ant1, ant2, phiSector, kDeltaPhiSect)){
	  combosToUseGlobal[phiSector].push_back(combo);
	}
      }
    }
  }
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Linearly interpolates between upsampled correlation points.
 *
 * @param pol is polarization of the antenna pair.
 * @param combo is the index of the antenna pair.
 * @param deltaT is the time to interpolate the correlation values at.
 * @returns the interpolated upsampled correlation value.
 */
Double_t CrossCorrelator::getInterpolatedUpsampledCorrelationValue(AnitaPol::AnitaPol_t pol,
								   Int_t combo, Double_t deltaT){

  Int_t offsetLow = floor(deltaT/correlationDeltaT);

  Double_t dt1 = offsetLow*correlationDeltaT;
  // Double_t dt2 = offsetHigh*ptr->correlationDeltaT;

  const Int_t offset = numSamplesUpsampled/2;
  offsetLow += offset;
  Int_t offsetHigh = offsetLow+1;

  // offsetLow = offsetLow < 0 ? offsetLow + numSamplesUpsampled : offsetLow;
  // offsetHigh = offsetHigh < 0 ? offsetHigh + numSamplesUpsampled : offsetHigh;

  Double_t c1 = crossCorrelationsUpsampled[pol][combo][offsetLow];
  Double_t c2 = crossCorrelationsUpsampled[pol][combo][offsetHigh];

  Double_t cInterp = (deltaT - dt1)*(c2 - c1)/(correlationDeltaT) + c1;

  return cInterp;

}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Gets an actual histogram of the zoomed in map.
 *
 * @param pol is the polarization
 * @param peakValue is the bin content of the image peak.
 * @param peakPhiDeg is the phi coordinate in degrees.
 * @param peakThetaDeg is the theta coordinate in degrees.
 * @param l3TrigPattern is a bit mask of the phi-sectors to use in reconstruct, default value is ALL_PHI_TRIGS (=0xffff).
 * @returns the TH2D histogram
 */
TH2D* CrossCorrelator::getMap(AnitaPol::AnitaPol_t pol, Double_t& peakValue,
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
  // Double_t thetaMin = -THETA_RANGE/2 + double(THETA_RANGE)/NUM_BINS_THETA;
  // Double_t thetaMax = THETA_RANGE/2 + double(THETA_RANGE)/NUM_BINS_THETA;
  // Double_t thetaMin = -THETA_RANGE/2;;
  // Double_t thetaMax = THETA_RANGE/2;
  Double_t thetaMin = MIN_THETA;
  Double_t thetaMax = MAX_THETA;

  TH2D* hImage = new TH2D(name, title,
			  NUM_BINS_PHI*NUM_PHI, phiMin, phiMax,
			  NUM_BINS_THETA, thetaMin, thetaMax);
  hImage->GetXaxis()->SetTitle("Azimuth (Degrees)");
  hImage->GetYaxis()->SetTitle("Elevation (Degrees)");

  // peakValue = -DBL_MAX;
  // peakPhiDeg = -DBL_MAX;
  // peakThetaDeg = -DBL_MAX;
  peakValue = -2;
  peakPhiDeg = -9999;
  peakThetaDeg = -9999;

  for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA; thetaBin++){
    // Int_t invertedThetaBin = NUM_BINS_THETA - thetaBin - 1;
    for(Int_t phiSector=0; phiSector<NUM_PHI; phiSector++){
      Int_t doPhiSector = RootTools::getBit(phiSector, l3TrigPattern);
      if(doPhiSector){
	for(Int_t phiBin = phiSector*NUM_BINS_PHI; phiBin < NUM_BINS_PHI*(phiSector+1); phiBin++){
	  hImage->SetBinContent(phiBin+1, thetaBin+1, coarseMap[pol][phiBin][thetaBin]);
	  // hImage->SetBinContent(phiBin+1, thetaBin+1, coarseMap[pol][phiBin][invertedThetaBin]);
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





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Gets an actual histogram of the zoomed in map.
 *
 * @param pol is the polarization
 * @param peakInd picks out which maxima to get, 0 is the largest, 1 is the second largest etc.
 * @returns the TH2D histogram
 */
TH2D* CrossCorrelator::getZoomMap(AnitaPol::AnitaPol_t pol, Int_t peakInd){

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
			  // NUM_BINS_THETA_ZOOM,
			  // -1*(zoomThetaMin[pol] + THETA_RANGE_ZOOM) + ZOOM_BIN_SIZE_THETA,
			  // -1*zoomThetaMin[pol] + ZOOM_BIN_SIZE_THETA);

  hImage->GetXaxis()->SetTitle("Azimuth (Degrees)");
  hImage->GetYaxis()->SetTitle("Elevation (Degrees)");

  for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA_ZOOM; thetaBin++){
    // Int_t invertedThetaBin = NUM_BINS_THETA_ZOOM - thetaBin - 1;
    for(Int_t phiBin = 0; phiBin < NUM_BINS_PHI_ZOOM; phiBin++){
      hImage->SetBinContent(phiBin+1, thetaBin+1, fineMap[pol][peakInd][thetaBin][phiBin]);
      // hImage->SetBinContent(phiBin+1, thetaBin+1, fineMap[pol][invertedThetaBin][phiBin]);
    }
  }

  return hImage;
}







//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates an interferometric map using plane wave deltaTs and antenna pairs from all phi-sectors.
 *
 * @param pol is the polarization.
 * @param imagePeak is the maximum value.
 * @param peakPhiDeg is azimuth of the maximum value (Degrees) relative to the ADU5 aft-fore.
 * @param peakThetaDeg is the elevation of the maximum value (Degrees).
 * @returns the interferometric map.
 */
TH2D* CrossCorrelator::makeGlobalImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
				       Double_t& peakThetaDeg){

  return getMap(pol, imagePeak, peakPhiDeg, peakThetaDeg);
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates an interferometric map using plane wave deltaTs and antenna pairs from all phi-sectors.
 *
 * @param pol is the polarization.
 * @returns the interferometric map.
 */
TH2D* CrossCorrelator::makeGlobalImage(AnitaPol::AnitaPol_t pol){
  Double_t imagePeak, peakPhiDeg, peakThetaDeg;
  return makeGlobalImage(pol, imagePeak, peakPhiDeg, peakThetaDeg);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates an interferometric map using plane wave deltaTs and l3Triggered pairs from all phi-sectors.
 *
 * @param pol is the polarization.
 * @param imagePeak is the maximum value.
 * @param peakPhiDeg is azimuth of the maximum value (Degrees) relative to the ADU5 aft-fore.
 * @param peakThetaDeg is the elevation of the maximum value (Degrees).
 * @param l3TrigPattern is the L3 trigger pattern, determines what antenna pairs to use in the recontruction.
 * @returns the interferometric map.
 */
TH2D* CrossCorrelator::makeTriggeredImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak,
					  Double_t& peakPhiDeg, Double_t& peakThetaDeg,
					  UShort_t l3TrigPattern){
  return getMap(pol, imagePeak, peakPhiDeg, peakThetaDeg, l3TrigPattern);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates an interferometric map using plane wave deltaTs centered around a particular phi/theta.
 *
 * @param pol is the polarization.
 * @param imagePeak is the maximum value.
 * @param peakPhiDeg is azimuth of the maximum value (Degrees) relative to the ADU5 aft-fore.
 * @param peakThetaDeg is the elevation of the maximum value (Degrees).
 * @param zoomCenterPhiDeg is azimuth to center the zoomed image on (Degrees) relative to the ADU5 aft-fore.
 * @param zoomCenterThetaDeg is the elevation to center the zoomed image on (Degrees).
 * @returns the interferometric map.
 */
TH2D* CrossCorrelator::makeZoomedImage(AnitaPol::AnitaPol_t pol,
				       Double_t& imagePeak, Double_t& peakPhiDeg, Double_t& peakThetaDeg,
				       Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg){

  reconstructZoom(pol, imagePeak, peakPhiDeg, peakThetaDeg,
		  zoomCenterPhiDeg, zoomCenterThetaDeg);
  return getZoomMap(pol);

}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates an interferometric map using plane wave deltaTs centered around a particular phi/theta.
 *
 * @param pol is the polarization.
 * @param imagePeak is the maximum value.
 * @param peakPhiDeg is azimuth of the maximum value (Degrees) relative to the ADU5 aft-fore.
 * @param peakThetaDeg is the elevation of the maximum value (Degrees).
 * @param l3TrigPattern is the L3 trig pattern, determines what antenna pairs to use in the recontruction.
 * @param zoomCenterPhiDeg is azimuth to center the zoomed image on (Degrees) relative to the ADU5 aft-fore.
 * @param zoomCenterThetaDeg is the elevation to center the zoomed image on (Degrees).
 * @returns the interferometric map.
 */
TH2D* CrossCorrelator::makeZoomedImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
				       Double_t& peakThetaDeg, UShort_t l3TrigPattern,
				       Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg){
  reconstructZoom(pol, imagePeak, peakPhiDeg, peakThetaDeg,
		  zoomCenterPhiDeg, zoomCenterThetaDeg);
  return getZoomMap(pol);

}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates an interferometric map using plane wave deltaTs centered around a particular phi/theta.
 *
 * @param pol is the polarization.
 * @param l3TrigPattern is the L3 trig pattern, determines what antenna pairs to use in the recontruction.
 * @param zoomCenterPhiDeg is azimuth to center the zoomed image on (Degrees) relative to the ADU5 aft-fore.
 * @param zoomCenterThetaDeg is the elevation to center the zoomed image on (Degrees).
 * @returns the interferometric map.
 */
TH2D* CrossCorrelator::makeZoomedImage(AnitaPol::AnitaPol_t pol, UShort_t l3TrigPattern,
				       Double_t zoomCenterPhiDeg,Double_t zoomCenterThetaDeg){

  Double_t imagePeak, peakPhiDeg, peakThetaDeg;
  // return makeZoomedImage(pol, imagePeak, peakPhiDeg, peakThetaDeg,
  // 			 l3TrigPattern, zoomCenterPhiDeg, zoomCenterThetaDeg);
  reconstructZoom(pol, imagePeak, peakPhiDeg, peakThetaDeg,
		  zoomCenterPhiDeg, zoomCenterThetaDeg);
  return getZoomMap(pol);
}










//---------------------------------------------------------------------------------------------------------
/**
 * @brief Wrapper function which launches the threaded functions, which fill the interferometric maps.
 *
 * @param pol is the polarization.
 * @param imagePeak is the maximum value.
 * @param peakPhiDeg is azimuth of the maximum value (Degrees) relative to the ADU5 aft-fore.
 * @param peakThetaDeg is the elevation of the maximum value (Degrees).
 *
 * This function is also responsible for merging the peak finding results of each of the threads.
 */
void CrossCorrelator::reconstruct(AnitaPol::AnitaPol_t pol, Double_t& imagePeak,
				  Double_t& peakPhiDeg, Double_t& peakThetaDeg){

  threadPol = pol;

  // LAUNCH THREADS HERE
  if(numThreads > 1){
    std::vector<TThread*> mapThreads(numThreads, NULL);
    for(long threadInd=0; threadInd<numThreads; threadInd++){
      TString name = TString::Format("threadMap%ld", threadInd);
      CrossCorrelator::threadArgs threadArgVals;
      threadArgVals.threadInd = threadInd;
      threadArgVals.ptr = this;
      mapThreads.at(threadInd) = new TThread(name.Data(),
					     CrossCorrelator::makeSomeOfImageThreaded,
					     // (void*)&threadArgsVec.at(threadInd))
					     (void*)&threadArgVals);
    }

    for(long threadInd=0; threadInd<numThreads; threadInd++){
      mapThreads.at(threadInd)->Run();
    }

    for(long threadInd=0; threadInd<numThreads; threadInd++){
      mapThreads.at(threadInd)->Join();
    }

    for(long threadInd=0; threadInd<numThreads; threadInd++){
      delete mapThreads.at(threadInd);
    }
    mapThreads.clear();
  }
  else{// call function directly
    CrossCorrelator::threadArgs threadArgVals;
    threadArgVals.threadInd = 0;
    threadArgVals.ptr = this;
    CrossCorrelator::makeSomeOfImageThreaded((void*)&threadArgVals);
  }

  // Combine peak search results from each thread.
  imagePeak = -DBL_MAX;
  for(Long_t threadInd=0; threadInd<numThreads; threadInd++){
    if(threadImagePeak[threadInd] > imagePeak){
      imagePeak = threadImagePeak[threadInd];
      peakPhiDeg = threadPeakPhiDeg[threadInd];
      peakThetaDeg = threadPeakThetaDeg[threadInd];
    }
  }

}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Static member function which fills the interferometric maps.
 *
 * @param voidPtrArgs contains a pointer to a CrossCorrelator::threadArgs struct
 *
 * This function contains the meat and bones of this class.
 * I've really tried to optimize this for speed, which means it's not very readable, sorry.
 */

void* CrossCorrelator::makeSomeOfImageThreaded(void* voidPtrArgs){

  // std::cout << "here 0" << std::endl;

  // Disgusting hacks to get ROOT threading to compile inside a class.
  CrossCorrelator::threadArgs* args = (CrossCorrelator::threadArgs*) voidPtrArgs;
  Long_t threadInd = args->threadInd;
  CrossCorrelator* ptr = args->ptr;
  AnitaPol::AnitaPol_t pol = ptr->threadPol;
  Int_t numThreads = ptr->numThreads;
  // const Int_t numThetaBinsPerThread = hImage->GetNbinsY()/numThreads;
  // const Int_t startThetaBin = threadInd*numThetaBinsPerThread;
  // const Int_t endThetaBin = startThetaBin+numThetaBinsPerThread;
  // const Int_t startPhiBin = 0;
  // const Int_t endPhiBin = hImage->GetNbinsX();


  // const Int_t numPhiBinsPerThread = (NUM_BINS_PHI*NUM_PHI)/numThreads;
  const Int_t numPhiSectorsPerThread = NUM_PHI/numThreads;
  const Int_t startThetaBin = 0;
  const Int_t endThetaBin = NUM_BINS_THETA;
  const Int_t startPhiSector = threadInd*numPhiSectorsPerThread;
  const Int_t endPhiSector = startPhiSector+numPhiSectorsPerThread;


  // std::cout << numPhiSectorsPerThread << "\t" << startThetaBin << "\t" << endThetaBin << "\t" << startPhiSector << "\t" << endPhiSector << std::endl;


  ptr->threadImagePeak[threadInd] = -DBL_MAX;
  Int_t peakPhiBin = -1;
  Int_t peakThetaBin = -1;

  // zero internal map
  for(Int_t phiSector = startPhiSector; phiSector < endPhiSector; phiSector++){
    Int_t startPhiBin = phiSector*NUM_BINS_PHI;
    Int_t endPhiBin = (phiSector+1)*NUM_BINS_PHI;
    // std::cout << startPhiBin << "\t" << endPhiBin << std::endl;
    for(Int_t phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){
      for(Int_t thetaBin = startThetaBin; thetaBin < endThetaBin; thetaBin++){
	ptr->coarseMap[pol][phiBin][thetaBin] = 0;
      }
    }
  }

  // std::cout << "here 1" << std::endl;

  std::vector<Int_t>* combosToUse = NULL;
  for(Int_t phiSector = startPhiSector; phiSector < endPhiSector; phiSector++){
    combosToUse = &ptr->combosToUseGlobal[phiSector];

    Int_t startPhiBin = phiSector*NUM_BINS_PHI;
    Int_t endPhiBin = (phiSector+1)*NUM_BINS_PHI;
    for(UInt_t comboInd=0; comboInd<combosToUse->size(); comboInd++){
      Int_t combo = combosToUse->at(comboInd);
      if(ptr->kOnlyThisCombo >= 0 && combo!=ptr->kOnlyThisCombo){
	continue;
      }
      for(Int_t phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){

	for(Int_t thetaBin = startThetaBin; thetaBin < endThetaBin; thetaBin++){
	  Int_t offsetLow = ptr->offsetLows[pol][combo][phiBin][thetaBin];
	  // Int_t offsetLow = 0;
	  // Int_t offsetLow = ptr->offsetLows[pol][phiBin][combo][thetaBin];
	  Double_t c1 = ptr->crossCorrelations[pol][combo][offsetLow];
	  Double_t c2 = ptr->crossCorrelations[pol][combo][offsetLow+1];
	  Double_t cInterp = ptr->interpPreFactors[pol][combo][phiBin][thetaBin]*(c2 - c1) + c1;

	  ptr->coarseMap[pol][phiBin][thetaBin] += cInterp;
	}
      }
    }
  }

  // std::cout << "here 2" << std::endl;

  for(Int_t phiSector = startPhiSector; phiSector < endPhiSector; phiSector++){
    combosToUse = &ptr->combosToUseGlobal[phiSector];

    Double_t normFactor = ptr->kOnlyThisCombo < 0 && combosToUse->size() > 0 ? combosToUse->size() : 1;
    // absorb the removed inverse FFT normalization
    normFactor*=(ptr->numSamples*ptr->numSamples);

    Int_t startPhiBin = phiSector*NUM_BINS_PHI;
    Int_t endPhiBin = (phiSector+1)*NUM_BINS_PHI;
    for(Int_t phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){
      for(Int_t thetaBin = startThetaBin; thetaBin < endThetaBin; thetaBin++){

	// if(combosToUse->size()>0 && ptr->kOnlyThisCombo < 0){
	//   ptr->coarseMap[pol][phiBin][thetaBin] /= combosToUse->size();
	// }
	ptr->coarseMap[pol][phiBin][thetaBin]/=normFactor;
	if(ptr->coarseMap[pol][phiBin][thetaBin] > ptr->threadImagePeak[threadInd]){
	  ptr->threadImagePeak[threadInd] = ptr->coarseMap[pol][phiBin][thetaBin];
	  peakPhiBin = phiBin;
	  peakThetaBin = thetaBin;
	}
      }
    }
  }

  // std::cout << "here 3" << std::endl;

  ptr->threadPeakPhiDeg[threadInd] = ptr->phiWaveLookup[peakPhiBin]*TMath::RadToDeg();
  ptr->threadPeakThetaDeg[threadInd] = ptr->thetaWaves[peakThetaBin]*TMath::RadToDeg();

  return NULL;
}






//---------------------------------------------------------------------------------------------------------
/**
 * @brief Wrapper function which launches the threaded functions, which fill the zoomed  maps.
 *
 * @param pol is the polarization.
 * @param imagePeak is the maximum value.
 * @param peakPhiDeg is azimuth of the maximum value (Degrees) relative to the ADU5 aft-fore.
 * @param peakThetaDeg is the elevation of the maximum value (Degrees).
 * @param zoomCenterPhiDeg should be the phi (deg) of the coarse map to reconstruct around.
 * @param zoomCenterThetaDeg should be the theta (deg) of the coarse map to reconstruct around.
 *
 * This function is also responsible for merging the peak finding results of each of the threads.
 */
void CrossCorrelator::reconstructZoom(AnitaPol::AnitaPol_t pol, Double_t& imagePeak,
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

  threadPol = pol;
  threadPeakIndex = peakIndex; // used to indentify the map
  Double_t deltaPhiDegPhi0 = RootTools::getDeltaAngleDeg(zoomCenterPhiDeg, getBin0PhiDeg());
  deltaPhiDegPhi0 = deltaPhiDegPhi0 < 0 ? deltaPhiDegPhi0 + DEGREES_IN_CIRCLE : deltaPhiDegPhi0;

  Int_t phiSector = floor(deltaPhiDegPhi0/PHI_RANGE);
  doUpsampledCrossCorrelationsThreaded(pol, phiSector); // sets threadPhiSector
  // akimaUpsampleCrossCorrelations(pol, phiSector); // sets threadPhiSector

  zoomCenterPhiDeg = (TMath::Nint(zoomCenterPhiDeg/ZOOM_BIN_SIZE_PHI))*ZOOM_BIN_SIZE_PHI;
  zoomCenterThetaDeg = (TMath::Nint(zoomCenterThetaDeg/ZOOM_BIN_SIZE_THETA))*ZOOM_BIN_SIZE_THETA;

  zoomPhiMin[pol] = zoomCenterPhiDeg - PHI_RANGE_ZOOM/2;
  zoomThetaMin[pol] = zoomCenterThetaDeg - THETA_RANGE_ZOOM/2;

  // LAUNCH THREADS HERE
  if(numThreads > 1){
    std::vector<TThread*> mapThreads(numThreads, NULL);
    for(long threadInd=0; threadInd<numThreads; threadInd++){
      TString name = TString::Format("threadZoomMap%ld", threadInd);
      CrossCorrelator::threadArgs threadArgVals;
      threadArgVals.threadInd = threadInd;
      threadArgVals.ptr = this;
      mapThreads.at(threadInd) = new TThread(name.Data(),
					     CrossCorrelator::makeSomeOfZoomImageThreaded,
					     (void*)&threadArgVals);
    }

    for(long threadInd=0; threadInd<numThreads; threadInd++){
      mapThreads.at(threadInd)->Run();
    }

    for(long threadInd=0; threadInd<numThreads; threadInd++){
      mapThreads.at(threadInd)->Join();
    }

    for(long threadInd=0; threadInd<numThreads; threadInd++){
      delete mapThreads.at(threadInd);
    }
    mapThreads.clear();
  }
  else{
    CrossCorrelator::threadArgs threadArgVals;
    threadArgVals.threadInd = 0;
    threadArgVals.ptr = this;
    CrossCorrelator::makeSomeOfZoomImageThreaded((void*)&threadArgVals);
  }


  // Combine peak search results from each thread.
  imagePeak = -DBL_MAX;
  Int_t peakPhiBin, peakThetaBin;
  for(Long_t threadInd=0; threadInd<numThreads; threadInd++){
    if(threadImagePeakZoom[threadInd] > imagePeak){
      imagePeak = threadImagePeakZoom[threadInd];
      peakPhiDeg = threadPeakPhiDegZoom[threadInd];
      peakThetaDeg = threadPeakThetaDegZoom[threadInd];
      peakPhiBin = threadPeakPhiBinZoom[threadInd];
      peakThetaBin = threadPeakThetaBinZoom[threadInd];
    }
  }
  // const int numAkimaBins = 5;

  // Int_t startThetaBin = peakThetaBin - numAkimaBins/2;
  // if(startThetaBin < 0){
  //   startThetaBin = 0;
  // }
  // Int_t endThetaBin = startThetaBin + numAkimaBins;
  // if(endThetaBin >= NUM_BINS_THETA_ZOOM){
  //   endThetaBin = NUM_BINS_THETA_ZOOM - 1;
  //   startThetaBin = endThetaBin - numAkimaBins;
  // }
  // std::vector<Double_t> thetaVecTheta(numAkimaBins, 0);
  // std::vector<Double_t> thetaVecVals(numAkimaBins, 0);
  // Int_t thetaInd = 0;
  // for(int thetaBin = startThetaBin; thetaBin < endThetaBin; thetaBin++){
  //   thetaVecTheta.at(thetaInd) = zoomThetaMin[pol] + thetaBin*ZOOM_BIN_SIZE_THETA;
  //   thetaVecVals.at(thetaInd) = fineMap[pol][peakIndex][thetaBin][peakPhiBin];
  //   thetaInd++;
  // }
  // ROOT::Math::Interpolator thetaInterp(thetaVecTheta,thetaVecVals,ROOT::Math::Interpolation::kAKIMA);

  // // peakThetaDeg = ;

  // double thetaDeg = thetaVecTheta.at(0);
  // double peakTheta = -1;
  // while(thetaDeg <= thetaVecTheta.at(numAkimaBins-1)){
  //   Double_t val = thetaInterp.Eval(thetaDeg);
  //   if(val > imagePeak){
  //     peakTheta = val;
  //     peakThetaDeg = thetaDeg;
  //   }
  //   thetaDeg += AKIMA_PEAK_INTERP_DELTA_THETA_DEG;
  // }


  // Int_t startPhiBin = peakPhiBin - numAkimaBins/2;
  // if(startPhiBin < 0){
  //   startPhiBin = 0;
  // }
  // Int_t endPhiBin = startPhiBin + numAkimaBins;
  // if(endPhiBin >= NUM_BINS_PHI_ZOOM){
  //   endPhiBin = NUM_BINS_PHI_ZOOM - 1;
  //   startPhiBin = endPhiBin - numAkimaBins;
  // }
  // std::vector<Double_t> phiVecPhi(numAkimaBins, 0);
  // std::vector<Double_t> phiVecVals(numAkimaBins, 0);
  // Int_t phiInd = 0;
  // for(int phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){
  //   phiVecPhi.at(phiInd) = zoomPhiMin[pol] + phiBin*ZOOM_BIN_SIZE_PHI;
  //   phiVecVals.at(phiInd) = fineMap[pol][peakIndex][peakThetaBin][phiBin];
  //   phiInd++;
  // }

  // ROOT::Math::Interpolator phiInterp(phiVecPhi,phiVecVals,ROOT::Math::Interpolation::kAKIMA);

  // double phiDeg = phiVecPhi.at(0);
  // double peakPhi = -1;
  // while(phiDeg <= phiVecPhi.at(numAkimaBins-1)){
  //   Double_t val = phiInterp.Eval(phiDeg);
  //   if(val > peakPhi){
  //     peakPhi = val;
  //     peakPhiDeg = phiDeg;
  //   }
  //   phiDeg += AKIMA_PEAK_INTERP_DELTA_PHI_DEG;
  // }

  // imagePeak = peakTheta > peakPhi ?  peakTheta : peakPhi;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Static member function which fills the interferometric maps.
 *
 * @param voidPtrArgs contains a pointer to a CrossCorrelator::threadArgs struct
 *
 * This function contains the meat and bones of this class.
 * I've really tried to optimize this for speed, which means it's not very readable, sorry.
 */
void* CrossCorrelator::makeSomeOfZoomImageThreaded(void* voidPtrArgs){
  // Disgusting hacks to get ROOT threading to compile inside a class.

  CrossCorrelator::threadArgs* args = (CrossCorrelator::threadArgs*) voidPtrArgs;
  Long_t threadInd = args->threadInd;
  CrossCorrelator* ptr = args->ptr;
  AnitaPol::AnitaPol_t pol = ptr->threadPol;
  Int_t peakIndex = ptr->threadPeakIndex;
  Int_t numThreads = ptr->numThreads;

  const Int_t numThetaBinsPerThread = NUM_BINS_THETA_ZOOM/numThreads;
  const Int_t startThetaBin = threadInd*numThetaBinsPerThread;
  const Int_t endThetaBin = startThetaBin+numThetaBinsPerThread;
  const Int_t startPhiBin = 0;
  const Int_t endPhiBin = NUM_BINS_PHI_ZOOM;


  // const Int_t numPhiBinsPerThread = NUM_BINS_PHI_ZOOM/numThreads;
  // const Int_t startThetaBin = 0;
  // const Int_t endThetaBin = NUM_BINS_THETA_ZOOM;
  // const Int_t startPhiBin = threadInd*numPhiBinsPerThread;
  // const Int_t endPhiBin = startPhiBin+numPhiBinsPerThread;

  // TThread::Lock();
  // std::cerr << threadInd << "\t" << startPhiBin << "\t" << endPhiBin << std::endl;
  // TThread::UnLock();

  ptr->threadImagePeakZoom[threadInd] = -DBL_MAX;
  Int_t peakPhiBin = -1;
  Int_t peakThetaBin = -1;

  Int_t kUseOffAxisDelay = ptr->kUseOffAxisDelay;

  Int_t phiSector = ptr->threadPhiSector;
  std::vector<Int_t>* combosToUse = &ptr->combosToUseGlobal[phiSector];

  Int_t phiZoomBase = TMath::Nint((ptr->zoomPhiMin[pol] - ptr->minPhiDegZoom)/ZOOM_BIN_SIZE_PHI);
  Int_t thetaZoomBase = TMath::Nint((ptr->zoomThetaMin[pol] - ptr->minThetaDegZoom)/ZOOM_BIN_SIZE_THETA);

  for(Int_t thetaBin = startThetaBin; thetaBin < endThetaBin; thetaBin++){
    for(Int_t phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){
      ptr->fineMap[pol][peakIndex][thetaBin][phiBin]=0;
    }
  }

  const Int_t offset = ptr->numSamplesUpsampled/2;
  for(UInt_t comboInd=0; comboInd<combosToUse->size(); comboInd++){
    Int_t combo = combosToUse->at(comboInd);
    if(ptr->kOnlyThisCombo >= 0 && combo!=ptr->kOnlyThisCombo){
      continue;
    }
    for(Int_t thetaBin = startThetaBin; thetaBin < endThetaBin; thetaBin++){
      Int_t zoomThetaInd = thetaZoomBase + thetaBin;
      Double_t partBA = ptr->partBAsZoom[pol][combo][zoomThetaInd];
      // Double_t dtFactor = ptr->zoomedCosThetaWaves[zoomThetaInd]/SPEED_OF_LIGHT_NS;
      Double_t dtFactor = ptr->dtFactors[zoomThetaInd];
      for(Int_t phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){
	Int_t zoomPhiInd = phiZoomBase + phiBin;
	Double_t offsetLowDouble = dtFactor*(partBA - ptr->part21sZoom[pol][combo][zoomPhiInd]);
	offsetLowDouble += kUseOffAxisDelay > 0 ? ptr->offAxisDelaysDivided[pol][combo][zoomPhiInd] : 0;
	// hack for floor()
	Int_t offsetLow = (int) offsetLowDouble - (offsetLowDouble < (int) offsetLowDouble);

	Double_t deltaT = (offsetLowDouble - offsetLow);
	offsetLow += offset;
	Double_t c1 = ptr->crossCorrelationsUpsampled[pol][combo][offsetLow];
	Double_t c2 = ptr->crossCorrelationsUpsampled[pol][combo][offsetLow+1];
	Double_t cInterp = deltaT*(c2 - c1) + c1;

	ptr->fineMap[pol][peakIndex][thetaBin][phiBin] += cInterp;
      }
    }
  }

  Double_t normFactor = ptr->kOnlyThisCombo < 0 && combosToUse->size() > 0 ? combosToUse->size() : 1;
  // absorb the removed inverse FFT normalization
  normFactor*=(ptr->numSamples*ptr->numSamples);

  for(Int_t thetaBin = startThetaBin; thetaBin < endThetaBin; thetaBin++){
    for(Int_t phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){
      ptr->fineMap[pol][peakIndex][thetaBin][phiBin] /= normFactor;
      // if(combosToUse->size()>0 && ptr->kOnlyThisCombo < 0){
      // 	ptr->fineMap[pol][peakIndex][thetaBin][phiBin] /= combosToUse->size();
      // }


      if(ptr->fineMap[pol][peakIndex][thetaBin][phiBin] > ptr->threadImagePeakZoom[threadInd]){
	ptr->threadImagePeakZoom[threadInd] = ptr->fineMap[pol][peakIndex][thetaBin][phiBin];
	peakPhiBin = phiBin;
	peakThetaBin = thetaBin;
      }
    }
  }

  // ptr->threadPeakPhiDeg[threadInd] = hImage->GetXaxis()->GetBinLowEdge(peakPhiBin+1);
  // ptr->threadPeakThetaDeg[threadInd] = hImage->GetYaxis()->GetBinLowEdge(peakThetaBin+1);


  ptr->threadPeakPhiDegZoom[threadInd] = ptr->zoomPhiMin[pol] + peakPhiBin*ZOOM_BIN_SIZE_PHI;
  ptr->threadPeakThetaDegZoom[threadInd] = ptr->zoomThetaMin[pol] + peakThetaBin*ZOOM_BIN_SIZE_THETA;
  ptr->threadPeakPhiBinZoom[threadInd] = peakPhiBin;
  ptr->threadPeakThetaBinZoom[threadInd] = peakThetaBin;

  return 0;

}









//---------------------------------------------------------------------------------------------------------
/**
 * @brief Deletes the waveform TGraphs in memory and removes dangling pointers.
 *
 * @param pol is the polarization to delete.
 */
void CrossCorrelator::deleteAllWaveforms(AnitaPol::AnitaPol_t pol){
  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
    if(grs[pol][ant]){
      delete grs[pol][ant];
      grs[pol][ant] = NULL;
    }
    if(grsResampled[pol][ant]){
      delete grsResampled[pol][ant];
      grsResampled[pol][ant] = NULL;
    }
  }
}








//---------------------------------------------------------------------------------------------------------
/**
 * @brief Used to insert phase center geometry files from Linda.
 *
 * @param pathToLindasFile is the relative path to Linda's file.
 * @param pol is the polarization of the channels under test.
 * @returns 0 if everything went without a problem, 1 if there was a problem.
*/
Int_t CrossCorrelator::directlyInsertGeometry(TString pathToLindasFile, AnitaPol::AnitaPol_t pol){

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








//---------------------------------------------------------------------------------------------------------
/**
 * @brief Called by wrapper functions. Turns internal correlation arrays into TGraphs to be seen by humans.
 *
 * @param numSamps tells this function whether to use the interpolated or un-interpolated correlations.
 * @param pol is the polarization.
 * @param ant1 is the first antenna.
 * @param ant2 is the second antenna.
 * @returns the cross-correlation TGraph.
*/
TGraph* CrossCorrelator::getCrossCorrelationGraphWorker(Int_t numSamps, AnitaPol::AnitaPol_t pol,
							Int_t ant1, Int_t ant2){
  // Primarily for debugging, put cross correlations in a TGraph

  Int_t combo = comboIndices[ant1][ant2];
  if(combo < 0){
    return NULL;
  }

  Double_t graphDt = correlationDeltaT;
  Double_t* corrPtr = crossCorrelationsUpsampled[pol][combo];
  if(numSamps != numSamplesUpsampled){
    numSamps = numSamples;
    graphDt = nominalSamplingDeltaT;
    corrPtr = crossCorrelations[pol][combo];
  }

  // Could actually perform the calculation here... but I'll leave it NULL for now
  if(corrPtr==NULL){
    return NULL;
  }

  std::vector<Double_t> offsets = std::vector<Double_t>(numSamps, 0);
  std::vector<Double_t> corrs = std::vector<Double_t>(numSamps, 0);

  // put the correlations in, in time order in memory
  for(Int_t i=0; i<numSamps; i++){
    Int_t offset = (i - numSamps/2);
    offsets.at(i) = offset*graphDt;
    corrs.at(i) = corrPtr[i];
  }


  TGraph* gr = new TGraph(numSamps,  &offsets[0], &corrs[0]);
  if(numSamps != numSamplesUpsampled){
    gr->SetName(TString::Format("grCorr_%d_%d", ant1, ant2));
    gr->SetTitle(TString::Format("Cross Correlation ant1 = %d, ant2 = %d", ant1, ant2));
  }
  else{
    gr->SetName(TString::Format("grCorrUpsampled_%d_%d", ant1, ant2));
    gr->SetTitle(TString::Format("Upsampled Cross Correlation ant1 = %d, ant2 = %d", ant1, ant2));
  }

  return gr;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Turns internal correlation arrays into TGraphs to be seen by humans.
 *
 * @param pol is the polarization.
 * @param ant1 is the first antenna.
 * @param ant2 is the second antenna.
 * @returns the cross-correlation TGraph.
*/
TGraph* CrossCorrelator::getCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2){
  // Primarily for debugging, put cross correlations in a TGraph
  return getCrossCorrelationGraphWorker(numSamples, pol, ant1, ant2);
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Turns internal upsampled correlation arrays into TGraphs to be seen by humans.
 *
 * @param pol is the polarization.
 * @param ant1 is the first antenna.
 * @param ant2 is the second antenna.
 * @returns the cross-correlation TGraph.
*/
TGraph* CrossCorrelator::getUpsampledCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2){
  return getCrossCorrelationGraphWorker(numSamplesUpsampled, pol, ant1, ant2);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Finds the phi-sector closest to a particular phi direction.
 *
 * @param pol is the polarization.
 * @param phiDeg is the azimuthal direction (Degrees) relative to ADU5 aft-fore.
 * @returns the phiSector (0->15)
*/
Int_t CrossCorrelator::getPhiSectorOfAntennaClosestToPhiDeg(AnitaPol::AnitaPol_t pol, Double_t phiDeg){
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



//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates the coherently summed waveform from the zero padded FFTs held in memory
 * @param pol is the polarization
 * @param phiDeg is the incoming phi direction (Degrees) relative to the ADU5 aft-fore.
 * @param thetaDeg is the incoming theta direction (Degrees).
 * @param maxDeltaPhiSect is the number of phi-sectors to contribute either side of the incoming phi-direction.
 * @param snr is an estimate of the signal-to-noise ratio of the coherently summed waveform, using the local max-to-min of the coherent waveform and the rms of the first few ns of the contributing uninterpolated waveforms.
 * @return the upsampled coherently summed waveform made from the zero padded ffts.
*/
TGraph* CrossCorrelator::makeUpsampledCoherentlySummedWaveform(AnitaPol::AnitaPol_t pol, Double_t phiDeg,
							       Double_t thetaDeg, Int_t maxDeltaPhiSect,
							       Double_t& snr){
  return makeCoherentWorker(pol, phiDeg, thetaDeg, maxDeltaPhiSect, snr, numSamplesUpsampled);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates the coherently summed waveform from the FFTs held in memory
 * @param pol is the polarization
 * @param phiDeg is the incoming phi direction (Degrees) relative to the ADU5 aft-fore.
 * @param thetaDeg is the incoming theta direction (Degrees).
 * @param maxDeltaPhiSect is the number of phi-sectors to contribute either side of the incoming phi-direction.
 * @param snr is an estimate of the signal-to-noise ratio of the coherently summed waveform, using the local max-to-min of the coherent waveform and the rms of the first few ns of the contributing uninterpolated waveforms.
 * @return the coherently summed waveform made from the (non-padded) ffts.
*/
TGraph* CrossCorrelator::makeCoherentlySummedWaveform(AnitaPol::AnitaPol_t pol, Double_t phiDeg,
						      Double_t thetaDeg, Int_t maxDeltaPhiSect,
						      Double_t& snr){
  return makeCoherentWorker(pol, phiDeg, thetaDeg, maxDeltaPhiSect, snr, numSamples);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Worker function to create the coherently summed waveform from either the regular ffts or the padded ffts.
 * @param pol is the polarization
 * @param phiDeg is the incoming phi direction (Degrees) relative to the ADU5 aft-fore.
 * @param thetaDeg is the incoming theta direction (Degrees).
 * @param maxDeltaPhiSect is the number of phi-sectors to contribute either side of the incoming phi-direction.
 * @param snr is an estimate of the signal-to-noise ratio of the coherently summed waveform, using the local max-to-min of the coherent waveform and the rms of the first few ns of the contributing uninterpolated waveforms.
 * @param nSamp is the number of samples in the time domain of the ffts (numSamples or numSamplesUpsampled)
 * @return the coherently summed waveform.
*/
TGraph* CrossCorrelator::makeCoherentWorker(AnitaPol::AnitaPol_t pol, Double_t phiDeg,
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

  Double_t theDeltaT = nSamp == numSamples ? nominalSamplingDeltaT : correlationDeltaT;

  Double_t phiRad = phiDeg*TMath::DegToRad();
  Double_t thetaRad = thetaDeg*TMath::DegToRad();
  Int_t numAnts = 0;

  Int_t centerPhiSector = getPhiSectorOfAntennaClosestToPhiDeg(pol, phiDeg);

  const Int_t firstAnt = centerPhiSector;


  std::pair<Int_t, Int_t> key(nSamp, 0);
  // vArray is actually internal memory managed by FancyFFTs... don't delete this!!!
  Double_t* vArray = FancyFFTs::getRealArray(key);
  std::complex<double>* cArray = FancyFFTs::getComplexArray(key);
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
  Double_t t0 = grsResampled[pol][firstAnt]->GetX()[0];
  for(Int_t samp=0; samp<nSamp; samp++){
    double time = t0 + samp*theDeltaT;
    tArray.at(samp) = time;
    // vArray[samp] *= interpRMS[pol][firstAnt]*lengthFftNorm; // Undo the normalization.
    vArray[samp] = grsResampled[pol][firstAnt]->Eval(time);
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

	Double_t t0 = grsResampled[pol][ant]->GetX()[0];

	for(Int_t samp=0; samp<grCoherent->GetN(); samp++){
	  double time = t0 + samp*theDeltaT + dt;

	  // grCoherent->GetY()[samp] += grsResampled[pol][ant]->Eval(time - deltaT); //interpRMS[pol][ant];
	  // grCoherent->GetY()[samp] += grsResampled[pol][ant]->Eval(time + deltaT); //interpRMS[pol][ant];
	  grCoherent->GetY()[samp] += grsResampled[pol][ant]->Eval(time + dt); //interpRMS[pol][ant];

	}
      }

      // Here we look at the RMS of the first few ns of the uninterpolated waveforms
      const Double_t timeToEvalRms = 10; // ns
      // for(Int_t samp3=0; samp3 < grs[pol][ant]->GetN(); samp3++){
      // start time of the interpolated waveform... but this could be front padded with zero...
      Double_t t0 = grsResampled[pol][ant]->GetX()[0];
      // this is when we've got the good stuff..., we really want to start counting from here
      Double_t t0Good = grs[pol][ant]->GetX()[0];
      for(Int_t samp3=0; samp3 < nSamp; samp3++){
	Double_t t = t0 + samp3*theDeltaT;
	if(t >= t0Good && t < t0Good + timeToEvalRms){
	  // rms += grs[pol][ant]->GetY()[samp3]*grs[pol][ant]->GetY()[samp3];
	  double V = grsResampled[pol][ant]->GetY()[samp3];
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
    TString name = nSamp == numSamples ? "grCoherent" : "grInterpCoherent";
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





// https://www.codeproject.com/Articles/3515/How-To-Organize-Template-Source-Code
// http://stackoverflow.com/questions/115703/storing-c-template-function-definitions-in-a-cpp-file
// and now I see why people hate template meta-programming...
template void CrossCorrelator::getNormalizedInterpolatedTGraphs<UsefulAnitaEvent>(UsefulAnitaEvent* realEvent, AnitaPol::AnitaPol_t pol);
template void CrossCorrelator::getNormalizedInterpolatedTGraphs<FilteredAnitaEvent>(FilteredAnitaEvent* realEvent, AnitaPol::AnitaPol_t pol);


template void CrossCorrelator::correlateEvent<UsefulAnitaEvent>(UsefulAnitaEvent* realEvent);
template void CrossCorrelator::correlateEvent<FilteredAnitaEvent>(FilteredAnitaEvent* realEvent);


template void CrossCorrelator::correlateEvent<UsefulAnitaEvent>(UsefulAnitaEvent* realEvent, AnitaPol::AnitaPol_t pol);
template void CrossCorrelator::correlateEvent<FilteredAnitaEvent>(FilteredAnitaEvent* realEvent, AnitaPol::AnitaPol_t pol);

template void CrossCorrelator::reconstructEvent<UsefulAnitaEvent>(UsefulAnitaEvent* usefulEvent, Int_t numFinePeaks=MAX_NUM_PEAKS, Int_t numCoarsePeaks=MAX_NUM_PEAKS);
template void CrossCorrelator::reconstructEvent<FilteredAnitaEvent>(FilteredAnitaEvent* usefulEvent, Int_t numFinePeaks=MAX_NUM_PEAKS, Int_t numCoarsePeaks=MAX_NUM_PEAKS);

template AnitaPol::AnitaPol_t CrossCorrelator::reconstructEventPeakPol<UsefulAnitaEvent>(UsefulAnitaEvent* usefulEvent, Int_t numFinePeaks=MAX_NUM_PEAKS, Int_t numCoarsePeaks=MAX_NUM_PEAKS);
template AnitaPol::AnitaPol_t CrossCorrelator::reconstructEventPeakPol<FilteredAnitaEvent>(FilteredAnitaEvent* usefulEvent, Int_t numFinePeaks=MAX_NUM_PEAKS, Int_t numCoarsePeaks=MAX_NUM_PEAKS);


template void CrossCorrelator::reconstructEvent<UsefulAnitaEvent>(UsefulAnitaEvent* usefulEvent, UsefulAdu5Pat& usefulPat, AnitaEventSummary* eventSummary);
template void CrossCorrelator::reconstructEvent<FilteredAnitaEvent>(FilteredAnitaEvent* usefulEvent, UsefulAdu5Pat& usefulPat, AnitaEventSummary* eventSummary);
