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
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Workhorse function to set internal variables.
 */
void CrossCorrelator::initializeVariables(){

  kDeltaPhiSect = 2;
  // Initialize with NULL otherwise very bad things will happen with gcc 
  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      grs[pol][ant] = NULL;
      grsResampled[pol][ant] = NULL;
      interpRMS[pol][ant] = 0;
    }
    lastEventNormalized[pol] = 0;
    eventNumber[pol] = 0;
  }

  maxDPhiDeg = 0;
  kOnlyThisCombo = -1;
  kUseOffAxisDelay = 1;
  numSamples = 2*NUM_SAMPLES; // Factor of two for padding. Turns circular xcor into linear xcor.
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

  fillDeltaTLookup();

  mapModeNames[kGlobal] = "Global";
  mapModeNames[kTriggered] = "Triggered";
  zoomModeNames[kZoomedOut] = "";
  zoomModeNames[kZoomedIn] = "Zoom";

  threadImage = NULL;
  threadPol = AnitaPol::kHorizontal;
  threadMapMode = kGlobal;
  threadZoomMode = kZoomedOut;
  threadRWave = 0;
  threadL3TrigPattern = 0;

  
  for(Long_t threadInd=0; threadInd<NUM_THREADS; threadInd++){
    CrossCorrelator::threadArgs threadArgVals;
    threadArgVals.threadInd = threadInd;
    threadArgVals.ptr = this;
    threadArgsVec.push_back(threadArgVals);

    // Creation of fftw plans is not thread safe, 
    // so we need to create plans before doing any fftw
    // plans inside a thread.
    FancyFFTs::makeNewPlanIfNeeded(numSamples, threadInd);
    FancyFFTs::makeNewPlanIfNeeded(numSamplesUpsampled, threadInd);
  }
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Prints some summary information about the CrossCorrelator to stdout.
 *
 */
void CrossCorrelator::printInfo(){
  std::cerr << __PRETTY_FUNCTION__ << std::endl;
  std::cerr << "\tupsample factor = " << UPSAMPLE_FACTOR << std::endl;
  std::cerr << "\tBin size theta (deg) = " << Double_t(THETA_RANGE)/NUM_BINS_THETA << std::endl;
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
  // Double_t phi0 = phiArrayDeg[0].at(0);
  Double_t phi0 = -45; //phiArrayDeg[0].at(0);  
  if(phi0 < -180){
    phi0+=DEGREES_IN_CIRCLE;
  }
  else if(phi0 >= 180){
    phi0-=DEGREES_IN_CIRCLE;
  }
  return phi0 - PHI_RANGE/2;
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Makes evenly re-sampled, normalized waveform graphs from the UsefulAnitaEvent.
 * 
 * @param usefulEvent points to the UsefulAnitaEvent of interest.
 * @param pol is the polarization of interest.
 */
void CrossCorrelator::getNormalizedInterpolatedTGraphs(UsefulAnitaEvent* usefulEvent,
						       AnitaPol::AnitaPol_t pol){

  // Delete any old waveforms (at start rather than end to leave in memory to be examined if need be)
  deleteAllWaveforms(pol);

  // Find the start time of all waveforms 
  std::vector<Double_t> earliestStart(NUM_POL, 100); // ns
  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
    grs[pol][ant] = usefulEvent->getGraph(ant, (AnitaPol::AnitaPol_t)pol);
      
    if(grs[pol][ant]->GetX()[0]<earliestStart.at(pol)){
      earliestStart.at(pol) = grs[pol][ant]->GetX()[0];
    }
  }

  // Interpolate with earliest start time 
  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
    grsResampled[pol][ant] = RootTools::interpolateWithStartTime(grs[pol][ant], earliestStart.at(pol),
								 nominalSamplingDeltaT, numSamples);
    Double_t mean; // don't care enough about this to store it anywhere.
    RootTools::normalize(grsResampled[pol][ant], mean, interpRMS[pol][ant]);
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
    
  lastEventNormalized[pol] = usefulEvent->eventNumber;
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Takes FFTs of the normalized, evenly resampled waveforms and puts them in memory
 * @param pol is which polarization to process
 *
 * Now also creates the zero padded FFTs used for upsampled cross-correlations.
 */
void CrossCorrelator::doFFTs(AnitaPol::AnitaPol_t pol){

  Int_t numFreqs = FancyFFTs::getNumFreqs(numSamples);
  Int_t numFreqsPadded = FancyFFTs::getNumFreqs(numSamplesUpsampled);

  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
    FancyFFTs::doFFT(numSamples, grsResampled[pol][ant]->GetY(), ffts[pol][ant]);
    FancyFFTs::zeroPadFFT(ffts[pol][ant], fftsPadded[pol][ant], numFreqs, numFreqsPadded);    
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

void CrossCorrelator::correlateEvent(UsefulAnitaEvent* usefulEvent){

  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){  
    correlateEvent(usefulEvent, (AnitaPol::AnitaPol_t)pol);
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
 * @param pol tells CrossCorrelator to only do this polarization.
 */
void CrossCorrelator::correlateEvent(UsefulAnitaEvent* usefulEvent, AnitaPol::AnitaPol_t pol){

  // Read TGraphs from events into memory (also deletes old TGraphs)
  getNormalizedInterpolatedTGraphs(usefulEvent, pol);

  // Generate set of ffts for cross correlation (each waveform only needs to be done once)
  doFFTs(pol);

  // Now cross correlate those already FFT'd waveforms
  doAllCrossCorrelationsThreaded(pol);

  // Safety check to make sure we don't do any hard work twice.
  eventNumber[pol] = usefulEvent->eventNumber;
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
  
  for(long threadInd=0; threadInd<NUM_THREADS; threadInd++){
    TString name = TString::Format("threadCorr%ld", threadInd);
    corrThreads.push_back(new TThread(name.Data(),
				      CrossCorrelator::doSomeCrossCorrelationsThreaded,
				      (void*)&threadArgsVec.at(threadInd))
			  );
  }
  
  for(long threadInd=0; threadInd<NUM_THREADS; threadInd++){
    corrThreads.at(threadInd)->Run();
  }

  for(long threadInd=0; threadInd<NUM_THREADS; threadInd++){
    corrThreads.at(threadInd)->Join();
  }

  for(long threadInd=0; threadInd<NUM_THREADS; threadInd++){
    delete corrThreads.at(threadInd);
  }

  corrThreads.clear();
  
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
  
  std::pair<UInt_t, Int_t> key(ptr->threadL3TrigPattern, ptr->kDeltaPhiSect);
  std::vector<Int_t> combosToUse = ptr->combosToUseTriggered[key];
  Int_t numCombosAllThreads = combosToUse.size();
  Int_t numCorrPerThread = numCombosAllThreads/NUM_THREADS;
  Int_t numRemainder = numCombosAllThreads%NUM_THREADS;

  Int_t startComboInd = threadInd*numCorrPerThread;

  Double_t stash[NUM_SAMPLES*UPSAMPLE_FACTOR*2];


  for(int comboInd=startComboInd; comboInd<startComboInd+numCorrPerThread; comboInd++){
    Int_t combo = combosToUse.at(comboInd);
    Int_t ant1 = ptr->comboToAnt1s.at(combo);
    Int_t ant2 = ptr->comboToAnt2s.at(combo);
    
    FancyFFTs::crossCorrelate(ptr->numSamplesUpsampled,
			      ptr->fftsPadded[pol][ant2],
			      ptr->fftsPadded[pol][ant1],
			      ptr->crossCorrelationsUpsampled[pol][combo],
			      threadInd);

    // experiment, copy negative times to behind positive times...
    for(Int_t samp=0; samp < ptr->numSamplesUpsampled; samp++){
      stash[samp] = ptr->crossCorrelationsUpsampled[pol][combo][samp];
    }
    const int offset = ptr->numSamplesUpsampled/2;

    // copies first half of original array (times >= 0) into second half of internal storage
    for(Int_t samp=0; samp < ptr->numSamplesUpsampled/2; samp++){
      ptr->crossCorrelationsUpsampled[pol][combo][samp+offset] = stash[samp];
    }
    // copies second half of original array (times < 0) into first half of internal storage
    for(Int_t samp=ptr->numSamplesUpsampled/2; samp < ptr->numSamplesUpsampled; samp++){
      ptr->crossCorrelationsUpsampled[pol][combo][samp-offset] = stash[samp];
    }    
  }
  
  if(threadInd < numRemainder){
    Int_t numDoneInAllThreads = NUM_THREADS*numCorrPerThread;
    Int_t comboInd = numDoneInAllThreads + threadInd;
    Int_t combo = combosToUse.at(comboInd);
    Int_t ant1 = ptr->comboToAnt1s.at(combo);
    Int_t ant2 = ptr->comboToAnt2s.at(combo);

    FancyFFTs::crossCorrelate(ptr->numSamplesUpsampled,
			      ptr->fftsPadded[pol][ant2],
			      ptr->fftsPadded[pol][ant1],
			      ptr->crossCorrelationsUpsampled[pol][combo],
			      threadInd);

    // experiment, copy negative times to behind positive times...
    for(Int_t samp=0; samp < ptr->numSamplesUpsampled; samp++){
      stash[samp] = ptr->crossCorrelationsUpsampled[pol][combo][samp];
    }
    const int offset = ptr->numSamplesUpsampled/2;

    // copies first half of original array (times >= 0) into second half of internal storage
    for(Int_t samp=0; samp < ptr->numSamplesUpsampled/2; samp++){
      ptr->crossCorrelationsUpsampled[pol][combo][samp+offset] = stash[samp];
    }
    // copies second half of original array (times < 0) into first half of internal storage
    for(Int_t samp=ptr->numSamplesUpsampled/2; samp < ptr->numSamplesUpsampled; samp++){
      ptr->crossCorrelationsUpsampled[pol][combo][samp-offset] = stash[samp];
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
  AnitaPol::AnitaPol_t pol = ptr->threadPol;

  Int_t numCorrPerThread = NUM_COMBOS/NUM_THREADS;
  Int_t startCombo = threadInd*numCorrPerThread;

  Double_t stash[NUM_SAMPLES*2];
  
  for(int combo=startCombo; combo<startCombo+numCorrPerThread; combo++){
    Int_t ant1 = ptr->comboToAnt1s.at(combo);
    Int_t ant2 = ptr->comboToAnt2s.at(combo);
    FancyFFTs::crossCorrelate(ptr->numSamples,
			      ptr->ffts[pol][ant2],
			      ptr->ffts[pol][ant1],
			      ptr->crossCorrelations[pol][combo],
			      threadInd);

    // experiment, copy negative times to behind positive times...
    for(Int_t samp=0; samp < ptr->numSamples; samp++){
      stash[samp] = ptr->crossCorrelations[pol][combo][samp];
    }
    const int offset = ptr->numSamples/2;

    // copies first half of original array (times >= 0) into second half of internal storage
    for(Int_t samp=0; samp < ptr->numSamples/2; samp++){
      ptr->crossCorrelations[pol][combo][samp+offset] = stash[samp];
    }
    // copies second half of original array (times < 0) into first half of internal storage
    for(Int_t samp=ptr->numSamples/2; samp < ptr->numSamples; samp++){
      ptr->crossCorrelations[pol][combo][samp-offset] = stash[samp];
    }
  }  
  return 0;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Static member function which generates the finely binned set of cross correlations from the FFTs held in memory.
 *
 * @param pol tells CrossCorrelator to only do this polarization.
 * @param l3TrigPattern is used to figure out which finely binned cross correlations are required.
 */
void CrossCorrelator::doUpsampledCrossCorrelationsThreaded(AnitaPol::AnitaPol_t pol, UShort_t l3TrigPattern){

  threadL3TrigPattern = l3TrigPattern;
  threadPol = pol;

  for(long threadInd=0; threadInd<NUM_THREADS; threadInd++){
    TString name = TString::Format("threadUpsampledCorr%ld", threadInd);
    upsampledCorrThreads.push_back(new TThread(name.Data(),
    					       CrossCorrelator::doSomeUpsampledCrossCorrelationsThreaded,
    					       (void*)&threadArgsVec.at(threadInd))
    				   );
  }
  
  for(long threadInd=0; threadInd<NUM_THREADS; threadInd++){
    upsampledCorrThreads.at(threadInd)->Run();
  }

  for(long threadInd=0; threadInd<NUM_THREADS; threadInd++){
    upsampledCorrThreads.at(threadInd)->Join();
  }

  for(long threadInd=0; threadInd<NUM_THREADS; threadInd++){
    delete upsampledCorrThreads.at(threadInd);
  }
  upsampledCorrThreads.clear();
  
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
  Int_t offset = maxInd >= numSamples/2 ? maxInd - numSamples : maxInd;
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
  Int_t offset = maxInd >= numSamplesUpsampled/2 ? maxInd - numSamplesUpsampled : maxInd;
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
 * @brief Get the off axis delay between a pair of antennas for an incoming plane wave
 *
 * @param pol is the polarization of the antennas
 * @param ant1 is the first antenna 
 * @param ant2 is the second antenna 
 * @param phiWave is the phi direction of the incoming wave in radians (relative to the ADU5)
 * @returns the off-axis delay in ns
 */
inline Double_t CrossCorrelator::getOffAxisDelay(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
					  Double_t phiWave){
  // From Linda's fits
  // -8.41435e-06 (quadratic term)
  // 1.15582e-08 (quartic term)

  Double_t phiWaveDeg = phiWave*TMath::RadToDeg();
  Double_t deltaPhiDeg1 = RootTools::getDeltaAngleDeg(phiArrayDeg[pol].at(ant1), phiWaveDeg);
  Double_t deltaPhiDeg2 = RootTools::getDeltaAngleDeg(phiArrayDeg[pol].at(ant2), phiWaveDeg);

  // const Double_t maxDeltaPhiDeg = 45;
  // // if(TMath::Abs(deltaPhiDeg1) > maxDeltaPhiDeg || TMath::Abs(deltaPhiDeg2) > maxDeltaPhiDeg){
  // //   return 0;
  // // }
  // deltaPhiDeg1 = TMath::Abs(deltaPhiDeg1) > maxDeltaPhiDeg ? maxDeltaPhiDeg : deltaPhiDeg1;
  // deltaPhiDeg2 = TMath::Abs(deltaPhiDeg2) > maxDeltaPhiDeg ? maxDeltaPhiDeg : deltaPhiDeg2;
  // // std::cout << deltaPhiDeg1 << "\t" << deltaPhiDeg2 << std::endl;

  // const Double_t quadraticTerm = -8.41435e-06;
  // const Double_t quarticTerm = 1.15582e-08;

  // Double_t delay1 = quadraticTerm*pow(deltaPhiDeg1,2) + quarticTerm*pow(deltaPhiDeg1, 4);
  // Double_t delay2 = quadraticTerm*pow(deltaPhiDeg2,2) + quarticTerm*pow(deltaPhiDeg2, 4);

  
  const int nPowsOf2 = 6;
  Double_t params[nPowsOf2] = {-1.68751e-05,
			       2.77815e-08,
			       -8.29351e-12,
			       1.15064e-15,
			       -7.71170e-20,
			       1.99661e-24};
  
  Double_t delay1 = 0;
  Double_t delay2 = 0;

  for(Int_t powInd=0; powInd < nPowsOf2; powInd++){
    delay1 += params[powInd]*pow(deltaPhiDeg1, 2*powInd);
    delay2 += params[powInd]*pow(deltaPhiDeg2, 2*powInd);    
  }
  
  return delay2 - delay1;
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
inline Double_t CrossCorrelator::getDeltaTExpected(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t phiWave, Double_t thetaWave){
  
  Double_t tanThetaW = tan(thetaWave);
  Double_t part1 = zArray[pol].at(ant1)*tanThetaW - rArray[pol].at(ant1)*cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol].at(ant1));
  Double_t part2 = zArray[pol].at(ant2)*tanThetaW - rArray[pol].at(ant2)*cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol].at(ant2));
  Double_t tdiff = 1e9*((cos(thetaWave) * (part2 - part1))/SPEED_OF_LIGHT); // Returns time in ns

  // tdiff += getOffAxisDelay(pol, ant1, ant2, phiWave);
  
  return tdiff;
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Function to calculate the time taken by a plane wave to cross the payload for the zoomed in map
 *
 * @param pol is the polarization of the antennas
 * @param ant1 is the first antenna 
 * @param ant2 is the second antenna 
 * @param phiIndex is the phi index of the zoomed phi bin, used to access cached values.
 * @param thetaIndex is the theta index of the zoomed phi bin, used to access cached values.
 * @returns time at ant2 - time at ant1 in ns.
 *
 * This function attempts to be faster by using cached sin/cos/tan values for each bin calculated in initializeVairables.
 */
inline Double_t CrossCorrelator::getDeltaTExpectedFast(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
						       Int_t phiIndex, Int_t thetaIndex){
  
  Double_t tanThetaW = zoomedTanThetaWaves[thetaIndex];  
  Double_t part1 = zArray[pol].at(ant1)*tanThetaW - zoomedCosPartLookup[phiIndex][pol][ant1];
  Double_t part2 = zArray[pol].at(ant2)*tanThetaW - zoomedCosPartLookup[phiIndex][pol][ant2];
  // Double_t tdiff = 1e9*((cos(thetaWave) * (part2 - part1))/SPEED_OF_LIGHT); // Returns time in ns 
  Double_t tdiff = 1e9*((zoomedCosThetaWaves[thetaIndex] * (part2 - part1))/SPEED_OF_LIGHT); // Returns time in ns

  // tdiff += getOffAxisDelay(pol, ant1, ant2, zoomedPhiWaveLookup[phiIndex]);
  // tdiff += offAxisDelays[pol][phiIndex][combo];
  
  return tdiff;
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Inserts the photogrammetry geometry into the CrossCorrelator internals *AND* the AnitaEventCalibrator internals
 */
void CrossCorrelator::insertPhotogrammetryGeometry(){
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->useKurtAnita3Numbers(1);
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
  geom->useKurtAnita3Numbers(0);
  
}



//---------------------------------------------------------------------------------------------------------
/**
 * @brief Function to calculate the time taken for a spherical wave to pass between antennas on the payload.
 *
 * @param pol is the polarization of the antennas
 * @param ant1 is the first antenna 
 * @param ant2 is the second antenna 
 * @param phiWave is the phi direction of the incoming wave in radians (relative to the ADU5)
 * @param thetaWave is the theta direction of the incoming wave in radians (relative to the ADU5)
 * @param rWave is the distance from the payload in meters
 * @returns time at ant2 - time at ant1 in ns.
 *
 * I originally made this function to look for the source of the mysterious self triggered blasts.
 * But I never really got anywhere, so I don't really know how well this works.
 */
Int_t CrossCorrelator::getDeltaTExpectedSpherical(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t phiWave, Double_t thetaWave, Double_t rWave){

  Double_t phi1 = TMath::DegToRad()*phiArrayDeg[pol].at(ant1);
  Double_t x1 = rArray[pol].at(ant1)*TMath::Cos(phi1);
  Double_t y1 = rArray[pol].at(ant1)*TMath::Sin(phi1);
  Double_t z1 = zArray[pol].at(ant1);

  Double_t phi2 = TMath::DegToRad()*phiArrayDeg[pol].at(ant2);
  Double_t x2 = rArray[pol].at(ant2)*TMath::Cos(phi2);
  Double_t y2 = rArray[pol].at(ant2)*TMath::Sin(phi2);
  Double_t z2 = zArray[pol].at(ant2);

  Double_t x0 = rWave*TMath::Cos(phiWave)*TMath::Cos(thetaWave);
  Double_t y0 = rWave*TMath::Sin(phiWave)*TMath::Cos(thetaWave);
  Double_t z0 = rWave*TMath::Sin(thetaWave);

  Double_t part1 = TMath::Sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0));
  Double_t part2 = TMath::Sqrt((x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0));

  Double_t tdiff = 1e9*(part2 - part1)/SPEED_OF_LIGHT; // 1e9 for seconds to nanoseconds 

  tdiff /= nominalSamplingDeltaT;
  return TMath::Nint(tdiff);
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
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Fills in an array of cached deltaTs between antenna pairs as a function of arrival direction
 */
void CrossCorrelator::fillDeltaTLookup(){

  Double_t phi0 = getBin0PhiDeg();
  for(Int_t polInd=0; polInd<NUM_POL; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(Int_t phiSector = 0; phiSector<NUM_PHI; phiSector++){
      for(Int_t phiInd = 0; phiInd < NUM_BINS_PHI; phiInd++){
  	Int_t phiBin = phiSector*NUM_BINS_PHI + phiInd;

  	Double_t phiDeg = phi0 + phiBin*Double_t(PHI_RANGE)/NUM_BINS_PHI;
	
  	Double_t phiWave = TMath::DegToRad()*phiDeg;
  	for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA; thetaBin++){
  	  Double_t thetaDeg = THETA_RANGE*((Double_t)thetaBin/NUM_BINS_THETA - 0.5);
  	  Double_t thetaWave = TMath::DegToRad()*thetaDeg;
  	  for(Int_t combo=0; combo<numCombos; combo++){
  	    Int_t ant1 = comboToAnt1s.at(combo);
  	    Int_t ant2 = comboToAnt2s.at(combo);
  	    Double_t deltaT = getDeltaTExpected(pol, ant1, ant2, phiWave, thetaWave);
  	    // Int_t offset = TMath::Nint(deltaT/nominalSamplingDeltaT);
  	    // deltaTs[pol][phiBin][thetaBin][combo] = deltaT; //offset;
  	    // deltaTs[pol][combo][thetaBin][phiBin] = deltaT; //offset;
	    Int_t offsetLow = floor(deltaT/nominalSamplingDeltaT);
	    offsetLows[pol][combo][thetaBin][phiBin] = offsetLow;
	    Double_t dt1 = offsetLow*nominalSamplingDeltaT;
	    interpPreFactors[pol][combo][thetaBin][phiBin] = (deltaT - dt1)/nominalSamplingDeltaT;

	    // Here we account for the fact that we are now time ordering the correlations
	    offsetLows[pol][combo][thetaBin][phiBin]+=numSamples/2;
	  }
  	}
      }
    }
  }

  // const Double_t thetaBinSize = (Double_t(THETA_RANGE)/NUM_BINS_THETA);
  // for(Int_t thetaIndex=0; thetaIndex < NUM_BINS_THETA; thetaIndex++){
  //   Double_t thetaWaveDeg = (thetaIndex-NUM_BINS_THETA/2)*thetaBinSize;
  //   Double_t thetaWave = thetaWaveDeg*TMath::DegToRad();
  //   thetaWaves[thetaIndex] = thetaWave;
  //   tanThetaWaves[thetaIndex] = tan(thetaWave);
  //   cosThetaWaves[thetaIndex] = cos(thetaWave);
  // }
  // const Double_t phiBinSize = Double_t(PHI_RANGE)/NUM_BINS_PHI;
  // for(Int_t phiIndex=0; phiIndex < NUM_BINS_PHI*NUM_PHI; phiIndex++){
  //   Double_t phiDeg = phi0 + phiIndex*phiBinSize;
  //   Double_t phiWave = TMath::DegToRad()*phiDeg;

  //   phiWaveLookup[phiIndex] = phiWave;
  //   for(Int_t pol=0; pol<AnitaPol::kNotAPol; pol++){

  //     for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
  // 	cosPartLookup[phiIndex][pol][ant] = rArray[pol].at(ant)*cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol].at(ant));
  //     }
  //   }
  // }
  minThetaDegZoom = -78.5;
  minPhiDegZoom = -59.25;

  for(Int_t thetaIndex=0; thetaIndex < NUM_BINS_THETA_ZOOM_TOTAL; thetaIndex++){
    // Double_t thetaWaveDeg = (thetaIndex-NUM_BINS_THETA_ZOOM_TOTAL/2)*ZOOM_BIN_SIZE_THETA;
    Double_t thetaWaveDeg = minThetaDegZoom + thetaIndex*ZOOM_BIN_SIZE_THETA;
    Double_t thetaWave = thetaWaveDeg*TMath::DegToRad();
    zoomedThetaWaves[thetaIndex] = thetaWave;
    zoomedTanThetaWaves[thetaIndex] = tan(thetaWave);
    zoomedCosThetaWaves[thetaIndex] = cos(thetaWave);

    for(Int_t pol=0; pol<AnitaPol::kNotAPol; pol++){
      for(Int_t combo=0; combo < NUM_COMBOS; combo++){
	Int_t ant1 = comboToAnt1s.at(combo);
	Int_t ant2 = comboToAnt2s.at(combo);
	partBAsZoom[pol][combo][thetaIndex] = zoomedTanThetaWaves[thetaIndex]*(zArray[pol].at(ant2)-zArray[pol].at(ant1));
      }
    }
  }

  for(Int_t phiIndex=0; phiIndex < NUM_BINS_PHI_ZOOM_TOTAL; phiIndex++){
    // Double_t phiWave = TMath::DegToRad()*phiIndex*ZOOM_BIN_SIZE_PHI;
    Double_t phiDeg = minPhiDegZoom + phiIndex*ZOOM_BIN_SIZE_PHI;
    Double_t phiWave = phiDeg*TMath::DegToRad();
    zoomedPhiWaveLookup[phiIndex] = phiIndex*ZOOM_BIN_SIZE_PHI;

    for(Int_t pol=0; pol<AnitaPol::kNotAPol; pol++){

      for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
	// zoomedCosPartLookup[phiIndex][pol][ant] = rArray[pol].at(ant)*cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol].at(ant));
	zoomedCosPartLookup[pol][ant][phiIndex] = rArray[pol].at(ant)*cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol].at(ant));	
      }    

      
      for(Int_t combo=0; combo < NUM_COMBOS; combo++){
	Int_t ant1 = comboToAnt1s.at(combo);
	Int_t ant2 = comboToAnt2s.at(combo);

	Double_t offAxisDelay = getOffAxisDelay((AnitaPol::AnitaPol_t)pol, ant1, ant2, phiWave);
	offAxisDelays[pol][phiIndex][combo] = offAxisDelay;
	
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
inline Bool_t CrossCorrelator::useCombo(Int_t ant1, Int_t ant2, Int_t phiSector, Int_t deltaPhiSect){

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


  // See rather large rant above.
  Bool_t extraCondition = true;
  if(deltaPhiSect < 0){
    extraCondition = (phiSectorOfAnt1==phiSector || phiSectorOfAnt2==phiSector);
  }
  
  return (ant1InRange && ant2InRange && extraCondition);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates a vector of antenna combo indices and puts them in combosToUseGlobal or combosToUseTriggered.
 *
 * @param mapMode is the type of reconstruction being used.
 * @param l3TrigPattern is the l3TrigPattern(H) being used in the reconstruction.
 */
void CrossCorrelator::fillCombosToUseIfNeeded(mapMode_t mapMode, UShort_t l3TrigPattern){
  
  std::pair<UInt_t, Int_t> key(l3TrigPattern, kDeltaPhiSect);
  
  if(mapMode==kTriggered){
    std::map<std::pair<UInt_t, Int_t>,std::vector<Int_t> >::iterator it = combosToUseTriggered.find(key);    
    if(it==combosToUseTriggered.end()){
      combosToUseTriggered[key] = std::vector<Int_t>();
      for(Int_t phiSector = 0; phiSector<NUM_PHI; phiSector++){
	UInt_t doPhiSector = ((l3TrigPattern >> phiSector) & 1);
	if(doPhiSector){
	  for(Int_t combo=0; combo<numCombos; combo++){
	    Int_t ant1 = comboToAnt1s.at(combo);
	    Int_t ant2 = comboToAnt2s.at(combo);

	    if(useCombo(ant1, ant2, phiSector, kDeltaPhiSect)){	      
	      if(std::find(combosToUseTriggered[key].begin(),
			   combosToUseTriggered[key].end(),
			   combo) == combosToUseTriggered[key].end()){
		combosToUseTriggered[key].push_back(combo);
	      }
	    }
	  }
	}
      }
    }
  }
  else if(mapMode==kGlobal){
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
 * @brief Linearly interpolates between correlation points.
 *
 * @param pol is polarization of the antenna pair.
 * @param combo is the index of the antenna pair.
 * @param deltaT is the time to interpolate the correlation values at.
 * @returns the interpolated correlation value.
 */
inline Double_t CrossCorrelator::getInterpolatedCorrelationValue(AnitaPol::AnitaPol_t pol, Int_t combo, Double_t deltaT){

  Int_t offsetLow = floor(deltaT/nominalSamplingDeltaT);
  Int_t offsetHigh = ceil(deltaT/nominalSamplingDeltaT);

  Double_t dt1 = offsetLow*nominalSamplingDeltaT;
  // Double_t dt2 = offsetHigh*ptr->correlationDeltaT;
	  
  offsetLow = offsetLow < 0 ? offsetLow + numSamples : offsetLow;
  offsetHigh = offsetHigh < 0 ? offsetHigh + numSamples : offsetHigh;
	  
  Double_t c1 = crossCorrelations[pol][combo][offsetLow];
  Double_t c2 = crossCorrelations[pol][combo][offsetHigh];

  Double_t cInterp = (deltaT - dt1)*(c2 - c1)/(nominalSamplingDeltaT) + c1;

  return cInterp;
  
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates a blank histogram for cross correlations.
 *
 * @param name is the histogram name
 * @param title is the histogram title
 * @returns an empty historgram, to be filled as an interferometric map.
 */
TH2D* CrossCorrelator::makeBlankImage(TString name, TString title){

  Double_t phiMin = getBin0PhiDeg();
  Double_t phiMax = phiMin + DEGREES_IN_CIRCLE;
  Double_t thetaMin = -THETA_RANGE/2;
  Double_t thetaMax = THETA_RANGE/2;

  TH2D* hImage = new TH2D(name, title,
			  NUM_BINS_PHI*NUM_PHI, phiMin, phiMax,
			  NUM_BINS_THETA, thetaMin, thetaMax);
  hImage->GetXaxis()->SetTitle("Azimuth (Degrees)");
  hImage->GetYaxis()->SetTitle("Elevation (Degrees)");
  hImage->Sumw2();
  return hImage;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates a blank histogram with a small range in azimuth and elevation.
 *
 * @param name is the histogram name
 * @param title is the histogram title
 * @param zoomCenterPhiDeg is the centre of the histogram is azimuth (Degrees) relative to the ADU5 aft-fore.
 * @param zoomCenterThetaDeg is the centre of the histogram is elevation (Degrees).
 * @returns an empty historgram, to be filled as an interferometric map.
 */
TH2D* CrossCorrelator::makeBlankZoomedImage(TString name, TString title,
					    Double_t zoomCenterPhiDeg,
					    Double_t zoomCenterThetaDeg){

  zoomCenterPhiDeg = (TMath::Nint(zoomCenterPhiDeg/ZOOM_BIN_SIZE_PHI))*ZOOM_BIN_SIZE_PHI;
  zoomCenterThetaDeg = (TMath::Nint(zoomCenterThetaDeg/ZOOM_BIN_SIZE_THETA))*ZOOM_BIN_SIZE_THETA;
  
  Double_t phiMin = zoomCenterPhiDeg - PHI_RANGE_ZOOM/2;
  Double_t phiMax = zoomCenterPhiDeg + PHI_RANGE_ZOOM/2;
  
  Double_t thetaMin = zoomCenterThetaDeg - THETA_RANGE_ZOOM/2;
  Double_t thetaMax = zoomCenterThetaDeg + THETA_RANGE_ZOOM/2;

  TH2D* hImage = new TH2D(name, title,
			  NUM_BINS_PHI_ZOOM, phiMin, phiMax,
			  NUM_BINS_THETA_ZOOM, thetaMin, thetaMax);
  hImage->GetXaxis()->SetTitle("Azimuth (Degrees)");
  hImage->GetYaxis()->SetTitle("Elevation (Degrees)");
  hImage->Sumw2();
  
  return hImage;

}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates an interferometric map using spherical wave deltaTs and antenna pairs from all phi-sectors.
 *
 * @param pol is the polarization.
 * @param rWave is the source distance in meters.
 * @returns the interferometric map.
 */
TH2D* CrossCorrelator::makeGlobalSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave){
  Double_t imagePeak, peakPhiDeg, peakThetaDeg;

  return makeGlobalSphericalImage(pol, rWave, imagePeak, peakPhiDeg, peakThetaDeg);  

}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates an interferometric map using spherical wave deltaTs and antenna pairs from all phi-sectors.
 *
 * @param pol is the polarization.
 * @param rWave is the source distance in meters.
 * @param imagePeak is the maximum value.
 * @param peakPhiDeg is azimuth of the maximum value (Degrees) relative to the ADU5 aft-fore.
 * @param peakThetaDeg is the elevation of the maximum value (Degrees).
 * @returns the interferometric map.
 */
TH2D* CrossCorrelator::makeGlobalSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave,
						Double_t& imagePeak, Double_t& peakPhiDeg,
						Double_t& peakThetaDeg){

  return makeImageThreaded(pol, rWave, imagePeak, peakPhiDeg, peakThetaDeg, ALL_PHI_TRIGS,
  			   kGlobal, kZoomedOut);
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates an interferometric map using spherical wave deltaTs and antenna pairs from l3 triggered phi-sectors.
 *
 * @param pol is the polarization.
 * @param rWave is the source distance in meters.
 * @param l3TrigPattern is the L3 trigger pattern, determines what antenna pairs to use in the recontruction.
 * @returns the interferometric map.
 */
TH2D* CrossCorrelator::makeTriggeredSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave, UShort_t l3TrigPattern){
  Double_t imagePeak, peakPhiDeg, peakThetaDeg;
  return makeTriggeredSphericalImage(pol, rWave, imagePeak, peakPhiDeg, peakThetaDeg, l3TrigPattern);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates an interferometric map using spherical wave deltaTs and antenna pairs from l3 triggered phi-sectors.
 *
 * @param pol is the polarization.
 * @param rWave is the source distance in meters.
 * @param imagePeak is the maximum value.
 * @param peakPhiDeg is azimuth of the maximum value (Degrees) relative to the ADU5 aft-fore.
 * @param peakThetaDeg is the elevation of the maximum value (Degrees).
 * @param l3TrigPattern is the L3 trigger pattern, determines what antenna pairs to use in the recontruction.
 * @returns the interferometric map.
 */
TH2D* CrossCorrelator::makeTriggeredSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave,
						   Double_t& imagePeak, Double_t& peakPhiDeg,
						   Double_t& peakThetaDeg, UShort_t l3TrigPattern){

  return makeImageThreaded(pol, rWave, imagePeak, peakPhiDeg, peakThetaDeg,
  			   l3TrigPattern, kTriggered, kZoomedOut);
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

  return makeImageThreaded(pol, 0, imagePeak, peakPhiDeg, peakThetaDeg, ALL_PHI_TRIGS,
  			   kGlobal, kZoomedOut);
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

  return makeImageThreaded(pol, 0, imagePeak, peakPhiDeg, peakThetaDeg,
  			   l3TrigPattern, kTriggered, kZoomedOut);  
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

  return makeImageThreaded(pol, 0, imagePeak, peakPhiDeg, peakThetaDeg, 0,
  			   kTriggered, kZoomedIn, zoomCenterPhiDeg, zoomCenterThetaDeg);

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

  return makeImageThreaded(pol, 0, imagePeak, peakPhiDeg, peakThetaDeg, l3TrigPattern,
  			   kTriggered, kZoomedIn, zoomCenterPhiDeg, zoomCenterThetaDeg);

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
  return makeZoomedImage(pol, imagePeak, peakPhiDeg, peakThetaDeg,
  			 l3TrigPattern, zoomCenterPhiDeg, zoomCenterThetaDeg);

}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Makes the histogram name and title, based on the exact type of map requested.
 *
 * @param name is a reference variable to the histogram name.
 * @param title is a reference variable to the histogram title.
 * @param mapMode is the type of reconstruction being used.
 * @param zoomMode specifies whether or not the histogram is zoomed in around a particular phi/theta.
 * @param rWave is the distance of the spherical wave source. If rWave = 0, then assume a planeWave.
 * @param pol is the polarization.
 */
void CrossCorrelator::createImageNameAndTitle(TString& name, TString& title, mapMode_t mapMode,
					      zoomMode_t zoomMode, Double_t rWave, AnitaPol::AnitaPol_t pol){

  // Generate unique name and title from eventNumber, polarization, spherical/plane wave, triggered vs global.
  name = "h" + mapModeNames[mapMode];
  name += rWave == 0 ? "" : "Spherical";
  name += zoomModeNames[zoomMode];
  name += pol == AnitaPol::kVertical ? "ImageV" : "ImageH";
  name += TString::Format("%u", eventNumber[pol]);

  title = TString::Format("Event %u ", eventNumber[pol]);
  title += (pol == AnitaPol::kVertical ? "VPOL" : "HPOL");
  title += rWave == 0 ? "" : TString::Format(" Spherical r=%4.2lf", rWave);
  title += " " + mapModeNames[mapMode];
  if(zoomMode==kZoomedIn){
    title += " Zoomed In";
  }
  title += " Map";
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates a blank image and sets internal variables for threading.
 *
 * @param pol is the polarization.
 * @param rWave is the distance of the spherical wave source. If rWave = 0, then assume a planeWave.
 * @param l3TrigPattern is the L3 trig pattern, determines what antenna pairs to use in the recontruction.
 * @param mapMode is the type of reconstruction being used.
 * @param zoomMode specifies whether or not the histogram is zoomed in around a particular phi/theta.
 * @param zoomCenterPhiDeg is azimuth to center the zoomed image on (Degrees) relative to the ADU5 aft-fore.
 * @param zoomCenterThetaDeg is the elevation to center the zoomed image on (Degrees). 
 * @returns an empty historgram, to be filled as an interferometric map.
 */
TH2D* CrossCorrelator::prepareForImageMaking(AnitaPol::AnitaPol_t pol, Double_t rWave, UShort_t l3TrigPattern,
					     mapMode_t mapMode, zoomMode_t zoomMode,
					     Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg){


  TString name, title;
  createImageNameAndTitle(name, title, mapMode, zoomMode, rWave, pol);
  
  TH2D* hImage = NULL;
  if(zoomMode == kZoomedOut){
    hImage = makeBlankImage(name, title);
  }
  else if(zoomMode == kZoomedIn){
    hImage = makeBlankZoomedImage(name, title, zoomCenterPhiDeg, zoomCenterThetaDeg);

    if(l3TrigPattern==0){
      Int_t phiSectorOfPeak = -1;
      Double_t bestDeltaPhiOfPeakToAnt = 360;
      for(int ant=0; ant < NUM_SEAVEYS; ant++){
	Double_t phiOfAnt = phiArrayDeg[pol].at(ant);
	Double_t deltaPhiOfPeakToAnt = TMath::Abs(RootTools::getDeltaAngleDeg(phiOfAnt, zoomCenterPhiDeg));
	if(deltaPhiOfPeakToAnt < bestDeltaPhiOfPeakToAnt){
	  bestDeltaPhiOfPeakToAnt = deltaPhiOfPeakToAnt;
	  phiSectorOfPeak = (ant % NUM_PHI);
	}
      }
      l3TrigPattern = (1 << phiSectorOfPeak);
    }
  }

  threadImage = hImage;
  threadPol = pol;
  threadMapMode = mapMode;
  threadZoomMode = zoomMode;
  threadRWave = rWave;
  threadL3TrigPattern = l3TrigPattern;

  fillCombosToUseIfNeeded(mapMode, l3TrigPattern);
  std::pair<UInt_t, Int_t> key(l3TrigPattern, kDeltaPhiSect);
  threadCombosToUse = combosToUseTriggered[key];
  
  if(zoomMode == kZoomedIn){
    doUpsampledCrossCorrelationsThreaded(pol, l3TrigPattern);
  }
  return hImage;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Wrapper function which launches the threaded functions, which fill the interferometric maps.
 *
 * @param pol is the polarization.
 * @param rWave is the distance of the spherical wave source. If rWave = 0, then assume a planeWave.
 * @param imagePeak is the maximum value.
 * @param peakPhiDeg is azimuth of the maximum value (Degrees) relative to the ADU5 aft-fore.
 * @param peakThetaDeg is the elevation of the maximum value (Degrees).
 * @param l3TrigPattern is the L3 trig pattern, determines what antenna pairs to use in the recontruction.
 * @param mapMode is the type of reconstruction being used.
 * @param zoomMode specifies whether or not the histogram is zoomed in around a particular phi/theta.
 * @param zoomCenterPhiDeg is azimuth to center the zoomed image on (Degrees) relative to the ADU5 aft-fore.
 * @param zoomCenterThetaDeg is the elevation to center the zoomed image on (Degrees). 
 * @returns the interferometric map.
 *
 * This function is also responsible for merging the peak finding results of each of the threads.
 */
TH2D* CrossCorrelator::makeImageThreaded(AnitaPol::AnitaPol_t pol, Double_t rWave, Double_t& imagePeak,
					 Double_t& peakPhiDeg, Double_t& peakThetaDeg, UShort_t l3TrigPattern,
					 mapMode_t mapMode, zoomMode_t zoomMode, Double_t zoomCenterPhiDeg,
					 Double_t zoomCenterThetaDeg){

  TH2D* hImage = prepareForImageMaking(pol, rWave, l3TrigPattern,
				       mapMode, zoomMode,  zoomCenterPhiDeg,
				       zoomCenterThetaDeg);

  // LAUNCH THREADS HERE
  for(long threadInd=0; threadInd<NUM_THREADS; threadInd++){
    TString name = TString::Format("threadMap%ld", threadInd);
    mapThreads.push_back(new TThread(name.Data(),
				     CrossCorrelator::makeSomeOfImageThreaded,
				     (void*)&threadArgsVec.at(threadInd))
			 );
  }

  for(long threadInd=0; threadInd<NUM_THREADS; threadInd++){
    mapThreads.at(threadInd)->Run();
  }

  for(long threadInd=0; threadInd<NUM_THREADS; threadInd++){
    mapThreads.at(threadInd)->Join();
  }

  for(long threadInd=0; threadInd<NUM_THREADS; threadInd++){
    delete mapThreads.at(threadInd);
  }
  mapThreads.clear();
  
  // Combine peak search results from each thread.
  imagePeak = -DBL_MAX;
  for(Long_t threadInd=0; threadInd<NUM_THREADS; threadInd++){
    if(threadImagePeak[threadInd] > imagePeak){
      imagePeak = threadImagePeak[threadInd];
      peakPhiDeg = threadPeakPhiDeg[threadInd];
      peakThetaDeg = threadPeakThetaDeg[threadInd];
    }
  }

  // Remove global pointer so we can't accidentally delete it
  threadImage = NULL;
  
  // Now it's yours
  return hImage;
}






//---------------------------------------------------------------------------------------------------------
/**
 * @brief Static member function which fills the interferometric maps.
 *
 * @param voidPtrArgs contains a pointer to a CrossCorrelator::threadArgs struct
 *
 * This function contains the meat and bones of this class.
 */
void* CrossCorrelator::makeSomeOfImageThreaded(void* voidPtrArgs){
  // Disgusting hacks to get ROOT threading to compile inside a class.
  CrossCorrelator::threadArgs* args = (CrossCorrelator::threadArgs*) voidPtrArgs;
  Long_t threadInd = args->threadInd;
  CrossCorrelator* ptr = args->ptr;
  AnitaPol::AnitaPol_t pol = ptr->threadPol;
  mapMode_t mapMode = ptr->threadMapMode;
  zoomMode_t zoomMode = ptr->threadZoomMode;
  Double_t rWave = ptr->threadRWave;
  TH2D* hImage = ptr->threadImage;

  // Int_t numThetaBinsPerThread = hImage->GetNbinsY()/NUM_THREADS;
  // const Int_t startThetaBin = threadInd*numThetaBinsPerThread;
  // const Int_t endThetaBin = startThetaBin+numThetaBinsPerThread;
  // const Int_t startPhiBin = 0;
  // const Int_t endPhiBin = hImage->GetNbinsX();

  Int_t numPhiBinsPerThread = hImage->GetNbinsX()/NUM_THREADS;
  const Int_t startThetaBin = 0;
  const Int_t endThetaBin = hImage->GetNbinsY();
  const Int_t startPhiBin = threadInd*numPhiBinsPerThread;
  const Int_t endPhiBin = startPhiBin+numPhiBinsPerThread;
  

  ptr->threadImagePeak[threadInd] = -DBL_MAX;
  Int_t peakPhiBin = -1;
  Int_t peakThetaBin = -1;
  std::vector<Int_t>* combosToUse = &ptr->threadCombosToUse;
  
  if(zoomMode==kZoomedOut && rWave==0){    

    // zero internal map
    for(Int_t thetaBin = startThetaBin; thetaBin < endThetaBin; thetaBin++){
      for(Int_t phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){
	ptr->coarseMap[thetaBin][phiBin] = 0;
      }
    }

    for(UInt_t comboInd=0; comboInd<combosToUse->size(); comboInd++){
      Int_t combo = combosToUse->at(comboInd);
      if(ptr->kOnlyThisCombo >= 0 && combo!=ptr->kOnlyThisCombo){
	continue;
      }

      // for(Int_t thetaBin = 0; thetaBin < hImage->GetNbinsY(); thetaBin++){            
      // 	for(Int_t phiBin = startPhiBin; phiBin < startPhiBin+numPhiBinsThread; phiBin++){
      // for(Int_t thetaBin = startThetaBin; thetaBin < startThetaBin+numThetaBinsPerThread; thetaBin++){	
      // 	for(Int_t phiBin = 0; phiBin < hImage->GetNbinsX(); phiBin++){
      for(Int_t thetaBin = startThetaBin; thetaBin < endThetaBin; thetaBin++){
	for(Int_t phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){	  
	  // Int_t phiSector = phiBin/NUM_BINS_PHI;
	  // if(mapMode==kGlobal){
	  //   combosToUse = &ptr->combosToUseGlobal[phiSector];
	  // }
	  // Double_t deltaT = ptr->deltaTs[pol][phiBin][thetaBin][combo];
	  // Double_t deltaT = ptr->deltaTs[pol][combo][thetaBin][phiBin];	  
	  // Int_t offsetLow  = floor(deltaT/ptr->nominalSamplingDeltaT);
	  // Double_t dt1 = offsetLow*ptr->nominalSamplingDeltaT;
	  // offsetLow += offset;
	  // Int_t offsetHigh = offsetLow+1;

	  Int_t offsetLow = ptr->offsetLows[pol][combo][thetaBin][phiBin];	  
	  Double_t c1 = ptr->crossCorrelations[pol][combo][offsetLow];
	  Double_t c2 = ptr->crossCorrelations[pol][combo][offsetLow+1];

	  // Double_t cInterp = (deltaT - dt1)*(c2 - c1)/(ptr->nominalSamplingDeltaT) + c1;

	  Double_t cInterp = ptr->interpPreFactors[pol][combo][thetaBin][phiBin]*(c2 - c1) + c1;
	  
	  // Double_t cInterp = ptr->getInterpolatedCorrelationValue(pol, combo, deltaT);
	  ptr->coarseMap[thetaBin][phiBin] += cInterp;
	  // ptr->coarseMap[thetaBin][phiBin] += ptr->crossCorrelations[pol][combo][offsetLow];	  
	}
      }
    }


    // for(Int_t phiBin = startPhiBin; phiBin < startPhiBin+numPhiBinsThread; phiBin++){
    //   for(Int_t thetaBin = 0; thetaBin < hImage->GetNbinsY(); thetaBin++){
    // for(Int_t thetaBin = startThetaBin; thetaBin < startThetaBin+numThetaBinsPerThread; thetaBin++){		
    //   for(Int_t phiBin = 0; phiBin < hImage->GetNbinsX(); phiBin++){
    for(Int_t thetaBin = startThetaBin; thetaBin < endThetaBin; thetaBin++){
      for(Int_t phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){	          
	if(combosToUse->size()>0 && ptr->kOnlyThisCombo < 0){
	  ptr->coarseMap[thetaBin][phiBin] /= combosToUse->size();
	}
	hImage->SetBinContent(phiBin+1, thetaBin+1, ptr->coarseMap[thetaBin][phiBin]);
	if(ptr->coarseMap[thetaBin][phiBin] > ptr->threadImagePeak[threadInd]){
	  ptr->threadImagePeak[threadInd] = ptr->coarseMap[thetaBin][phiBin];
	  peakPhiBin = phiBin;
	  peakThetaBin = thetaBin;
	}
      }
    }
  }

  // If we are in zoomed in & plane wave mode then calculate
  // Now using Fast function, this has some intermediate caching of the deltaT calc.
  else if(zoomMode==kZoomedIn && rWave==0){
    
    Int_t phiZoomBase = TMath::Nint((hImage->GetXaxis()->GetBinLowEdge(1) - ptr->minPhiDegZoom)/ZOOM_BIN_SIZE_PHI);
    Int_t thetaZoomBase = TMath::Nint((hImage->GetYaxis()->GetBinLowEdge(1) - ptr->minThetaDegZoom)/ZOOM_BIN_SIZE_THETA);

    // for(Int_t thetaBin = startThetaBin; thetaBin < startThetaBin+numThetaBinsPerThread; thetaBin++){		
    //   for(Int_t phiBin = 0; phiBin < hImage->GetNbinsX(); phiBin++){	  	    
    // for(Int_t thetaBin = 0; thetaBin < hImage->GetNbinsY(); thetaBin++){
    //   for(Int_t phiBin = startPhiBin; phiBin < startPhiBin+numPhiBinsThread; phiBin++){
    for(Int_t thetaBin = startThetaBin; thetaBin < endThetaBin; thetaBin++){
      for(Int_t phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){	          	
	ptr->fineMap[thetaBin][phiBin]=0;
      }
    }

    // const Int_t offset = ptr->numSamplesUpsampled/2;
    // for(UInt_t comboInd=0; comboInd<combosToUse->size(); comboInd++){
    //   Int_t combo = combosToUse->at(comboInd);
    //   if(ptr->kOnlyThisCombo >= 0 && combo!=ptr->kOnlyThisCombo){
    // 	continue;
    //   }
    //   Int_t ant1 = ptr->comboToAnt1s.at(combo);
    //   Int_t ant2 = ptr->comboToAnt2s.at(combo);
      
    //   for(Int_t thetaBin = startThetaBin; thetaBin < startThetaBin+numThetaBinsPerThread; thetaBin++){		
    // 	Int_t zoomThetaInd = thetaZoomBase + thetaBin;

    // 	Double_t partA = ptr->zArray[pol].at(ant1)*ptr->zoomedTanThetaWaves[zoomThetaInd];
    // 	Double_t partB = ptr->zArray[pol].at(ant2)*ptr->zoomedTanThetaWaves[zoomThetaInd];
    // 	Double_t dtFactor = ptr->zoomedCosThetaWaves[zoomThetaInd]/SPEED_OF_LIGHT_NS;
	
    // 	for(Int_t phiBin = 0; phiBin < hImage->GetNbinsX(); phiBin++){
    // 	  Int_t zoomPhiInd = phiZoomBase + phiBin;
    // 	  Double_t part1 = partA - ptr->zoomedCosPartLookup[pol][ant1][zoomPhiInd];
    // 	  Double_t part2 = partB - ptr->zoomedCosPartLookup[pol][ant2][zoomPhiInd];
    // 	  ptr->interpPreFactorsZoom[combo][thetaBin][phiBin] = dtFactor*(part2 - part1);
    // 	  ptr->offsetLowsZoom[combo][thetaBin][phiBin] = floor(ptr->interpPreFactorsZoom[combo][thetaBin][phiBin]/ptr->correlationDeltaT);
    // 	  // Double_t dt1 = offsetLow*ptr->correlationDeltaT;
    // 	  ptr->interpPreFactorsZoom[combo][thetaBin][phiBin] -= ptr->offsetLowsZoom[combo][thetaBin][phiBin]*ptr->correlationDeltaT;
    // 	  ptr->interpPreFactorsZoom[combo][thetaBin][phiBin] /= ptr->correlationDeltaT;
	  
    // 	  ptr->offsetLowsZoom[combo][thetaBin][phiBin] += offset;
    // 	}
    //   }
    // }

    // for(UInt_t comboInd=0; comboInd<combosToUse->size(); comboInd++){
    //   Int_t combo = combosToUse->at(comboInd);
    //   if(ptr->kOnlyThisCombo >= 0 && combo!=ptr->kOnlyThisCombo){
    // 	continue;
    //   }
    //   for(Int_t thetaBin = startThetaBin; thetaBin < startThetaBin+numThetaBinsPerThread; thetaBin++){		
    // 	for(Int_t phiBin = 0; phiBin < hImage->GetNbinsX(); phiBin++){	  	    

    // 	  Int_t offsetLow = ptr->offsetLowsZoom[combo][thetaBin][phiBin];
    // 	  Double_t c1 = ptr->crossCorrelationsUpsampled[pol][combo][offsetLow];
    // 	  Double_t c2 = ptr->crossCorrelationsUpsampled[pol][combo][offsetLow+1];

    // 	  // Double_t cInterp = (deltaT - dt1)*(c2 - c1)/(ptr->correlationDeltaT) + c1;
    // 	  Double_t cInterp = ptr->interpPreFactorsZoom[combo][thetaBin][phiBin]*(c2 - c1) + c1;

    // 	  ptr->fineMap[thetaBin][phiBin] += cInterp;
    // 	}
    //   }
    // }

    const Int_t offset = ptr->numSamplesUpsampled/2;
    for(UInt_t comboInd=0; comboInd<combosToUse->size(); comboInd++){
      Int_t combo = combosToUse->at(comboInd);
      if(ptr->kOnlyThisCombo >= 0 && combo!=ptr->kOnlyThisCombo){
    	continue;
      }
      // Int_t ant1 = ptr->comboToAnt1s.at(combo);
      // Int_t ant2 = ptr->comboToAnt2s.at(combo);
      
      for(Int_t thetaBin = startThetaBin; thetaBin < endThetaBin; thetaBin++){
	// for(Int_t thetaBin = startThetaBin; thetaBin < startThetaBin+numThetaBinsPerThread; thetaBin++){
	// for(Int_t thetaBin = 0; thetaBin < hImage->GetNbinsY(); thetaBin++){
    	Int_t zoomThetaInd = thetaZoomBase + thetaBin;

    	// Double_t partA = ptr->zArray[pol].at(ant1)*ptr->zoomedTanThetaWaves[zoomThetaInd];
    	// Double_t partB = ptr->zArray[pol].at(ant2)*ptr->zoomedTanThetaWaves[zoomThetaInd];
    	// Double_t partA = ptr->zArray[pol].at(ant1)*ptr->zoomedTanThetaWaves[zoomThetaInd];
	// Double_t partBA = partB - partA;	
    	// Double_t partBA = ptr->zoomedTanThetaWaves[zoomThetaInd]*(ptr->zArray[pol].at(ant2)-ptr->zArray[pol].at(ant1));
    	Double_t partBA = ptr->partBAsZoom[pol][combo][zoomThetaInd];
    	Double_t dtFactor = ptr->zoomedCosThetaWaves[zoomThetaInd]/SPEED_OF_LIGHT_NS;
	
    	// for(Int_t phiBin = startPhiBin; phiBin < startPhiBin+numPhiBinsThread; phiBin++){
    	// for(Int_t phiBin = 0; phiBin < hImage->GetNbinsX(); phiBin++){
	for(Int_t phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){
    	  Int_t zoomPhiInd = phiZoomBase + phiBin;
    	  // Double_t part1 = partA - ptr->zoomedCosPartLookup[pol][ant1][zoomPhiInd];
    	  // Double_t part2 = partB - ptr->zoomedCosPartLookup[pol][ant2][zoomPhiInd];
    	  // Double_t deltaT = dtFactor*(part2 - part1);
    	  // Double_t deltaT = dtFactor*(partB - ptr->zoomedCosPartLookup[pol][ant2][zoomPhiInd] - partA + ptr->zoomedCosPartLookup[pol][ant1][zoomPhiInd]);
    	  // Double_t deltaT = dtFactor*(partB - partA - ptr->zoomedCosPartLookup[pol][ant2][zoomPhiInd] + ptr->zoomedCosPartLookup[pol][ant1][zoomPhiInd]);
    	  // Double_t deltaT = dtFactor*(partBA - ptr->zoomedCosPartLookup[pol][ant2][zoomPhiInd] + ptr->zoomedCosPartLookup[pol][ant1][zoomPhiInd]);

	  // Double_t part21 = ptr->zoomedCosPartLookup[pol][ant2][zoomPhiInd] - ptr->zoomedCosPartLookup[pol][ant1][zoomPhiInd];
    	  // Double_t deltaT = dtFactor*(partBA - part21);

    	  Double_t deltaT = dtFactor*(partBA - ptr->part21sZoom[pol][combo][zoomPhiInd]);  
	  
    	  Int_t offsetLow = floor(deltaT/ptr->correlationDeltaT);
    	  // Double_t dt1 = offsetLow*ptr->correlationDeltaT;
    	  deltaT -= offsetLow*ptr->correlationDeltaT;
    	  deltaT /= ptr->correlationDeltaT;
    	  offsetLow += offset;
    	  Double_t c1 = ptr->crossCorrelationsUpsampled[pol][combo][offsetLow];
    	  Double_t c2 = ptr->crossCorrelationsUpsampled[pol][combo][offsetLow+1];

    	  // Double_t cInterp = (deltaT - dt1)*(c2 - c1)/(ptr->correlationDeltaT) + c1;
    	  Double_t cInterp = deltaT*(c2 - c1) + c1;

    	  ptr->fineMap[thetaBin][phiBin] += cInterp;
    	}
      }
    }    

    // for(Int_t thetaBin = startThetaBin; thetaBin < startThetaBin+numThetaBinsPerThread; thetaBin++){		
    //   for(Int_t phiBin = 0; phiBin < hImage->GetNbinsX(); phiBin++){	  	    
    // for(Int_t thetaBin = 0; thetaBin < hImage->GetNbinsY(); thetaBin++){    
    //   for(Int_t phiBin = startPhiBin; phiBin < startPhiBin+numPhiBinsThread; phiBin++){
    for(Int_t thetaBin = startThetaBin; thetaBin < endThetaBin; thetaBin++){
      for(Int_t phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){	          		
	if(combosToUse->size()>0 && ptr->kOnlyThisCombo < 0){
	  ptr->fineMap[thetaBin][phiBin] /= combosToUse->size();
	}

	if(ptr->fineMap[thetaBin][phiBin] > ptr->threadImagePeak[threadInd]){
	  ptr->threadImagePeak[threadInd] = ptr->fineMap[thetaBin][phiBin];
	  peakPhiBin = phiBin;
	  peakThetaBin = thetaBin;
	}

	hImage->SetBinContent(phiBin + 1, thetaBin + 1, ptr->fineMap[thetaBin][phiBin]);
	// hImage->SetBinContent(phiBin + 1, thetaBin + 1, threadInd+1);	
      }
    }
  }

	
  // If we are in zoomed out & spherical wave mode then calculate
  else{
    // for(Int_t phiBin = startPhiBin; phiBin < startPhiBin+numPhiBinsThread; phiBin++){
    //   Int_t phiSector = phiBin/NUM_BINS_PHI;
    //   if(mapMode==kGlobal){
    // 	combosToUse = &ptr->combosToUseGlobal[phiSector];
    //   }
    
    //   Double_t phiWave = hImage->GetXaxis()->GetBinLowEdge(phiBin+1)*TMath::DegToRad();

    //   for(Int_t thetaBin = 0; thetaBin < hImage->GetNbinsY(); thetaBin++){
    // 	Double_t thetaWave = hImage->GetYaxis()->GetBinLowEdge(thetaBin+1)*TMath::DegToRad();

    // 	Double_t correlations = 0;
    // 	for(UInt_t comboInd=0; comboInd<combosToUse->size(); comboInd++){
    // 	  Int_t combo = combosToUse->at(comboInd);
    // 	  if(ptr->kOnlyThisCombo >= 0 && combo!=ptr->kOnlyThisCombo){
    // 	    continue;
    // 	  }
	  
    // 	  Int_t offset = 0;
    
    // 	  Int_t ant1 = ptr->comboToAnt1s.at(combo);
    // 	  Int_t ant2 = ptr->comboToAnt2s.at(combo);
    // 	  offset = ptr->getDeltaTExpectedSpherical(pol, ant1, ant2, phiWave, thetaWave, rWave);
    // 	  offset = offset < 0 ? offset + ptr->numSamples : offset;
    // 	  correlations += ptr->crossCorrelations[pol][combo][offset];
    // 	}
    // 	if(combosToUse->size()>0 && ptr->kOnlyThisCombo < 0){
    // 	  correlations /= combosToUse->size();
    // 	}
    // 	hImage->SetBinContent(phiBin + 1, thetaBin + 1, correlations);
    // 	if(correlations > ptr->threadImagePeak[threadInd]){
    // 	  ptr->threadImagePeak[threadInd] = correlations;
    // 	  peakPhiBin = phiBin;
    // 	  peakThetaBin = thetaBin;
    // 	}	
    //   }
    // }
  }

  ptr->threadPeakPhiDeg[threadInd] = hImage->GetXaxis()->GetBinLowEdge(peakPhiBin+1);
  ptr->threadPeakThetaDeg[threadInd] = hImage->GetYaxis()->GetBinLowEdge(peakThetaBin+1);

  return 0;
  
}








//---------------------------------------------------------------------------------------------------------
/**
 * @brief Makes a TGraph which indicates the L3 triggered regions of an interferometric map.
 *
 * @param name is the TGraph name.
 * @param l3TrigPattern is the L3 trig pattern, determines what portion of the map is highlighted.
 * @param fillStyle specifies the ROOT fill style, should be something "translucent".
 * @param col is the color to fill the TGraph in.
 * @returns the TGraph that highlights the L3 triggered phi-sectors
 *
 * Several people have indicated that they don't like this being in the MagicDisplay. 
 * So I will probably turn it off.
 */
TGraph* CrossCorrelator::makeTrigPatternGraph(TString name, UShort_t l3TrigPattern, Color_t col, Int_t fillStyle){

  // Something pretty for MagicDisplay integration.
  Double_t phiMin = getBin0PhiDeg();
  Double_t thetaMin = -THETA_RANGE/2;
  Double_t thetaMax = THETA_RANGE/2;
  
  TGraph* gr = new TGraph();
  gr->SetName(name);

  Int_t numPoints = gr->GetN();
  gr->SetPoint(numPoints, phiMin - PHI_RANGE, thetaMin);
  numPoints++;

  // Loop over extended range -1 to NUM_PHI+1
  for(Int_t phiSect=-1; phiSect<NUM_PHI+1; phiSect++){

    Int_t phiSectBit = phiSect;
    if(phiSectBit< 0){
      phiSectBit += NUM_PHI;
    }
    else if(phiSect >= NUM_PHI){
      phiSectBit -= NUM_PHI;
    }
    
    if(RootTools::getBit(phiSectBit, l3TrigPattern)){
      gr->SetPoint(numPoints, phiMin+PHI_RANGE*phiSect, thetaMax);
      numPoints++;
      gr->SetPoint(numPoints, phiMin+PHI_RANGE*(phiSect+1), thetaMax);
      numPoints++;
    }
    else{
      gr->SetPoint(numPoints, phiMin+PHI_RANGE*phiSect, thetaMin);
      numPoints++;
      gr->SetPoint(numPoints, phiMin+PHI_RANGE*(phiSect+1), thetaMin);
      numPoints++;
    }
  }
  gr->SetPoint(numPoints, phiMin+PHI_RANGE*(NUM_PHI+1), thetaMin);
  numPoints++;

  gr->SetLineColor(col);
  gr->SetFillColor(col);
  gr->SetFillStyle(fillStyle);
  return gr;
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
 * @brief Creates 2D hist showing the interpolated correlation values at used at a interferometric map position
 *
 * @param pol is the polarization
 * @param l3TrigPattern is the L3 trigger pattern used to figure out the combinatorics
 * @param phiDeg is the phi direction (Degrees) relative to the ADU5 aft-fore.
 * @param thetaDeg is the theta direction (Degrees).
 * @returns the correlation summary histogram for the given direction.
 */
TH2D* CrossCorrelator::makeCorrelationSummaryHistogram(AnitaPol::AnitaPol_t pol, UShort_t l3TrigPattern,
						       Double_t phiDeg, Double_t thetaDeg)
{


  // Get combinatorics
  fillCombosToUseIfNeeded(CrossCorrelator::kTriggered, l3TrigPattern);
  std::pair<UInt_t, Int_t> key(l3TrigPattern, kDeltaPhiSect);
  const std::vector<Int_t>& combos = combosToUseTriggered[key];
    

  // Try to center the correlation summary histogram on the l3Triggered phi-sectors  
  Int_t firstPhi = -1;
  for(Int_t phi=0; phi<NUM_PHI; phi++){
    Int_t trig = RootTools::getBit(phi, l3TrigPattern);
    if(trig > 0){
      firstPhi = phi;
      break;
    }
  }
  Int_t binShift = (NUM_PHI/2 - firstPhi)*NUM_RING;
  binShift = binShift < 0 ? binShift + NUM_SEAVEYS : binShift;


  // Book histogram and set x, y axis labels to antenna names
  TString name = "hCorrelationSummary";
  TString polLetter = pol == AnitaPol::kHorizontal ? "H" : "V";  
  name += polLetter;
  name += TString::Format("%u", eventNumber[pol]);

  phiDeg = phiDeg < 0 ? phiDeg + DEGREES_IN_CIRCLE : phiDeg;
  phiDeg = phiDeg >= 0 ? phiDeg - DEGREES_IN_CIRCLE : phiDeg;
  TString title = TString::Format("Correlation values at #theta=%4.2lf, #phi=%4.2lf for event %u ",
				  thetaDeg, phiDeg, eventNumber[pol]);

  title += pol==AnitaPol::kHorizontal ? "HPOL" : "VPOL";
  title += "; ; ; Correlation coefficient";
  TH2D* h = new TH2D(name, title, NUM_SEAVEYS, 0, NUM_SEAVEYS, NUM_SEAVEYS, 0, NUM_SEAVEYS);
  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    TString lab = TString::Format("%d", 1 + (ant/NUM_RING));
    if(ant%NUM_RING==0){
      lab += "T"+polLetter;
    }
    else if((ant%NUM_RING)==1){
      lab += "M"+polLetter;
    }
    else{
      lab += "B"+polLetter;
    }
    int binInd = ant + binShift;
    binInd = binInd >= NUM_SEAVEYS ? binInd - NUM_SEAVEYS : binInd;    
    h->GetXaxis()->SetBinLabel(binInd+1, lab.Data());
    h->GetYaxis()->SetBinLabel(binInd+1, lab.Data());
  }

  // Convert location to radians
  Double_t thetaWave = thetaDeg*TMath::DegToRad();
  Double_t phiWave = phiDeg*TMath::DegToRad();

  // Loop through antenna pairs, get deltaT and correlations and fill
  for(UInt_t comboInd=0; comboInd < combos.size(); comboInd++){
    const Int_t combo = combos.at(comboInd);      

    Int_t ant1 = comboToAnt1s.at(combo);
    Int_t ant2 = comboToAnt2s.at(combo);

    Double_t dt = getDeltaTExpected(pol, ant1, ant2, phiWave, thetaWave);
    Double_t cInterp = getInterpolatedUpsampledCorrelationValue(pol, combo, dt);

    Int_t phi1 = ant1%NUM_PHI;
    Int_t phi2 = ant2%NUM_PHI;

    Int_t ring1 = ant1/NUM_PHI;
    Int_t ring2 = ant2/NUM_PHI;

    Int_t theBinX = phi1*NUM_RING + ring1 + 1;
    Int_t theBinY = phi2*NUM_RING + ring2 + 1;

    theBinX += binShift;
    theBinY += binShift;

    theBinX = theBinX >= NUM_SEAVEYS ? theBinX - NUM_SEAVEYS : theBinX;
    theBinY = theBinY >= NUM_SEAVEYS ? theBinY - NUM_SEAVEYS : theBinY;      
      
    h->SetBinContent(theBinX, theBinY, cInterp);
    h->SetBinContent(theBinY, theBinX, cInterp);      
  }
  
  return h;  
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates 2D hist showing the times of interpolated correlation maxima for many antenna pairs.
 *
 * @param pol is the polarization
 * @param l3TrigPattern is the L3 trigger pattern used to figure out the combinatorics
 * @param phiDeg is the phi direction (Degrees) relative to the ADU5 aft-fore.
 * @param thetaDeg is the theta direction (Degrees).
 * @param corThresh is a threshold value, that the maximum interpolated correlation needs to be above.
 * @returns the correlation summary histogram for the given direction.
 */

TH2D* CrossCorrelator::makeDeltaTSummaryHistogram(AnitaPol::AnitaPol_t pol, UShort_t l3TrigPattern,
						  Double_t phiDeg, Double_t thetaDeg, Double_t corThresh)
{
  /* 
     This one should compare the time of the maximum vs. the time expected at at phiDeg/thetaDeg 
     Finds the difference between the nearest local maximum in correlation vs. the expected time offset given by the arrival direction thetaDeg, phiDeg.
  */
  

  // Get combinatorics
  fillCombosToUseIfNeeded(CrossCorrelator::kTriggered, l3TrigPattern);
  std::pair<UInt_t, Int_t> key(l3TrigPattern, kDeltaPhiSect);
  const std::vector<Int_t>& combos = combosToUseTriggered[key];



  // Try to center the correlation summary histogram on the l3Triggered phi-sectors  
  Int_t firstPhi = -1;
  for(Int_t phi=0; phi<NUM_PHI; phi++){
    Int_t trig = RootTools::getBit(phi, l3TrigPattern);
    if(trig > 0){
      firstPhi = phi;
      break;
    }
  }
  Int_t binShift = (NUM_PHI/2 - firstPhi)*NUM_RING;
  binShift = binShift < 0 ? binShift + NUM_SEAVEYS : binShift;


  // Book histogram and set x, y axis labels to antenna names
  TString name = "hCorrelationSummary";
  TString polLetter = pol == AnitaPol::kHorizontal ? "H" : "V"; 
  name += polLetter;
  name += TString::Format("%u", eventNumber[pol]);

  phiDeg = phiDeg < 0 ? phiDeg + DEGREES_IN_CIRCLE : phiDeg;
  phiDeg = phiDeg >= 0 ? phiDeg - DEGREES_IN_CIRCLE : phiDeg;
  TString title = TString::Format("Correlation values at #theta=%4.2lf, #phi=%4.2lf for event %u ",
				  thetaDeg, phiDeg, eventNumber[pol]);

  title += pol==AnitaPol::kHorizontal ? "HPOL" : "VPOL";
  title += "; ; ; Correlation coefficient";
  TH2D* h = new TH2D(name, title, NUM_SEAVEYS, 0, NUM_SEAVEYS, NUM_SEAVEYS, 0, NUM_SEAVEYS);
  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    TString lab = TString::Format("%d", 1 + (ant/NUM_RING));
    if(ant%NUM_RING==0){
      lab += "T"+polLetter;
    }
    else if((ant%NUM_RING)==1){
      lab += "M"+polLetter;
    }
    else{
      lab += "B"+polLetter;
    }
    int binInd = ant + binShift;
    binInd = binInd >= NUM_SEAVEYS ? binInd - NUM_SEAVEYS : binInd;    
    h->GetXaxis()->SetBinLabel(binInd+1, lab.Data());
    h->GetYaxis()->SetBinLabel(binInd+1, lab.Data());
  }

  // Convert location to radians
  Double_t thetaWave = thetaDeg*TMath::DegToRad();
  Double_t phiWave = phiDeg*TMath::DegToRad();

  // Loop through antenna pairs, get deltaT and correlations and fill
  for(UInt_t comboInd=0; comboInd < combos.size(); comboInd++){
    const Int_t combo = combos.at(comboInd);      

    Int_t ant1 = comboToAnt1s.at(combo);
    Int_t ant2 = comboToAnt2s.at(combo);

    Double_t dtExp = getDeltaTExpected(pol, ant1, ant2, phiWave, thetaWave);
    
    
    Double_t tToExp = 1e9;
    Double_t corVal = 0;
    Int_t corInd = 0;
    Double_t theDtMaxima;
    for(int samp=0; samp < numSamplesUpsampled; samp++){
      Int_t samp0 = samp - 1;
      Int_t samp1 = samp;
      Int_t samp2 = samp + 1;

      samp0 = samp0 < 0 ? samp0 + numSamplesUpsampled : samp0;
      samp2 = samp2 >= numSamplesUpsampled ? samp2 - numSamplesUpsampled : samp2;
      
      if(crossCorrelationsUpsampled[pol][combo][samp0] < crossCorrelationsUpsampled[pol][combo][samp1]
	 && crossCorrelationsUpsampled[pol][combo][samp2] < crossCorrelationsUpsampled[pol][combo][samp1]){

	Double_t dtMaxima = samp1 < numSamplesUpsampled/2 ? samp1*correlationDeltaT : (samp1-numSamplesUpsampled)*correlationDeltaT;

	if(TMath::Abs(dtMaxima - dtExp) < tToExp && corVal > corThresh){
	  tToExp = dtMaxima - dtExp;
	  corVal = crossCorrelationsUpsampled[pol][combo][samp1];
	  corInd = samp1;
	  theDtMaxima = dtMaxima;
	}
      }
    }    


    Int_t phi1 = ant1%NUM_PHI;
    Int_t phi2 = ant2%NUM_PHI;

    Int_t ring1 = ant1/NUM_PHI;
    Int_t ring2 = ant2/NUM_PHI;

    Int_t theBinX = phi1*NUM_RING + ring1 + 1;
    Int_t theBinY = phi2*NUM_RING + ring2 + 1;

    theBinX += binShift;
    theBinY += binShift;

    theBinX = theBinX >= NUM_SEAVEYS ? theBinX - NUM_SEAVEYS : theBinX;
    theBinY = theBinY >= NUM_SEAVEYS ? theBinY - NUM_SEAVEYS : theBinY;      
      
    h->SetBinContent(theBinX, theBinY, tToExp);
  }
  
  return h;  
}






//---------------------------------------------------------------------------------------------------------
/**
 * @brief Used to insert phase center geometry files from Linda.
 * 
 * @param pathToLindasFile is the relative path to Linda's file.
 * @param pol is the polarization of the channels under test.
 * @returns 0 if everything went without a problem, 1 if the text file could not be opened.
*/
Int_t CrossCorrelator::directlyInsertGeometry(TString pathToLindasFile, AnitaPol::AnitaPol_t pol){
  
  // Since I am simulataneously testing many of Linda's geometries on lots of different files
  // I need the help of a machine to check I'm testing the geometry I think I'm testing.
  
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->useKurtAnita3Numbers(0); // i.e. definitely use the numbers I am inserting.
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

    Double_t newR = geom->rPhaseCentreFromVerticalHornKurtAnita3[ant][pol] + dr;
    Double_t newPhi = geom->azPhaseCentreFromVerticalHornKurtAnita3[ant][pol] + dPhiRad;
    Double_t newZ = geom->zPhaseCentreFromVerticalHornKurtAnita3[ant][pol] + dz;
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

  }

  return 0;
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Creates a test event with a pulse from a known source location.
 * 
 * @param phiDegSource is the source phi direction (Degrees) relative to ADU5 aft-fore.
 * @param thetaDegSource is the source theta direction (Degrees).
 * @param rSource is the distance to the source of spherical waves in meters.
 *
 * Currently this just handles spherical plane waves. 
*/
void CrossCorrelator::correlateEventTest(Double_t phiDegSource, Double_t thetaDegSource, Double_t rSource){
  // Generates a set of delta function like waveforms, correlates them 

  // Delete any old waveforms (at start rather than end to leave in memory to be examined if need be)

  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){  
    deleteAllWaveforms((AnitaPol::AnitaPol_t)pol);
  }


  for(Int_t polInd = AnitaPol::kHorizontal; polInd < NUM_POL; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t)polInd;
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      std::vector<Double_t> newVs = std::vector<Double_t>(numSamples, 0);
      std::vector<Double_t> newTimes = std::vector<Double_t>(numSamples, 0);
      Int_t dt = numSamples/2 + getDeltaTExpectedSpherical(pol, 0, ant,
							    phiDegSource*TMath::DegToRad(), 
							    thetaDegSource*TMath::DegToRad(), 
							    rSource);
      for(Int_t samp=0; samp<numSamples; samp++){
	newTimes[samp] = nominalSamplingDeltaT*samp;
	if(samp==dt){
	  newVs[samp] = numSamples;
	}
      }
      grsResampled[pol][ant] = new TGraph(numSamples, &newTimes[0], &newVs[0]);
      RootTools::normalize(grsResampled[pol][ant]);
    }
    doFFTs((AnitaPol::AnitaPol_t)pol);
    doAllCrossCorrelationsThreaded((AnitaPol::AnitaPol_t)pol);
  }
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

  for(Int_t i=0; i<numSamps; i++){
    Int_t offset = (i - numSamps/2);
    offsets.at(i) = offset*graphDt;
    Int_t j = offset < 0 ? offset + numSamps : offset;
    corrs.at(i) = corrPtr[j];
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


  Double_t theDeltaT = nSamp == numSamples ? nominalSamplingDeltaT : correlationDeltaT;
  
  Double_t phiRad = phiDeg*TMath::DegToRad();
  Double_t thetaRad = thetaDeg*TMath::DegToRad();
  Int_t numAnts = 0;

  Int_t centerPhiSector = getPhiSectorOfAntennaClosestToPhiDeg(pol, phiDeg);

  const Int_t firstAnt = centerPhiSector;

  std::pair<Int_t, Int_t> key(nSamp, 0);
  // vArray is actually internal memory managed by FancyFFTs... don't delete this!!!
  Double_t* vArray = FancyFFTs::getRealArray(key);

  if(nSamp==numSamples){
    FancyFFTs::doInvFFT(nSamp, ffts[pol][firstAnt], false);
  }
  else{
    FancyFFTs::doInvFFT(nSamp, fftsPadded[pol][firstAnt], false);    
  }

  std::vector<Double_t> tArray(nSamp, 0);
  Double_t t0 = grsResampled[pol][firstAnt]->GetX()[0];
  for(Int_t samp=0; samp<nSamp; samp++){
    tArray.at(samp) = t0 + samp*theDeltaT;
    vArray[samp] *= interpRMS[pol][firstAnt]; // Undo the normalization.
  }

  // sum of rms of first few ns of each waveform
  Double_t rms = 0;
  
  // Factor of two here to drop the zero padding at the back of the waveform
  // which was used during correlations.
  TGraph* grCoherent = new TGraph(nSamp/2, &tArray[0], &vArray[0]);	
  
  for(Int_t deltaPhiSect=-maxDeltaPhiSect; deltaPhiSect<=maxDeltaPhiSect; deltaPhiSect++){

    Int_t phiSector = deltaPhiSect + centerPhiSector;
    phiSector = phiSector < 0 ? phiSector + NUM_PHI : phiSector;
    phiSector = phiSector >= NUM_PHI ? phiSector - NUM_PHI : phiSector;

    for(Int_t ring=0; ring<NUM_RING; ring++){
      Int_t ant= phiSector + ring*NUM_PHI;

      // Here we do the inverse FFT on the padded FFTs in memory.
      // The output is now in the vArray
      if(nSamp==numSamples){
	FancyFFTs::doInvFFT(nSamp, ffts[pol][ant], false);
      }
      else{
	FancyFFTs::doInvFFT(nSamp, fftsPadded[pol][ant], false);      
      }
	
      if(firstAnt!=ant){ // Don't do the first antenna twice
	Double_t deltaT = getDeltaTExpected(pol, firstAnt, ant, phiRad, thetaRad);
	Int_t offset1 = floor(deltaT/theDeltaT);
	Int_t offset2 = ceil(deltaT/theDeltaT);

	// How far between the samples we need to interpolate.
	Double_t ddt = deltaT - theDeltaT*offset1;
	
	for(Int_t samp=0; samp<grCoherent->GetN(); samp++){
	  Int_t samp1 = samp + offset1;
	  Int_t samp2 = samp + offset2;
	  if(samp1 >= 0 && samp2 < grCoherent->GetN()){

	    Double_t v1 = vArray[samp1];
	    Double_t v2 = vArray[samp2];

	    Double_t vInterp = (ddt)*(v2 - v1)/theDeltaT + v1;

	    grCoherent->GetY()[samp] += vInterp*interpRMS[pol][ant];
	  }
	}
      }

      // Here we look at the RMS of the first few ns of the uninterpolated waveforms
      const Double_t timeToEvalRms = 10; // ns
      for(Int_t samp3=0; samp3 < grs[pol][ant]->GetN(); samp3++){
	if(grs[pol][ant]->GetX()[samp3] - grs[pol][ant]->GetX()[0] < timeToEvalRms){
	  rms += grs[pol][ant]->GetY()[samp3]*grs[pol][ant]->GetY()[samp3];
	}
      }
      numAnts++;
    }
  }

  // Normalize
  if(numAnts > 0){
    TString name;
    TString title;
    for(Int_t samp=0; samp<grCoherent->GetN(); samp++){
      grCoherent->GetY()[samp]/=numAnts;
    }

    if(pol==AnitaPol::kHorizontal){
      name = TString::Format("grCoherentH_%u", eventNumber[pol]);
      title = "HPOL ";
    }
    else{
      name = TString::Format("grCoherentV_%u", eventNumber[pol]);
      title = "VPOL ";      
    }

    title += TString::Format("Coherently Summed Waveform for arrival direction elevation %4.2lf (Deg) and azimuth %4.2lf (Deg); Time (ns); Voltage (mV)", thetaDeg, phiDeg);

    grCoherent->SetName(name);
    grCoherent->SetTitle(title);
    
    rms/=numAnts;
    rms = TMath::Sqrt(rms);

    Double_t maxY = 0;
    Double_t maxX = 0;
    Double_t minY = 0;
    Double_t minX = 0;
    RootTools::getLocalMaxToMin(grCoherent, maxY, maxX, minY, minX);
    snr = (maxY - minY)/(2*rms);   
  }  
  return grCoherent;  
}
