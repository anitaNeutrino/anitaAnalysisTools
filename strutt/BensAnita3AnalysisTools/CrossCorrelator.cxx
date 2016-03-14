/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             A Cross Correlator to interact with the ROOTified ANITA-3 data and do some interferometry. 
*************************************************************************************************************** */

#include "CrossCorrelator.h"

/************************************************************************************************************
Constructor and destructor functions
************************************************************************************************************/

/*!
  \brief Constructor
*/
CrossCorrelator::CrossCorrelator(){
  initializeVariables();
}

/*!
  \brief Destructor
*/
CrossCorrelator::~CrossCorrelator(){
  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    deleteAllWaveforms((AnitaPol::AnitaPol_t)pol);
  }
}



/*!
  \brief Workhorse function to set internal variables.
*/
void CrossCorrelator::initializeVariables(){

  kDebug = false;
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

  kZeroChannel16BH = false;
  numSamples = 2*NUM_SAMPLES; // Factor of two for padding 
  numSamplesUpsampled = numSamples*UPSAMPLE_FACTOR; // For upsampling

  nominalSamplingDeltaT = NOMINAL_SAMPLING_DELTAT;
  correlationDeltaT = nominalSamplingDeltaT/UPSAMPLE_FACTOR;

  deltaTMax = 0;
  deltaTMin = 0;

  do5PhiSectorCombinatorics();
  
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  // geom->useKurtAnitaIIINumbers(1);
  // AnitaEventCalibrator* cal = AnitaEventCalibrator::Instance();
  // for(int surf=0; surf<12; surf++){
  //   for(int chan=0; chan < 9; chan++){
  //     cal->relativePhaseCenterToAmpaDelays[surf][chan] = 0;
  //   }
  // }
  // std::cout << this << std::endl;
  // AnitaEventCalibrator* cal = AnitaEventCalibrator::Instance();

  // std::cout << "pol\tant\trArray\tzArray\tphiArray" << std::endl;
  for(Int_t pol=0; pol < NUM_POL; pol++){
    for(int ant=0; ant<NUM_SEAVEYS; ant++){
      rArray[pol].push_back(geom->getAntR(ant, AnitaPol::AnitaPol_t(pol)));
      zArray[pol].push_back(geom->getAntZ(ant, AnitaPol::AnitaPol_t(pol)));
      phiArrayDeg[pol].push_back(geom->getAntPhiPositionRelToAftFore(ant, AnitaPol::AnitaPol_t(pol))*TMath::RadToDeg());

      // Int_t surf, chan, ant2;
      // geom->getSurfChanAntFromRingPhiPol(AnitaRing::AnitaRing_t (ant/NUM_PHI), ant%NUM_PHI, AnitaPol::AnitaPol_t (pol), surf, chan, ant2);
      // std::cout << pol << "\t" << ant << "\t" << rArray[pol].at(ant) << "\t" << zArray[pol].at(ant) << "\t" << phiArrayDeg[pol].at(ant) << "\t" << cal->relativePhaseCenterToAmpaDelays[surf][chan] << std::endl;
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
    usefulPats.push_back(new UsefulAdu5Pat());
  }
}


void CrossCorrelator::printInfo(){
  std::cerr << __PRETTY_FUNCTION__ << std::endl;
  std::cerr << "\tupsample factor = " << UPSAMPLE_FACTOR << std::endl;
  std::cerr << "\tdeltaT max = " << deltaTMax << std::endl;
  std::cerr << "\tdeltaT min = " << deltaTMin << std::endl;
  std::cerr << "\tBin size theta (deg) = " << Double_t(THETA_RANGE)/NUM_BINS_THETA << std::endl;
  std::cerr << "\tBin size phi (deg) = " << Double_t(PHI_RANGE)/NUM_BINS_PHI << std::endl;
  std::cerr << "\tdeltaTs array size = "
	    << sizeof(dtIndex_t)*NUM_POL*numCombos*NUM_PHI*NUM_BINS_PHI*NUM_BINS_THETA
	    << " bytes" << std::endl;
}






/*!
  \brief Single function to get the angle of the first bin of the interferometric histogram.
  \returns position of antenna 0 in ADU5Pat coordinates, offset by half a phi-sector.
*/

Double_t CrossCorrelator::getBin0PhiDeg(){
  // Double_t phi0 = phiArrayDeg[0].at(0);
  Double_t phi0 = -45; //phiArrayDeg[0].at(0);  
  if(phi0 < -180){
    phi0+=360;
  }
  else if(phi0 >= 180){
    phi0-=360;
  }
  return phi0 - PHI_RANGE/2;
}









/************************************************************************************************************
Waveform manipulation functions
************************************************************************************************************/


/*!
  \brief Loops through all waveform graphs in the UsefulAnitaEvent and makes an evenly re-sampled, normalized copy of each one.
  \param usefulEvent points to the UsefulAnitaEvent of interest
*/

void CrossCorrelator::getNormalizedInterpolatedTGraphs(UsefulAnitaEvent* usefulEvent,
						       AnitaPol::AnitaPol_t pol){
  // Potentially needed in a few places, so it gets its own function 

  // Pretty much just for profiling 
  // if(usefulEvent->eventNumber!=lastEventNormalized[pol]){

    // Delete any old waveforms (at start rather than end to leave in memory to be examined if need be)
    deleteAllWaveforms(pol);

    // Find the start time of all waveforms 
    std::vector<Double_t> earliestStart(NUM_POL, 100); // ns
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      grs[pol][ant] = usefulEvent->getGraph(ant, (AnitaPol::AnitaPol_t)pol);
      

      if(kZeroChannel16BH==true){      
	//// DELIBERATE ZEROING OF 16BH. THIS IS THE CHANNEL WHERE WE DON'T HAVE THE CABLE DELAY
	//// I WANT TO FIT THE PITCH/ROLL OFFSET WITH REASONABLE CABLE DELAYS!!!!
	//// DOES THIS HELP THE POINTING??????
	if(ant==47 && pol==0){
	  for(Int_t samp=0; samp < grs[pol][ant]->GetN(); samp++){
	    grs[pol][ant]->GetY()[samp] = 0;
	  }
	}
      }
      
      if(grs[pol][ant]->GetX()[0]<earliestStart.at(pol)){
	earliestStart.at(pol) = grs[pol][ant]->GetX()[0];
      }
    }

    // Interpolate with earliest start time 
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      grsResampled[pol][ant] = interpolateWithStartTime(grs[pol][ant], earliestStart.at(pol));
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
  // }
}


void CrossCorrelator::doFFTs(AnitaPol::AnitaPol_t pol){

  // deleteAllFFTs(pol); // First delete any FFTs that might still be knocking around

  // Generate a new set of FFTs.
  // for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
  //   ffts[pol][ant] = FancyFFTs::doFFT(numSamples, grsResampled[pol][ant]->GetY(), true);    
  // }
  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
    FancyFFTs::doFFT(numSamples, grsResampled[pol][ant]->GetY(), ffts[pol][ant]);    
  }
}



TGraph* CrossCorrelator::interpolateWithStartTime(TGraph* grIn, Double_t startTime){
  return RootTools::interpolateWithStartTime(grIn, startTime, nominalSamplingDeltaT, numSamples);
}

















/************************************************************************************************************
All correlation functions
************************************************************************************************************/


/*!
  \brief Correlate the event

  First step in interferometry, pass pointer to UsefulAnitaEvent in here, which populates the internal arrays of normalized, interpolated waveforms and set of cross correlations. These internal values are then used as look up tables for generating interferometic images.
*/
void CrossCorrelator::correlateEvent(UsefulAnitaEvent* usefulEvent){

  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){  
    correlateEvent(usefulEvent, (AnitaPol::AnitaPol_t)pol);
  }
}


/*!
  \brief Correlate the event

  First step in interferometry, pass pointer to UsefulAnitaEvent in here, which populates the internal arrays of normalized, interpolated waveforms and set of cross correlations. These internal values are then used as look up tables for generating interferometic images.
*/

void CrossCorrelator::correlateEvent(UsefulAnitaEvent* usefulEvent, AnitaPol::AnitaPol_t pol){

  // Read TGraphs from events into memory (also deletes old TGraphs)
  getNormalizedInterpolatedTGraphs(usefulEvent, pol);

  // Generate set of ffts for cross correlation (each waveform only needs to be done once)
  doFFTs(pol);

  // Now cross correlate those already FFT'd waveforms
  // doAllCrossCorrelations(pol);
  doAllCrossCorrelationsThreaded(pol);

  // Safety check to make sure we don't do any hard work twice.
  eventNumber[pol] = usefulEvent->eventNumber;
}



/*!
  \brief Loop over both polarizations and all combinations, and generate set of cross correlations.
*/
void CrossCorrelator::doAllCrossCorrelations(AnitaPol::AnitaPol_t pol){

  // Delete old cross correlations first 
  // deleteCrossCorrelations(pol);

  // Loop over combinations generating set of cross correlations.
  for(Int_t combo=0; combo<numCombos; combo++){
    Int_t ant1 = comboToAnt1s.at(combo);
    Int_t ant2 = comboToAnt2s.at(combo);
    FancyFFTs::crossCorrelate(numSamples, ffts[pol][ant2], ffts[pol][ant1], crossCorrelations[pol][combo]);        
  }
}


void CrossCorrelator::doAllCrossCorrelationsThreaded(AnitaPol::AnitaPol_t pol){

  // Delete old cross correlations first 
  // deleteCrossCorrelations(pol);

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

void* CrossCorrelator::doSomeCrossCorrelationsThreaded(void* voidPtrArgs){

  // Disgusting hacks to get ROOT threading to compile inside a class.
  CrossCorrelator::threadArgs* args = (CrossCorrelator::threadArgs*) voidPtrArgs;
  Long_t threadInd = args->threadInd;
  CrossCorrelator* ptr = args->ptr;
  AnitaPol::AnitaPol_t pol = ptr->threadPol;

  Int_t numCorrPerThread = NUM_COMBOS/NUM_THREADS;
  Int_t startCombo = threadInd*numCorrPerThread;

  for(int combo=startCombo; combo<startCombo+numCorrPerThread; combo++){
    Int_t ant1 = ptr->comboToAnt1s.at(combo);
    Int_t ant2 = ptr->comboToAnt2s.at(combo);
    FancyFFTs::crossCorrelate(ptr->numSamples,
			      ptr->ffts[pol][ant2],
			      ptr->ffts[pol][ant1],
			      ptr->crossCorrelations[pol][combo],
			      threadInd);
  }  
  return 0;
}









void CrossCorrelator::doUpsampledCrossCorrelations(AnitaPol::AnitaPol_t pol, UShort_t l3TrigPattern){

  // deleteAllPaddedFFTs(pol);
  // deleteUpsampledCrossCorrelations(pol);
  std::vector<Int_t> doneCombos(numCombos, 0);
  std::vector<Int_t> combosToUse = combosToUseTriggered[l3TrigPattern];
  std::vector<Int_t> donePadding(NUM_SEAVEYS, 0);
  
  Int_t numFreqs = FancyFFTs::getNumFreqs(numSamples);
  Int_t numFreqsPadded = FancyFFTs::getNumFreqs(numSamplesUpsampled);

  for(int phiSector=0; phiSector<NUM_PHI; phiSector++){
    Int_t numComboInd = combosToUse.size();
    for(Int_t comboInd=0; comboInd < numComboInd; comboInd++){
      Int_t combo=combosToUse.at(comboInd);

      // if(crossCorrelationsUpsampled[pol][combo]==NULL){	
      if(doneCombos.at(combo)==0){

	Int_t ant1 = comboToAnt1s.at(combo);
	Int_t ant2 = comboToAnt2s.at(combo);

	if(donePadding.at(ant1)==0){
	  FancyFFTs::zeroPadFFT(ffts[pol][ant1], fftsPadded[pol][ant1], numFreqs, numFreqsPadded);
	  donePadding.at(ant1)=1;
	}
	if(donePadding.at(ant2)==0){	
	  FancyFFTs::zeroPadFFT(ffts[pol][ant2], fftsPadded[pol][ant1], numFreqs, numFreqsPadded);
	  donePadding.at(ant2)=1;	  
	}
	
	FancyFFTs::crossCorrelate(numSamplesUpsampled,
				  fftsPadded[pol][ant2],
				  fftsPadded[pol][ant1],
				  crossCorrelationsUpsampled[pol][combo]);
	
	doneCombos.at(combo) = 1;
      }
    }
  }
}


void CrossCorrelator::doUpsampledCrossCorrelationsThreaded(AnitaPol::AnitaPol_t pol, UShort_t l3TrigPattern){

  threadL3TrigPattern = l3TrigPattern;
  threadPol = pol;

  std::vector<Int_t> combosToUse = combosToUseTriggered[l3TrigPattern];
  std::vector<Int_t> donePadding(NUM_SEAVEYS, 0);
  
  Int_t numFreqs = FancyFFTs::getNumFreqs(numSamples);
  Int_t numFreqsPadded = FancyFFTs::getNumFreqs(numSamplesUpsampled);

  Int_t numComboInd = combosToUse.size();
  for(Int_t comboInd=0; comboInd < numComboInd; comboInd++){
    Int_t combo=combosToUse.at(comboInd);

    Int_t ant1 = comboToAnt1s.at(combo);
    Int_t ant2 = comboToAnt2s.at(combo);

    if(donePadding.at(ant1)==0){
      FancyFFTs::zeroPadFFT(ffts[pol][ant1], fftsPadded[pol][ant1], numFreqs, numFreqsPadded);
      donePadding.at(ant1)=1;
    }
    if(donePadding.at(ant2)==0){	
      FancyFFTs::zeroPadFFT(ffts[pol][ant2], fftsPadded[pol][ant2], numFreqs, numFreqsPadded);
      donePadding.at(ant2)=1;	  
    }
  }

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

void* CrossCorrelator::doSomeUpsampledCrossCorrelationsThreaded(void* voidPtrArgs){

  // Disgusting hacks to get ROOT threading to compile inside a class.
  CrossCorrelator::threadArgs* args = (CrossCorrelator::threadArgs*) voidPtrArgs;
  Long_t threadInd = args->threadInd;
  CrossCorrelator* ptr = args->ptr;
  AnitaPol::AnitaPol_t pol = ptr->threadPol;
  
  std::vector<Int_t> combosToUse = ptr->combosToUseTriggered[ptr->threadL3TrigPattern];
  Int_t numCombosAllThreads = combosToUse.size();
  Int_t numCorrPerThread = numCombosAllThreads/NUM_THREADS;
  Int_t numRemainder = numCombosAllThreads%NUM_THREADS;

  Int_t startComboInd = threadInd*numCorrPerThread;
  
  for(int comboInd=startComboInd; comboInd<startComboInd+numCorrPerThread; comboInd++){
    Int_t combo = combosToUse.at(comboInd);
    Int_t ant1 = ptr->comboToAnt1s.at(combo);
    Int_t ant2 = ptr->comboToAnt2s.at(combo);
    
    FancyFFTs::crossCorrelate(ptr->numSamplesUpsampled,
			      ptr->fftsPadded[pol][ant2],
			      ptr->fftsPadded[pol][ant1],
			      ptr->crossCorrelationsUpsampled[pol][combo],
			      threadInd);
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
    
  }
  return 0;
}










/*!
  \brief Interface to FFT library to generate cross correlations.
*/
Double_t* CrossCorrelator::crossCorrelateFourier(TGraph* gr1, TGraph* gr2){
  return FancyFFTs::crossCorrelate(gr1->GetN(), gr1->GetY(), gr2->GetY()) ;  
}

/*!
  \brief Interface to FFT library to generate cross correlations.
*/
Double_t* CrossCorrelator::crossCorrelateFourier(Int_t numSamplesInTimeDomain,
						 std::complex<Double_t>* fft1,
						 std::complex<Double_t>* fft2,
						 Int_t threadInd){  
  // return FancyFFTs::crossCorrelate(numSamplesInTimeDomain, fft1, fft2, threadInd);
  return FancyFFTs::crossCorrelate(numSamplesInTimeDomain, fft1, fft2, threadInd);    
}








/************************************************************************************************************
Get correlation summary information
************************************************************************************************************/



void CrossCorrelator::getMaxCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t combo,
						 Double_t& time, Double_t& value){

  Int_t maxInd = RootTools::getIndexOfMaximum(numSamples, crossCorrelations[pol][combo]);
  Int_t offset = maxInd >= numSamples/2 ? maxInd - numSamples : maxInd;
  value = crossCorrelations[pol][combo][maxInd];
  time = offset*nominalSamplingDeltaT;
}


void CrossCorrelator::getMaxCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
						 Double_t& time, Double_t& value){
  Int_t combo = comboIndices[ant1][ant2];
  getMaxCorrelationTimeValue(pol, combo, time, value);
}

void CrossCorrelator::getMaxUpsampledCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t combo,
							  Double_t& time, Double_t& value){

 Int_t maxInd = RootTools::getIndexOfMaximum(numSamplesUpsampled, crossCorrelationsUpsampled[pol][combo]);
  Int_t offset = maxInd >= numSamplesUpsampled/2 ? maxInd - numSamplesUpsampled : maxInd;
  value = crossCorrelationsUpsampled[pol][combo][maxInd];
  time = offset*correlationDeltaT;
}


void CrossCorrelator::getMaxUpsampledCorrelationTimeValue(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
							  Double_t& time, Double_t& value){
  Int_t combo = comboIndices[ant1][ant2];
  getMaxUpsampledCorrelationTimeValue(pol, combo, time, value);
}




























/************************************************************************************************************
Calculate deltaT between two antennas (for a plane wave unless function name says otherwise)
************************************************************************************************************/

Double_t CrossCorrelator::getDeltaTExpected(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t phiWave, Double_t thetaWave){
  
  Double_t tanThetaW = tan(-thetaWave);

  // if(kDebug==true){
  //   std::cerr << "deltaTExpected tanThetaW = " << tanThetaW << std::endl;
  //   std::cerr << "deltaTExpected thetaWave = " << thetaWave << std::endl;    
  // }
  Double_t part1 = zArray[pol].at(ant1)*tanThetaW - rArray[pol].at(ant1)*cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol].at(ant1));
  Double_t part2 = zArray[pol].at(ant2)*tanThetaW - rArray[pol].at(ant2)*cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol].at(ant2));
  Double_t tdiff = 1e9*((cos(-thetaWave) * (part2 - part1))/SPEED_OF_LIGHT); // Returns time in ns
  return tdiff;
}


Double_t CrossCorrelator::getDeltaTExpectedPat(Int_t ant1, Int_t ant2,Double_t phiWave, Double_t thetaWave, Int_t print){
  AnitaGeomTool* fUPGeomTool = AnitaGeomTool::Instance();
  
  Double_t phi1=fUPGeomTool->getAntPhiPositionRelToAftFore(ant1);
  Double_t r1=fUPGeomTool->getAntR(ant1);
  Double_t z1=fUPGeomTool->getAntZ(ant1);

  Double_t phi2=fUPGeomTool->getAntPhiPositionRelToAftFore(ant2);
  Double_t r2=fUPGeomTool->getAntR(ant2);
  Double_t z2=fUPGeomTool->getAntZ(ant2);

  if(print){
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr << phiWave << "\t" << thetaWave << std::endl;
    std::cerr << TMath::Tan(thetaWave) << "\t" << tan(thetaWave) << std::endl;
    std::cerr << phi1 << "\t" << r1 << "\t" << z1 << std::endl;
    std::cerr << phi2 << "\t" << r2 << "\t" << z2 << std::endl;
  }

  Double_t tanThetaW=TMath::Tan(thetaWave);
  Double_t part1=z1*tanThetaW - r1 * TMath::Cos(phiWave-phi1);
  Double_t part2=z2*tanThetaW - r2 * TMath::Cos(phiWave-phi2);

  if(print){
    std::cerr << tanThetaW << "\t" << part1 << "\t" << part2 << std::endl;
    std::cerr << std::endl;
  }
  
  return  1e9*((TMath::Cos(thetaWave) * (part1 - part2))/C_LIGHT);    //returns time in ns
}


Double_t CrossCorrelator::getDeltaTExpectedFast(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
						Int_t phiIndex, Int_t thetaIndex){
  
  Double_t tanThetaW = zoomedTanThetaWaves[thetaIndex];
  // if(kDebug==true){
  //   Double_t thetaWave = zoomedThetaWaves[thetaIndex];
  //   TThread::Lock();
  //   std::cerr << "deltaTExpectedFast tanThetaW = " << tanThetaW << "\t" << std::endl;
  //   std::cerr << "deltaTExpectedFast thetaWave = " << thetaWave << "\t" << std::endl;    
  //   TThread::UnLock();
  // }
  
  Double_t part1 = zArray[pol].at(ant1)*tanThetaW - zoomedCosPartLookup[phiIndex][pol][ant1];
  Double_t part2 = zArray[pol].at(ant2)*tanThetaW - zoomedCosPartLookup[phiIndex][pol][ant2];
  // Double_t tdiff = 1e9*((cos(-thetaWave) * (part2 - part1))/SPEED_OF_LIGHT); // Returns time in ns 
  Double_t tdiff = 1e9*((zoomedCosThetaWaves[thetaIndex] * (part2 - part1))/SPEED_OF_LIGHT); // Returns time in ns 
  return tdiff;
}


// Moving parts...
// tan(-thetaW) NUM_BINS_THETA
// Part1/2 this is adding up two numbers
//    1st bit... zArray[pol][ant1]*tanThetaW  NUM_POL*NUM_SEAVEYS*NUM_BINS_THETA
//    2nd bit... rArray[pol][ant1]*cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol][ant]) NUM_POL*NUM_SEAVEYS*NUM_BINS_PHI*NUM_PHI
// cos(-thetaW) NUM_BINS_THETA

// So going from an array of size ...
// [NUM_POL][NUM_PHI*NUM_BINS_PHI][NUM_BINS_THETA][NUM_COMBOS]; sizeof(dtIndex) (char)
// 2 * 16 * 25 * 150 * 336 = 40320000 = 4e7 * 1 byte = 47 MB.

// To
// + 150
// + 2 * 48 * 150
// + 2 * 48 * 25*16
// + 150
// =  150   +   2 * 48 * 150   +   2 * 48 * 25*16  +   150 = 53100 = 5.31e5 * sizeof(dtIndex) * 
// #define NUM_COMBOS 336
// #define NUM_THREADS 4
 
// // Typical number of samples in waveform
// #define NUM_SAMPLES 256
// #define UPSAMPLE_FACTOR 40

// // // Image definitions
// #define NUM_BINS_THETA 150
// #define NUM_BINS_PHI 25


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
  fillDeltaTLookup();
  geom->useKurtAnita3Numbers(0);
  
}




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









/************************************************************************************************************
Precalculate DeltaTs during initialization where appropriate
************************************************************************************************************/

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

void CrossCorrelator::fillDeltaTLookup(){
  
  Double_t phi0 = getBin0PhiDeg();
  for(Int_t polInd=0; polInd<NUM_POL; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(Int_t phiSector = 0; phiSector<NUM_PHI; phiSector++){
      for(Int_t phiInd = 0; phiInd < NUM_BINS_PHI; phiInd++){
	Int_t phiBin = phiSector*NUM_BINS_PHI + phiInd;

	Double_t phiDeg = phi0 + phiBin*Double_t(PHI_RANGE)/NUM_BINS_PHI;

	// if(polInd==0){
	//   std::cerr.precision(8);
	//   std::cerr << phi0 << "\t" << phiBin << "\t" << phiDeg << std::endl;
	// }
	
	Double_t phiWave = TMath::DegToRad()*phiDeg;
	for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA; thetaBin++){
	  Double_t thetaDeg = THETA_RANGE*((Double_t)thetaBin/NUM_BINS_THETA - 0.5);
	  Double_t thetaWave = TMath::DegToRad()*thetaDeg;
	  for(Int_t combo=0; combo<numCombos; combo++){
	    Int_t ant1 = comboToAnt1s.at(combo);
	    Int_t ant2 = comboToAnt2s.at(combo);
	    Double_t deltaT = getDeltaTExpected(pol, ant1, ant2, phiWave, thetaWave);
	    Int_t offset = TMath::Nint(deltaT/nominalSamplingDeltaT);
	    if(offset > deltaTMax){
	      deltaTMax = offset;
	    }
	    if(offset < deltaTMin){
	      deltaTMin = offset;
	    }
	    deltaTs[pol][phiBin][thetaBin][combo] = offset;
	  }
	}
      }
    }
  }

  for(Int_t thetaIndex=0; thetaIndex < NUM_BINS_THETA_ZOOM_TOTAL; thetaIndex++){
    Double_t thetaWaveDeg = (thetaIndex-NUM_BINS_THETA_ZOOM_TOTAL/2)*ZOOM_BIN_SIZE_THETA;
    Double_t thetaWave = thetaWaveDeg*TMath::DegToRad();
    zoomedThetaWaves[thetaIndex] = thetaWave;
    // kDebug=true;
    // if(kDebug){
    //   std::cerr << thetaIndex << "\t" << thetaWave << "\t" << thetaWave*TMath::RadToDeg() << std::endl;
    // }
    zoomedTanThetaWaves[thetaIndex] = tan(-thetaWave);
    zoomedCosThetaWaves[thetaIndex] = cos(-thetaWave);
  }
  // throw 0;

  for(Int_t phiIndex=0; phiIndex < NUM_BINS_PHI_ZOOM_TOTAL; phiIndex++){
    Double_t phiWave = TMath::DegToRad()*phiIndex*ZOOM_BIN_SIZE_PHI;
    zoomedPhiWaveLookup[phiIndex] = phiIndex*ZOOM_BIN_SIZE_PHI;
    for(Int_t pol=0; pol<AnitaPol::kNotAPol; pol++){
      for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
	zoomedCosPartLookup[phiIndex][pol][ant] = rArray[pol].at(ant)*cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol].at(ant));
	 // = cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol][ant]);
	
      }
    }
  }
}


Bool_t CrossCorrelator::useCombo(Int_t ant1, Int_t ant2, Int_t phiSector){

  Bool_t ant1InRange = TMath::Abs(phiSector - (ant1%NUM_PHI))<=DELTA_PHI_SECT;
  ant1InRange = ant1InRange || TMath::Abs(phiSector - (ant1%NUM_PHI))>=(NUM_PHI-DELTA_PHI_SECT);

  Bool_t ant2InRange = TMath::Abs(phiSector - (ant2%NUM_PHI))<=DELTA_PHI_SECT;
  ant2InRange = ant2InRange || TMath::Abs(phiSector - (ant2%NUM_PHI))>=(NUM_PHI-DELTA_PHI_SECT);  

  return (ant1InRange && ant2InRange);
}
    
void CrossCorrelator::fillCombosToUseIfNeeded(mapMode_t mapMode, UShort_t l3TrigPattern){
  
  if(mapMode==kTriggered){
    std::map<UInt_t,std::vector<Int_t> >::iterator it = combosToUseTriggered.find(l3TrigPattern);
    if(it==combosToUseTriggered.end()){
      combosToUseTriggered[l3TrigPattern] = std::vector<Int_t>();
      for(Int_t phiSector = 0; phiSector<NUM_PHI; phiSector++){
	UInt_t doPhiSector = ((l3TrigPattern >> phiSector) & 1);
	if(doPhiSector){
	  for(Int_t combo=0; combo<numCombos; combo++){
	    Int_t ant1 = comboToAnt1s.at(combo);
	    Int_t ant2 = comboToAnt2s.at(combo);
	    if(useCombo(ant1, ant2, phiSector)){
	      if(std::find(combosToUseTriggered[l3TrigPattern].begin(),
			   combosToUseTriggered[l3TrigPattern].end(),
			   combo) == combosToUseTriggered[l3TrigPattern].end()){
		combosToUseTriggered[l3TrigPattern].push_back(combo);
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
	  if(useCombo(ant1, ant2, phiSector)){
	    combosToUseGlobal[phiSector].push_back(combo);
	  }
	}
      }
    }
  }
}





























/************************************************************************************************************
Image generation functions.
************************************************************************************************************/


TH2D* CrossCorrelator::makeBlankImage(TString name, TString title){

  Double_t phiMin = getBin0PhiDeg();
  Double_t phiMax = phiMin + 360;
  Double_t thetaMin = -THETA_RANGE/2;
  Double_t thetaMax = THETA_RANGE/2;

  TH2D* hImage = new TH2D(name, title,
			  NUM_BINS_PHI*NUM_PHI, phiMin, phiMax,
			  NUM_BINS_THETA, thetaMin, thetaMax);
  hImage->GetXaxis()->SetTitle("Azimuth (Degrees)");
  hImage->GetYaxis()->SetTitle("Elevation (Degrees)");
  return hImage;
}


TH2D* CrossCorrelator::makeBlankZoomedImage(TString name, TString title,
					    Double_t zoomCenterPhiDeg,
					    Double_t zoomCenterThetaDeg){

  // std::cout << zoomCenterPhiDeg << "\t" << zoomCenterThetaDeg << "\t";
  zoomCenterPhiDeg = (TMath::Nint(zoomCenterPhiDeg/ZOOM_BIN_SIZE_PHI))*ZOOM_BIN_SIZE_PHI;
  zoomCenterThetaDeg = (TMath::Nint(zoomCenterThetaDeg/ZOOM_BIN_SIZE_THETA))*ZOOM_BIN_SIZE_THETA;

  // if(kDebug){
  //   std::cerr << "zoomCenterPhiDeg = " << zoomCenterPhiDeg << "\nzoomCenterThetaDeg = "
  // 	      << zoomCenterThetaDeg << std::endl;
  // }
  
  // std::cout << zoomCenterPhiDeg << "\t" << zoomCenterThetaDeg << "\t" << std::endl;
  
  Double_t phiMin = zoomCenterPhiDeg - PHI_RANGE_ZOOM/2;
  Double_t phiMax = zoomCenterPhiDeg + PHI_RANGE_ZOOM/2;
  // Double_t thetaMin = zoomCenterThetaDeg - THETA_RANGE_ZOOM/2;
  // Double_t thetaMax = zoomCenterThetaDeg + THETA_RANGE_ZOOM/2;
  
  Double_t thetaMin = zoomCenterThetaDeg - THETA_RANGE_ZOOM/2;
  Double_t thetaMax = zoomCenterThetaDeg + THETA_RANGE_ZOOM/2;

  TH2D* hImage = new TH2D(name, title,
			  NUM_BINS_PHI_ZOOM, phiMin, phiMax,
			  NUM_BINS_THETA_ZOOM, thetaMin, thetaMax);
  hImage->GetXaxis()->SetTitle("Azimuth (Degrees)");
  hImage->GetYaxis()->SetTitle("Elevation (Degrees)");

  // if(kDebug){
  //   std::cerr << "In making blank zoomed in histogram" << std::endl;
  //   std::cerr << "zoomCenterPhiDeg = " << zoomCenterPhiDeg << "\n"
  // 	      << "hImage->GetXaxis()->GetBinLowEdge(1) = " << hImage->GetXaxis()->GetBinLowEdge(1) << std::endl;
  // }
  
  return hImage;

}


TH2D* CrossCorrelator::makeGlobalSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave){
  Double_t imagePeak, peakPhiDeg, peakThetaDeg;

  return makeGlobalSphericalImage(pol, rWave, imagePeak, peakPhiDeg, peakThetaDeg);  

}

TH2D* CrossCorrelator::makeGlobalSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave,
						Double_t& imagePeak, Double_t& peakPhiDeg,
						Double_t& peakThetaDeg){

  return makeImageThreaded(pol, rWave, imagePeak, peakPhiDeg, peakThetaDeg, ALL_PHI_TRIGS,
  			   kGlobal, kZoomedOut);
  // return makeImage(pol, rWave, imagePeak, peakPhiDeg, peakThetaDeg, ALL_PHI_TRIGS,
  // 		   kGlobal, kZoomedOut);
}


TH2D* CrossCorrelator::makeTriggeredSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave, UShort_t l3TrigPattern){
  Double_t imagePeak, peakPhiDeg, peakThetaDeg;
  return makeTriggeredSphericalImage(pol, rWave, imagePeak, peakPhiDeg, peakThetaDeg, l3TrigPattern);
}

TH2D* CrossCorrelator::makeTriggeredSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave,
						   Double_t& imagePeak, Double_t& peakPhiDeg,
						   Double_t& peakThetaDeg, UShort_t l3TrigPattern){

  return makeImageThreaded(pol, rWave, imagePeak, peakPhiDeg, peakThetaDeg,
  			   l3TrigPattern, kTriggered, kZoomedOut);
  // return makeImage(pol, rWave, imagePeak, peakPhiDeg, peakThetaDeg,
  // 		   l3TrigPattern, kTriggered, kZoomedOut);    
}


TH2D* CrossCorrelator::makeGlobalImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
				       Double_t& peakThetaDeg){

  return makeImageThreaded(pol, 0, imagePeak, peakPhiDeg, peakThetaDeg, ALL_PHI_TRIGS,
  			   kGlobal, kZoomedOut);
  // return makeImage(pol, 0, imagePeak, peakPhiDeg, peakThetaDeg, ALL_PHI_TRIGS,
  // 		   kGlobal, kZoomedOut);  
}

TH2D* CrossCorrelator::makeGlobalImage(AnitaPol::AnitaPol_t pol){
  Double_t imagePeak, peakPhiDeg, peakThetaDeg;
  return makeGlobalImage(pol, imagePeak, peakPhiDeg, peakThetaDeg);
}


TH2D* CrossCorrelator::makeTriggeredImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak,
					  Double_t& peakPhiDeg, Double_t& peakThetaDeg,
					  UShort_t l3TrigPattern){

  return makeImageThreaded(pol, 0, imagePeak, peakPhiDeg, peakThetaDeg,
  			   l3TrigPattern, kTriggered, kZoomedOut);  
  // return makeImage(pol, 0, imagePeak, peakPhiDeg, peakThetaDeg,
  // 		   l3TrigPattern, kTriggered, kZoomedOut);  
}


TH2D* CrossCorrelator::makeZoomedImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
				       Double_t& peakThetaDeg, UShort_t l3TrigPattern,
				       Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg){

  return makeImageThreaded(pol, 0, imagePeak, peakPhiDeg, peakThetaDeg, l3TrigPattern,
  			   kTriggered, kZoomedIn, zoomCenterPhiDeg, zoomCenterThetaDeg);
  // return makeImage(pol, 0, imagePeak, peakPhiDeg, peakThetaDeg, l3TrigPattern,
  // 		   kTriggered, kZoomedIn, zoomCenterPhiDeg, zoomCenterThetaDeg);

}

TH2D* CrossCorrelator::makeZoomedImage(AnitaPol::AnitaPol_t pol, UShort_t l3TrigPattern,
				       Double_t zoomCenterPhiDeg,Double_t zoomCenterThetaDeg){

  Double_t imagePeak, peakPhiDeg, peakThetaDeg;
  return makeZoomedImage(pol, imagePeak, peakPhiDeg, peakThetaDeg,
  			 l3TrigPattern, zoomCenterPhiDeg, zoomCenterThetaDeg);

}

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


TH2D* CrossCorrelator::prepareForImageMaking(AnitaPol::AnitaPol_t pol, Double_t rWave, UShort_t l3TrigPattern,
					     mapMode_t mapMode, zoomMode_t zoomMode,
					     Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg){
  fillCombosToUseIfNeeded(mapMode, l3TrigPattern);

  TString name, title;
  createImageNameAndTitle(name, title, mapMode, zoomMode, rWave, pol);
  
  TH2D* hImage = NULL;
  if(zoomMode == kZoomedOut){
    hImage = makeBlankImage(name, title);
  }
  else if(zoomMode == kZoomedIn){
    hImage = makeBlankZoomedImage(name, title, zoomCenterPhiDeg, zoomCenterThetaDeg);
  }

  threadImage = hImage;
  threadPol = pol;
  threadMapMode = mapMode;
  threadZoomMode = zoomMode;
  threadRWave = rWave;
  threadL3TrigPattern = l3TrigPattern;
  
  if(zoomMode == kZoomedIn){
    doUpsampledCrossCorrelationsThreaded(pol, l3TrigPattern);
  }
  return hImage;
}



TH2D* CrossCorrelator::makeImage(AnitaPol::AnitaPol_t pol, Double_t rWave, Double_t& imagePeak,
				 Double_t& peakPhiDeg, Double_t& peakThetaDeg, UShort_t l3TrigPattern,
				 mapMode_t mapMode, zoomMode_t zoomMode, Double_t zoomCenterPhiDeg,
				 Double_t zoomCenterThetaDeg){

  TH2D* hImage = prepareForImageMaking(pol, rWave, l3TrigPattern,
				       mapMode, zoomMode, zoomCenterPhiDeg,
				       zoomCenterThetaDeg);

  imagePeak = -DBL_MAX;
  Int_t peakPhiBin = -1;
  Int_t peakThetaBin = -1;
  std::vector<Int_t> combosToUse;
  Int_t lastPhiSector = 0;
  if(mapMode==kTriggered){
    combosToUse = combosToUseTriggered[l3TrigPattern];
  }
  else{
    combosToUse = combosToUseGlobal[0];
  }
  
  Int_t phiZoomBase = TMath::Nint(hImage->GetXaxis()->GetBinLowEdge(1)/ZOOM_BIN_SIZE_PHI);
  Int_t thetaZoomBase = TMath::Nint(hImage->GetYaxis()->GetBinLowEdge(1)/ZOOM_BIN_SIZE_THETA + NUM_BINS_THETA_ZOOM_TOTAL/2);
  
  for(Int_t phiBin = 0; phiBin < hImage->GetNbinsX(); phiBin++){
    Int_t phiSector = zoomMode==kZoomedIn ? 0 : phiBin/NUM_BINS_PHI;
    if(phiSector!=lastPhiSector && mapMode==kGlobal){
      combosToUse = combosToUseGlobal[phiSector];
      lastPhiSector = phiSector;
    }

    Int_t zoomPhiInd = phiZoomBase + phiBin;
    zoomPhiInd = zoomPhiInd < 0 ? zoomPhiInd + NUM_BINS_PHI_ZOOM_TOTAL : zoomPhiInd;
    zoomPhiInd = zoomPhiInd >= NUM_BINS_PHI_ZOOM_TOTAL ? zoomPhiInd - NUM_BINS_PHI_ZOOM_TOTAL : zoomPhiInd;  
    
    Double_t phiWave = hImage->GetXaxis()->GetBinLowEdge(phiBin+1)*TMath::DegToRad();

    for(Int_t thetaBin = 0; thetaBin < hImage->GetNbinsY(); thetaBin++){
      Double_t thetaWave = hImage->GetYaxis()->GetBinLowEdge(thetaBin+1)*TMath::DegToRad();
      Int_t zoomThetaInd = thetaZoomBase + thetaBin;

      Double_t correlations = 0;

      for(UInt_t comboInd=0; comboInd<combosToUse.size(); comboInd++){
	Int_t combo = combosToUse.at(comboInd);
	Int_t offset = 0;

	// If we are in zoomed out & plane wave mode then use the lookup table for a big speedup
	if(zoomMode==kZoomedOut && rWave==0){
          offset = deltaTs[pol][phiBin][thetaBin][combo];
	  offset = offset < 0 ? offset + numSamples : offset;
	  correlations += crossCorrelations[pol][combo][offset];
	}
	// If we are in zoomed in & plane wave mode then calculate
	else if(zoomMode==kZoomedIn && rWave==0){
	  Int_t ant1 = comboToAnt1s.at(combo);
	  Int_t ant2 = comboToAnt2s.at(combo);


	  Double_t deltaT = 0;
	  if(zoomThetaInd >= 0 && zoomThetaInd < NUM_BINS_THETA_ZOOM_TOTAL){
	    deltaT = getDeltaTExpectedFast(pol, ant1, ant2,
					   zoomPhiInd,
					   zoomThetaInd);
	  }
	  else{
	    // Double_t thetaDeg = hImage->GetYaxis()->GetBinLowEdge(thetaBin+1);
	    // Double_t phiDeg = hImage->GetXaxis()->GetBinLowEdge(phiBin+1);
	    // Double_t thetaWave = thetaDeg*TMath::DegToRad();
	    // Double_t phiWave = phiDeg*TMath::DegToRad();
	    deltaT = getDeltaTExpected(pol, ant1, ant2, phiWave, thetaWave);
	  }

	  
	  // Double_t deltaT = getDeltaTExpectedFast(pol, ant1, ant2,
	  // 					  phiBin + phiZoomBase,
	  // 					  thetaBin + thetaZoomBase);

	  if(kDebug){
	      Double_t deltaT_slow = getDeltaTExpected(pol, ant1, ant2, phiWave, thetaWave);
	      
	      if(TMath::Abs(deltaT_slow - deltaT) > 0.0001){
	      std::cerr << ant1 << "\t" << ant2 << "\t" << phiWave << "\t" << thetaWave << "\t" << deltaT_slow << "\t" << deltaT << std::endl;
	    }
	  }
	  
	  offset = TMath::Nint(deltaT/correlationDeltaT);	  
	  offset = offset < 0 ? offset + numSamplesUpsampled : offset;
	  correlations += crossCorrelationsUpsampled[pol][combo][offset];
	}
	// If we are in zoomed out & spherical wave mode then calculate
	else{
	  Int_t ant1 = comboToAnt1s.at(combo);
	  Int_t ant2 = comboToAnt2s.at(combo);
	  offset = getDeltaTExpectedSpherical(pol, ant1, ant2,
					      phiWave, thetaWave, rWave);
	  offset = offset < 0 ? offset + numSamples : offset;
	  correlations += crossCorrelations[pol][combo][offset];
	}
      }
      if(combosToUse.size()>0){
	correlations /= combosToUse.size();
      }
      hImage->SetBinContent(phiBin + 1, thetaBin + 1, correlations);
      if(correlations > imagePeak){
	imagePeak = correlations;
	peakPhiBin = phiBin;
	peakThetaBin = thetaBin;
      }
    }
  }

  peakPhiDeg = hImage->GetXaxis()->GetBinLowEdge(peakPhiBin+1);
  peakThetaDeg = hImage->GetYaxis()->GetBinLowEdge(peakThetaBin+1);

  return hImage;
}

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

  Int_t numPhiBinsThread = hImage->GetNbinsX()/NUM_THREADS;
  Int_t startPhiBin = threadInd*numPhiBinsThread;

  ptr->threadImagePeak[threadInd] = -DBL_MAX;
  Int_t peakPhiBin = -1;
  Int_t peakThetaBin = -1;
  std::vector<Int_t> combosToUse;
  if(mapMode==kTriggered){
    combosToUse = ptr->combosToUseTriggered[ptr->threadL3TrigPattern];
  }

  Int_t phiZoomBase = TMath::Nint(hImage->GetXaxis()->GetBinLowEdge(1)/ZOOM_BIN_SIZE_PHI);
  Int_t thetaZoomBase = TMath::Nint(hImage->GetYaxis()->GetBinLowEdge(1)/ZOOM_BIN_SIZE_THETA + NUM_BINS_THETA_ZOOM_TOTAL/2);

  // if(ptr->kDebug){
  //   if(zoomMode==kZoomedIn){
  //     TThread::Lock();
  //     std::cerr << "BASIC NUMBERS:" << std::endl;
  //     std::cerr << "threadInd = " << threadInd << "\tphiZoomBase = " << phiZoomBase
  // 		<< "\tthetaZoomBase = " << thetaZoomBase
  // 		<< "\thImage->GetXaxis()->GetBinLowEdge(1) = " << hImage->GetXaxis()->GetBinLowEdge(1)
  // 		<< "\thImage->GetYaxis()->GetBinLowEdge(1) = " << hImage->GetYaxis()->GetBinLowEdge(1)
  // 		<< "\tNUM_BINS_PHI_ZOOM_TOTAL = " << NUM_BINS_PHI_ZOOM_TOTAL
  // 		<< "\tNUM_PHI = " << NUM_PHI
  // 		<< "\tPHI_RANGE = " << PHI_RANGE
  // 		<< "\tNUM_BINS_THETA_ZOOM_TOTAL = " << NUM_BINS_THETA_ZOOM_TOTAL
  // 		<< std::endl;
  //     TThread::UnLock();
  //   }
  // }

  if(zoomMode==kZoomedOut && rWave==0){
    for(Int_t phiBin = startPhiBin; phiBin < startPhiBin+numPhiBinsThread; phiBin++){
      Int_t phiSector = zoomMode==kZoomedIn ? 0 : phiBin/NUM_BINS_PHI;
      if(mapMode==kGlobal){
	combosToUse = ptr->combosToUseGlobal[phiSector];
      }
      for(Int_t thetaBin = 0; thetaBin < hImage->GetNbinsY(); thetaBin++){      
	Double_t correlations = 0;
	for(UInt_t comboInd=0; comboInd<combosToUse.size(); comboInd++){
	  Int_t combo = combosToUse.at(comboInd);
	  Int_t offset = 0;

	  offset = ptr->deltaTs[pol][phiBin][thetaBin][combo];
	  offset = offset < 0 ? offset + ptr->numSamples : offset;
	  correlations += ptr->crossCorrelations[pol][combo][offset];
	}
	if(combosToUse.size()>0){
	  correlations /= combosToUse.size();
	}
	hImage->SetBinContent(phiBin + 1, thetaBin + 1, correlations);
	if(correlations > ptr->threadImagePeak[threadInd]){
	  ptr->threadImagePeak[threadInd] = correlations;
	  peakPhiBin = phiBin;
	  peakThetaBin = thetaBin;
	}
      }
    }
  }

  // If we are in zoomed in & plane wave mode then calculate
  // Now using Fast function, this has some intermediate caching of the deltaT calc.
  else if(zoomMode==kZoomedIn && rWave==0){
    
    for(Int_t phiBin = startPhiBin; phiBin < startPhiBin+numPhiBinsThread; phiBin++){
      Int_t phiSector = zoomMode==kZoomedIn ? 0 : phiBin/NUM_BINS_PHI;
      if(mapMode==kGlobal){
	combosToUse = ptr->combosToUseGlobal[phiSector];
      }

      // Take account of wrapping in phi...
      Int_t zoomPhiInd = phiZoomBase + phiBin;
      zoomPhiInd = zoomPhiInd < 0 ? zoomPhiInd + NUM_BINS_PHI_ZOOM_TOTAL : zoomPhiInd;
      zoomPhiInd = zoomPhiInd >= NUM_BINS_PHI_ZOOM_TOTAL ? zoomPhiInd - NUM_BINS_PHI_ZOOM_TOTAL : zoomPhiInd;

      for(Int_t thetaBin = 0; thetaBin < hImage->GetNbinsY(); thetaBin++){
	Int_t zoomThetaInd = thetaZoomBase + thetaBin;
      
	Double_t correlations = 0;
	for(UInt_t comboInd=0; comboInd<combosToUse.size(); comboInd++){
	  Int_t combo = combosToUse.at(comboInd);
	  // Int_t offset = 0;

	  Int_t ant1 = ptr->comboToAnt1s.at(combo);
	  Int_t ant2 = ptr->comboToAnt2s.at(combo);

	  Double_t deltaT = 0;
	  if(ptr->kDebug){

	    if(zoomThetaInd >= 0 && zoomThetaInd < NUM_BINS_THETA_ZOOM_TOTAL){
	      deltaT = ptr->getDeltaTExpectedFast(pol, ant1, ant2,
						  zoomPhiInd,
						  zoomThetaInd);

	      Double_t phiDeg = hImage->GetXaxis()->GetBinLowEdge(phiBin+1);

	      phiDeg = phiDeg < 0 ? phiDeg + 360 : phiDeg;
	      if(TMath::Abs(phiDeg - ptr->zoomedPhiWaveLookup[zoomPhiInd]) > 0.0001){
		TThread::Lock();
		std::cerr.precision(10);
		std::cerr << phiDeg << "\t" << ptr->zoomedPhiWaveLookup[zoomPhiInd] << std::endl;
		TThread::UnLock();
	      }
	    }
	    else{
	      Double_t thetaDeg = hImage->GetYaxis()->GetBinLowEdge(thetaBin+1);
	      Double_t phiDeg = hImage->GetXaxis()->GetBinLowEdge(phiBin+1);
	      Double_t thetaWave = thetaDeg*TMath::DegToRad();
	      Double_t phiWave = phiDeg*TMath::DegToRad();
	      deltaT = ptr->getDeltaTExpected(pol, ant1, ant2, phiWave, thetaWave);
	    }

	    Double_t thetaDeg = hImage->GetYaxis()->GetBinLowEdge(thetaBin+1);
	    Double_t phiDeg = hImage->GetXaxis()->GetBinLowEdge(phiBin+1);

	    // Double_t deltaPhiDeg = phiDeg - hImage->GetXaxis()->GetBinLowEdge(phiBin);
	    Double_t thetaWave = thetaDeg*TMath::DegToRad();
	    Double_t phiWave = phiDeg*TMath::DegToRad();

	    Double_t deltaT_slow = ptr->getDeltaTExpected(pol, ant1, ant2, phiWave, thetaWave);
	    Double_t deltaT_slow2 = ptr->usefulPats[threadInd]->getDeltaTExpected(ant2, ant1, phiWave, -1*thetaWave);
	    Double_t deltaT_slow3 = ptr->getDeltaTExpectedPat(ant2, ant1, phiWave, -1*thetaWave);
	    
	    // if(TMath::Abs(deltaT_slow - deltaT) > 0.0001){
	    if(TMath::Abs(deltaT_slow - deltaT) > 0.0001 || TMath::Abs(deltaT_slow2 - deltaT) > 0.0001
	       || TMath::Abs(deltaT_slow3 - deltaT) > 0.0001){
	      TThread::Lock();	      
	      // std::cerr << threadInd << "\t" << ant1 << "\t" << ant2 << "\t"
	      // 		<< phiWave << "\t" << thetaWave << "\t" << deltaT << "\t" 
	      // 		<< deltaT_slow << "\t" << deltaT_slow2 << "\t" << deltaT_slow3 << std::endl;


	      
	      Double_t deltaT_slow2b = ptr->usefulPats[threadInd]->getDeltaTExpected(ant2, ant1, phiWave, -1*thetaWave, 1);
	      Double_t deltaT_slow3b = ptr->getDeltaTExpectedPat(ant2, ant1, phiWave, -1*thetaWave, 1);

	      std::cerr << threadInd << "\t" << ant1 << "\t" << ant2 << "\t"
			<< phiWave << "\t" << thetaWave << std::endl;
	      std::cerr << "deltaT       = " << deltaT << std::endl;
	      std::cerr << "deltaT_slow  = " << deltaT_slow << std::endl;
	      std::cerr << "deltaT_slow2 = " << deltaT_slow2 << std::endl;
	      std::cerr << "deltaT_slow3 = " << deltaT_slow3 << std::endl;
	      std::cerr << "deltaT_slow2b = " << deltaT_slow2b << std::endl;
	      std::cerr << "deltaT_slow3b = " << deltaT_slow3b << std::endl;
	      std::cerr << std::endl;


	      // std::cerr << "thetaDeg = " << thetaDeg      
	      // 		<< "\tzoomThetaInd = " << zoomThetaInd << std::endl;
	      // std::cerr << "phiDeg = " << phiDeg
	      // 		<< "\tzoomPhiInd = " << zoomPhiInd << std::endl;
	      // std::cerr << "zoomedThetaWaves[zoomThetaInd] " << ptr->zoomedThetaWaves[zoomThetaInd] << std::endl;
	      // std::cerr << "TMath::RadToDeg()*zoomedThetaWaves[zoomThetaInd] "
	      // 		<< TMath::RadToDeg()*ptr->zoomedThetaWaves[zoomThetaInd] << std::endl;
	      // std::cerr << deltaT_slow2 << std::endl;
	      TThread::UnLock();	      
	    }
	  }
	  else{
	    if(zoomThetaInd >= 0 && zoomThetaInd < NUM_BINS_THETA_ZOOM_TOTAL){
	      deltaT = ptr->getDeltaTExpectedFast(pol, ant1, ant2,
						  zoomPhiInd,
						  zoomThetaInd);
	    }
	    else{
	      Double_t thetaDeg = hImage->GetYaxis()->GetBinLowEdge(thetaBin+1);
	      Double_t phiDeg = hImage->GetXaxis()->GetBinLowEdge(phiBin+1);
	      Double_t thetaWave = thetaDeg*TMath::DegToRad();
	      Double_t phiWave = phiDeg*TMath::DegToRad();
	      deltaT = ptr->getDeltaTExpected(pol, ant1, ant2, phiWave, thetaWave);
	    }
	  }
	  // offset = TMath::Nint(deltaT/ptr->correlationDeltaT);

	  Int_t offsetLow = floor(deltaT/ptr->correlationDeltaT);
	  Int_t offsetHigh = ceil(deltaT/ptr->correlationDeltaT);

	  Double_t dt1 = offsetLow*ptr->correlationDeltaT;
	  // Double_t dt2 = offsetHigh*ptr->correlationDeltaT;
	  
	  offsetLow = offsetLow < 0 ? offsetLow + ptr->numSamplesUpsampled : offsetLow;
	  offsetHigh = offsetHigh < 0 ? offsetHigh + ptr->numSamplesUpsampled : offsetHigh;	  
	  
	  Double_t c1 = ptr->crossCorrelationsUpsampled[pol][combo][offsetLow];
	  Double_t c2 = ptr->crossCorrelationsUpsampled[pol][combo][offsetHigh];

	  Double_t cInterp = (deltaT - dt1)*(c2 - c1)/(ptr->correlationDeltaT) + c1;

	  // TThread::Lock();	      	  
	  // std::cerr << offsetLow << "\t" << offsetHigh << std::endl;
	  // std::cerr << dt1 << "\t" << dt2 << "\t" << deltaT << std::endl;
	  // std::cerr << c1 << "\t" << c2 << "\t" << cInterp << std::endl;	  
	  // TThread::UnLock();	      

	  // offset = offset < 0 ? offset + ptr->numSamplesUpsampled : offset;

	  // correlations += ptr->crossCorrelationsUpsampled[pol][combo][offset];
	  correlations += cInterp;
	}
	if(combosToUse.size()>0){
	  correlations /= combosToUse.size();
	}
	hImage->SetBinContent(phiBin + 1, thetaBin + 1, correlations);
	if(correlations > ptr->threadImagePeak[threadInd]){
	  ptr->threadImagePeak[threadInd] = correlations;
	  peakPhiBin = phiBin;
	  peakThetaBin = thetaBin;
	}
      }
    }
  }

	
  // If we are in zoomed out & spherical wave mode then calculate
  else{
    for(Int_t phiBin = startPhiBin; phiBin < startPhiBin+numPhiBinsThread; phiBin++){
      Int_t phiSector = phiBin/NUM_BINS_PHI;
      if(mapMode==kGlobal){
	combosToUse = ptr->combosToUseGlobal[phiSector];
      }
    
      Double_t phiWave = hImage->GetXaxis()->GetBinLowEdge(phiBin+1)*TMath::DegToRad();

      for(Int_t thetaBin = 0; thetaBin < hImage->GetNbinsY(); thetaBin++){
	Double_t thetaWave = hImage->GetYaxis()->GetBinLowEdge(thetaBin+1)*TMath::DegToRad();

	Double_t correlations = 0;
	for(UInt_t comboInd=0; comboInd<combosToUse.size(); comboInd++){
	  Int_t combo = combosToUse.at(comboInd);
	  Int_t offset = 0;
    
	  Int_t ant1 = ptr->comboToAnt1s.at(combo);
	  Int_t ant2 = ptr->comboToAnt2s.at(combo);
	  offset = ptr->getDeltaTExpectedSpherical(pol, ant1, ant2, phiWave, thetaWave, rWave);
	  offset = offset < 0 ? offset + ptr->numSamples : offset;
	  correlations += ptr->crossCorrelations[pol][combo][offset];
	}
	if(combosToUse.size()>0){
	  correlations /= combosToUse.size();
	}
	hImage->SetBinContent(phiBin + 1, thetaBin + 1, correlations);
	if(correlations > ptr->threadImagePeak[threadInd]){
	  ptr->threadImagePeak[threadInd] = correlations;
	  peakPhiBin = phiBin;
	  peakThetaBin = thetaBin;
	}	
      }
    }
  }

  ptr->threadPeakPhiDeg[threadInd] = hImage->GetXaxis()->GetBinLowEdge(peakPhiBin+1);
  ptr->threadPeakThetaDeg[threadInd] = hImage->GetYaxis()->GetBinLowEdge(peakThetaBin+1);

  return 0;
  
}




Double_t CrossCorrelator::findImagePeak(TH2D* hist, Double_t& imagePeak,
					Double_t& imagePeakTheta, Double_t& imagePeakPhi){

  Int_t nx = hist->GetNbinsX();
  Int_t ny = hist->GetNbinsY();

  imagePeak = -DBL_MAX;
  for(Int_t by = 1; by<=ny; by++){
    for(Int_t bx = 1; bx<=nx; bx++){
      Double_t val = hist->GetBinContent(bx, by);
      if(val > imagePeak){
	imagePeak = val;
	imagePeakPhi = hist->GetXaxis()->GetBinLowEdge(bx);
	imagePeakTheta = hist->GetYaxis()->GetBinLowEdge(by);
      }
    }
  }

  return imagePeak;
  
}





/************************************************************************************************************
Functions to delete pointers to internal variables
************************************************************************************************************/
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















/************************************************************************************************************
Functions for debugging or testing
************************************************************************************************************/


/*!
  \brief Used to validate that I am testing the geometry I think I am. Returns 0 on success and 1 on failure.
  \param pathToLindasFile is the relative path to Linda's file.
  \param pol is the polarization of the channels under test.
*/
Int_t CrossCorrelator::validateGeometry(TString pathToLindasFile, AnitaPol::AnitaPol_t pol){
  
  // Since I am simulataneously testing many of Linda's geometries on lots of different files
  // I need the help of a machine to check I'm testing the geometry I think I'm testing.
  
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  AnitaEventCalibrator* cal = AnitaEventCalibrator::Instance();


  std::ifstream lindasNums(pathToLindasFile.Data());

  Int_t ant;
  Double_t dr, dPhiRad, dz, dt;

  Double_t sumOfErrors = 0;
  while(lindasNums >> ant >> dr >> dz >> dPhiRad >> dt){
    
    Int_t surf, chan, ant2;
    geom->getSurfChanAntFromRingPhiPol(AnitaRing::AnitaRing_t (ant/NUM_PHI), ant%NUM_PHI, pol,
				       surf, chan, ant2);

    Double_t newR = geom->rPhaseCentreFromVerticalHornKurtAnita3[ant][pol] + dr;
    Double_t newPhi = geom->azPhaseCentreFromVerticalHornKurtAnita3[ant][pol] + dPhiRad;
    Double_t newZ = geom->zPhaseCentreFromVerticalHornKurtAnita3[ant][pol] + dz;
    Double_t newT = dt;
    
    Double_t geomR = geom->rPhaseCentreFromVerticalHorn[ant][pol];
    Double_t geomPhi = geom->azPhaseCentreFromVerticalHorn[ant][pol];
    Double_t geomZ = geom->zPhaseCentreFromVerticalHorn[ant][pol];
    Double_t calT = cal->relativePhaseCenterToAmpaDelays[surf][chan];

    Double_t ddr = geomR - newR;
    Double_t ddPhi = geomPhi - newPhi;
    if(ddPhi >= TMath::Pi()){
      ddPhi -= TMath::TwoPi();
    }
    else if(ddPhi < -TMath::Pi()){
      ddPhi += TMath::TwoPi();
    }
    Double_t ddZ = geomZ - newZ;
    Double_t ddT = calT - newT;
           
    // std::cerr << ant << "\t";
    // std::cerr << ddr << "\t";
    // std::cerr << ddPhi << "\t";
    // std::cerr << ddZ << "\t";
    // std::cerr << ddT << "\t";
    // std::cerr << std::endl;

    sumOfErrors += TMath::Abs(ddr);
    sumOfErrors += TMath::Abs(ddPhi);
    sumOfErrors += TMath::Abs(ddZ);
    sumOfErrors += TMath::Abs(ddT);
    
  }

  const Double_t errorTolerance = 1e-5;
  if(sumOfErrors < errorTolerance){
    return 0;
  }
  else{
    return 1;
  }
}


/*!
  \brief Used to validate that I am testing the geometry I think I am. Returns 1 if could not open file.
  \param pathToLindasFile is the relative path to Linda's file.
  \param pol is the polarization of the channels under test.
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

  Double_t sumPhiRad = 0;
  
  while(lindasNums >> ant >> dr >> dz >> dPhiRad >> dt){

    sumPhiRad += dPhiRad;
    
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

    // if(ant==47){
    //   std::cout << "inside CC" << "\t";
    //   std::cout << geom->azPhaseCentreFromVerticalHornKurtAnita3[ant][pol] << "\t"
    // 		<< dPhiRad << "\t" << geom->azPhaseCentreFromVerticalHorn[ant][pol]
    // 		<< "\tNOT " << geom->azPhaseCentreFromVerticalHornKurtAnita3[ant][pol] - dPhiRad
    // 		<< std::endl;
    //   std::cout << "geom = " << geom << std::endl;
    // }
  }

  std::cerr << "Average phi-shift = " << TMath::RadToDeg()*(sumPhiRad / NUM_SEAVEYS) << std::endl;
  
  return 0;
}



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

TGraph* CrossCorrelator::getCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2){
  // Primarily for debugging, put cross correlations in a TGraph 
  return getCrossCorrelationGraphWorker(numSamples, pol, ant1, ant2);
}




TGraph* CrossCorrelator::getUpsampledCrossCorrelationGraph(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2){
  // Primarily for debugging, put cross correlations in a TGraph 

  return getCrossCorrelationGraphWorker(numSamplesUpsampled, pol, ant1, ant2);
}


TGraph* CrossCorrelator::makeCoherentlySummedWaveform(AnitaPol::AnitaPol_t pol, Double_t phiDeg, Double_t thetaDeg, UInt_t l3Trigger){
  // Sums all l3 triggered phi-sector waveforms together to make the coherently summed waveform

  Double_t phiRad = phiDeg*TMath::DegToRad();
  Double_t thetaRad = thetaDeg*TMath::DegToRad();
  Int_t numAnts = 0;
  Int_t firstAnt = -1;

  std::pair<Int_t, Int_t> key(numSamplesUpsampled, 0);
  TGraph* grCoherent = NULL;
  Double_t* vArray = FancyFFTs::fReals[key];
  
  for(Int_t phiSector=0; phiSector<NUM_PHI; phiSector++){
    Int_t doPhiSector = RootTools::getBit(phiSector, l3Trigger);
    if(doPhiSector > 0){
      for(Int_t ring=0; ring<NUM_RING; ring++){
	Int_t ant= phiSector + ring*NUM_PHI;
	FancyFFTs::doInvFFT(numSamplesUpsampled, fftsPadded[pol][ant], false);

	if(firstAnt==-1){
	  std::vector<Double_t> tArray(numSamplesUpsampled, 0);
	  Double_t t0 = grsResampled[pol][ant]->GetX()[0];
	  for(Int_t samp=0; samp<numSamplesUpsampled; samp++){
	    tArray.at(samp) = t0 + samp*correlationDeltaT;
	    vArray[samp] *= interpRMS[pol][ant];
	  }
	  firstAnt = ant;
	  grCoherent = new TGraph(numSamplesUpsampled, &tArray[0], &vArray[0]);
	}
	else{
	  Double_t deltaT = getDeltaTExpected(pol, firstAnt, ant, phiRad, thetaRad);
	  Int_t offset = TMath::Nint(deltaT/correlationDeltaT);
	  for(Int_t samp=0; samp<numSamplesUpsampled; samp++){
	    Int_t samp2 = samp + offset;
	    if(samp2 >= 0 && samp2 < numSamplesUpsampled){
	      grCoherent->GetY()[samp] += vArray[samp2]*interpRMS[pol][ant];
	    }
	  }
	}
	numAnts++;
      }
    }
  }

  // Normalize
  if(numAnts > 0){
    TString name;
    TString title;
    for(Int_t samp=0; samp<numSamplesUpsampled; samp++){
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

    title += TString::Format("Coherently Summed Waveform for arrival direction elevation %lf (Deg) and azimuth %lf (Deg); Time (ns); Voltage (mV)", thetaDeg, phiDeg);

    grCoherent->SetName(name);
    grCoherent->SetTitle(title);    
  }
  
  return grCoherent;
  
}
