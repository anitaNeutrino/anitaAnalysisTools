/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             A Cross Correlator to interact with the ROOTified ANITA-3 data and do some interferometry. 
*************************************************************************************************************** */

#include "CrossCorrelator.h"

ClassImp(CrossCorrelator)


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
    deleteCrossCorrelations((AnitaPol::AnitaPol_t)pol);
    deleteAllFFTs((AnitaPol::AnitaPol_t)pol);    
  }
}


/*!
  \brief Workhorse function to set internal variables.
*/
void CrossCorrelator::initializeVariables(){

  // Initialize with NULL otherwise very bad things will happen with gcc 
  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      grs[pol][ant] = NULL;
      grsInterp[pol][ant] = NULL;
      interpRMS[pol][ant] = 0;
      ffts[pol][ant] = NULL;
      fftsPadded[pol][ant] = NULL;      
    }
    for(int combo=0; combo<NUM_COMBOS; combo++){
      crossCorrelations[pol][combo] = NULL;
      crossCorrelationsUpsampled[pol][combo] = NULL;      
    }
    lastEventNormalized[pol] = 0;
    eventNumber[pol] = 0;    
  }


  upsampleFactor = 40; // 0.38ns -> ~10ps

  numSamples = 2*NUM_SAMPLES; // Factor of two for padding 
  numSamplesUpsampled = numSamples*upsampleFactor; // For upsampling
  
  nominalSamplingDeltaT = 1./2.6;
  correlationDeltaT = nominalSamplingDeltaT/upsampleFactor;


  deltaTMax = 0;
  deltaTMin = 0;

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  for(int ant=0; ant<NUM_SEAVEYS; ant++){
    rArray.push_back(geom->getAntR(ant));
    zArray.push_back(geom->getAntZ(ant));
    phiArrayDeg.push_back(geom->getAntPhiPositionRelToAftFore(ant)*TMath::RadToDeg());
  }

  do5PhiSectorCombinatorics();
  fillDeltaTLookup();
  printInfo();

  mapModeNames[kGlobal] = "Global";
  mapModeNames[kTriggered] = "Triggered";
  zoomModeNames[kZoomedOut] = "";
  zoomModeNames[kZoomedIn] = "Zoom";

  for(Long_t threadInd=0; threadInd<NUM_THREADS; threadInd++){
    CrossCorrelator::threadArgs threadArgVals;
    threadArgVals.threadInd = threadInd;
    threadArgVals.ptr = this;
    threadArgsVec.push_back(threadArgVals);

    FancyFFTs::makeNewPlanIfNeeded(numSamples, threadInd);
    FancyFFTs::makeNewPlanIfNeeded(numSamplesUpsampled, threadInd);
  }

  for(Long_t threadInd=0; threadInd<NUM_THREADS; threadInd++){
    
    TString name = TString::Format("threadCorr%ld", threadInd);
    corrThreads.push_back(new TThread(name.Data(),
				      CrossCorrelator::doSomeCrossCorrelationsThreaded,
				      (void*)&threadArgsVec.at(threadInd))
			  );

    // name = TString::Format("threadUpsampledCorr%ld", threadInd);
    // upsampledCorrThreads.push_back(new TThread(name.Data(),
    // 					       CrossCorrelator::doSomeUpsampledCrossCorrelationsThreaded,
    // 					       (void*)&threadArgsVec.at(threadInd))
    // 				   );
    
    name = TString::Format("threadMap%ld", threadInd);
    mapThreads.push_back(new TThread(name.Data(),
				     CrossCorrelator::makeSomeOfImageThreaded,
				     (void*)&threadArgsVec.at(threadInd))
			 );
  }

}



void CrossCorrelator::printInfo(){
  std::cerr << __PRETTY_FUNCTION__ << std::endl;
  std::cerr << "\tupsample factor = " << upsampleFactor << std::endl;
  std::cerr << "\tdeltaT max = " << deltaTMax << std::endl;
  std::cerr << "\tdeltaT min = " << deltaTMin << std::endl;
  std::cerr << "\tBin size theta (deg) = " << Double_t(THETA_RANGE)/NUM_BINS_THETA << std::endl;
  std::cerr << "\tBin size phi (deg) = " << Double_t(PHI_RANGE)/NUM_BINS_PHI << std::endl;
  std::cerr << "\tdeltaTs array size = " << sizeof(dtIndex_t)*numCombos*NUM_PHI*NUM_BINS_PHI*NUM_BINS_THETA << " bytes" << std::endl;
}






/*!
  \brief Single function to get the angle of the first bin of the interferometric histogram.
  \returns position of antenna 0 in ADU5Pat coordinates, offset by half a phi-sector.
*/

Double_t CrossCorrelator::getBin0PhiDeg(){
  Double_t phi0 = phiArrayDeg.at(0);
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
void CrossCorrelator::getNormalizedInterpolatedTGraphs(UsefulAnitaEvent* usefulEvent, AnitaPol::AnitaPol_t pol){
  // Potentially needed in a few places, so it gets its own function 

  // Pretty much just for profiling 
  // if(usefulEvent->eventNumber!=lastEventNormalized[pol]){

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
      grsInterp[pol][ant] = interpolateWithStartTime(grs[pol][ant], earliestStart.at(pol));
      Double_t mean; // don't care enough about this to store it anywhere.
      RootTools::normalize(grsInterp[pol][ant], mean, interpRMS[pol][ant]);
    }
    lastEventNormalized[pol] = usefulEvent->eventNumber;
  // }
}


void CrossCorrelator::doFFTs(AnitaPol::AnitaPol_t pol){

  deleteAllFFTs(pol); // First delete any FFTs that might still be knocking around

  // Generate a new set of FFTs.
  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
    ffts[pol][ant] = FancyFFTs::doFFT(numSamples, grsInterp[pol][ant]->GetY(), true);    
  }
}

/*!
  \brief Wrapper function for ROOT's interpolator, can zero pad the front to start from a particular time.
  \param grIn points to the TGraph containing the waveform to interpolate
  \param startTime is the start time for interpolation: zero pads if this is earlier than the TGraph start time.
*/
TGraph* CrossCorrelator::interpolateWithStartTime(TGraph* grIn, Double_t startTime){

  // std::vector<Double_t> newTimes = std::vector<Double_t>(numSamplesUpsampled, 0);
  // std::vector<Double_t> newVolts = std::vector<Double_t>(numSamplesUpsampled, 0);
  std::vector<Double_t> newTimes = std::vector<Double_t>(numSamples, 0);
  std::vector<Double_t> newVolts = std::vector<Double_t>(numSamples, 0);    
  Double_t thisStartTime = grIn->GetX()[0];
  Double_t lastTime = grIn->GetX()[grIn->GetN()-1];


  // Quantizes the start and end times so data poInt_ts lie at Int_teger multiples of nominal sampling 
  // startTime = correlationDeltaT*TMath::Nint(startTime/correlationDeltaT + 0.5);
  // lastTime = correlationDeltaT*TMath::Nint(lastTime/correlationDeltaT - 0.5);
  startTime = nominalSamplingDeltaT*TMath::Nint(startTime/nominalSamplingDeltaT + 0.5);
  lastTime = nominalSamplingDeltaT*TMath::Nint(lastTime/nominalSamplingDeltaT - 0.5);  

   //ROOT Int_terpolator object constructor takes std::vector objects
  std::vector<Double_t> tVec(grIn->GetX(), grIn->GetX() + grIn->GetN());
  std::vector<Double_t> vVec(grIn->GetY(), grIn->GetY() + grIn->GetN());
  
  // This is ROOT's Int_terpolator object
  ROOT::Math::Interpolator chanInterp(tVec,vVec,ROOT::Math::Interpolation::kAKIMA);
  
  // Put new data Int_to arrays
  Double_t time = startTime;
  // for(Int_t samp = 0; samp < numSamplesUpsampled; samp++){
  for(Int_t samp = 0; samp < numSamples; samp++){    
    newTimes.at(samp) = time;
    if(time >= thisStartTime && time <= lastTime){
      newVolts.at(samp) = chanInterp.Eval(time);
    }
    else{
      newVolts.at(samp) = 0;
    }
    time += nominalSamplingDeltaT;
  }

  // return new TGraph(numSamplesUpsampled, &newTimes[0], &newVolts[0]);
  return new TGraph(numSamples, &newTimes[0], &newVolts[0]);  

}





















/************************************************************************************************************
All correlation functions
************************************************************************************************************/


/*!
  \brief Loops through the set of cross correlations and returns a vector of vectors containing the maximum times for the specified polarization.
  \returns A vector of length numCombos containing the set of times of maximum correlation.
*/
std::vector<Double_t> CrossCorrelator::getMaxCorrelationTimes(AnitaPol::AnitaPol_t pol){
  
  std::vector<Double_t> peakTimes(NUM_COMBOS, 0);

  for(Int_t combo=0; combo<NUM_COMBOS; combo++){
    Int_t maxInd = RootTools::getIndexOfMaximum(numSamples, crossCorrelations[pol][combo]);
    Int_t offset = maxInd >= numSamples/2 ? maxInd - numSamples : maxInd;
    peakTimes.at(combo) = offset*nominalSamplingDeltaT;      
  }
  return peakTimes;
}

/*!
  \brief Loops through the set of cross correlations and returns a vector of vectors containing the maximum times.
  \returns A vector length 2 of vector length numCombos containing the set of times of maximum correlation.
*/
std::vector<std::vector<Double_t> > CrossCorrelator::getMaxCorrelationTimes(){
  
  std::vector<std::vector<Double_t> > peakTimes;
  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    peakTimes.push_back(getMaxCorrelationTimes((AnitaPol::AnitaPol_t)pol));
  }
  return peakTimes;
}


/*!
  \brief Loops through the set of cross correlations and returns a vector of vectors containing the maximum correlation values for the specified polarization.
  \returns A vector of length numCombos containing the set of maximum correlation values.
*/
std::vector<Double_t> CrossCorrelator::getMaxCorrelationValues(AnitaPol::AnitaPol_t pol){
  
  std::vector<Double_t> peakVals(NUM_COMBOS, 0);
  for(Int_t combo=0; combo<NUM_COMBOS; combo++){
    Int_t maxInd = RootTools::getIndexOfMaximum(numSamples, crossCorrelations[pol][combo]);
    peakVals.at(combo) = crossCorrelations[pol][combo][maxInd];
  }
  return peakVals;
}


/*!
  \brief Loops through the set of cross correlations and returns a vector of vectors containing the maximum correlation values.
  \returns A vector length 2 of vector length numCombos containing the set of maximum correlation values.
*/
std::vector<std::vector<Double_t> > CrossCorrelator::getMaxCorrelationValues(){
  
  std::vector<std::vector<Double_t> > peakVals;
  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    peakVals.push_back(getMaxCorrelationValues((AnitaPol::AnitaPol_t)pol));
  }
  return peakVals;
}



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
  deleteCrossCorrelations(pol);

  // Loop over combinations generating set of cross correlations.
  for(Int_t combo=0; combo<numCombos; combo++){
    Int_t ant1 = comboToAnt1s.at(combo);
    Int_t ant2 = comboToAnt2s.at(combo);
    // crossCorrelations[pol][combo] = crossCorrelateFourier(grsInterp[pol][ant2], grsInterp[pol][ant1]);
    crossCorrelations[pol][combo] = crossCorrelateFourier(numSamples, ffts[pol][ant2], ffts[pol][ant1]);
  }
}


/*!
  \brief Loop over both polarizations and all combinations, and generate set of cross correlations.
*/

void CrossCorrelator::doAllCrossCorrelationsThreaded(AnitaPol::AnitaPol_t pol){

  // Delete old cross correlations first 
  deleteCrossCorrelations(pol);

  // Set variable for use in threads
  threadPol = pol;
  
  for(long threadInd=0; threadInd<NUM_THREADS; threadInd++){
    corrThreads.at(threadInd)->Run();
  }

  for(long threadInd=0; threadInd<NUM_THREADS; threadInd++){
    corrThreads.at(threadInd)->Join();
  }
}

void* CrossCorrelator::doSomeCrossCorrelationsThreaded(void* voidPtrArgs){

  // Disgusting hacks to get ROOT threading to compile inside a class.
  CrossCorrelator::threadArgs* args = (CrossCorrelator::threadArgs*) voidPtrArgs;
  Long_t threadInd = args->threadInd;
  CrossCorrelator* ptr = args->ptr;
  AnitaPol::AnitaPol_t pol = ptr->threadPol;


  
  Int_t numCorrPerThread = NUM_COMBOS/NUM_THREADS;
  Int_t startCombo = threadInd*numCorrPerThread;
  // TThread::Lock();
  // std::cout << threadInd << "\t" << ptr << "\t" << pol << "\t" << startCombo << std::endl;
  // TThread::UnLock();

  
  for(int combo=startCombo; combo<startCombo+numCorrPerThread; combo++){
    Int_t ant1 = ptr->comboToAnt1s.at(combo);
    Int_t ant2 = ptr->comboToAnt2s.at(combo);
    
    ptr->crossCorrelations[pol][combo] = ptr->crossCorrelateFourier(ptr->numSamples,
								    ptr->ffts[pol][ant2],
								    ptr->ffts[pol][ant1],
								    threadInd);
    // TThread::Lock();
    // std::cout << "threadInd = " << threadInd << "\tcombo = " << combo
    // 	      << "\tant1 = " << ant1 << "\tant2 = " << ant2 << std::endl;
    // for(int i=0; i<5; i++){
    //   std::cout << ptr->crossCorrelations[pol][combo][i] << ", ";
    // }
    // std::cout << std::endl << std::endl;
    // TThread::UnLock();
    
  }
  
  return 0;
}

void CrossCorrelator::correlateForZoomedImage(AnitaPol::AnitaPol_t pol, UInt_t l3TrigPattern){
  // std::cerr << __PRETTY_FUNCTION__ << std::endl;

  deleteUpsampledCrossCorrelations(pol);
  deleteAllPaddedFFTs(pol);

  std::vector<Int_t> combosToUse = combosToUseTriggered[l3TrigPattern];
  
  Int_t numFreqs = FancyFFTs::getNumFreqs(numSamples);
  Int_t numFreqsPadded = FancyFFTs::getNumFreqs(numSamplesUpsampled);

  for(int phiSector=0; phiSector<NUM_PHI; phiSector++){
    Int_t numComboInd = combosToUse.size();
    for(Int_t comboInd=0; comboInd < numComboInd; comboInd++){
      Int_t combo=combosToUse.at(comboInd);

      if(crossCorrelationsUpsampled[pol][combo]==NULL){

	Int_t ant1 = comboToAnt1s.at(combo);
	Int_t ant2 = comboToAnt2s.at(combo);

	// std::cerr << combo << "\t" << ant1 << "\t" << fftsPadded[pol][ant1] << "\t"
	// 	  << ant2 << "\t" << fftsPadded[pol][ant2] << std::endl;
      
	if(fftsPadded[pol][ant1]==NULL){
	  // std::cerr << pol << "\t" << ant1 << std::endl;
	  fftsPadded[pol][ant1] = FancyFFTs::zeroPadFFT(ffts[pol][ant1], numFreqs, numFreqsPadded);
	  // std::cerr << pol << "\t" << ant1 << "\t" << fftsPadded[pol][ant1] << std::endl;
	}
	if(fftsPadded[pol][ant2]==NULL){
	  // std::cerr << pol << "\t" << ant2 << std::endl;
	  fftsPadded[pol][ant2] = FancyFFTs::zeroPadFFT(ffts[pol][ant2], numFreqs, numFreqsPadded);
	  // std::cerr << pol << "\t" << ant2 << "\t" << fftsPadded[pol][ant2] << std::endl;
	}

	crossCorrelationsUpsampled[pol][combo] = crossCorrelateFourier(numSamplesUpsampled,
								       fftsPadded[pol][ant2],
								       fftsPadded[pol][ant1]);
      }
    }
  }
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
// Double_t* CrossCorrelator::crossCorrelateFourier(FFTWComplex* fft1, FFTWComplex* fft2){
Double_t* CrossCorrelator::crossCorrelateFourier(Int_t numSamplesInTimeDomain,
						 std::complex<Double_t>* fft1,
						 std::complex<Double_t>* fft2,
						 Int_t threadInd){  
  return FancyFFTs::crossCorrelate(numSamplesInTimeDomain, fft1, fft2, threadInd);  
}











/************************************************************************************************************
Calculate deltaT between two antennas (for a plane wave unless function name says otherwise)
************************************************************************************************************/

Double_t CrossCorrelator::getDeltaTExpected(Int_t ant1, Int_t ant2, Double_t phiWave, Double_t thetaWave){
  Double_t tanThetaW = tan(-thetaWave);
  Double_t part1 = zArray.at(ant1)*tanThetaW - rArray.at(ant1)*cos(phiWave-TMath::DegToRad()*phiArrayDeg.at(ant1));
  Double_t part2 = zArray.at(ant2)*tanThetaW - rArray.at(ant2)*cos(phiWave-TMath::DegToRad()*phiArrayDeg.at(ant2));
  Double_t tdiff = 1e9*((cos(-thetaWave) * (part2 - part1))/SPEED_OF_LIGHT); // Returns time in ns
  return tdiff;
}


Int_t CrossCorrelator::getDeltaTExpectedSpherical(Int_t ant1, Int_t ant2, Double_t phiWave, Double_t thetaWave, Double_t rWave){

  Double_t phi1 = TMath::DegToRad()*phiArrayDeg.at(ant1);
  Double_t x1 = rArray.at(ant1)*TMath::Cos(phi1);
  Double_t y1 = rArray.at(ant1)*TMath::Sin(phi1);
  Double_t z1 = zArray.at(ant1);

  Double_t phi2 = TMath::DegToRad()*phiArrayDeg.at(ant2);
  Double_t x2 = rArray.at(ant2)*TMath::Cos(phi2);
  Double_t y2 = rArray.at(ant2)*TMath::Sin(phi2);
  Double_t z2 = zArray.at(ant2);

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

void CrossCorrelator::fillDeltaTLookupZoomed(Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg, UInt_t l3TrigPattern){

  std::vector<Int_t> combosToUse = combosToUseTriggered[l3TrigPattern];
  Double_t phiMin = zoomCenterPhiDeg - PHI_RANGE_ZOOM/2;
  Double_t thetaMin = zoomCenterThetaDeg - THETA_RANGE_ZOOM/2;

  for(Int_t phiBin = 0; phiBin < NUM_BINS_PHI_ZOOM; phiBin++){
    Double_t phiDeg = phiMin + phiBin*Double_t(PHI_RANGE_ZOOM)/NUM_BINS_PHI_ZOOM;
    Double_t phiWave = TMath::DegToRad()*phiDeg;
    for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA_ZOOM; thetaBin++){
      Double_t thetaDeg = thetaMin + THETA_RANGE_ZOOM*((Double_t)thetaBin/NUM_BINS_THETA_ZOOM);
      Double_t thetaWave = TMath::DegToRad()*thetaDeg;
      for(UInt_t comboInd=0; comboInd<combosToUse.size(); comboInd++){
	Int_t combo = combosToUse.at(comboInd);
	Int_t ant1 = comboToAnt1s.at(combo);
	Int_t ant2 = comboToAnt2s.at(combo);
	Double_t deltaT = getDeltaTExpected(ant1, ant2, phiWave, thetaWave);
	Int_t offset = TMath::Nint(deltaT/correlationDeltaT);
	deltaTsZoom[phiBin][thetaBin][combo] = offset;
      }
    }
  }
}

void CrossCorrelator::fillDeltaTLookup(){
  
  Double_t phi0 = getBin0PhiDeg();
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
	  Double_t deltaT = getDeltaTExpected(ant1, ant2, phiWave, thetaWave);
	  Int_t offset = TMath::Nint(deltaT/nominalSamplingDeltaT);
	  if(offset > deltaTMax){
	    deltaTMax = offset;
	  }
	  if(offset < deltaTMin){
	    deltaTMin = offset;
	  }
	  deltaTs[phiBin][thetaBin][combo] = offset;
	}
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
    
void CrossCorrelator::fillCombosToUseIfNeeded(mapMode_t mapMode, UInt_t l3TrigPattern){
  
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

  Double_t phiMin = zoomCenterPhiDeg - PHI_RANGE_ZOOM/2;
  Double_t phiMax = zoomCenterPhiDeg + PHI_RANGE_ZOOM/2;
  Double_t thetaMin = zoomCenterThetaDeg - THETA_RANGE_ZOOM/2;
  Double_t thetaMax = zoomCenterThetaDeg + THETA_RANGE_ZOOM/2;

  TH2D* hImage = new TH2D(name, title,
			  NUM_BINS_PHI_ZOOM, phiMin, phiMax,
			  NUM_BINS_THETA_ZOOM, thetaMin, thetaMax);
  hImage->GetXaxis()->SetTitle("Azimuth (Degrees)");
  hImage->GetYaxis()->SetTitle("Elevation (Degrees)");
  return hImage;

}


TH2D* CrossCorrelator::makeGlobalSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave){
  Double_t imagePeak, peakPhiDeg, peakThetaDeg;
  return makeGlobalSphericalImage(pol, rWave, imagePeak, peakPhiDeg, peakThetaDeg);
}

TH2D* CrossCorrelator::makeGlobalSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave, Double_t& imagePeak,
						Double_t& peakPhiDeg, Double_t& peakThetaDeg){
  return makeImage(pol, rWave, imagePeak, peakPhiDeg, peakThetaDeg, ALL_PHI_TRIGS,
		   kGlobal, kZoomedOut);
}


TH2D* CrossCorrelator::makeTriggeredSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave, UInt_t l3TrigPattern){
  Double_t imagePeak, peakPhiDeg, peakThetaDeg;
  return makeTriggeredSphericalImage(pol, rWave, imagePeak, peakPhiDeg, peakThetaDeg, l3TrigPattern);
}

TH2D* CrossCorrelator::makeTriggeredSphericalImage(AnitaPol::AnitaPol_t pol, Double_t rWave,
						   Double_t& imagePeak, Double_t& peakPhiDeg,
						   Double_t& peakThetaDeg, UInt_t l3TrigPattern){
  return makeImage(pol, rWave, imagePeak, peakPhiDeg, peakThetaDeg, l3TrigPattern, kTriggered, kZoomedOut);
}


TH2D* CrossCorrelator::makeGlobalImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
				       Double_t& peakThetaDeg){
  return makeImage(pol, 0, imagePeak, peakPhiDeg, peakThetaDeg, ALL_PHI_TRIGS,
		   kGlobal, kZoomedOut);
}


TH2D* CrossCorrelator::makeTriggeredImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
					  Double_t& peakThetaDeg, UInt_t l3TrigPattern){
  return makeImage(pol, 0, imagePeak, peakPhiDeg, peakThetaDeg, l3TrigPattern, kTriggered, kZoomedOut);
}


TH2D* CrossCorrelator::makeZoomedImage(AnitaPol::AnitaPol_t pol, Double_t& imagePeak, Double_t& peakPhiDeg,
				       Double_t& peakThetaDeg, UInt_t l3TrigPattern,
				       Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg){

  return makeImage(pol, 0, imagePeak, peakPhiDeg, peakThetaDeg, l3TrigPattern,
		   kTriggered, kZoomedIn, zoomCenterPhiDeg, zoomCenterThetaDeg);

}

TH2D* CrossCorrelator::makeZoomedImage(AnitaPol::AnitaPol_t pol, UInt_t l3TrigPattern,
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

TH2D* CrossCorrelator::makeImage(AnitaPol::AnitaPol_t pol, Double_t rWave, Double_t& imagePeak,
				 Double_t& peakPhiDeg, Double_t& peakThetaDeg, UInt_t l3TrigPattern,
				 mapMode_t mapMode, zoomMode_t zoomMode, Double_t zoomCenterPhiDeg,
				 Double_t zoomCenterThetaDeg){


  fillCombosToUseIfNeeded(mapMode, l3TrigPattern);

  TString name, title;
  createImageNameAndTitle(name, title, mapMode, zoomMode, rWave, pol);
  
  TH2D* hImage = NULL;
  if(zoomMode == kZoomedOut){
    hImage = makeBlankImage(name, title);    
  }
  else if(zoomMode == kZoomedIn){
    hImage = makeBlankZoomedImage(name, title, zoomCenterPhiDeg, zoomCenterThetaDeg);
    correlateForZoomedImage(pol, l3TrigPattern); // Pads FFTs for more finely grained correlation.
    fillDeltaTLookupZoomed(zoomCenterPhiDeg, zoomCenterThetaDeg, l3TrigPattern);
  }

  // std::cout << zoomMode << "\t" << zoomModeNames[zoomMode] << "\t" << l3TrigPattern << std::endl;
  // for(int phiSect=0; phiSect<NUM_PHI; phiSect++){
  //   std::cout << phiSect << "\t" << combosToUse[phiSect].size() << std::endl;
  // }
  // std::cout << std::endl << std::endl;
  imagePeak = -DBL_MAX;
  Int_t peakPhiBin = -1;
  Int_t peakThetaBin = -1;
  for(Int_t phiInd = 0; phiInd < hImage->GetNbinsX(); phiInd++){
    Int_t phiBin = phiInd;
    Int_t phiSector = zoomMode==kZoomedIn ? 0 : phiBin/NUM_BINS_PHI;
    std::vector<Int_t> combosToUse;
    if(mapMode==kTriggered){
      combosToUse = combosToUseTriggered[l3TrigPattern];
    }
    else{
      combosToUse = combosToUseGlobal[phiSector];
    }
    
    Double_t phiWave = hImage->GetXaxis()->GetBinLowEdge(phiBin+1)*TMath::DegToRad();
    for(Int_t thetaBin = 0; thetaBin < hImage->GetNbinsY(); thetaBin++){
      Double_t thetaWave = hImage->GetYaxis()->GetBinLowEdge(thetaBin+1)*TMath::DegToRad();

      Double_t correlations = 0;
      for(UInt_t comboInd=0; comboInd<combosToUse.size(); comboInd++){
	Int_t combo = combosToUse.at(comboInd);
	Int_t offset = 0;

	// If we are in zoomed out & plane wave mode then use the lookup table for a big speedup
	if(zoomMode==kZoomedOut && rWave==0){
	  offset = deltaTs[phiBin][thetaBin][combo];
	  offset = offset < 0 ? offset + numSamples : offset;
	  correlations += crossCorrelations[pol][combo][offset];
	}
	// If we are in zoomed in & plane wave mode then calculate
	else if(zoomMode==kZoomedIn && rWave==0){
	  offset = deltaTsZoom[phiBin][thetaBin][combo];
	  offset = offset < 0 ? offset + numSamplesUpsampled : offset;
	  correlations += crossCorrelationsUpsampled[pol][combo][offset];
	}
	// If we are in zoomed out & spherical wave mode then calculate
	else{
	  Int_t ant1 = comboToAnt1s.at(combo);
	  Int_t ant2 = comboToAnt2s.at(combo);
	  offset = getDeltaTExpectedSpherical(ant1, ant2, phiWave, thetaWave, rWave);
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


void* CrossCorrelator::makeSomeOfImageThreaded(void* voidPtrArgs){
  // Disgusting hacks to get ROOT threading to compile inside a class.
  CrossCorrelator::threadArgs* args = (CrossCorrelator::threadArgs*) voidPtrArgs;
  Long_t threadInd = args->threadInd;
  CrossCorrelator* ptr = args->ptr;
  AnitaPol::AnitaPol_t pol = ptr->threadPol;

  threadInd++;
  ptr++;
  pol=AnitaPol::kHorizontal;
  
  return 0;
  
}




Double_t CrossCorrelator::findImagePeak(TH2D* hist, Double_t& imagePeakTheta, Double_t& imagePeakPhi){

  Int_t nx = hist->GetNbinsX();
  Int_t ny = hist->GetNbinsY();

  Double_t maxVal = -2;
  for(Int_t by = 1; by<=ny; by++){
    for(Int_t bx = 1; bx<=nx; bx++){
      Double_t val = hist->GetBinContent(bx, by);
      if(val > maxVal){
	maxVal = val;
	imagePeakPhi = hist->GetXaxis()->GetBinLowEdge(bx);
	imagePeakTheta = hist->GetYaxis()->GetBinLowEdge(by);
      }
    }
  }

  return maxVal;
  
}





/************************************************************************************************************
Functions to delete pointers to internal variables
************************************************************************************************************/
void CrossCorrelator::deleteCrossCorrelations(AnitaPol::AnitaPol_t pol){
  for(Int_t comboInd=0; comboInd<NUM_COMBOS; comboInd++){
    if(crossCorrelations[pol][comboInd] != NULL){
      delete [] crossCorrelations[pol][comboInd];
      crossCorrelations[pol][comboInd] = NULL;
    }
  }
}

void CrossCorrelator::deleteUpsampledCrossCorrelations(AnitaPol::AnitaPol_t pol){
  for(Int_t comboInd=0; comboInd<NUM_COMBOS; comboInd++){
    if(crossCorrelationsUpsampled[pol][comboInd] != NULL){
      delete [] crossCorrelationsUpsampled[pol][comboInd];
      crossCorrelationsUpsampled[pol][comboInd] = NULL;
    }
  }
}


void CrossCorrelator::deleteAllWaveforms(AnitaPol::AnitaPol_t pol){
  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
    if(grs[pol][ant]){
      delete grs[pol][ant];
      grs[pol][ant] = NULL;
    }
    if(grsInterp[pol][ant]){
      delete grsInterp[pol][ant];
      grsInterp[pol][ant] = NULL;
    }
  }
}

void CrossCorrelator::deleteAllFFTs(AnitaPol::AnitaPol_t pol){
  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
    if(ffts[pol][ant]){
      delete [] ffts[pol][ant];
      ffts[pol][ant] = NULL;
    }
  }
}


void CrossCorrelator::deleteAllPaddedFFTs(AnitaPol::AnitaPol_t pol){
  for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
    if(fftsPadded[pol][ant]){
      delete [] fftsPadded[pol][ant];
      fftsPadded[pol][ant] = NULL;
    }
  }
}


















/************************************************************************************************************
Functions for debugging or testing
************************************************************************************************************/

void CrossCorrelator::correlateEventTest(Double_t phiDegSource, Double_t thetaDegSource, Double_t rSource){
  // Generates a set of delta function like waveforms, correlates them 

  // Delete any old waveforms (at start rather than end to leave in memory to be examined if need be)

  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){  
    deleteAllWaveforms((AnitaPol::AnitaPol_t)pol);
  }


  for(Int_t pol = AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      std::vector<Double_t> newVs = std::vector<Double_t>(numSamples, 0);
      std::vector<Double_t> newTimes = std::vector<Double_t>(numSamples, 0);
      Int_t dt = numSamples/2 + getDeltaTExpectedSpherical(0, ant,
							    phiDegSource*TMath::DegToRad(), 
							    thetaDegSource*TMath::DegToRad(), 
							    rSource);
      for(int samp=0; samp<numSamples; samp++){
	newTimes[samp] = nominalSamplingDeltaT*samp;
	if(samp==dt){
	  newVs[samp] = numSamples;
	}
      }
      grsInterp[pol][ant] = new TGraph(numSamples, &newTimes[0], &newVs[0]);    
      RootTools::normalize(grsInterp[pol][ant]);
    }
    doAllCrossCorrelations((AnitaPol::AnitaPol_t)pol);
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


