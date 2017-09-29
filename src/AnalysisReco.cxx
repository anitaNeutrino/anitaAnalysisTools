#include "AnalysisReco.h"
#include "InterferometricMap.h"
#include "InterferometryCache.h"
#include "AcclaimFilters.h"
#include "ResponseManager.h"
#include "NoiseMonitor.h"
#include <numeric>
#include "RootTools.h"
#include "TGraphInteractive.h"
#include "ImpulsivityMeasure.h"
#include "TPaveText.h"

ClassImp(Acclaim::AnalysisReco);



/** 
 * Default constructor
 */
Acclaim::AnalysisReco::AnalysisReco(){
  initializeInternals();
}


/** 
 * Default destructor
 */
Acclaim::AnalysisReco::~AnalysisReco(){

  if(spawnedCrossCorrelator && fCrossCorr){
    delete fCrossCorr;
    fCrossCorr = NULL;
  }


  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    if(heatMaps[polInd]){
      delete heatMaps[polInd];
      heatMaps[polInd] = NULL;
    }
  }

  nicelyDeleteInternalFilteredEvents();  
}


void Acclaim::AnalysisReco::nicelyDeleteInternalFilteredEvents(){
  if(fEvMin != NULL){
    delete fEvMin;
    fEvMin = NULL;
  }
  if(fEvMinDeco != NULL){
    delete fEvMinDeco;
    fEvMinDeco = NULL;
  }
  if(fEvDeco != NULL){
    delete fEvDeco;
    fEvDeco = NULL;  
  }  
}




/**
 * Calculate relevant quantities and fill the waveform info in the AnitaEventSummary
 *
 * This section may change rapidly as the analysis continues...
 *
 * @param pol is the polarisation
 * @param info is the selected WaveformInfo in the ANITA event summary
 * @param fEv is the FilteredAnitaEvent from which we are coherently summing
 * @param waveStore is the internal storage in the AnalysisReco class (used to save a copy for MagicDisplay)
 * @param h is the InterferometricMap which contains the peak direction in which we want to coherently sum
 */
void Acclaim::AnalysisReco::fillWaveformInfo(AnitaPol::AnitaPol_t pol,
                                             AnitaEventSummary::WaveformInfo& info,
                                             const FilteredAnitaEvent* fEv,
                                             AnalysisWaveform** waveStore,
                                             InterferometricMap* h,
                                             NoiseMonitor* noiseMonitor) {

  // This section will fill the WaveformInfo as rapidly as I can keep up...
  // Some variables are marked with TODO, because they haven't been implemented yet.



  ////////////////////////////////////////////////////////////////
  // Double_t snr; ///Signal to Noise of waveform
  // Double_t peakHilbert; /// peak of hilbert envelope
  // Double_t peakVal;  /// peak value
  // Double_t xPolPeakVal;  // Peak of xpol trace
  // Double_t xPolPeakHilbert;  // Peak of xpol hilbert Envelope

  // Double_t I,Q,U,V;  // Stokes Parameters

  // Double_t totalPower;  ///Total power in waveform
  // Double_t totalPowerXpol;  ///Total power in xPol waveform

  Int_t peakPhiSector = h->getPeakPhiSector();
  const std::vector<Int_t>& theAnts = phiSectorToCoherentAnts(peakPhiSector);

  Double_t peak, phiDeg, thetaDeg;
  h->getPeakInfo(peak, phiDeg, thetaDeg);
  Double_t largestPeakToPeak = 0;
  AnalysisWaveform* coherentWave = coherentlySum(fEv, pol, theAnts, phiDeg, thetaDeg, &largestPeakToPeak);

  // Internal storage
  if(waveStore[0]){
    delete waveStore[0];
  }
  waveStore[0] = coherentWave;
  
  info.snr = 0;
  if(noiseMonitor){  
    double noise = 0;
    for(unsigned antInd = 0; antInd < theAnts.size(); antInd++){
      int ant = theAnts[antInd];
      double thisRMS = noiseMonitor->getRMS(pol, ant, fEv->getHeader()->realTime);

      noise += thisRMS;
      
    }
    noise /= theAnts.size();
    
    info.snr = largestPeakToPeak/(2*noise);
  }

  
  const TGraphAligned* grHilbert = coherentWave->hilbertEnvelope();
  info.peakHilbert = TMath::MaxElement(grHilbert->GetN(), grHilbert->GetY());

  const TGraphAligned* gr = coherentWave->even();
  info.peakVal = TMath::MaxElement(gr->GetN(), gr->GetY());

  Double_t localMaxVolts, localMinVolts, localMaxTime, localMinTime;
  RootTools::getLocalMaxToMin(gr, localMaxVolts, localMaxTime, localMinVolts, localMinTime);
  info.localMaxToMin = localMaxVolts - localMinVolts;
  info.localMaxToMinTime = localMaxTime - localMinTime;

  Int_t gMaxInd = TMath::LocMax(gr->GetN(), gr->GetY());
  Int_t gMinInd = TMath::LocMin(gr->GetN(), gr->GetY());
  info.globalMaxToMin = gr->GetY()[gMaxInd] - gr->GetY()[gMinInd];
  info.globalMaxToMinTime = gr->GetX()[gMaxInd] - gr->GetX()[gMinInd];
  

  AnitaPol::AnitaPol_t xPol = RootTools::swapPol(pol);
  Double_t largestPeakToPeakXPol = 0;
  AnalysisWaveform* xPolCoherentWave = coherentlySum(fEv, xPol, theAnts, phiDeg, thetaDeg, &largestPeakToPeakXPol); // cross polarisation

  // Internal storage
  if(waveStore[1]){
    delete waveStore[1];
  }
  waveStore[1] = xPolCoherentWave;

  const TGraphAligned* grHilbertX = xPolCoherentWave->hilbertEnvelope();
  info.xPolPeakHilbert = TMath::MaxElement(grHilbertX->GetN(), grHilbertX->GetY());

  const TGraphAligned* grX = xPolCoherentWave->even();
  info.xPolPeakVal = TMath::MaxElement(grX->GetN(), grX->GetY());

  AnalysisWaveform* hWave = pol == AnitaPol::kHorizontal ? coherentWave : xPolCoherentWave;
  AnalysisWaveform* vWave = pol == AnitaPol::kHorizontal ? xPolCoherentWave : coherentWave;  

  const TGraphAligned* grV = vWave->even();
  const TGraphAligned* grH = hWave->even();

  const TGraphAligned* grV_hilbert = vWave->hilbertTransform()->even();
  const TGraphAligned* grH_hilbert = hWave->hilbertTransform()->even();
  
  // Stokes params (I,Q,U,V) - from the wikipedia https://en.wikipedia.org/wiki/Stokes_parameters#Examples
  // (1,  1,  0,  0) HPol linearly polarised
  // (1, -1,  0,  0) VPol linearly polarised
  // (1,  0,  1,  0) +45 deg linearly polarised
  // (1,  0, -1,  0) -45 deg linearly polarised
  // (1,  0,  0,  1) Right hand circularaly polarised
  // (1,  0,  0, -1) Left hand circularaly polarised
  // (1,  0,  0,  0) Unpolarised

  FFTtools::stokesParameters(grH->GetN(),
                             grH->GetY(),  grH_hilbert->GetY(),
                             grV->GetY(), grV_hilbert->GetY(),
                             &(info.I), &(info.Q), &(info.U), &(info.V));

  info.totalPower = RootTools::getTimeIntegratedPower(coherentWave->even());
  info.totalPowerXpol = RootTools::getTimeIntegratedPower(xPolCoherentWave->even());

  // //Shape parameters, computed using hilbert envelope
  // // This should probably taken out into its own class
  // Double_t riseTime_10_90;  /// Rise time of hilbert env from 10% to 90% of peak
  // Double_t riseTime_10_50;  /// Rise time of hilbert env from 10% to 50% of peak
  // Double_t fallTime_90_10;  /// Fall time of hilbert env from 90% to 10% of peak
  // Double_t fallTime_50_10;  /// Fall time of hilbert env from 50% to 10% of peak
  // Double_t width_50_50;   /// Width from first envelope crossing of 50 percent of peak to last
  // Double_t width_10_10;  /// Width from first envelope crossing of 10 percent of peak to last
  // Double_t power_10_10;  /// Power enclosed within 10_10 width
  // Double_t power_50_50;  /// Power enclosed within 50_50 width
  // Double_t peakTime;  // Time that peak hilbert env occurs
  // Double_t peakMoments[5];  // moments about Peak  (1st - 5th moments)

  int maxIndex = std::find(grHilbert->GetY(), grHilbert->GetY() + grHilbert->GetN(), info.peakHilbert) - grHilbert->GetY();
  if(maxIndex >= grHilbert->GetN()){
    std::cerr << "Error in " << __PRETTY_FUNCTION__
              << ", unable to find time of Hilbert envelope peak using std::find()."
              << "peakHilbert = " << info.peakHilbert << ", numPoints = " << grHilbert->GetN()
              << ", maxIndex = " << maxIndex << std::endl;
  }
  info.peakTime = grHilbert->GetX()[maxIndex];

  grHilbert->getMoments(sizeof(info.peakMoments)/sizeof(double), info.peakTime, info.peakMoments);

  // set impulsivity measure defined in AnitaAnalysisSummary
  info.impulsivityMeasure = impulsivity::impulsivityMeasure(coherentWave, NULL, maxIndex, true);

  for(int w=0; w < AnitaEventSummary::numFracPowerWindows; w++){
    const double powerFraction = (w+1)*0.1;
    std::pair<double, double> window = RootTools::findSmallestWindowContainingFracOfPower(grHilbert, powerFraction);
    // info.narrowestWidths[w] = window.second - window.first;
    info.fracPowerWindowBegins[w] = window.first;
    info.fracPowerWindowEnds[w] = window.second;
  }

  ///////////////////////////////////////////////////////////////////////////////////
  // //spectrum info
  // Double_t bandwidth[peaksPerSpectrum];  /// bandwidth of each peak (implementation defined, may not be comparable between analyses) 
  // Double_t peakFrequency[peaksPerSpectrum]; //peak frequency of power spectrum 
  // Double_t peakPower[peaksPerSpectrum]; //power within +/- bandwidth of each peak 
  // Double_t spectrumSlope;  ///  Slope of line fit to spectrum (in log-space, so this is spectral-index) 
  // Double_t spectrumIntercept; /// Intercept of line fit to spectrum (in log-space) 

  // we have interpolated the coherently summed waveform,  probably significantly to extract useful information
  // from the time domain. Therefore there's a bunch of trailing high frequency crap that I want to ignore
  // The maximum useful frequency is the nominal sampling

  const double maxFreqIfIHadNotIntepolated = FancyFFTs::getNumFreqs(NUM_SAMP)*(1./(NOMINAL_SAMPLING_DELTAT*NUM_SAMP));
  const double df = coherentWave->deltaF(); // the actual deltaF
  int numGoodFreqBins = floor(maxFreqIfIHadNotIntepolated/df);

  // std::cout << maxFreq << std::endl;
  const TGraphAligned* grPow = coherentWave->power();

  if(numGoodFreqBins > grPow->GetN()){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", expected at least " << numGoodFreqBins << " frequency bins, only got " << grPow->GetN()
              << ". Did you downsample the coherently summed waveform?" << std::endl;
    numGoodFreqBins = grPow->GetN();
  }

  TGraph grMaxima;
  for(int i=1; i < numGoodFreqBins-1; i++){
    if(grPow->GetY()[i] > grPow->GetY()[i-i]  && grPow->GetY()[i] > grPow->GetY()[i+i] ){
      grMaxima.SetPoint(grMaxima.GetN(),  grPow->GetX()[i], grPow->GetY()[i]);
    }

    // This check should be redundent now, but won't hurt to leave it in...
    if(grPow->GetX()[i] >= maxFreqIfIHadNotIntepolated){
      break;
    }
  }

  // Sort my collection of local maxima,
  // The second argument, false, means sort the local maxima in decreasing order (i.e. biggest first)
  grMaxima.Sort(TGraph::CompareY, false);

  // Put the N largest local maxima into the AnitaEventSummary::WaveformInfo
  for(int i=0; i < AnitaEventSummary::peaksPerSpectrum; i++){
    if(i < grMaxima.GetN()){
      info.bandwidth[i] = 0; // for now, let's just try the peak value and frequency
      info.peakPower[i] = grMaxima.GetY()[i];
      info.peakFrequency[i] = grMaxima.GetX()[i];
    }
    else{
      info.bandwidth[i] = 0;
      info.peakPower[i] = 0;
      info.peakFrequency[i] = 0;
    }
  }
  
  const TGraphAligned* grPow_db = coherentWave->powerdB();

  const int firstSlopeBin = TMath::Nint(fSlopeFitStartFreqGHz/df);
  const int lastSlopeBin  = TMath::Nint(fSlopeFitEndFreqGHz/df);
  const int numSlopeBins = lastSlopeBin - firstSlopeBin;

  const double mean_f      = TMath::Mean(numSlopeBins, &grPow_db->GetX()[firstSlopeBin]);
  const double mean_pow_db = TMath::Mean(numSlopeBins, &grPow_db->GetY()[firstSlopeBin]);

  double gradNumerator = 0;
  double gradDenominator = 0;
  for(int i=firstSlopeBin; i < lastSlopeBin; i++){
    double dx = grPow_db->GetX()[i] - mean_f;
    double dy = grPow_db->GetY()[i] - mean_pow_db;

    gradNumerator   += dy*dx;
    gradDenominator += dx*dx;
  }

  if(gradDenominator <= 0){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", got gradDenominator = " << gradDenominator << " will get NaN or Inf..." << std::endl;
  }
  info.spectrumSlope = gradNumerator/gradDenominator;
  info.spectrumIntercept = mean_pow_db - info.spectrumSlope*mean_f;

}






/**
 * Creates the histogram with the easting/northing dimensions for Antarctica
 *
 * @param name is the histogram name
 * @param title is the histogram title
 *
 * @return a pointer to the newly created histogram
 */
TProfile2D* Acclaim::AnalysisReco::makeHeatMap(const TString& name, const TString& title){

  const int coarseness = 20;
  Int_t nx, ny;
  RampdemReader::getNumXY(nx, ny);
  nx/=coarseness;
  ny/=coarseness;
  Double_t xMin, xMax, yMin, yMax;
  RampdemReader::getMapCoordinates(xMin, yMin, xMax, yMax);

  return new TProfile2D(name, title, nx, xMin, xMax, ny, yMin, yMax);    
}






/** 
 * Fill Peng's ChannelInfo object, it might come in handy...
 * 
 * @param fEv is the filtered event
 * @param sum is the event summary to fill
 */
void Acclaim::AnalysisReco::fillChannelInfo(const FilteredAnitaEvent* fEv, AnitaEventSummary* sum){
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(int ant=0; ant<NUM_SEAVEYS; ant++){
      const AnalysisWaveform* wf = fEv->getFilteredGraph(ant, pol);
      const TGraphAligned* gr = wf->even();
      const TGraphAligned* he = wf->hilbertEnvelope();

      sum->channels[polInd][ant].rms = gr->GetRMS();
      sum->channels[polInd][ant].avgPower = RootTools::getTimeIntegratedPower(gr);
      Double_t maxY, maxX, minY, minX;
      RootTools::getLocalMaxToMin(gr, maxY, maxX, minY, minX);
      sum->channels[polInd][ant].snr = (maxY - minY)/sum->channels[polInd][ant].rms;
      sum->channels[polInd][ant].peakHilbert = TMath::MaxElement(he->GetN(), he->GetY());
    }
  }
}


/** 
 * Reconstructs the FilteredAnitaEvent and fills the pointing and coherent bits of the AnitaEventSummary
 * 
 * @param fEv is the FilteredAnitaEvent to reconstruct
 * @param usefulPat is the useful version of the ANITA gps object, can trace direction to the ground
 * @param eventSummary is the AnitaEventSummary to fill during the recontruction
 * @param noiseMonitor is an optional parameter, which points to a NoiseMonitor, which tracks the RMS of MinBias events (used in SNR calculation)
 */

void Acclaim::AnalysisReco::process(const FilteredAnitaEvent * fEv, AnitaEventSummary * eventSummary, NoiseMonitor* noiseMonitor, TruthAnitaEvent* truth) {


  // spawn a CrossCorrelator if we need one
  if(!fCrossCorr){
    fCrossCorr = new CrossCorrelator();
    spawnedCrossCorrelator = true;
    dtCache.init(fCrossCorr, this);
  }

  // FilteredAnitaEvent fEvMin(fEv->getUsefulAnitaEvent(), fMinFilter, fEv->getGPS(), fEv->getHeader(), false); // just the minimum filter
  // FilteredAnitaEvent fEvMinDeco(fEv->getUsefulAnitaEvent(), fMinDecoFilter, fEv->getGPS(), fEv->getHeader(), false); // minimum + deconvolution
  // FilteredAnitaEvent fEvDeco(fEv, fMinDecoFilter); // extra deconvolution
  nicelyDeleteInternalFilteredEvents();
  fEvMin = new FilteredAnitaEvent(fEv->getUsefulAnitaEvent(), fMinFilter, fEv->getGPS(), fEv->getHeader(), false); // just the minimum filter
  fEvMinDeco = new FilteredAnitaEvent(fEv->getUsefulAnitaEvent(), fMinDecoFilter, fEv->getGPS(), fEv->getHeader(), false); // minimum + deconvolution
  fEvDeco = new FilteredAnitaEvent(fEv, fMinDecoFilter); // extra deconvolution

  fillChannelInfo(fEv, eventSummary);

  eventSummary->eventNumber = fEv->getHeader()->eventNumber;

  chooseAntennasForCoherentlySumming(fCoherentDeltaPhi);

  for(Int_t polInd=0; polInd < AnitaPol::kNotAPol; polInd++){

    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;

    fCrossCorr->correlateEvent(fEv, pol);

    // do the coarsely binned reconstruction
    reconstruct(pol);

    std::vector<Double_t> coarseMapPeakValues;
    std::vector<Double_t> coarseMapPeakPhiDegs;
    std::vector<Double_t> coarseMapPeakThetaDegs;
    coarseMaps[pol]->findPeakValues(fNumPeaks, coarseMapPeakValues, coarseMapPeakPhiDegs, coarseMapPeakThetaDegs);
    coarseMaps[pol]->addGpsInfo(fEv->getGPS());
    coarseMaps[pol]->addTruthInfo(truth);

    if(fDoHeatMap > 0){
      if(!heatMaps[pol]){
        TString name = pol == AnitaPol::kHorizontal ? "heatMapH" : "heatMapV";
        heatMaps[pol] = makeHeatMap(name, name);
      }
      coarseMaps[pol]->project(heatMaps[pol], fHeatMapHorizonKm);
    }

    eventSummary->nPeaks[pol] = fNumPeaks;

    for(Int_t peakInd=0; peakInd < fNumPeaks; peakInd++){
      reconstructZoom(pol, peakInd, coarseMapPeakPhiDegs.at(peakInd), coarseMapPeakThetaDegs.at(peakInd));

      InterferometricMap* h = fineMaps[pol][peakInd];
      if(h){

        // for plotting
        h->addGpsInfo(fEv->getGPS());
        h->addTruthInfo(truth);


        h->getPeakInfo(eventSummary->peak[pol][peakInd].value,
                       eventSummary->peak[pol][peakInd].phi,
                       eventSummary->peak[pol][peakInd].theta,
                       eventSummary->peak[pol][peakInd].chisq);

        // fill in difference between rough and fine
        eventSummary->peak[pol][peakInd].dphi_rough = eventSummary->peak[pol][peakInd].phi - coarseMapPeakPhiDegs.at(peakInd);
        eventSummary->peak[pol][peakInd].dtheta_rough = eventSummary->peak[pol][peakInd].theta - coarseMapPeakThetaDegs.at(peakInd);

        // based on Cosmin's comments in AnitaAnalysisSummary.h
        eventSummary->peak[pol][peakInd].phi_separation = RootTools::getDeltaAngleDeg(eventSummary->peak[pol][peakInd].phi, eventSummary->peak[pol][0].phi);

        // hwAngle, is the angle between the peak and the nearest L3 phi-sector trigger of that polarisation        
        setTriggerInfoFromPeakPhi(fEv->getHeader(), pol, h->getPeakPhiSector(), eventSummary->peak[pol][peakInd]);

        // coherent
        // coherent_filtered
        // deconvolved
        // deconvolved_filtered
      
        fillWaveformInfo(pol, eventSummary->coherent_filtered[pol][peakInd],    fEv,        wfCoherentFiltered[pol][peakInd],    h, noiseMonitor);
        fillWaveformInfo(pol, eventSummary->coherent[pol][peakInd]         ,    fEvMin,     wfCoherent[pol][peakInd],            h, noiseMonitor);
        fillWaveformInfo(pol, eventSummary->deconvolved_filtered[pol][peakInd], fEvDeco,    wfDeconvolvedFiltered[pol][peakInd], h, noiseMonitor);
        fillWaveformInfo(pol, eventSummary->deconvolved[pol][peakInd],          fEvMinDeco, wfDeconvolved[pol][peakInd],         h, noiseMonitor);

        
        if(fEv->getGPS() != NULL){
          UsefulAdu5Pat usefulPat(fEv->getGPS());
          int success = 0;
          if(eventSummary->peak[pol][peakInd].theta < 0){ // work around for bug in traceBackToContinent
            Double_t phiWave = TMath::DegToRad()*eventSummary->peak[pol][peakInd].phi;
            Double_t thetaWave = -1*TMath::DegToRad()*eventSummary->peak[pol][peakInd].theta;

            // *   Returns 0 if never hits the ground, even with maximum adjustment
            // *   Returns 1 if hits the ground with no adjustment
            // *   Returns 2 if it hits the ground with adjustment
            success = usefulPat.traceBackToContinent(phiWave, thetaWave,
                                                     &eventSummary->peak[pol][peakInd].longitude,
                                                     &eventSummary->peak[pol][peakInd].latitude,
                                                     &eventSummary->peak[pol][peakInd].altitude,
                                                     &eventSummary->peak[pol][peakInd].theta_adjustment_needed);
          }
          if(success==0){
            eventSummary->peak[pol][peakInd].longitude = -9999;
            eventSummary->peak[pol][peakInd].latitude = -9999;
            eventSummary->peak[pol][peakInd].altitude = -9999;
            eventSummary->peak[pol][peakInd].theta_adjustment_needed = -9999;
            eventSummary->peak[pol][peakInd].distanceToSource = -9999;
          }
          else{
            eventSummary->peak[pol][peakInd].distanceToSource = SPEED_OF_LIGHT_NS*usefulPat.getTriggerTimeNsFromSource(eventSummary->peak[pol][peakInd].latitude,
                                                                                                                       eventSummary->peak[pol][peakInd].longitude,
                                                                                                                       eventSummary->peak[pol][peakInd].altitude);
          }
        }
      }
    }
  }

  // keep internal copy of summary
  summary = (*eventSummary);
}





/** 
 * Sets the phi depended trigger stuff in the peak
 * 
 * @param header is the RawAnitaHeader
 * @param h is the interferometric map
 * @param peak is a reference
 */
void Acclaim::AnalysisReco::setTriggerInfoFromPeakPhi(const RawAnitaHeader* header,
                                                      AnitaPol::AnitaPol_t pol,
                                                      Int_t peakPhiSector,
                                                      AnitaEventSummary::PointingHypothesis& peak){


  AnitaPol::AnitaPol_t xPol = RootTools::swapPol(pol);
  // const RawAnitaHeader* header = fEv->getHeader();
  peak.hwAngle = -9999;
  peak.hwAngleXPol = -9999;

  for(int phi=0; phi < NUM_PHI; phi++){
    int phiWasTriggered = header->isInL3Pattern(phi, pol);
    if(phiWasTriggered > 0){
      double phiSectorPhi = InterferometricMap::getPhiSectorCenterPhiDeg(phi);
      double dPhi = RootTools::getDeltaAngleDeg(peak.phi, phiSectorPhi);
      if(TMath::Abs(dPhi) < TMath::Abs(peak.hwAngle)){
        peak.hwAngle = dPhi;
      }
    }
    int phiWasTriggeredXPol = header->isInL3Pattern(phi, xPol);
    if(phiWasTriggeredXPol > 0){
      double phiSectorPhi = InterferometricMap::getPhiSectorCenterPhiDeg(phi);
      double dPhi = RootTools::getDeltaAngleDeg(peak.phi, phiSectorPhi);
      if(TMath::Abs(dPhi) < TMath::Abs(peak.hwAngleXPol)){
        peak.hwAngleXPol = dPhi;
      }
    }
  }

  peak.triggered = header->isInL3Pattern(peakPhiSector, pol);
  peak.triggered_xpol = header->isInL3Pattern(peakPhiSector, xPol);

  peak.masked = header->isInPhiMaskOffline(peakPhiSector, pol) || header->isInL1MaskOffline(peakPhiSector, pol);
  peak.masked_xpol = header->isInPhiMaskOffline(peakPhiSector, xPol) || header->isInL1MaskOffline(peakPhiSector, xPol);

  return;
}


/**
 * Set default values and zero pointers for dynamically initialised heap dwelling members
 *
 */
void Acclaim::AnalysisReco::initializeInternals(){

  fCrossCorr = NULL;
  spawnedCrossCorrelator = false;
  fWhichResponseDir = 0;
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  for(Int_t pol=0; pol < AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant<NUM_SEAVEYS; ant++){
      rArray[pol].push_back(geom->getAntR(ant, AnitaPol::AnitaPol_t(pol)));
      zArray[pol].push_back(geom->getAntZ(ant, AnitaPol::AnitaPol_t(pol)));
      phiArrayDeg[pol].push_back(geom->getAntPhiPositionRelToAftFore(ant, AnitaPol::AnitaPol_t(pol))*TMath::RadToDeg());
    }
  }

  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    coarseMaps[polInd] = NULL;
    heatMaps[polInd] = NULL;
    for(int peakInd=0; peakInd < AnitaEventSummary::maxDirectionsPerPol; peakInd++){
      fineMaps[polInd][peakInd] = NULL;
      for(int xPol=0; xPol < AnitaPol::kNotAPol; xPol++){
        wfCoherentFiltered[polInd][peakInd][xPol] = NULL;
        wfCoherent[polInd][peakInd][xPol] = NULL;
        wfDeconvolved[polInd][peakInd][xPol] = NULL;
        wfDeconvolvedFiltered[polInd][peakInd][xPol] = NULL;
      }
    }
  }

  fUseOffAxisDelay = 1;
  fCoherentDeltaPhi = 2;
  fLastCoherentDeltaPhi = -1;
  fNumPeaks = 3;
  fDoHeatMap = 0;
  fHeatMapHorizonKm = 800;
  fCoherentDtNs = 0.01;
  fSlopeFitStartFreqGHz = 0.18;
  fSlopeFitEndFreqGHz = 1.3;
  
  const TString minFiltName = "Minimum";
  fMinFilter = Filters::findDefaultStrategy(minFiltName);
  if(!fMinFilter){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", unable to find default filter " << minFiltName << std::endl;
  }  

  const TString minDecoFiltName = "Deconvolve";
  fMinDecoFilter = Filters::findDefaultStrategy(minDecoFiltName);
  if(!fMinDecoFilter){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", unable to find default filter " << minDecoFiltName << std::endl;
  }

  fEvMin = NULL;
  fEvMinDeco = NULL;
  fEvDeco = NULL;
  
}




/** 
 * Get the off boresight delay for a single antenna
 * 
 * @param deltaPhiDeg is the phi angle from boresight
 * 
 * @return the delay
 */
Double_t Acclaim::AnalysisReco::singleAntennaOffAxisDelay(Double_t deltaPhiDeg) const {


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
  static const Double_t params[numPowers] = {0.00000e+00, 0.00000e+00, -1.68751e-05, 0.00000e+00,
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



/** 
 * Get the off axis delay between two antennas for a given phi angle
 * 
 * @param pol is the polarisation
 * @param ant1 is the first antenna
 * @param ant2 is the second antennas
 * @param phiDeg is the position in degrees in payload coordinates relative to ADU5 aft-fore
 * 
 * @return the relative off axis delay
 */
Double_t Acclaim::AnalysisReco::relativeOffAxisDelay(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
                                                     Double_t phiDeg) const {

  Double_t deltaPhiDeg1 = RootTools::getDeltaAngleDeg(phiArrayDeg[pol].at(ant1), phiDeg);
  Double_t deltaPhiDeg2 = RootTools::getDeltaAngleDeg(phiArrayDeg[pol].at(ant2), phiDeg);
  Double_t delay1 = singleAntennaOffAxisDelay(deltaPhiDeg1);
  Double_t delay2 = singleAntennaOffAxisDelay(deltaPhiDeg2);
  return (delay1-delay2);
}




/** 
 * Get the expected delay between antenna pairs for a given direction
 * 
 * @param pol is the polarisation
 * @param ant1 is the first antenna
 * @param ant2 is the second antenna
 * @param phiWave is the incoming plane wave direction in radians in payload coordinate relative to ADU5 aft-fore
 * @param thetaWave is the incoming plane wave direction in radians (theta=0 is horizontal, +ve theta is up)
 * 
 * @return the expected time difference in nano-seconds
 */
Double_t Acclaim::AnalysisReco::getDeltaTExpected(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t phiWave, Double_t thetaWave) const {

  // Double_t tanThetaW = tan(thetaWave);
  Double_t tanThetaW = tan(-1*thetaWave);
  Double_t part1 = zArray[pol].at(ant1)*tanThetaW - rArray[pol].at(ant1)*cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol].at(ant1));
  Double_t part2 = zArray[pol].at(ant2)*tanThetaW - rArray[pol].at(ant2)*cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol].at(ant2));
  Double_t tdiff = 1e9*((cos(thetaWave) * (part2 - part1))/SPEED_OF_LIGHT); // Returns time in ns

  return tdiff;
}




/** 
 * Inserts the photogrammetry geometry from AnitaGeomTool into this classes copy of the antenna position. 
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !WARNING! This class overwrites the "extra cable delays" in AnitaEventCalibrator with zero! 
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * Probably don't use this, or if you do, use with caution.
 * 
 */
void Acclaim::AnalysisReco::insertPhotogrammetryGeometry(){
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

  if(!fCrossCorr){
    fCrossCorr = new CrossCorrelator();
    spawnedCrossCorrelator = true;
  }
  
  dtCache.init(fCrossCorr, this, true);
  geom->usePhotogrammetryNumbers(0);

}






/** 
 * Get the coarse interferometric map in class memory. It is the user's responsibility to delete!
 * 
 * @param pol is the polarisation of the map to get.
 * 
 * @return the internal copy of the interferometric map
 */
Acclaim::InterferometricMap* Acclaim::AnalysisReco::getMap(AnitaPol::AnitaPol_t pol){
  InterferometricMap* h = coarseMaps[pol];
  coarseMaps[pol] = NULL;
  return h;
}






/** 
 * Get the fine interferometric map in class memory from the chosen peak. It is the user's responsibility to delete!
 * 
 * @param pol is the polarisation of the map to get.
 * 
 * @return the internal copy of the interferometric map
 */
Acclaim::InterferometricMap* Acclaim::AnalysisReco::getZoomMap(AnitaPol::AnitaPol_t pol, UInt_t peakInd){

  InterferometricMap* h = fineMaps[pol][peakInd];
  if(h==NULL){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unable to find fineMap with pol " << pol 
              << " for peakInd = " << peakInd << " about to return NULL." << std::endl;
  }
  fineMaps[pol][peakInd] = NULL;

  return h;
}







/** 
 * Do the coarsely binned reconstruction.
 * Overwrites the last map in memory, if there was one, creates one otherwise.
 * 
 * @param pol is the polarisation to reconstruct
 */
void Acclaim::AnalysisReco::reconstruct(AnitaPol::AnitaPol_t pol, const Adu5Pat* pat) {

  if(coarseMaps[pol]==NULL)
  {
    // std::cerr << "new coarse map " << pol << std::endl;
    // I don't own the map or there isn't one, so I'll make a new one
    coarseMaps[pol] = new InterferometricMap();
  }  
  else{// I own the map so I can overwrite it to avoid allocating memory
    // std::cerr << "old coarse map " << pol << std::endl;
    coarseMaps[pol]->Reset();
  }

  if(pat){
    coarseMaps[pol]->addGpsInfo(pat);
  }

  coarseMaps[pol]->Fill(pol, fCrossCorr, &dtCache);
}




/** 
 * 
 * 
 * @param pol 
 * @param peakIndex 
 * @param zoomCenterPhiDeg 
 * @param zoomCenterThetaDeg 
 */
void Acclaim::AnalysisReco::reconstructZoom(AnitaPol::AnitaPol_t pol, Int_t peakIndex, Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg, const Adu5Pat* pat) {

  // Some kind of sanity check here due to the unterminating while loop inside RootTools::getDeltaAngleDeg
  if(zoomCenterPhiDeg < -500 || zoomCenterThetaDeg < -500 ||
     zoomCenterPhiDeg >= 500 || zoomCenterThetaDeg >= 500){

    std::cerr << "Warning in "<< __PRETTY_FUNCTION__ << " in " << __FILE__
              << ". zoomCenterPhiDeg = " << zoomCenterPhiDeg
              << " and zoomCenterThetaDeg = " << zoomCenterThetaDeg
              << " ...these values look suspicious so I'm skipping this reconstruction."
              << " eventNumber = " << fCrossCorr->eventNumber[pol] << std::endl;
    return;
  }

  Double_t deltaPhiDegPhi0 = RootTools::getDeltaAngleDeg(zoomCenterPhiDeg, InterferometricMap::getBin0PhiDeg());
  deltaPhiDegPhi0 = deltaPhiDegPhi0 < 0 ? deltaPhiDegPhi0 + DEGREES_IN_CIRCLE : deltaPhiDegPhi0;

  Int_t phiSector = floor(deltaPhiDegPhi0/PHI_RANGE);

  InterferometricMap* h = new InterferometricMap(peakIndex, phiSector, zoomCenterPhiDeg, PHI_RANGE_ZOOM, zoomCenterThetaDeg, THETA_RANGE_ZOOM);
  if(pat){
    h->addGpsInfo(pat);
  }
  h->Fill(pol, fCrossCorr, &dtCache);  

  // std::cout << h->GetName() << std::endl;
  // std::cout << "in " << __PRETTY_FUNCTION__ << ": " << fineMaps[pol][peakIndex] << std::endl;
  if(fineMaps[pol][peakIndex]){
    delete fineMaps[pol][peakIndex];
  }
  fineMaps[pol][peakIndex] = h;
}







/** 
 * Directly insert some geometry generated by Linda.
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * !WARNING! This class overwrites the "extra cable delays" in AnitaEventCalibrator with zero! 
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * Probably don't use this, or if you do, use with caution.
 * 
 * @param pathToLindasFile points to the file containing the antenna positions/cable delays.
 * @param pol is the polarisation to overwrite.
 * 
 * @return 1 on error, 0 if successful.
 */
Int_t Acclaim::AnalysisReco::directlyInsertGeometry(TString pathToLindasFile, AnitaPol::AnitaPol_t pol){

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







/** 
 * Generate a set of antennas to use to generate the coherent waveform.
 * 
 * @param coherentDeltaPhi is the number of neighbouring phi-sectors to use.
 */
void Acclaim::AnalysisReco::chooseAntennasForCoherentlySumming(int coherentDeltaPhi){

  // std::cerr << __PRETTY_FUNCTION__ << std::endl;

  if(coherentDeltaPhi!=fLastCoherentDeltaPhi){

    for(int peakPhiSector = 0; peakPhiSector < NUM_PHI; peakPhiSector++){
      phiSectorToAnts[peakPhiSector] = std::vector<Int_t>();

      // phiSectorToAnts[peakPhiSector].clear();
      for(int deltaPhiSect=-fCoherentDeltaPhi; deltaPhiSect <= fCoherentDeltaPhi; deltaPhiSect++){

        Int_t phiSector = peakPhiSector + deltaPhiSect;
        phiSector = phiSector < 0        ? phiSector + NUM_PHI : phiSector;
        phiSector = phiSector >= NUM_PHI ? phiSector - NUM_PHI : phiSector;

        for(int ring = 0; ring < AnitaRing::kNotARing; ring++){
          int ant = AnitaGeomTool::getAntFromPhiRing(phiSector, AnitaRing::AnitaRing_t(ring));
          phiSectorToAnts[peakPhiSector].push_back(ant);
        }
      }
    }
    fLastCoherentDeltaPhi = fCoherentDeltaPhi;
  }
}



/** 
 * Modularising the coherently summed waveforms
 * 
 * @param theAnts are the antenna indices
 * @param pol is the polarisation
 * @param peakPhiDeg is the phi angle relative to ADU5 aft-fore (degrees)
 * @param peakThetaDeg is the theta angle (0 = horizontal, +ve theta is up)
 * @param dts is a vector which gets filled with the time offsets
 */
void Acclaim::AnalysisReco::directionAndAntennasToDeltaTs(const std::vector<Int_t>& theAnts, AnitaPol::AnitaPol_t pol,
                                                          Double_t peakPhiDeg, Double_t peakThetaDeg, std::vector<double>& dts) {

  dts.clear(); // first empty the vector
  Double_t phiWave = peakPhiDeg*TMath::DegToRad();
  Double_t thetaWave = peakThetaDeg*TMath::DegToRad();
  for(unsigned antInd=0; antInd < theAnts.size(); antInd++){
    int ant = theAnts.at(antInd);
    Double_t dt = getDeltaTExpected(pol, theAnts.at(0), ant, phiWave, thetaWave);
    dts.push_back(dt);
  }
}




/** 
 * Takes a polarisation, set of antennas, and coherently sums the waveforms in the FilteredAnitaEvents with the time delays implied by the incoming direction
 * 
 * @param fEv is the FilteredAnitaEvent from which to draw the waveforms to coherently sum
 * @param pol is the polarisation of the waveforms in the FilteredAnitaEvent 
 * @param theAnts is the set of antennas to use in the coherent sum
 * @param peakPhiDeg is the incoming phi-direction in degrees in payload coordinates relative to ADU5 aft-fore
 * @param peakThetaDeg is the elevation in degrees (theta=0 is horizontal, +ve theta is up) in payload coordinates relative to ADU5 aft-fore
 * 
 * @return the coherently summed waveform
 */
AnalysisWaveform* Acclaim::AnalysisReco::coherentlySum(const FilteredAnitaEvent* fEv, AnitaPol::AnitaPol_t pol,
                                                       const std::vector<Int_t>& theAnts,
                                                       Double_t peakPhiDeg, Double_t peakThetaDeg, Double_t* biggestPeakToPeak) {

  Int_t biggest = -1;
  Double_t largestPeakToPeak = 0;
  for(unsigned antInd=0; antInd < theAnts.size(); antInd++){
    int ant = theAnts.at(antInd);

    const AnalysisWaveform* wf = fEv->getFilteredGraph(ant, pol);
    const TGraphAligned* gr = wf->even();

    Double_t vMax, vMin, tMax, tMin;
    RootTools::getLocalMaxToMin((TGraph *)gr, vMax, tMax, vMin, tMin);

    if(vMax - vMin > largestPeakToPeak){
      largestPeakToPeak = vMax - vMin;
      biggest = ant;
    }
    else if(largestPeakToPeak <= 0){
      std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", the waveform max/min aren't sensible?" << std::endl;
      std::cerr << vMax << "\t" << vMin << "\t" << tMax << "\t" << tMin << std::endl;
    }
  }

  if(biggestPeakToPeak){
    (*biggestPeakToPeak) = largestPeakToPeak;
  }

  if(biggest < 0 || biggest >= NUM_SEAVEYS){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", I couldn't find a waveform where vMax - vMin > 0. "
              << "Something's wrong, and I'm probably about to vomit a stack trace all over your terminal..." << std::endl;
  }

  // now we've found the channel with the biggest peak-to-peak
  // prepare to coherently sum the waves with this channel first
  std::vector<const AnalysisWaveform*> waves;
  waves.push_back(fEv->getFilteredGraph(biggest, pol));
  std::vector<int> orderedAnts(1, biggest);
  for(unsigned antInd=0; antInd < theAnts.size(); antInd++){
    int ant = theAnts.at(antInd);
    if(ant != biggest){
      waves.push_back(fEv->getFilteredGraph(ant, pol));
      orderedAnts.push_back(ant);
    }
  }

  std::vector<double> dts;
  directionAndAntennasToDeltaTs(orderedAnts, pol, peakPhiDeg, peakThetaDeg, dts);
  
  return coherentlySum(waves, dts);
}







/** 
 * Check vector lengths are correct, trims/pads the dts if required
 * 
 * @param waves the waveforms to sum
 * @param dts the time offsets to apply
 * 
 * @return the size of both vectors
 */
size_t Acclaim::AnalysisReco::checkWavesAndDtsMatch(std::vector<const AnalysisWaveform*>& waves, std::vector<Double_t>& dts) {
  size_t s = waves.size();
  if(s < 1){    
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", nothing to sum." << std::endl;
  }  
  else if(s != dts.size()){    
    const char* action = dts.size() < waves.size() ? "padding" : "trimming";
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unequal vectors (waves.size() = " << waves.size() << ", dts.size() = " << dts.size() << ")\n"
              << action << " dts..." << std::endl;
    while(s != dts.size()){
      if(dts.size() < s){
        dts.push_back(0);
      }
      else{
        dts.pop_back();	
      }
    }
  }
  return s;
}





/** 
 * Coherently sum waveforms
 * 
 * @param waves is a set of waveforms to coherently sum
 * @param dts is the set of time offsets (ns) to apply to those waveforms when summing
 * 
 * @return a newly created AnalysisWaveform, created by coherently summing the input waves
 */
AnalysisWaveform* Acclaim::AnalysisReco::coherentlySum(std::vector<const AnalysisWaveform*>& waves, std::vector<Double_t>& dts) {

  if(checkWavesAndDtsMatch(waves, dts)==0){
    std::cerr << "Nothing to sum.. about to return NULL" << std::endl;
    return NULL;
  }  

  const TGraph* gr0 = waves[0]->even();
  const double totalTime = gr0->GetX()[gr0->GetN()-1] - gr0->GetX()[0];
  const int interpN = floor(totalTime/fCoherentDtNs);
  
  std::vector<double> zeros(interpN, 0);
  AnalysisWaveform* coherentWave = new AnalysisWaveform(interpN, &zeros[0], fCoherentDtNs, gr0->GetX()[0]);
  TGraphAligned* grCoherent = coherentWave->updateEven();

  for(UInt_t i=0; i < waves.size(); i++){
    const TGraphAligned* grU = waves[i]->even();
    // std::cout << i << " here" << std::endl;
    
    const double t0 = grU->GetX()[0];
    const double tN = grU->GetX()[grU->GetN()-1];
    
    for(int samp=0; samp < grCoherent->GetN(); samp++){
      double t = grCoherent->GetX()[samp];
      double tPlusDt = t + dts[i];
      if(tPlusDt > t0 && tPlusDt < tN){
        double y = waves[i]->evalEven(tPlusDt, AnalysisWaveform::EVAL_AKIMA);
        if(TMath::IsNaN(y)){ // Hopefully this is a bit redundent
          std::cerr << "Warning in " << __PRETTY_FUNCTION__ << " interpolation returned " << y << " ignoring sample" << std::endl;
        }
        else{
          grCoherent->GetY()[samp] += y;
        }
      }
    }
  }

  for(int samp=0; samp < grCoherent->GetN(); samp++){
    grCoherent->GetY()[samp]/=waves.size();
  };    

  return coherentWave;
}




/** 
 * Generate the waveforms that do into the coherent sum
 * 
 * @param waves is a set of waveforms to coherently sum
 * @param dts is the set of time offsets (ns) to apply to those waveforms when summing
 * 
 * @return a newly created AnalysisWaveform, created by coherently summing the input waves
 */
void Acclaim::AnalysisReco::wavesInCoherent(std::vector<const AnalysisWaveform*>& waves, std::vector<Double_t>& dts, std::vector<TGraphAligned*>& grs){

  for(unsigned w=0; w < waves.size(); w++){
    const TGraphAligned* grOld = waves[w]->even();
    TGraphAligned* gr = new TGraphAligned(grOld->GetN(), grOld->GetX(), grOld->GetY());

    for(int i=0; i < gr->GetN(); i++){
      gr->GetX()[i] -= dts[w];
    }
    
    grs.push_back(gr);
  }
}













/** 
 * Function for MagicDisplay
 * 
 * @param summaryPad is the pad if it already exists (makes a new Canvas if passed NULL)
 * @param pol is the polarization to draw
 */
void Acclaim::AnalysisReco::drawSummary(TPad* wholePad, AnitaPol::AnitaPol_t pol){

  const int numColsForNow = 3;
  EColor peakColors[AnitaPol::kNotAPol][numColsForNow] = {{kBlack, EColor(kMagenta+2), EColor(kViolet+2)},
                                                          {kBlue,  EColor(kSpring+4),  EColor(kPink + 10)}}; 

  if(wholePad==NULL){
    UInt_t eventNumber = fCrossCorr->eventNumber[pol];
    TString canName = TString::Format("can%u", eventNumber);
    TString polSuffix = pol == AnitaPol::kHorizontal ? "HPol" : "VPol";
    canName += polSuffix;
    TString canTitle = TString::Format("Event %u - ", eventNumber) + polSuffix;
    wholePad = new TCanvas(canName);
  }
  wholePad->Clear();
  

  TPad* wholeTitlePad = RootTools::makeSubPad(wholePad, 0, 0.95, 1, 1, TString::Format("%d_title", (int)pol));
  TPaveText *wholeTitle = new TPaveText(0, 0, 1, 1);  
  TString wholeTitleText; // = TString::Format(");
  wholeTitleText += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";
  wholeTitleText += " Reconstruction";
  wholeTitle->AddText(wholeTitleText);
  wholeTitle->SetBit(kCanDelete, true);
  wholeTitle->SetLineWidth(0);
  wholeTitle->Draw();
  
  TPad* coarseMapPad = RootTools::makeSubPad(wholePad, 0, 0.75, 1, 0.95, TString::Format("%d_coarse", (int)pol));
  
  InterferometricMap* hCoarse = coarseMaps[pol];  

  TPad* finePeaksAndCoherent = RootTools::makeSubPad(wholePad, 0, 0.35, 1, 0.75, "peaks");

  std::list<InterferometricMap*> drawnFineMaps;
  Double_t coherentMax = -1e9, coherentMin = 1e9;
  
  const int nFine = 3; // maybe discover this dynamically?
  for(int peakInd = 0; peakInd < nFine; peakInd++){
    double yUp = 1 - double(peakInd)/nFine;
    double yLow = yUp - double(1)/nFine;    
    TPad* finePeak = RootTools::makeSubPad(finePeaksAndCoherent, 0, yLow, 0.2, yUp, "fine");

    
    InterferometricMap* h = fineMaps[pol][peakInd];
    if(h){
      h->SetTitleSize(1);
      h->GetXaxis()->SetTitleSize(0.01);
      h->GetYaxis()->SetTitleSize(0.01);
      h->SetLineColor(peakColors[pol][peakInd]);
      h->SetLineWidth(3);

      h->DrawGroup("col");
      drawnFineMaps.push_back(h);      
      
      coarseMapPad->cd();

      const TGraphInteractive* grEdge = h->getEdgeBoxGraph();
      hCoarse->addGuiChild(*grEdge, "l");
    }

    TPad* coherentPad    = RootTools::makeSubPad(finePeaksAndCoherent, 0.2, yLow, 0.6, yUp, "coherent");
    TPad* coherentFftPad = RootTools::makeSubPad(finePeaksAndCoherent, 0.6, yLow,   1, yUp, "coherentFFT");

    // TODO, make this selectable?
    AnalysisWaveform* coherentWave = wfCoherentFiltered[pol][peakInd][0];
    AnitaEventSummary::WaveformInfo &coherentInfo = summary.coherent_filtered[pol][peakInd];
    
    AnalysisWaveform* coherentUnfilteredWave = wfCoherent[pol][peakInd][0];
    AnitaEventSummary::WaveformInfo &coherentUnfilteredInfo = summary.coherent[pol][peakInd];
    
    // AnalysisWaveform* coherentUnfilteredWave = wfDeconvolved[pol][peakInd][0];    
    // AnalysisWaveform* coherentUnfilteredWave = wfDeconvolvedFiltered[pol][peakInd][0];
    if(coherentWave && coherentUnfilteredWave){

      const char* opt = "al";

      // don't want to be able to edit it by accident so copy it...
      TGraphAligned* gr2 = const_cast<TGraphAligned*>(coherentWave->even());
      TGraphAligned* gr4 = const_cast<TGraphAligned*>(coherentUnfilteredWave->even());

      gr2->SetFillColor(0);
      gr4->SetFillColor(0);
      gr2->SetLineColor(peakColors[pol][peakInd]);
      gr4->SetLineColor(kRed);      
      gr4->SetLineStyle(3);

      TGraphInteractive* gr = new TGraphInteractive(gr2, "l");
      gr->SetBit(kCanDelete); // Let ROOT track and handle deletion
      
      TString title = "Coherent Filtered;Time (ns); Amplitiude (mV)";
      gr->SetTitle(title);

      gr4->SetTitle("Coherent Unfiltered");
      gr->addGuiChild(*gr4, "l");
      
      coherentPad->cd();

      gr->DrawGroup(opt);

      TGraphAligned* grPower2 = const_cast<TGraphAligned*>(coherentWave->powerdB());
      TGraphAligned* grPower4 = const_cast<TGraphAligned*>(coherentUnfilteredWave->powerdB());

      grPower2->SetLineColor(peakColors[pol][peakInd]);
      grPower4->SetLineColor(kRed);

      TGraphInteractive* grPower = new TGraphInteractive(grPower2, "l");
      grPower->SetBit(kCanDelete); // Let ROOT track and handle deletion
      const double maxFreqGHz = 1.3;
      grPower->GetXaxis()->SetRangeUser(0, maxFreqGHz);

      title = "PSD Coherent Filtered;Frequency (GHz);Power Spectral Density (dBm/MHz)";
      grPower->SetTitle(title);

      TGraphInteractive* grSlope = new TGraphInteractive(0, NULL, NULL, "l");
      grSlope->SetPoint(0, fSlopeFitStartFreqGHz, coherentInfo.spectrumSlope*fSlopeFitStartFreqGHz + coherentInfo.spectrumIntercept);
      grSlope->SetPoint(1, fSlopeFitEndFreqGHz, coherentInfo.spectrumSlope*fSlopeFitEndFreqGHz + coherentInfo.spectrumIntercept);
      grSlope->SetLineColor(grPower->GetLineColor());
      
      grPower->addGuiChild(grSlope); // pass it a pointer and it takes ownership
      

      grPower4->SetTitle("PSD Coherent Unfiltered");
      grPower->addGuiChild(*grPower4, "l");
      

      TGraphInteractive* grSlope4 = new TGraphInteractive(0, NULL, NULL, "l");
      grSlope4->SetPoint(0, fSlopeFitStartFreqGHz, coherentUnfilteredInfo.spectrumSlope*fSlopeFitStartFreqGHz + coherentUnfilteredInfo.spectrumIntercept);
      grSlope4->SetPoint(1, fSlopeFitEndFreqGHz, coherentUnfilteredInfo.spectrumSlope*fSlopeFitEndFreqGHz + coherentUnfilteredInfo.spectrumIntercept);
      grSlope4->SetLineColor(grPower4->GetLineColor());
      grPower->addGuiChild(grSlope4);


      double df = coherentWave->deltaF();
      TGraphInteractive* grCoherentPeakMarker = new TGraphInteractive(0, NULL, NULL, "p");
      grCoherentPeakMarker->SetMarkerStyle(4);
      grCoherentPeakMarker->SetMarkerSize(0.75);
      grCoherentPeakMarker->SetMarkerColor(grPower->GetLineColor());
      
      for(int i=0; i < AnitaEventSummary::peaksPerSpectrum; i++){
        if(coherentInfo.peakFrequency[i] > 0){
          int powerBin = TMath::Nint(coherentInfo.peakFrequency[i]/df);
          grCoherentPeakMarker->SetPoint(grCoherentPeakMarker->GetN(), grPower4->GetX()[powerBin], grPower4->GetY()[powerBin]);
        }
      }
      grPower->addGuiChild(grCoherentPeakMarker);

      double df2 = coherentUnfilteredWave->deltaF();
      TGraphInteractive* grCoherentUnfilteredPeakMarker = new TGraphInteractive(0, NULL, NULL, "p");
      grCoherentUnfilteredPeakMarker->SetMarkerStyle(4);
      grCoherentUnfilteredPeakMarker->SetMarkerSize(0.75);
      grCoherentUnfilteredPeakMarker->SetMarkerColor(grPower4->GetLineColor());
      
      for(int i=0; i < AnitaEventSummary::peaksPerSpectrum; i++){
        if(coherentUnfilteredInfo.peakFrequency[i] > 0){
          int powerBin = TMath::Nint(coherentUnfilteredInfo.peakFrequency[i]/df2);
          grCoherentUnfilteredPeakMarker->SetPoint(grCoherentUnfilteredPeakMarker->GetN(), grPower4->GetX()[powerBin], grPower4->GetY()[powerBin]);
        }
      }
      grPower->addGuiChild(grCoherentUnfilteredPeakMarker);
          

      coherentFftPad->cd();
      grPower->DrawGroup(opt);
      
    }
    else{
      std::cerr << "missing coherent pointer(s)?\t" << pol << "\t" << peakInd << "\t" << coherentWave << "\t" << coherentUnfilteredWave << std::endl;
    }
  }

  coarseMapPad->cd();
  hCoarse->DrawGroup("colz");  

  if(drawnFineMaps.size() > 0){
    std::list<InterferometricMap*>::iterator it = drawnFineMaps.begin();
    Double_t polMax = -1e9;
    Double_t polMin = 1e9;    
    for(; it!=drawnFineMaps.end(); ++it){
      polMax = (*it)->GetMaximum() > polMax ? (*it)->GetMaximum() : polMax;
      polMin = (*it)->GetMinimum() < polMin ? (*it)->GetMinimum() : polMin;      
    }
    for(it=drawnFineMaps.begin(); it!=drawnFineMaps.end(); ++it){
      (*it)->SetMaximum(polMax);
      (*it)->SetMinimum(polMin);
    }
    if(hCoarse){
      hCoarse->SetMaximum(polMax);
      // hCoarse->SetMinimum(polMin);      
    }
    while(drawnFineMaps.size() > 0){
      drawnFineMaps.pop_back();
    }
    // std::cout << polMax << "\t" << polMin << std::endl;
  }

  TPad* textPad = RootTools::makeSubPad(wholePad, 0, 0, 1, 0.35, "text");
  
  for(int peakInd=0; peakInd < nFine; peakInd++){
    double xlow = double(peakInd)/nFine;
    double xup = xlow + double(1.)/nFine;
    
    TPaveText *title = new TPaveText(xlow, 0.9, xup, 1);
    // title->SetBorderSize(0);
    title->SetBit(kCanDelete, true);
    title->SetTextColor(peakColors[pol][peakInd]);
    title->SetLineColor(peakColors[pol][peakInd]);

    // TString titleText = pol == AnitaPol::kHorizontal ? "HPol" : "VPol";
    TString titleText;
    titleText += TString::Format(" Peak %d", peakInd);    
    title->AddText(titleText);
    
    title->Draw();
    
    TPaveText *pt = new TPaveText(xlow, 0, xup, 0.9);
    // pt->SetBorderSize(0);
    pt->SetBit(kCanDelete, true);
    pt->SetTextColor(peakColors[pol][peakInd]);
    pt->SetLineColor(peakColors[pol][peakInd]);
    
    pt->AddText(TString::Format("Image peak = %4.4lf", summary.peak[pol][peakInd].value));
    pt->AddText(TString::Format("#phi_{fine} = %4.2lf#circ", summary.peak[pol][peakInd].phi));
    pt->AddText(TString::Format("#theta_{fine} = %4.2lf#circ",summary.peak[pol][peakInd].theta));
    pt->AddText(TString::Format("Hilbert peak = %4.2lf mV, ", summary.coherent[pol][peakInd].peakHilbert));            

    pt->AddText(TString::Format("Latitude = %4.2lf #circ", summary.peak[pol][peakInd].latitude));
    pt->AddText(TString::Format("Longitude = %4.2lf #circ", summary.peak[pol][peakInd].longitude));
    pt->AddText(TString::Format("Altitude = %4.2lf #circ", summary.peak[pol][peakInd].altitude));                    
    
    pt->Draw();
  }
  
  wholePad->SetBorderSize(2);
    
}



/** 
 * @brief Coherently summed filtered (un-deconvolved) waveform accessor for external processes. 
 * Only works once per event processed as ownership is transferred to the function caller.
 * 
 * @param pol is the polarisation
 * @param peakInd corresponds to which peak in the interferometric map (0 is highest up to fNumPeaks)
 * @param xPol direct or cross polarisation compared to pol.
 * 
 * @return pointer to the AnalysisWaveform
 */
AnalysisWaveform* Acclaim::AnalysisReco::getCoherentFiltered(AnitaPol::AnitaPol_t pol, Int_t peakInd, bool xPol){
  AnalysisWaveform* wf = wfCoherentFiltered[pol][peakInd][xPol];
  wfCoherentFiltered[pol][peakInd][xPol] = NULL;
  return wf;  
}


/** 
 * @brief Coherently summed (un-filtered, un-deconvolved) waveform accessor for external processes. 
 * Only works once per event processed as ownership is transferred to the function caller.
 * 
 * @param pol is the polarisation
 * @param peakInd corresponds to which peak in the interferometric map (0 is highest up to fNumPeaks)
 * @param xPol direct or cross polarisation compared to pol.
 * 
 * @return pointer to the AnalysisWaveform
 */
AnalysisWaveform* Acclaim::AnalysisReco::getCoherent(AnitaPol::AnitaPol_t pol, Int_t peakInd, bool xPol){
  AnalysisWaveform* wf = wfCoherent[pol][peakInd][xPol];
  wfCoherent[pol][peakInd][xPol] = NULL;
  return wf;  
}



/** 
 * @brief Coherently summed (un-filtered) deconvolved waveform accessor for external processes. 
 * Only works once per event processed as ownership is transferred to the function caller.
 * 
 * @param pol is the polarisation
 * @param peakInd corresponds to which peak in the interferometric map (0 is highest up to fNumPeaks)
 * @param xPol direct or cross polarisation compared to pol.
 * 
 * @return pointer to the AnalysisWaveform
 */
AnalysisWaveform* Acclaim::AnalysisReco::getDeconvolved(AnitaPol::AnitaPol_t pol, Int_t peakInd, bool xPol){
  AnalysisWaveform* wf = wfDeconvolved[pol][peakInd][xPol];
  wfDeconvolved[pol][peakInd][xPol] = NULL;
  return wf;
}



/** 
 * @brief Coherently summed filtered deconvolved waveform accessor for external processes. 
 * Only works once per event processed as ownership is transferred to the function caller.
 * 
 * @param pol is the polarisation
 * @param peakInd corresponds to which peak in the interferometric map (0 is highest up to fNumPeaks)
 * @param xPol direct or cross polarisation compared to pol.
 * 
 * @return pointer to the AnalysisWaveform
 */

AnalysisWaveform* Acclaim::AnalysisReco::getDeconvolvedFiltered(AnitaPol::AnitaPol_t pol, Int_t peakInd, bool xPol){
  AnalysisWaveform* wf = wfDeconvolvedFiltered[pol][peakInd][xPol];
  wfDeconvolvedFiltered[pol][peakInd][xPol] = NULL;
  return wf;
}




/** 
 * Access for internally produce minimally filtered event
 * Only works once per event processed as ownership is transferred to the function caller.
 * 
 * @return The minimally filtered version of the event
 */
FilteredAnitaEvent* Acclaim::AnalysisReco::getEvMin(){
  FilteredAnitaEvent* f = fEvMin;
  fEvMin = NULL;
  return f;
}



/** 
 * Access for internally produced minimally filtered, deconvolved event
 * Only works once per event processed as ownership is transferred to the function caller.
 * 
 * @return The minimally filtered, deconvolved version of the event
 */
FilteredAnitaEvent* Acclaim::AnalysisReco::getEvMinDeco(){
  FilteredAnitaEvent* f = fEvMinDeco;
  fEvMinDeco = NULL;
  return f;
}



/** 
 * Access for internally produced filtered, deconvolved event
 * Only works once per event processed as ownership is transferred to the function caller.
 * 
 * @return The deconvolved version of the event
 */
FilteredAnitaEvent* Acclaim::AnalysisReco::getEvDeco(){
  FilteredAnitaEvent* f = fEvDeco;
  fEvDeco = NULL;
  return f;
}
