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

  if(fSpawnedCrossCorrelator && fCrossCorr){
    delete fCrossCorr;
    fCrossCorr = NULL;
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




void Acclaim::AnalysisReco::fillWaveformInfo(AnitaPol::AnitaPol_t pol,
                                             AnitaEventSummary::WaveformInfo& info,
                                             const FilteredAnitaEvent* fEv,
                                             AnalysisWaveform** waveStore,
                                             InterferometricMap* h,
                                             NoiseMonitor* noiseMonitor) {


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
  {
    double noise = 0;
    for(unsigned antInd = 0; antInd < theAnts.size(); antInd++){
      int ant = theAnts[antInd];

      double thisRMS = 0;
      if(noiseMonitor){
	thisRMS = noiseMonitor->getRMS(pol, ant, fEv->getHeader()->realTime);
      }
      else{
	static int warnCount = 0;
	const TGraphAligned* gr = fEv->getFilteredGraph(ant, pol)->even();
	thisRMS = gr->GetRMS(2);
	const int numWarnings = 10;
	if(warnCount < numWarnings){
	  std::cerr << "Warning ("<< warnCount + 1 << "/" << numWarnings << ") in "
		    << __PRETTY_FUNCTION__ << ", getting RMS from entire waveform rather than NoiseMonitor."
		    << " this many overestimate the noise... thisRMS = " << thisRMS << std::endl;
	  warnCount++;
	}
      }
      noise += thisRMS;
      
    }
    noise /= theAnts.size();
    
    info.snr = largestPeakToPeak/(2*noise);

    // std::cout << "SNR = " << info.snr << std::endl;
  }

  
  const TGraphAligned* grHilbert = coherentWave->hilbertEnvelope();
  info.peakHilbert = TMath::MaxElement(grHilbert->GetN(), grHilbert->GetY());

  const TGraphAligned* gr = coherentWave->even();
  info.peakVal = TMath::MaxElement(gr->GetN(), gr->GetY());

  // Double_t localMaxVolts, localMinVolts, localMaxTime, localMinTime;
  // RootTools::getLocalMaxToMin(gr, localMaxVolts, localMaxTime, localMinVolts, localMinTime);
  // info.localMaxToMin = localMaxVolts - localMinVolts;
  // info.localMaxToMinTime = localMaxTime - localMinTime;

  // Int_t gMaxInd = TMath::LocMax(gr->GetN(), gr->GetY());
  // Int_t gMinInd = TMath::LocMin(gr->GetN(), gr->GetY());
  // info.globalMaxToMin = gr->GetY()[gMaxInd] - gr->GetY()[gMinInd];
  // info.globalMaxToMinTime = gr->GetX()[gMaxInd] - gr->GetX()[gMinInd];
  

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

  if(fFillSpectrumInfo){
    const double maxFreqIfIHadNotIntepolated = FancyFFTs::getNumFreqs(NUM_SAMP)*(1./(NOMINAL_SAMPLING_DELTAT*NUM_SAMP));
    const double df = coherentWave->deltaF(); // the actual deltaF
    int numGoodFreqBins = floor(maxFreqIfIHadNotIntepolated/df);

    // std::cout << maxFreq << std::endl;
    const TGraphAligned* grPow = coherentWave->power();

    if(numGoodFreqBins > grPow->GetN()){
      std::cerr << "Warning in " << __PRETTY_FUNCTION__ << " for run " << fCurrentRun << ", eventNumber " << fCurrentEventNumber
		<< ", expected at least " << numGoodFreqBins << " frequency bins, only got " << grPow->GetN()
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
      std::cerr << "Warning in " << __PRETTY_FUNCTION__ << " for run "
		<< fCurrentRun << ", eventNumber " << fCurrentEventNumber
		<< ", got gradDenominator = " << gradDenominator << " will get NaN or Inf..." << std::endl;
    }
    info.spectrumSlope = gradNumerator/gradDenominator;
    info.spectrumIntercept = mean_pow_db - info.spectrumSlope*mean_f;

  }
}











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



void Acclaim::AnalysisReco::process(const FilteredAnitaEvent * fEv, AnitaEventSummary * sum, NoiseMonitor* noiseMonitor, TruthAnitaEvent* truth) {

  fCurrentEventNumber = fEv->getHeader()->eventNumber;
  fCurrentRun = fEv->getHeader()->run;

  // spawn a CrossCorrelator if we need one
  if(!fCrossCorr){
    fCrossCorr = new CrossCorrelator();
    fSpawnedCrossCorrelator = true;
    fDtCache.init(fCrossCorr, this);
  }

  // FilteredAnitaEvent fEvMin(fEv->getUsefulAnitaEvent(), fMinFilter, fEv->getGPS(), fEv->getHeader(), false); // just the minimum filter
  // FilteredAnitaEvent fEvMinDeco(fEv->getUsefulAnitaEvent(), fMinDecoFilter, fEv->getGPS(), fEv->getHeader(), false); // minimum + deconvolution
  // FilteredAnitaEvent fEvDeco(fEv, fMinDecoFilter); // extra deconvolution
  nicelyDeleteInternalFilteredEvents();
  fEvMin = new FilteredAnitaEvent(fEv->getUsefulAnitaEvent(), fMinFilter, fEv->getGPS(), fEv->getHeader(), false); // just the minimum filter
  fEvMinDeco = new FilteredAnitaEvent(fEv->getUsefulAnitaEvent(), fMinDecoFilter, fEv->getGPS(), fEv->getHeader(), false); // minimum + deconvolution
  fEvDeco = new FilteredAnitaEvent(fEv, fMinDecoFilter); // extra deconvolution

  if(fFillChannelInfo){
    fillChannelInfo(fEv, sum);
  }

  sum->eventNumber = fEv->getHeader()->eventNumber;

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

    sum->nPeaks[pol] = fNumPeaks;

    for(Int_t peakInd=0; peakInd < fNumPeaks; peakInd++){
      reconstructZoom(pol, peakInd, coarseMapPeakPhiDegs.at(peakInd), coarseMapPeakThetaDegs.at(peakInd));

      InterferometricMap* h = fineMaps[pol][peakInd];
      if(h){

        // for plotting
        h->addGpsInfo(fEv->getGPS());
        h->addTruthInfo(truth);


        h->getPeakInfo(sum->peak[pol][peakInd].value,
                       sum->peak[pol][peakInd].phi,
                       sum->peak[pol][peakInd].theta,
                       sum->peak[pol][peakInd].chisq);

        // fill in difference between rough and fine
        sum->peak[pol][peakInd].dphi_rough = sum->peak[pol][peakInd].phi - coarseMapPeakPhiDegs.at(peakInd);
        sum->peak[pol][peakInd].dtheta_rough = sum->peak[pol][peakInd].theta - coarseMapPeakThetaDegs.at(peakInd);

        // based on Cosmin's comments in AnitaAnalysisSummary.h
        sum->peak[pol][peakInd].phi_separation = RootTools::getDeltaAngleDeg(sum->peak[pol][peakInd].phi, sum->peak[pol][0].phi);

        // hwAngle, is the angle between the peak and the nearest L3 phi-sector trigger of that polarisation        
        setTriggerInfoFromPeakPhi(fEv->getHeader(), pol, h->getPeakPhiSector(), sum->peak[pol][peakInd]);

        // coherent
        // coherent_filtered
        // deconvolved
        // deconvolved_filtered
      
        fillWaveformInfo(pol, sum->coherent_filtered[pol][peakInd],    fEv,      fCoherentFiltered[pol][peakInd],    h, noiseMonitor);
        fillWaveformInfo(pol, sum->deconvolved_filtered[pol][peakInd], fEvDeco,  fDeconvolvedFiltered[pol][peakInd], h, noiseMonitor);

	h->setResolutionEstimateFromWaveformSNR(sum->deconvolved_filtered[pol][peakInd].snr);

	if(fFillUnfiltered){
	  fillWaveformInfo(pol, sum->coherent[pol][peakInd], fEvMin,        fCoherent[pol][peakInd],            h, noiseMonitor);	
	  fillWaveformInfo(pol, sum->deconvolved[pol][peakInd], fEvMinDeco, fDeconvolved[pol][peakInd],         h, noiseMonitor);
	}

	// if(sum->eventNumber==60840043){
	//   double x = sum->trainingDeconvolvedFiltered().fracPowerWindowGradient();
	//   for(int i=0; i < 5; i++){
	//     std::cout << sum->trainingDeconvolvedFiltered().fracPowerWindowBegins[i] << "\t" << sum->trainingDeconvolvedFiltered().fracPowerWindowEnds[i] << std::endl;
	//   }
	//   std::cout << x << std::endl;
	// }

        if(fEv->getGPS() != NULL){
          UsefulAdu5Pat usefulPat(fEv->getGPS());
          int success = 0;
          if(sum->peak[pol][peakInd].theta < 0){ // work around for bug in traceBackToContinent
            Double_t phiWave = TMath::DegToRad()*sum->peak[pol][peakInd].phi;
            Double_t thetaWave = -1*TMath::DegToRad()*sum->peak[pol][peakInd].theta;

            success = usefulPat.traceBackToContinent3(phiWave, thetaWave,
						      &sum->peak[pol][peakInd].longitude,
						      &sum->peak[pol][peakInd].latitude,
						      &sum->peak[pol][peakInd].altitude,
						      &sum->peak[pol][peakInd].theta_adjustment_needed);
          }
          if(success==0){
            sum->peak[pol][peakInd].longitude = -9999;
            sum->peak[pol][peakInd].latitude = -9999;
            sum->peak[pol][peakInd].altitude = -9999;
            sum->peak[pol][peakInd].theta_adjustment_needed = -9999;
            sum->peak[pol][peakInd].distanceToSource = -9999;
          }
          else{
            sum->peak[pol][peakInd].distanceToSource = SPEED_OF_LIGHT_NS*usefulPat.getTriggerTimeNsFromSource(sum->peak[pol][peakInd].latitude,
                                                                                                                       sum->peak[pol][peakInd].longitude,
                                                                                                                       sum->peak[pol][peakInd].altitude);
          }
        }
      }
    }
  }

  // keep internal copy of summary
  fSummary = (*sum);
}




void Acclaim::AnalysisReco::fillPowerFlags(const FilteredAnitaEvent* fEv, AnitaEventSummary::EventFlags& flags){

  const double bandsLowGHz[AnitaEventSummary::numBlastPowerBands] = {0.15-1e-10, 0};
  const double bandsHighGHz[AnitaEventSummary::numBlastPowerBands] = {0.25+1e-10, 999};

  // reset the values
  for(int band=0; band < AnitaEventSummary::numBlastPowerBands; band++){
    flags.middleOrBottomPower[band] = 0;
    flags.middleOrBottomAnt[band] = -1;
    flags.middleOrBottomPol[band] = 2;
    flags.topPower[band] = 0;
  }
  
  for(Int_t polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
    for(UInt_t phi=0; phi < NUM_PHI; phi++){
      for(UInt_t ring=AnitaRing::kMiddleRing; ring < AnitaRing::kNotARing; ring++){
  	Int_t ant = ring*NUM_PHI + phi;

	
	const AnalysisWaveform* wf = fEv->getRawGraph(ant, pol);
	const TGraphAligned* grPow = wf->power();
	const double df_GHz = grPow->GetX()[1] - grPow->GetX()[0];

	Double_t powThisAnt[AnitaEventSummary::numBlastPowerBands] = {0};

	for(int i=0; i < grPow->GetN(); i++){
	  const double f_GHz = grPow->GetX()[i];
	  for(int band=0; band < AnitaEventSummary::numBlastPowerBands; band++){
	    if(f_GHz >= bandsLowGHz[band] && f_GHz < bandsHighGHz[band]){
	      powThisAnt[band] +=grPow->GetY()[i]*df_GHz;
	    }
	  }
	}

	for(int band=0; band < AnitaEventSummary::numBlastPowerBands; band++){	
	  if(powThisAnt[band] > flags.middleOrBottomPower[band]){
	    flags.middleOrBottomPower[band] = powThisAnt[band];
	    flags.middleOrBottomAnt[band] = ant;
	    flags.middleOrBottomPol[band] = polInd;
	  }
	}
      }
    }
  }

  
  for(int band=0; band < AnitaEventSummary::numBlastPowerBands; band++){
    int ant = flags.middleOrBottomAnt[band] % NUM_PHI;
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) flags.middleOrBottomPol[band];
    const AnalysisWaveform* wf = fEv->getRawGraph(ant, pol);
    const TGraphAligned* grPow = wf->power();
    const double df_GHz = grPow->GetX()[1] - grPow->GetX()[0];
    
    for(int i=0; i < grPow->GetN(); i++){
      const double f_GHz = grPow->GetX()[i];
      for(int band=0; band < AnitaEventSummary::numBlastPowerBands; band++){
	if(f_GHz >= bandsLowGHz[band] && f_GHz < bandsHighGHz[band]){
	  flags.topPower[band] +=grPow->GetY()[i]*df_GHz;
	}
      }
    }
  }  
}




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


void Acclaim::AnalysisReco::initializeInternals(){

  fCrossCorr = NULL;
  fSpawnedCrossCorrelator = false;
  fWhichResponseDir = 0;
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  for(Int_t pol=0; pol < AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant<NUM_SEAVEYS; ant++){
      fRArray[pol].push_back(geom->getAntR(ant, AnitaPol::AnitaPol_t(pol)));
      fZArray[pol].push_back(geom->getAntZ(ant, AnitaPol::AnitaPol_t(pol)));
      fPhiArrayDeg[pol].push_back(geom->getAntPhiPositionRelToAftFore(ant, AnitaPol::AnitaPol_t(pol))*TMath::RadToDeg());
    }
  }

  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    coarseMaps[polInd] = NULL;
    for(int peakInd=0; peakInd < AnitaEventSummary::maxDirectionsPerPol; peakInd++){
      fineMaps[polInd][peakInd] = NULL;
      for(int xPol=0; xPol < AnitaPol::kNotAPol; xPol++){
        fCoherentFiltered[polInd][peakInd][xPol] = NULL;
        fCoherent[polInd][peakInd][xPol] = NULL;
        fDeconvolved[polInd][peakInd][xPol] = NULL;
        fDeconvolvedFiltered[polInd][peakInd][xPol] = NULL;
      }
    }
  }

  fDebug = 0;
  fUseOffAxisDelay = 1;
  fCoherentDeltaPhi = 2;
  fLastCoherentDeltaPhi = -1;
  fNumPeaks = 3;
  fCoherentDtNs = 0.01;
  fSlopeFitStartFreqGHz = 0.18;
  fSlopeFitEndFreqGHz = 1.3;
  fFillChannelInfo = 0;
  fFillSpectrumInfo = 0;
  fFillUnfiltered = 0;
  fMeanPowerFlagLowFreqGHz = 0;
  fMeanPowerFlagHighFreqGHz = 0;
  fCurrentEventNumber = 0;
  fCurrentRun = 0;

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

  chooseAntennasForCoherentlySumming(fCoherentDeltaPhi);
  
}




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



Double_t Acclaim::AnalysisReco::relativeOffAxisDelay(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
                                                     Double_t phiDeg) const {

  Double_t deltaPhiDeg1 = RootTools::getDeltaAngleDeg(fPhiArrayDeg[pol].at(ant1), phiDeg);
  Double_t deltaPhiDeg2 = RootTools::getDeltaAngleDeg(fPhiArrayDeg[pol].at(ant2), phiDeg);
  Double_t delay1 = singleAntennaOffAxisDelay(deltaPhiDeg1);
  Double_t delay2 = singleAntennaOffAxisDelay(deltaPhiDeg2);
  return (delay1-delay2);
}




Double_t Acclaim::AnalysisReco::getDeltaTExpected(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t phiWave, Double_t thetaWave) const {

  // Double_t tanThetaW = tan(thetaWave);
  Double_t tanThetaW = tan(-1*thetaWave);
  Double_t part1 = fZArray[pol].at(ant1)*tanThetaW - fRArray[pol].at(ant1)*cos(phiWave-TMath::DegToRad()*fPhiArrayDeg[pol].at(ant1));
  Double_t part2 = fZArray[pol].at(ant2)*tanThetaW - fRArray[pol].at(ant2)*cos(phiWave-TMath::DegToRad()*fPhiArrayDeg[pol].at(ant2));
  Double_t tdiff = 1e9*((cos(thetaWave) * (part2 - part1))/SPEED_OF_LIGHT); // Returns time in ns

  return tdiff;
}




void Acclaim::AnalysisReco::insertPhotogrammetryGeometry(){
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->usePhotogrammetryNumbers(1);
  for(Int_t pol=0; pol < AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant<NUM_SEAVEYS; ant++){
      fRArray[pol].at(ant) = geom->getAntR(ant, AnitaPol::AnitaPol_t(pol));
      fZArray[pol].at(ant) = geom->getAntZ(ant, AnitaPol::AnitaPol_t(pol));
      fPhiArrayDeg[pol].at(ant) = geom->getAntPhiPositionRelToAftFore(ant, AnitaPol::AnitaPol_t(pol))*TMath::RadToDeg();
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
    fSpawnedCrossCorrelator = true;
  }
  
  fDtCache.init(fCrossCorr, this, true);
  geom->usePhotogrammetryNumbers(0);

}






Acclaim::InterferometricMap* Acclaim::AnalysisReco::getMap(AnitaPol::AnitaPol_t pol){
  InterferometricMap* h = coarseMaps[pol];
  coarseMaps[pol] = NULL;
  return h;
}






Acclaim::InterferometricMap* Acclaim::AnalysisReco::getZoomMap(AnitaPol::AnitaPol_t pol, UInt_t peakInd){

  InterferometricMap* h = fineMaps[pol][peakInd];
  if(h==NULL){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << " for run "
	      << fCurrentRun << ", eventNumber " << fCurrentEventNumber
	      << ", unable to find fineMap with pol " << pol
              << " for peakInd = " << peakInd << " about to return NULL." << std::endl;
  }
  fineMaps[pol][peakInd] = NULL;

  return h;
}







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

  coarseMaps[pol]->Fill(pol, fCrossCorr, &fDtCache);
}




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
  h->Fill(pol, fCrossCorr, &fDtCache);  

  // std::cout << h->GetName() << std::endl;
  // std::cout << "in " << __PRETTY_FUNCTION__ << ": " << fineMaps[pol][peakIndex] << std::endl;
  if(fineMaps[pol][peakIndex]){
    delete fineMaps[pol][peakIndex];
  }
  fineMaps[pol][peakIndex] = h;
}


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

void Acclaim::AnalysisReco::chooseAntennasForCoherentlySumming(int coherentDeltaPhi){

  // std::cerr << __PRETTY_FUNCTION__ << std::endl;

  if(coherentDeltaPhi!=fLastCoherentDeltaPhi){

    for(int peakPhiSector = 0; peakPhiSector < NUM_PHI; peakPhiSector++){
      fPhiSectorToAnts[peakPhiSector] = std::vector<Int_t>();

      // fPhiSectorToAnts[peakPhiSector].clear();
      for(int deltaPhiSect=-fCoherentDeltaPhi; deltaPhiSect <= fCoherentDeltaPhi; deltaPhiSect++){

        Int_t phiSector = peakPhiSector + deltaPhiSect;
        phiSector = phiSector < 0        ? phiSector + NUM_PHI : phiSector;
        phiSector = phiSector >= NUM_PHI ? phiSector - NUM_PHI : phiSector;

        for(int ring = 0; ring < AnitaRing::kNotARing; ring++){
          int ant = AnitaGeomTool::getAntFromPhiRing(phiSector, AnitaRing::AnitaRing_t(ring));
          fPhiSectorToAnts[peakPhiSector].push_back(ant);
        }
      }
    }
    fLastCoherentDeltaPhi = fCoherentDeltaPhi;
  }
}

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
      std::cerr << "Warning in " << __PRETTY_FUNCTION__ << " for run "
		<< fCurrentRun << ", eventNumber " << fCurrentEventNumber << ", the waveform max/min aren't sensible?" << std::endl;
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

size_t Acclaim::AnalysisReco::checkWavesAndDtsMatch(std::vector<const AnalysisWaveform*>& waves, std::vector<Double_t>& dts) {
  size_t s = waves.size();
  if(s < 1){    
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << " for run " << fCurrentRun
	      << ", eventNumber " << fCurrentEventNumber << ", nothing to sum." << std::endl;
  }  
  else if(s != dts.size()){    
    const char* action = dts.size() < waves.size() ? "padding" : "trimming";
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << " for run " << fCurrentRun
	      << ", eventNumber " << fCurrentEventNumber << ", unequal vectors (waves.size() = "
	      << waves.size() << ", dts.size() = " << dts.size() << ")\n"
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

    if(fDebug){
      std::cerr << "Debug in " << __PRETTY_FUNCTION__ << ", waves[" << i << "] has even with " << grU->GetN() << " points. "; 
      const TGraphAligned* grA = waves[i]->uneven();

      std::cerr << "(uneven has " << grA->GetN() << " points)." << std::endl;
      for(int samp=1; samp < grU->GetN(); samp++){
        if(grU->GetX()[samp] - grU->GetX()[samp-1] <= 0){
          std::cerr << "Debug in " << __PRETTY_FUNCTION__ << ", x not monotonically increasing!" << std::endl;
	  std::cerr << "x[" << samp-1 << "] = " << grU->GetX()[samp-1] << std::endl;
	  std::cerr << "x[" << samp   << "] = " << grU->GetX()[samp]   << std::endl;
          break;
        }
      }
    }

    int numNan = 0;
    for(int samp=0; samp < grCoherent->GetN(); samp++){
      double t = grCoherent->GetX()[samp];

      double tPlusDt = t + dts[i];
      if(tPlusDt > t0 && tPlusDt < tN){
        double y = waves[i]->evalEven(tPlusDt, AnalysisWaveform::EVAL_AKIMA);
        if(TMath::IsNaN(y)){ // Hopefully this is a bit redundent
	  numNan++;
          // std::cerr << "Warning in " << __PRETTY_FUNCTION__ << " for run " << fCurrentRun << ", eventNumber " << fCurrentEventNumber << " interpolation returned " << y << " ignoring sample" << std::endl;
          // std::cerr << "Warning in " << __PRETTY_FUNCTION__ << " for run " << fCurrentRun << ", eventNumber " << fCurrentEventNumber << " interpolation returned " << y << " for time " << tPlusDt << " for sample " << samp << " of " << grCoherent->GetN() << std::endl;
        }
        else{
          grCoherent->GetY()[samp] += y;
        }
      }
    }

    if(numNan > 0){
      std::cerr << "Warning in " << __PRETTY_FUNCTION__ << " for run "
		<< fCurrentRun << ", eventNumber " << fCurrentEventNumber << " interpolation returned "
		<< numNan << " NaNs for wave " << i << std::endl;
    }
  }

  for(int samp=0; samp < grCoherent->GetN(); samp++){
    grCoherent->GetY()[samp]/=waves.size();
  };    

  return coherentWave;
}

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
    AnalysisWaveform* primaryWave = fCoherentFiltered[pol][peakInd][0];
    AnitaEventSummary::WaveformInfo &primaryInfo = fSummary.coherent_filtered[pol][peakInd];
    
    AnalysisWaveform* secondaryWave = fDeconvolvedFiltered[pol][peakInd][0];
    AnitaEventSummary::WaveformInfo &secondaryInfo = fSummary.deconvolved_filtered[pol][peakInd];

    // AnalysisWaveform* secondaryWave = fCoherent[pol][peakInd][0];
    // AnitaEventSummary::WaveformInfo &secondaryInfo = fSummary.coherent[pol][peakInd];
    
    // AnalysisWaveform* secondaryWave = fDeconvolved[pol][peakInd][0];
    // AnalysisWaveform* secondaryWave = fDeconvolvedFiltered[pol][peakInd][0];
    if(primaryWave && secondaryWave){

      const char* opt = "al";

      // don't want to be able to edit it by accident so copy it...
      TGraphAligned* gr2 = const_cast<TGraphAligned*>(primaryWave->even());
      TGraphAligned* gr4 = const_cast<TGraphAligned*>(secondaryWave->even());

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

      TGraphAligned* grPower2 = const_cast<TGraphAligned*>(primaryWave->powerdB());
      TGraphAligned* grPower4 = const_cast<TGraphAligned*>(secondaryWave->powerdB());

      grPower2->SetLineColor(peakColors[pol][peakInd]);
      grPower4->SetLineColor(kRed);

      TGraphInteractive* grPower = new TGraphInteractive(grPower2, "l");
      grPower->SetBit(kCanDelete); // Let ROOT track and handle deletion
      const double maxFreqGHz = 1.3;
      grPower->GetXaxis()->SetRangeUser(0, maxFreqGHz);

      title = "PSD Coherent Filtered;Frequency (GHz);Power Spectral Density (dBm/MHz)";
      grPower->SetTitle(title);

      TGraphInteractive* grSlope = new TGraphInteractive(0, NULL, NULL, "l");
      grSlope->SetPoint(0, fSlopeFitStartFreqGHz, primaryInfo.spectrumSlope*fSlopeFitStartFreqGHz + primaryInfo.spectrumIntercept);
      grSlope->SetPoint(1, fSlopeFitEndFreqGHz, primaryInfo.spectrumSlope*fSlopeFitEndFreqGHz + primaryInfo.spectrumIntercept);
      grSlope->SetLineColor(grPower->GetLineColor());
      
      grPower->addGuiChild(grSlope); // pass it a pointer and it takes ownership
      

      grPower4->SetTitle("PSD Coherent Unfiltered");
      grPower->addGuiChild(*grPower4, "l");
      

      TGraphInteractive* grSlope4 = new TGraphInteractive(0, NULL, NULL, "l");
      grSlope4->SetPoint(0, fSlopeFitStartFreqGHz, secondaryInfo.spectrumSlope*fSlopeFitStartFreqGHz + secondaryInfo.spectrumIntercept);
      grSlope4->SetPoint(1, fSlopeFitEndFreqGHz, secondaryInfo.spectrumSlope*fSlopeFitEndFreqGHz + secondaryInfo.spectrumIntercept);
      grSlope4->SetLineColor(grPower4->GetLineColor());
      grPower->addGuiChild(grSlope4);


      double df = primaryWave->deltaF();
      TGraphInteractive* grCoherentPeakMarker = new TGraphInteractive(0, NULL, NULL, "p");
      grCoherentPeakMarker->SetMarkerStyle(4);
      grCoherentPeakMarker->SetMarkerSize(0.75);
      grCoherentPeakMarker->SetMarkerColor(grPower->GetLineColor());
      
      for(int i=0; i < AnitaEventSummary::peaksPerSpectrum; i++){
        if(primaryInfo.peakFrequency[i] > 0){
          int powerBin = TMath::Nint(primaryInfo.peakFrequency[i]/df);
          grCoherentPeakMarker->SetPoint(grCoherentPeakMarker->GetN(), grPower4->GetX()[powerBin], grPower4->GetY()[powerBin]);
        }
      }
      grPower->addGuiChild(grCoherentPeakMarker);

      double df2 = secondaryWave->deltaF();
      TGraphInteractive* grCoherentUnfilteredPeakMarker = new TGraphInteractive(0, NULL, NULL, "p");
      grCoherentUnfilteredPeakMarker->SetMarkerStyle(4);
      grCoherentUnfilteredPeakMarker->SetMarkerSize(0.75);
      grCoherentUnfilteredPeakMarker->SetMarkerColor(grPower4->GetLineColor());
      
      for(int i=0; i < AnitaEventSummary::peaksPerSpectrum; i++){
        if(secondaryInfo.peakFrequency[i] > 0){
          int powerBin = TMath::Nint(secondaryInfo.peakFrequency[i]/df2);
          grCoherentUnfilteredPeakMarker->SetPoint(grCoherentUnfilteredPeakMarker->GetN(), grPower4->GetX()[powerBin], grPower4->GetY()[powerBin]);
        }
      }
      grPower->addGuiChild(grCoherentUnfilteredPeakMarker);
          

      coherentFftPad->cd();
      grPower->DrawGroup(opt);
      
    }
    else{
      std::cerr << "missing coherent pointer(s)?\t" << pol << "\t" << peakInd << "\t" << primaryWave << "\t" << secondaryWave << std::endl;
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
    
    pt->AddText(TString::Format("Image peak = %4.4lf", fSummary.peak[pol][peakInd].value));
    pt->AddText(TString::Format("#phi_{fine} = %4.2lf#circ", fSummary.peak[pol][peakInd].phi));
    pt->AddText(TString::Format("#theta_{fine} = %4.2lf#circ",fSummary.peak[pol][peakInd].theta));
    pt->AddText(TString::Format("Hilbert peak = %4.2lf mV, ", fSummary.coherent[pol][peakInd].peakHilbert));            

    pt->AddText(TString::Format("Latitude = %4.2lf #circ", fSummary.peak[pol][peakInd].latitude));
    pt->AddText(TString::Format("Longitude = %4.2lf #circ", fSummary.peak[pol][peakInd].longitude));
    pt->AddText(TString::Format("Altitude = %4.2lf #circ", fSummary.peak[pol][peakInd].altitude));                    
    
    pt->Draw();
  }
  
  wholePad->SetBorderSize(2);
    
}



AnalysisWaveform* Acclaim::AnalysisReco::getCoherentFiltered(AnitaPol::AnitaPol_t pol, Int_t peakInd, bool xPol){
  AnalysisWaveform* wf = fCoherentFiltered[pol][peakInd][xPol];
  fCoherentFiltered[pol][peakInd][xPol] = NULL;
  return wf;  
}


AnalysisWaveform* Acclaim::AnalysisReco::getCoherent(AnitaPol::AnitaPol_t pol, Int_t peakInd, bool xPol){
  AnalysisWaveform* wf = fCoherent[pol][peakInd][xPol];
  fCoherent[pol][peakInd][xPol] = NULL;
  return wf;  
}



AnalysisWaveform* Acclaim::AnalysisReco::getDeconvolved(AnitaPol::AnitaPol_t pol, Int_t peakInd, bool xPol){
  AnalysisWaveform* wf = fDeconvolved[pol][peakInd][xPol];
  fDeconvolved[pol][peakInd][xPol] = NULL;
  return wf;
}



AnalysisWaveform* Acclaim::AnalysisReco::getDeconvolvedFiltered(AnitaPol::AnitaPol_t pol, Int_t peakInd, bool xPol){
  AnalysisWaveform* wf = fDeconvolvedFiltered[pol][peakInd][xPol];
  fDeconvolvedFiltered[pol][peakInd][xPol] = NULL;
  return wf;
}




FilteredAnitaEvent* Acclaim::AnalysisReco::getEvMin(){
  FilteredAnitaEvent* f = fEvMin;
  fEvMin = NULL;
  return f;
}



FilteredAnitaEvent* Acclaim::AnalysisReco::getEvMinDeco(){
  FilteredAnitaEvent* f = fEvMinDeco;
  fEvMinDeco = NULL;
  return f;
}



FilteredAnitaEvent* Acclaim::AnalysisReco::getEvDeco(){
  FilteredAnitaEvent* f = fEvDeco;
  fEvDeco = NULL;
  return f;
}
