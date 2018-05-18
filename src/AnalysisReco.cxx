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

  if(fSpawnedCrossCorrelator && (fCrossCorr != nullptr)){
    delete fCrossCorr;
    fCrossCorr = nullptr;
  }  
  nicelyDeleteInternalFilteredEvents();  
}


void Acclaim::AnalysisReco::nicelyDeleteInternalFilteredEvents(){
  if(fEvMin != nullptr){
    delete fEvMin;
    fEvMin = nullptr;
  }
  if(fEvMinDeco != nullptr){
    delete fEvMinDeco;
    fEvMinDeco = nullptr;
  }
  if(fEvDeco != nullptr){
    delete fEvDeco;
    fEvDeco = nullptr;
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
  if(waveStore[0] != nullptr){
    delete waveStore[0];
  }
  waveStore[0] = coherentWave;
  
  info.snr = 0;
  {
    double noise = 0;
    for(int ant : theAnts){
      double thisRMS = 0;
      if(noiseMonitor != nullptr){
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
  AnalysisWaveform* xPolCoherentWave = coherentlySum(fEv, xPol, theAnts, phiDeg, thetaDeg, &largestPeakToPeakXPol, &gr->GetX()[0]); // cross polarisation

  // Internal storage
  if(waveStore[1] != nullptr){
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
                             grH->GetY(), grH_hilbert->GetY(),
                             grV->GetY(), grV_hilbert->GetY(),
                             &info.I,      &info.Q,     &info.U,     &info.V,
			     &info.max_dI, &info.max_dQ,&info.max_dU,&info.max_dV,
			     true);

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
  info.impulsivityMeasure = impulsivity::impulsivityMeasure(coherentWave, nullptr, maxIndex, true);

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

  if(fFillSpectrumInfo != 0){
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
    auto pol = (AnitaPol::AnitaPol_t) polInd;
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
  if(fCrossCorr == nullptr){
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

  if(fFillChannelInfo != 0){
    fillChannelInfo(fEv, sum);
  }

  sum->eventNumber = fEv->getHeader()->eventNumber;

  chooseAntennasForCoherentlySumming(fCoherentDeltaPhi);

  for(Int_t polInd=0; polInd < AnitaPol::kNotAPol; polInd++){

    auto pol = (AnitaPol::AnitaPol_t) polInd;

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
      if(h != nullptr){

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

	if(fFillUnfiltered != 0){
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

        if(fEv->getGPS() != nullptr){
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
    auto pol = (AnitaPol::AnitaPol_t) polInd;
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
    auto pol = (AnitaPol::AnitaPol_t) flags.middleOrBottomPol[band];
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

  peak.triggered = (header->isInL3Pattern(peakPhiSector, pol) != 0);
  peak.triggered_xpol = (header->isInL3Pattern(peakPhiSector, xPol) != 0);

  peak.masked = (header->isInPhiMaskOffline(peakPhiSector, pol) != 0) || (header->isInL1MaskOffline(peakPhiSector, pol) != 0);
  peak.masked_xpol = (header->isInPhiMaskOffline(peakPhiSector, xPol) != 0) || (header->isInL1MaskOffline(peakPhiSector, xPol) != 0);
}


void Acclaim::AnalysisReco::initializeInternals(){

  fCrossCorr = nullptr;
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
    coarseMaps[polInd] = nullptr;
    for(int peakInd=0; peakInd < AnitaEventSummary::maxDirectionsPerPol; peakInd++){
      fineMaps[polInd][peakInd] = nullptr;
      for(int xPol=0; xPol < AnitaPol::kNotAPol; xPol++){
        fCoherentFiltered[polInd][peakInd][xPol] = nullptr;
        fCoherent[polInd][peakInd][xPol] = nullptr;
        fDeconvolved[polInd][peakInd][xPol] = nullptr;
        fDeconvolvedFiltered[polInd][peakInd][xPol] = nullptr;
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

  fDrawNPeaks = 3;
  fDrawDomain = kTimeDomain;
  fDrawCoherent=1;
  fDrawDedispersed=1;
  fDrawXPol=1;
  fDrawXPolDedispersed=1;

  const TString minFiltName = "Minimum";
  fMinFilter = Filters::findDefaultStrategy(minFiltName);
  if(fMinFilter == nullptr){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", unable to find default filter " << minFiltName << std::endl;
  }  

  const TString minDecoFiltName = "Deconvolve";
  fMinDecoFilter = Filters::findDefaultStrategy(minDecoFiltName);
  if(fMinDecoFilter == nullptr){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", unable to find default filter " << minDecoFiltName << std::endl;
  }

  fEvMin = nullptr;
  fEvMinDeco = nullptr;
  fEvDeco = nullptr;

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


Double_t Acclaim::AnalysisReco::getDeltaTExpected(AnitaPol::AnitaPol_t pol, Int_t ant, Double_t phiWave, Double_t thetaWave) const {

  // Double_t tanThetaW = tan(thetaWave);
  Double_t tanThetaW = tan(-1*thetaWave);
  Double_t part = fZArray[pol].at(ant)*tanThetaW - fRArray[pol].at(ant)*cos(phiWave-TMath::DegToRad()*fPhiArrayDeg[pol].at(ant));
  Double_t tdiff = 1e9*((cos(thetaWave) * (part))/SPEED_OF_LIGHT); // Returns time in ns

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
  for(auto & delaysPerSurf : cal->relativePhaseCenterToAmpaDelays){
    for(double& delayThisChannel : delaysPerSurf){
      delayThisChannel = 0;
    }
  }

  if(fCrossCorr == nullptr){
    fCrossCorr = new CrossCorrelator();
    fSpawnedCrossCorrelator = true;
  }
  
  fDtCache.init(fCrossCorr, this, true);
  geom->usePhotogrammetryNumbers(0);

}






Acclaim::InterferometricMap* Acclaim::AnalysisReco::getMap(AnitaPol::AnitaPol_t pol){
  InterferometricMap* h = coarseMaps[pol];
  coarseMaps[pol] = nullptr;
  return h;
}






Acclaim::InterferometricMap* Acclaim::AnalysisReco::getZoomMap(AnitaPol::AnitaPol_t pol, UInt_t peakInd){

  InterferometricMap* h = fineMaps[pol][peakInd];
  if(h==nullptr){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << " for run "
	      << fCurrentRun << ", eventNumber " << fCurrentEventNumber
	      << ", unable to find fineMap with pol " << pol
              << " for peakInd = " << peakInd << " about to return NULL." << std::endl;
  }
  fineMaps[pol][peakInd] = nullptr;

  return h;
}







void Acclaim::AnalysisReco::reconstruct(AnitaPol::AnitaPol_t pol, const Adu5Pat* pat) {

  if(coarseMaps[pol]==nullptr)
  {
    // std::cerr << "new coarse map " << pol << std::endl;
    // I don't own the map or there isn't one, so I'll make a new one
    coarseMaps[pol] = new InterferometricMap();
  }  
  else{// I own the map so I can overwrite it to avoid allocating memory
    // std::cerr << "old coarse map " << pol << std::endl;
    coarseMaps[pol]->Reset();
  }

  if(pat != nullptr){
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

  auto* h = new InterferometricMap(peakIndex, phiSector, zoomCenterPhiDeg, PHI_RANGE_ZOOM, zoomCenterThetaDeg, THETA_RANGE_ZOOM);
  if(pat != nullptr){
    h->addGpsInfo(pat);
  }
  h->Fill(pol, fCrossCorr, &fDtCache);  

  // std::cout << h->GetName() << std::endl;
  // std::cout << "in " << __PRETTY_FUNCTION__ << ": " << fineMaps[pol][peakIndex] << std::endl;
  if(fineMaps[pol][peakIndex] != nullptr){
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
  if(static_cast<int>(lindasNums.is_open())==0){
    return 1; // This is an error
  }

  Int_t ant;
  Double_t dr, dPhiRad, dz, dt;

  while(lindasNums >> ant >> dr >> dz >> dPhiRad >> dt){

    Int_t surf, chan, ant2;
    AnitaGeomTool::getSurfChanAntFromRingPhiPol(AnitaRing::AnitaRing_t (ant/NUM_PHI), ant%NUM_PHI, pol,
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
  dts.reserve(theAnts.size());

  Double_t phiWave = peakPhiDeg*TMath::DegToRad();
  Double_t thetaWave = peakThetaDeg*TMath::DegToRad();
  for(int ant : theAnts){
    // Double_t dt = getDeltaTExpected(pol, theAnts.at(0), ant, phiWave, thetaWave);
    Double_t dt = getDeltaTExpected(pol, ant, phiWave, thetaWave);
    dts.push_back(dt);
  }
}

AnalysisWaveform* Acclaim::AnalysisReco::coherentlySum(const FilteredAnitaEvent* fEv, AnitaPol::AnitaPol_t pol,
                                                       const std::vector<Int_t>& theAnts,
                                                       Double_t peakPhiDeg, Double_t peakThetaDeg,
						       Double_t* biggestPeakToPeak, Double_t* forceT0) {

  Int_t biggest = -1;
  Double_t largestPeakToPeak = 0;
  std::vector<const AnalysisWaveform*> waves;
  waves.reserve(theAnts.size());

  for(int ant : theAnts){
    const AnalysisWaveform* wf = fEv->getFilteredGraph(ant, pol);
    waves.push_back(wf);
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

  if(biggestPeakToPeak != nullptr){
    (*biggestPeakToPeak) = largestPeakToPeak;
  }

  if(biggest < 0 || biggest >= NUM_SEAVEYS){
    std::cerr << "Error in " << __PRETTY_FUNCTION__
	      << ", I couldn't find a waveform where vMax - vMin > 0. "
              << "Something's wrong!" << std::endl;
  }

  std::vector<double> dts;
  directionAndAntennasToDeltaTs(theAnts, pol, peakPhiDeg, peakThetaDeg, dts);
  
  return coherentlySum(waves, dts, forceT0);
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

AnalysisWaveform* Acclaim::AnalysisReco::coherentlySum(std::vector<const AnalysisWaveform*>& waves, std::vector<Double_t>& dts, const Double_t* forceT0) {

  if(checkWavesAndDtsMatch(waves, dts)==0){
    std::cerr << "Nothing to sum.. about to return NULL" << std::endl;
    return nullptr;
  }

  const double extraTimeNs = 20; //2*(*std::max_element(dts.begin(), dts.end()));
  const TGraph* gr0 = waves[0]->even();
  const double totalTime = gr0->GetX()[gr0->GetN()-1] - gr0->GetX()[0] + extraTimeNs;
  const int interpN = floor(totalTime/fCoherentDtNs);
  
  std::vector<double> zeros(interpN, 0);
  double t0 = forceT0 != nullptr ? *forceT0 : gr0->GetX()[0] - 0.5*extraTimeNs;
  auto* coherentWave = new AnalysisWaveform(interpN, &zeros[0], fCoherentDtNs, t0);
  TGraphAligned* grCoherent = coherentWave->updateEven();

  for(UInt_t i=0; i < waves.size(); i++){
    const TGraphAligned* grU = waves[i]->even();
    // std::cout << i << " here" << std::endl;
    
    const double t0 = grU->GetX()[0];
    const double tN = grU->GetX()[grU->GetN()-1];

    if(fDebug != 0){
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
    auto* gr = new TGraphAligned(grOld->GetN(), grOld->GetX(), grOld->GetY());

    for(int i=0; i < gr->GetN(); i++){
      gr->GetX()[i] -= dts[w];
    }
    
    grs.push_back(gr);
  }
}


inline TString nameLargestPeak(Int_t peakInd){
  TString title;
  switch (peakInd){
  case 0:  title = "Largest Peak ";                                  break;
  case 1:  title = "2nd Largest Peak ";                              break;
  case 2:  title = "3rd Largest Peak ";                              break;
  default: title = TString::Format("%dth Largest Peak",  peakInd+1); break;
  }
  return title;
}


void Acclaim::AnalysisReco::drawSummary(TPad* wholePad, AnitaPol::AnitaPol_t pol){

  const int numColsForNow = 3;
  EColor peakColors[AnitaPol::kNotAPol][numColsForNow] = {{kBlack, EColor(kMagenta+2), EColor(kViolet+2)},
                                                          {kBlue,  EColor(kSpring+4),  EColor(kPink + 10)}};

  // will go from top to bottom, like this...
  // 0-1: vpol/hpol Reconstruction
  // 1-2: Coarse Map
  // 2-3: Fine map layers title
  // 3-4: Fine maps (divided into fDrawNPeaks internally)
  // 4-5: Detailed data tables
  std::vector<double> wholePadVLayers = {1, 0.95, 0.75, 0.72, 0.3,  0};
  Int_t wholePadLayer = 0;

  if(fDrawNPeaks==1){ // tweak layout as there's too much space for just 1 fine map
    wholePadVLayers = {1, 0.95, 0.65, 0.62, 0.35,  0};
  }

  if(wholePad==nullptr){
    UInt_t eventNumber = fCrossCorr->eventNumber[pol];
    TString canName = TString::Format("can%u", eventNumber);
    TString polSuffix = pol == AnitaPol::kHorizontal ? "HPol" : "VPol";
    canName += polSuffix;
    TString canTitle = TString::Format("Event %u - ", eventNumber) + polSuffix;
    wholePad = new TCanvas(canName);
  }
  wholePad->Clear();
  
  TPad* wholeTitlePad = RootTools::makeSubPad(wholePad,
					      0, wholePadVLayers.at(wholePadLayer+1),
					      1, wholePadVLayers.at(wholePadLayer),
					      TString::Format("%d_title", (int)pol));
  wholePadLayer++;

  (void) wholeTitlePad;
  auto *wholeTitle = new TPaveText(0, 0, 1, 1);
  TString wholeTitleText; // = TString::Format(");
  wholeTitleText += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";
  wholeTitleText += " Reconstruction";
  wholeTitle->AddText(wholeTitleText);
  wholeTitle->SetBit(kCanDelete, true);
  wholeTitle->SetLineWidth(0);
  wholeTitle->Draw();
  
  TPad* coarseMapPad = RootTools::makeSubPad(wholePad,
					     0, wholePadVLayers.at(wholePadLayer+1),
					     1, wholePadVLayers.at(wholePadLayer),
					     TString::Format("%d_coarse", (int)pol));
  wholePadLayer++;

  InterferometricMap* hCoarse = coarseMaps[pol];

  const double fracFinePeak = 0.2; // fraction of pad width for the fine peak, the rest is split between coherent/dedispsered
  const EColor xPolColor = kRed;
  const double xPolTitleBoxWidth = 0.1;

  TPad* tempPad = RootTools::makeSubPad(wholePad,
					0, wholePadVLayers.at(wholePadLayer+1),
					1, wholePadVLayers.at(wholePadLayer),
					"tempTitlePad");
  wholePadLayer++;
  tempPad->cd();

  std::vector<TPaveText*> paves;
  paves.push_back(new TPaveText(0, 0, fracFinePeak, 1));
  paves.back()->AddText("Fine Map");
  paves.back()->SetTextAlign(kVAlignCenter + kHAlignCenter);

  paves.push_back(new TPaveText(fracFinePeak, 0, fracFinePeak+0.5*(1-fracFinePeak)-0.1, 1));
  paves.back()->AddText("Coherent");
  paves.back()->SetTextAlign(kVAlignCenter + kHAlignCenter);
  if(fDrawXPol != 0){
    double x2 = paves.back()->GetX2();
    double x1 = x2 - xPolTitleBoxWidth;
    paves.back()->SetX2(x1);
    paves.back()->SetTextAlign(kVAlignCenter + kHAlignRight);

    paves.push_back(new TPaveText(x1, 0, x2, 1));
    paves.back()->AddText("(Cross-Pol)");
    paves.back()->SetTextAlign(kVAlignCenter + kHAlignLeft);
    ((TText*)paves.back()->GetListOfLines()->Last())->SetTextColor(xPolColor);
  }

  paves.push_back(new TPaveText(fracFinePeak+0.5*(1-fracFinePeak), 0, 1, 1));
  paves.back()->AddText("Dedispersed");
  paves.back()->SetTextAlign(kVAlignCenter + kHAlignCenter);

  if(fDrawXPolDedispersed != 0){
    double x2 = paves.back()->GetX2();
    double x1 = x2 - xPolTitleBoxWidth;
    paves.back()->SetX2(x1);
    paves.back()->SetTextAlign(kVAlignCenter + kHAlignRight);

    paves.push_back(new TPaveText(x1, 0, x2, 1));
    paves.back()->AddText("(Cross-Pol)");
    ((TText*)paves.back()->GetListOfLines()->Last())->SetTextColor(xPolColor);
    paves.back()->SetTextAlign(kVAlignCenter + kHAlignLeft);
  }

  for(auto & pave : paves){
    pave->SetBit(kCanDelete);
    pave->SetShadowColor(0);
    pave->SetLineWidth(0);
    pave->SetLineColor(0);
    pave->Draw();
  }

  TPad* finePeaksAndCoherent = RootTools::makeSubPad(wholePad,
						     0, wholePadVLayers.at(wholePadLayer+1),
						     1, wholePadVLayers.at(wholePadLayer),
						     "peaks");
  wholePadLayer++;

  std::list<InterferometricMap*> drawnFineMaps;

  const int nFine = TMath::Min(fNumPeaks, fDrawNPeaks);
  for(int peakInd = 0; peakInd < nFine; peakInd++){
    double yUp = 1 - double(peakInd)/nFine;
    double yLow = yUp - double(1)/nFine;
    TPad* finePeak = RootTools::makeSubPad(finePeaksAndCoherent, 0, yLow, fracFinePeak, yUp, "fine");
    (void) finePeak;
    
    InterferometricMap* h = fineMaps[pol][peakInd];
    if(h != nullptr){
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

    double x0 = fCoherentFiltered[pol][peakInd][0]->even()->GetX()[0] - 10;
    double xN = fCoherentFiltered[pol][peakInd][0]->even()->GetX()[fCoherentFiltered[pol][peakInd][0]->even()->GetN()-1] + 10;

    for(int j=0; j <= 1; j++){ // j loops over a coherent/dedispersed
      std::vector<AnalysisWaveform*> wavesToDraw;
      std::vector<const AnitaEventSummary::WaveformInfo*> waveInfo;
      std::vector<EColor> waveColors;
      std::vector<Style_t> lineStyles;
      std::vector<TString> legText;
      if(j==0) {
	if(fDrawCoherent > 0){
	  wavesToDraw.push_back(fCoherentFiltered[pol][peakInd][0]);
	  waveInfo.push_back(&fSummary.coherent_filtered[pol][peakInd]);
	  waveColors.push_back(peakColors[pol][peakInd]);
	  lineStyles.push_back(1);
	  legText.emplace_back("Co-pol");
	}
	if(fDrawXPol > 0){
	  wavesToDraw.push_back(fCoherentFiltered[pol][peakInd][1]);
	  waveInfo.push_back(nullptr);
	  // waveColors.push_back(peakColors[pol][peakInd]);
	  waveColors.push_back(xPolColor);
	  lineStyles.push_back(7);
	  legText.emplace_back("Cross-pol");
	}
      }
      else{
	if(fDrawDedispersed > 0){
	  wavesToDraw.push_back(fDeconvolvedFiltered[pol][peakInd][0]);
	  waveInfo.push_back(&fSummary.deconvolved_filtered[pol][peakInd]);
	  waveColors.push_back(peakColors[pol][peakInd]);
	  lineStyles.push_back(1);
	  legText.emplace_back("Co-pol");
	}
	if(fDrawXPolDedispersed > 0){
	  wavesToDraw.push_back(fDeconvolvedFiltered[pol][peakInd][1]);
	  waveInfo.push_back(nullptr);
	  waveColors.push_back(xPolColor); //peakColors[pol][peakInd]);
	  lineStyles.push_back(7);
	  legText.emplace_back("Cross-pol");
	}
      }

      // if(fDebug) std::cerr << "Draw flags:\t" << fDrawCoherent << "\t" << fDrawXPol << "\t" << fDrawDedispersed << "\t" << fDrawXPolDedispersed << std::endl;
      // if(fDebug) std::cerr << "This gives us " << wavesToDraw.size() << " waveforms" << std::endl;
      TGraphInteractive* gr = nullptr;
    
      for(int i=0; i < wavesToDraw.size(); i++){
	if(fDrawDomain==kTimeDomain){
	  const TGraphAligned* gr2 = wavesToDraw[i]->even();
      
	  auto* gr3 = new TGraphInteractive(gr2, "l");
	  gr3->SetFillColor(0);
	  gr3->SetLineColor(waveColors[i]);
	  gr3->SetMarkerColor(waveColors[i]);
	  gr3->SetLineStyle(lineStyles[i]);

	  gr3->GetXaxis()->SetTitleSize(1);
	  gr3->GetYaxis()->SetTitleSize(1);

	  gr3->SetBit(kCanDelete); // Let ROOT track and handle deletion

	  if(gr != nullptr){
	    gr->addGuiChild(gr3);
	  }
	  else {
	    gr = gr3;

	    TString title = nameLargestPeak(peakInd);
	    title += (j == 0 ? "Coherent Waveform" : "Dedispersed Waveform");
	    title += (pol == AnitaPol::kHorizontal ? " HPol" : " VPol");
	    title += ";Time (ns);Amplitude (mV)";

	    gr->SetTitle(title);
	    gr->GetXaxis()->SetRangeUser(x0, xN);
	  }
	}
	else{ // freq domain
	  const TGraphAligned* grPower2 = wavesToDraw[i]->powerdB();

	  auto* grPower = new TGraphInteractive(grPower2, "l");

	  grPower->SetLineColor(waveColors[i]);
	  grPower->SetMarkerColor(waveColors[i]);
	  grPower->SetLineStyle(lineStyles[i]);

	  grPower->SetBit(kCanDelete); // Let ROOT track and handle deletion

	  if(gr == nullptr){
	    gr = grPower;
	    const double maxFreqGHz = 1.3;
	    gr->GetXaxis()->SetRangeUser(0, maxFreqGHz);
	    TString title = nameLargestPeak(peakInd);
	    title += (j == 0 ? "Coherent Waveform" : "Dedispersed Waveform");
	    title += (pol == AnitaPol::kHorizontal ? " HPol" : " VPol");
	    title += ";Frequency (GHz);Power Spectral Density (dBm/MHz)";
	    grPower->SetTitle(title);
	  }
	  else{
	    gr->addGuiChild(grPower);
	  }

	  if((waveInfo[i] != nullptr) && (fFillSpectrumInfo != 0) && (gr != nullptr)){
	    double y0 = waveInfo[i]->spectrumSlope*fSlopeFitStartFreqGHz + waveInfo[i]->spectrumIntercept;
	    auto* grSlope = new TGraphInteractive(1, &fSlopeFitStartFreqGHz, &y0, "l");
	    grSlope->SetPoint(1, fSlopeFitEndFreqGHz, waveInfo[i]->spectrumSlope*fSlopeFitEndFreqGHz + waveInfo[i]->spectrumIntercept);
	    grSlope->SetLineColor(grPower->GetLineColor());
	    gr->addGuiChild(grSlope); // pass it a pointer and it takes ownership

	    double df = wavesToDraw[i]->deltaF();
	    auto* grCoherentPeakMarker = new TGraphInteractive(0, nullptr, nullptr, "p");
	    grCoherentPeakMarker->SetMarkerStyle(4);
	    grCoherentPeakMarker->SetMarkerSize(0.75);
	    grCoherentPeakMarker->SetMarkerColor(grPower->GetLineColor());

	    for(double k : waveInfo[i]->peakFrequency){
	      if(k > 0){
		int powerBin = TMath::Nint(k/df);
		grCoherentPeakMarker->SetPoint(grCoherentPeakMarker->GetN(), grPower->GetX()[powerBin], grPower->GetY()[powerBin]);
	      }
	    }
	    gr->addGuiChild(grCoherentPeakMarker);
	  }
	}
      }
      if(j==0){
	TPad* coherentPad = RootTools::makeSubPad(finePeaksAndCoherent, fracFinePeak, yLow, fracFinePeak + 0.5*(1 - fracFinePeak), yUp, "coherent");
	coherentPad->cd();
      }
      else{
	TPad* dedispersedPad = RootTools::makeSubPad(finePeaksAndCoherent, fracFinePeak + 0.5*(1 - fracFinePeak), yLow, 1, yUp, "dedispersed");
	dedispersedPad->cd();
      }
      if(gr != nullptr){
	gr->DrawGroup("al");
	// if(peakInd==0){
	//   gPad->Update(); // necessary to get root to paint it into existence
	//   TPaveText* titlePave = (TPaveText*)gPad->FindObject("title");
	//   if(titlePave){
	//     std::cout << titlePave->GetSize() << std::endl;
	//   }
	// }
      }
    }
  }

  coarseMapPad->cd();
  hCoarse->DrawGroup("colz");  

  if(!drawnFineMaps.empty()){
    auto it = drawnFineMaps.begin();
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
    if(hCoarse != nullptr){
      hCoarse->SetMaximum(polMax);
      // hCoarse->SetMinimum(polMin);      
    }
    while(!drawnFineMaps.empty()){
      drawnFineMaps.pop_back();
    }
    // std::cout << polMax << "\t" << polMin << std::endl;
  }

  const double epsilon = 0.001; // so there's enough room for line border on the edges of the pave text
  TPad* textPad = RootTools::makeSubPad(wholePad,
					0, wholePadVLayers.at(wholePadLayer+1),
					1, wholePadVLayers.at(wholePadLayer),
					"text");
  wholePadLayer++;
  (void) textPad;
  for(int peakInd=0; peakInd < nFine; peakInd++){
    double xlow = double(peakInd)/nFine;
    double xup = xlow + double(1.)/nFine - epsilon;
    
    auto *title = new TPaveText(xlow, 0.9, xup, 1);
    // title->SetBorderSize(0);
    title->SetBit(kCanDelete, true);
    title->SetTextColor(peakColors[pol][peakInd]);
    title->SetLineColor(peakColors[pol][peakInd]);

    // TString titleText = pol == AnitaPol::kHorizontal ? "HPol" : "VPol";
    TString titleText;
    titleText += TString::Format(" Peak %d", peakInd);    
    title->AddText(titleText);
    
    title->Draw();
    
    auto *pt = new TPaveText(xlow, epsilon, xup, 0.9);
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
  fCoherentFiltered[pol][peakInd][xPol] = nullptr;
  return wf;  
}


AnalysisWaveform* Acclaim::AnalysisReco::getCoherent(AnitaPol::AnitaPol_t pol, Int_t peakInd, bool xPol){
  AnalysisWaveform* wf = fCoherent[pol][peakInd][xPol];
  fCoherent[pol][peakInd][xPol] = nullptr;
  return wf;  
}



AnalysisWaveform* Acclaim::AnalysisReco::getDeconvolved(AnitaPol::AnitaPol_t pol, Int_t peakInd, bool xPol){
  AnalysisWaveform* wf = fDeconvolved[pol][peakInd][xPol];
  fDeconvolved[pol][peakInd][xPol] = nullptr;
  return wf;
}



AnalysisWaveform* Acclaim::AnalysisReco::getDeconvolvedFiltered(AnitaPol::AnitaPol_t pol, Int_t peakInd, bool xPol){
  AnalysisWaveform* wf = fDeconvolvedFiltered[pol][peakInd][xPol];
  fDeconvolvedFiltered[pol][peakInd][xPol] = nullptr;
  return wf;
}




FilteredAnitaEvent* Acclaim::AnalysisReco::getEvMin(){
  FilteredAnitaEvent* f = fEvMin;
  fEvMin = nullptr;
  return f;
}



FilteredAnitaEvent* Acclaim::AnalysisReco::getEvMinDeco(){
  FilteredAnitaEvent* f = fEvMinDeco;
  fEvMinDeco = nullptr;
  return f;
}



FilteredAnitaEvent* Acclaim::AnalysisReco::getEvDeco(){
  FilteredAnitaEvent* f = fEvDeco;
  fEvDeco = nullptr;
  return f;
}
