#include "AnalysisReco.h"
#include "InterferometricMap.h"
#include "InterferometryCache.h"


Acclaim::AnalysisReco::AnalysisReco(){
  initializeInternals();
}

Acclaim::AnalysisReco::~AnalysisReco(){

  if(spawnedCrossCorrelator && cc){
    delete cc;
  }
  
}





void Acclaim::AnalysisReco::process(const FilteredAnitaEvent * usefulEvent, UsefulAdu5Pat* usefulPat ,AnitaEventSummary * eventSummary) const{
  

  if(!cc){
    cc = new CrossCorrelator();
    spawnedCrossCorrelator = true;
    // std::cout << "here" << "\t" << spawnedCrossCorrelator << std::endl;
    dtCache.init(cc, this);
  }

  
  const int thisNumPeaks = 3;

  eventSummary->eventNumber = usefulEvent->eventNumber;

  for(Int_t polInd=0; polInd < AnitaPol::kNotAPol; polInd++){

    
    AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;

    cc->correlateEvent(usefulEvent, pol);

    // do the coarsely grained reconstruction
    reconstruct(pol);

    std::vector<Double_t> coarseMapPeakValues;
    std::vector<Double_t> coarseMapPeakPhiDegs;
    std::vector<Double_t> coarseMapPeakThetaDegs;    
    coarseMaps[pol]->findPeakValues(thisNumPeaks, coarseMapPeakValues, coarseMapPeakPhiDegs, coarseMapPeakThetaDegs);

    eventSummary->nPeaks[pol] = thisNumPeaks;
    
    
    for(Int_t peakInd=0; peakInd < thisNumPeaks; peakInd++){
      reconstructZoom(pol, peakInd, coarseMapPeakPhiDegs.at(peakInd), coarseMapPeakThetaDegs.at(peakInd));

      std::map<Int_t, InterferometricMap*>::iterator it = fineMaps[pol].find(peakInd);
      InterferometricMap* h = NULL;
      if(it!=fineMaps[pol].end()){
	h = it->second;
      }
      else{
	std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unable to find finely binned histogram for peak " << peakInd << std::endl;
      }

      h->getPeakInfo(eventSummary->peak[pol][peakInd].value, eventSummary->peak[pol][peakInd].phi, eventSummary->peak[pol][peakInd].theta);      

      // fill in difference between rough and fine
      eventSummary->peak[pol][peakInd].dphi_rough = eventSummary->peak[pol][peakInd].phi - coarseMapPeakPhiDegs.at(peakInd);
      eventSummary->peak[pol][peakInd].dtheta_rough = eventSummary->peak[pol][peakInd].theta - coarseMapPeakThetaDegs.at(peakInd);

      // based on Cosmin's comments in AnitaAnalysisSummary.h
      eventSummary->peak[pol][peakInd].phi_separation = peakInd == 0 ? 1000 : RootTools::getDeltaAngleDeg(eventSummary->peak[pol][peakInd].phi, eventSummary->peak[pol][0].phi);
      
      // AnalysisWaveform* coherent = coherentlySum(usefulEvent, h);


      AnalysisWaveform* coherentWave = coherentlySum(usefulEvent, h);
      
      std::map<Int_t, AnalysisWaveform*>::iterator it2 = coherent[pol].find(peakInd);
      if(it2!=coherent[pol].end()){
	if(it2->second != NULL){
	  delete it2->second;
	}
      }
      coherent[pol][peakInd] = coherentWave;

      const TGraphAligned* grHilbert = coherentWave->hilbertEnvelope();
      eventSummary->coherent[pol][peakInd].peakHilbert = TMath::MaxElement(grHilbert->GetN(), grHilbert->GetY());

      if(usefulPat != NULL){
      
	Double_t phiWave = TMath::DegToRad()*eventSummary->peak[pol][peakInd].phi;
	Double_t thetaWave = -1*TMath::DegToRad()*eventSummary->peak[pol][peakInd].theta;

	// *   Returns 0 if never hits the ground, even with maximum adjustment
	// *   Returns 1 if hits the ground with no adjustment
	// *   Returns 2 if it hits the ground with adjustment      
	int success = usefulPat->traceBackToContinent(phiWave, thetaWave, 
						      &eventSummary->peak[pol][peakInd].latitude,
						      &eventSummary->peak[pol][peakInd].longitude,
						      &eventSummary->peak[pol][peakInd].altitude,
						      &eventSummary->peak[pol][peakInd].theta_adjustment_needed);

	if(success==0){
	  eventSummary->peak[pol][peakInd].latitude = -9999;
	  eventSummary->peak[pol][peakInd].longitude = -9999;
	  eventSummary->peak[pol][peakInd].altitude = -9999;
	  eventSummary->peak[pol][peakInd].distanceToSource = -9999;
	}
	else{
	  eventSummary->peak[pol][peakInd].distanceToSource = SPEED_OF_LIGHT_NS*usefulPat->getTriggerTimeNsFromSource(eventSummary->peak[pol][peakInd].latitude,
														      eventSummary->peak[pol][peakInd].longitude,
														      eventSummary->peak[pol][peakInd].altitude);	
	}
      }

      summary = (*eventSummary);
      
    }
  }
}



void Acclaim::AnalysisReco::initializeInternals(){

  cc = NULL;
  spawnedCrossCorrelator = false;
  
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

  kUseOffAxisDelay = 1;  
  coherentDeltaPhi = 1; // +/- this many phi-sectors when coherently summing waves
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

Double_t Acclaim::AnalysisReco::relativeOffAxisDelay(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2,
					       Double_t phiDeg) const {

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


Double_t Acclaim::AnalysisReco::getDeltaTExpected(AnitaPol::AnitaPol_t pol, Int_t ant1, Int_t ant2, Double_t phiWave, Double_t thetaWave) const{

  // Double_t tanThetaW = tan(thetaWave);
  Double_t tanThetaW = tan(-1*thetaWave);
  Double_t part1 = zArray[pol].at(ant1)*tanThetaW - rArray[pol].at(ant1)*cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol].at(ant1));
  Double_t part2 = zArray[pol].at(ant2)*tanThetaW - rArray[pol].at(ant2)*cos(phiWave-TMath::DegToRad()*phiArrayDeg[pol].at(ant2));
  Double_t tdiff = 1e9*((cos(thetaWave) * (part2 - part1))/SPEED_OF_LIGHT); // Returns time in ns

  return tdiff;
}



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

  if(!cc){
    cc = new CrossCorrelator();
    spawnedCrossCorrelator = true;
  }
  dtCache.init(cc, this, true);
  geom->usePhotogrammetryNumbers(0);

}







Acclaim::InterferometricMap* Acclaim::AnalysisReco::getMap(AnitaPol::AnitaPol_t pol){
  InterferometricMap* h = coarseMaps[pol];
  coarseMaps[pol] = NULL;
  return h;
}





Acclaim::InterferometricMap* Acclaim::AnalysisReco::getZoomMap(AnitaPol::AnitaPol_t pol, UInt_t peakInd){

  InterferometricMap* h = NULL;

  std::map<Int_t, InterferometricMap*>::iterator it = fineMaps[pol].find((int)peakInd);
  if(it!=fineMaps[pol].end()){
    // could still be null though!
    h = it->second;
    it->second = NULL;
  }

  if(h==NULL){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unable to find fineMap with pol " << pol 
	      << " for peakInd = " << peakInd << " about to return NULL." << std::endl;
  }

  return h;
}








void Acclaim::AnalysisReco::reconstruct(AnitaPol::AnitaPol_t pol) const{

  if(coarseMaps[pol]==NULL)
  {
    // std::cerr << "new coarse map " << pol << std::endl;
    // I don't own the map or there isn't one, so I'll make a new one
    coarseMaps[pol] = new InterferometricMap();
  }  
  else{// I own the map so I can overwrite it to avoid allocating memory
    // std::cerr << "old coarse map " << pol << std::endl;

    for(int phiBin=1; phiBin<=coarseMaps[pol]->GetNbinsPhi(); phiBin++){
      for(int thetaBin=1; thetaBin<=coarseMaps[pol]->GetNbinsTheta(); thetaBin++){
	coarseMaps[pol]->SetBinContent(phiBin, thetaBin, 0);	
      }
    }
  }

  coarseMaps[pol]->Fill(pol, cc, &dtCache);
}


void Acclaim::AnalysisReco::reconstructZoom(AnitaPol::AnitaPol_t pol, Int_t peakIndex, Double_t zoomCenterPhiDeg, Double_t zoomCenterThetaDeg) const{

  // Some kind of sanity check here due to the unterminating while loop inside RootTools::getDeltaAngleDeg
  if(zoomCenterPhiDeg < -500 || zoomCenterThetaDeg < -500 ||
     zoomCenterPhiDeg >= 500 || zoomCenterThetaDeg >= 500){

    std::cerr << "Warning in "<< __PRETTY_FUNCTION__ << " in " << __FILE__
	      << ". zoomCenterPhiDeg = " << zoomCenterPhiDeg
	      << " and zoomCenterThetaDeg = " << zoomCenterThetaDeg
	      << " ...these values look suspicious so I'm skipping this reconstruction."
	      << " eventNumber = " << cc->eventNumber[pol] << std::endl;
    return;
  }

  Double_t deltaPhiDegPhi0 = RootTools::getDeltaAngleDeg(zoomCenterPhiDeg, InterferometricMap::getBin0PhiDeg());
  deltaPhiDegPhi0 = deltaPhiDegPhi0 < 0 ? deltaPhiDegPhi0 + DEGREES_IN_CIRCLE : deltaPhiDegPhi0;

  Int_t phiSector = floor(deltaPhiDegPhi0/PHI_RANGE);

  InterferometricMap* h = new InterferometricMap(peakIndex, phiSector, zoomCenterPhiDeg, PHI_RANGE_ZOOM, zoomCenterThetaDeg, THETA_RANGE_ZOOM);
  h->Fill(pol, cc, &dtCache);  

  // std::cout << h->GetName() << std::endl;
  
  std::map<Int_t, InterferometricMap*>::iterator it = fineMaps[pol].find(peakIndex);
  if(it!=fineMaps[pol].end() && it->second != NULL){
    // std::cerr << "trying to delete... " << it->first << "\t" << it->second << std::endl;
    delete it->second;
    fineMaps[pol].erase (it);
    // delete it->second;
    // it->second = NULL;
  }

  // for(it=fineMaps[pol].begin(); it!=fineMaps[pol].end(); ++it){
  //   std::cout << "the fineMaps[ " << pol << "] contains " << it->first << "\t" << it->second << std::endl;
  // }

  fineMaps[pol][peakIndex] = h;
  
  // for(it=fineMaps[pol].begin(); it!=fineMaps[pol].end(); ++it){
  //   std::cout << "the fineMaps[ " << pol << "] contains " << it->first << "\t" << it->second << std::endl;
  // }


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


AnalysisWaveform* Acclaim::AnalysisReco::coherentlySum(const FilteredAnitaEvent* fEv, const InterferometricMap* h) const{

  Int_t peakPhiSector = h->getPeakPhiSector();
  AnitaPol::AnitaPol_t pol = h->getPol();

  Int_t biggest = -1;
  Double_t largestPeakToPeak = 0;
  std::vector<Int_t> ants;
  for(int deltaPhiSect=-coherentDeltaPhi; deltaPhiSect <= coherentDeltaPhi; deltaPhiSect++){

    Int_t phiSector = peakPhiSector + deltaPhiSect;
    phiSector = phiSector < 0        ? phiSector + NUM_PHI : phiSector;
    phiSector = phiSector >= NUM_PHI ? phiSector - NUM_PHI : phiSector;
    
    for(int ring = 0; ring < AnitaRing::kNotARing; ring++){
      int ant = AnitaGeomTool::getAntFromPhiRing(phiSector, AnitaRing::AnitaRing_t(ring));
      ants.push_back(ant);
      
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
  }

  if(biggest < 0 || biggest >= NUM_SEAVEYS){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", I couldn't find a waveform where vMax - vMin > 0. "
	      << "Something's wrong, and I'm probably about to vomit a stack trace all over your terminal..." << std::endl;
  }

  // now we've found the channel with the biggest peak-to-peak
  // coherently sum the waves with this channel first
  std::vector<const AnalysisWaveform*> waves;
  waves.push_back(fEv->getFilteredGraph(biggest, pol));
  for(unsigned i=0; i < ants.size(); i++){
    if(ants[i] != biggest){
      waves.push_back(fEv->getFilteredGraph(ants[i], pol));
      // std::cout << waves.size() << "\t" << ants[i] << std::endl;
    }
  }

  Double_t ip, phiDeg, thetaDeg;
  h->getPeakInfo(ip, phiDeg, thetaDeg);
  Double_t phiWave = phiDeg*TMath::DegToRad();
  Double_t thetaWave = thetaDeg*TMath::DegToRad();  
    
  std::vector<Double_t> dts(1, 0); // 0 offset for biggest antenna

  for(unsigned i=0; i < ants.size(); i++){
    if(ants[i]!=biggest){
    
      Double_t dt = getDeltaTExpected(pol, biggest, ants[i], phiWave, thetaWave);      
      dts.push_back(dt);
      // std::cout << dts.size() << "\t" << ants[i] << std::endl;      
    }    
  }
					   
  // std::cerr << "return " << __PRETTY_FUNCTION__ << std::endl;
  return coherentlySum(waves, dts);
}



AnalysisWaveform* Acclaim::AnalysisReco::coherentlySum(std::vector<const AnalysisWaveform*>& waves, std::vector<Double_t>& dts) const{
  
  if(waves.size() < 1){    
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", nothing to sum.. about to return NULL" << std::endl;
    return NULL;
  }  
  else if(waves.size() != dts.size()){    
    const char* action = dts.size() < waves.size() ? "padding" : "trimming";
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unequal vectors (waves.size() = " << waves.size() << ", dts.size() = " << dts.size() << ")\n"
	      << action << " dts..." << std::endl;
    while(waves.size() != dts.size()){
      if(dts.size() < waves.size()){
	dts.push_back(0);
      }
      else{
	dts.pop_back();	
      }
    }
  }

  // std::cerr << "here 1\t"  << dts.size() << "\t" << waves.size() << std::endl;
  AnalysisWaveform* coherentWave = new AnalysisWaveform((*waves[0]));
  // std::cerr << "here 2\t"  << coherentWave << std::endl;  
  
  TGraphAligned* grCoherent = coherentWave->updateEven();
  for(UInt_t i=1; i < waves.size(); i++){
    for(int samp=0; samp < grCoherent->GetN(); samp++){
      double t = grCoherent->GetX()[samp];
      grCoherent->GetY()[samp] += waves[i]->evalEven(t + dts[i]);
    };    
  }

  // std::cerr << "here 3\t"  << grCoherent << std::endl;    

  for(int samp=0; samp < grCoherent->GetN(); samp++){
    grCoherent->GetY()[samp]/=waves.size();
  };    

  // std::cerr << "return " << __PRETTY_FUNCTION__ << "\t" << coherentWave << std::endl;  
  return coherentWave;
}












/** 
 * Helper function for drawSummary
 * draws a new TPad inside the parent TPad (cds into it too) with 
 * @param parentPad the pad to draw inside of
 * @param xlow is relative to parent pad
 * @param ylow is relative to parent pad
 * @param xup is relative to parent pad
 * @param yup is relative to parent pad
 * @param suffix is the suffix to append to the pad name
 * 
 * @return the new subPad
 */
TPad* Acclaim::AnalysisReco::makeSubPad(TPad* parentPad, double xlow, double ylow, double xup, double yup, TString suffix){

  parentPad->cd(); // go into parent pad (assume it's drawn, would that matter?)
  TString subPadName = TString::Format("%s_%s", parentPad->GetName(), suffix.Data());
  TPad* subPad = new TPad(subPadName, subPadName, xlow, ylow, xup, yup);
  subPad->Draw();
  subPad->cd();
  return subPad;
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
    UInt_t eventNumber = cc->eventNumber[pol];
    TString canName = TString::Format("can%u", eventNumber);
    TString polSuffix = pol == AnitaPol::kHorizontal ? "HPol" : "VPol";
    canName += polSuffix;
    TString canTitle = TString::Format("Event %u - ", eventNumber) + polSuffix;
    wholePad = new TCanvas(canName);
  }
  wholePad->Clear();
  

  TPad* wholeTitlePad = makeSubPad(wholePad, 0, 0.95, 1, 1, TString::Format("%d_title", (int)pol));
  TPaveText *wholeTitle = new TPaveText(0, 0, 1, 1);  
  TString wholeTitleText; // = TString::Format(");
  wholeTitleText += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";
  wholeTitleText += " Reconstruction";
  wholeTitle->AddText(wholeTitleText);
  wholeTitle->SetBit(kCanDelete, true);
  wholeTitle->SetLineWidth(0);
  wholeTitle->Draw();
  
  TPad* coarseMapPad = makeSubPad(wholePad, 0, 0.75, 1, 0.95, TString::Format("%d_coarse", (int)pol));
  InterferometricMap* hCoarse = coarseMaps[pol];
  if(hCoarse){
    hCoarse->Draw("colz");
  }

  // std::vector<TGraph> grPeaks;

  TPad* finePeaksAndCoherent = makeSubPad(wholePad, 0, 0.35, 1, 0.75, "peaks");

  std::list<InterferometricMap*> drawnFineMaps;
  std::vector<TGraphAligned*> drawnCoherent;
  Double_t coherentMax = -1e9, coherentMin = 1e9;
  
  const int nFine = fineMaps[pol].size();
  for(int peakInd = 0; peakInd < nFine; peakInd++){
    double yUp = 1 - double(peakInd)/nFine;
    double yLow = yUp - double(1)/nFine;    
    TPad* finePeak = makeSubPad(finePeaksAndCoherent, 0, yLow, 0.2, yUp, "fine");
    
    std::map<Int_t, InterferometricMap*>::iterator it = fineMaps[pol].find(peakInd);
    if(it!=fineMaps[pol].end()){
      InterferometricMap* h = it->second;
      if(h){
	h->SetTitleSize(1);
	h->GetXaxis()->SetTitleSize(0.01);
	h->GetYaxis()->SetTitleSize(0.01);        

	h->Draw("col");
	drawnFineMaps.push_back(h);

	
	TGraph& gr = h->getPeakPointGraph();
	gr.Draw("psame");

	TGraph& gr2 = h->getEdgeBoxGraph();

	gr2.Draw("lsame");
	gr.SetMarkerColor(peakColors[pol][peakInd]);
	gr.SetMarkerStyle(8); // dot
	gr2.SetLineColor(peakColors[pol][peakInd]);
	gr2.SetLineWidth(4);
	

	
	
	coarseMapPad->cd();
	TGraph* gr3 = (TGraph*) gr2.Clone();
	gr3->SetBit(kCanDelete, true);// leave to ROOT garbage collector
	
	gr3->SetLineWidth(2);
	// gr.Draw("psame");
	gr3->Draw("lsame");	
      }
    }

    TPad* coherentPad    = makeSubPad(finePeaksAndCoherent, 0.2, yLow, 0.6, yUp, "coherent");
    TPad* coherentFftPad = makeSubPad(finePeaksAndCoherent, 0.6, yLow,   1, yUp, "coherentFFT");

    std::map<Int_t, AnalysisWaveform*>::iterator it2 = coherent[pol].find(peakInd);
    if(it2!=coherent[pol].end()){
      AnalysisWaveform* coherentWave = it2->second;
      if(coherentWave){
	
	const char* opt = "al";

	// don't want to be able to edit it by accident so copy it...
	const TGraphAligned* gr = coherentWave->even();
	TGraphAligned* gr2 = (TGraphAligned*) gr->Clone();

	gr2->SetFillColor(0);
	gr2->SetBit(kCanDelete, true); // let ROOT's garbage collector worry about it
	// std::cout << gr2 << "\t" << pol << "\t" << peakInd << "\t" << gr2->IsOnHeap() << std::endl;
	TString title = TString::Format("Coherently summed wave - peak %d", peakInd);
	gr2->SetTitle(title);
	gr2->SetLineColor(peakColors[pol][peakInd]);

	coherentPad->cd();	
	gr2->Draw(opt);

	drawnCoherent.push_back(gr2);
	coherentMax = gr2->GetMaximum() > coherentMax ? gr2->GetMaximum() : coherentMax;
	coherentMin = gr2->GetMinimum() < coherentMin ? gr2->GetMinimum() : coherentMin;	
	
	const TGraphAligned* grPower = coherentWave->powerdB();
	TGraphAligned* grPower2 = (TGraphAligned*) grPower->Clone();
	grPower2->SetBit(kCanDelete, true); // let ROOT's garbage collector worry about it
	title = TString::Format("PSD coherently summed wave - peak %d", peakInd);
	grPower2->SetTitle(title);
	grPower2->SetLineColor(peakColors[pol][peakInd]);

	coherentFftPad->cd();
	grPower2->Draw(opt);

	// std::cout << "here \t" << summaryGraphs.size() << std::endl;
      }
      else{
	std::cerr << "missing coherent pointer?\t" << pol << "\t" << peakInd << std::endl;
      }
    }
    else{
      std::cerr << "missing coherent in map?\t" << pol << "\t" << peakInd << std::endl;
    }

  }

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

  if(drawnCoherent.size() > 0){
    for(unsigned i=0; i < drawnCoherent.size(); i++){
      drawnCoherent.at(i)->SetMaximum(coherentMax);
      drawnCoherent.at(i)->SetMinimum(coherentMin);      
    }
  }


  TPad* textPad = makeSubPad(wholePad, 0, 0, 1, 0.35, "text");
  
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
