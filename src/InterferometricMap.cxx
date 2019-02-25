#include "InterferometricMap.h"
#include "AnalysisReco.h" // for the geometric definitions
#include "InterferometryCache.h"

#include "TAxis.h"
#include "TMath.h"
#include "TMatrixD.h"

#include <iostream>
#include "TString.h"
#include "UsefulAdu5Pat.h"
#include "RampdemReader.h"
#include "TruthAnitaEvent.h"
#include "TGraphAntarctica.h"
#include "TArrowAntarctica.h"
#include "TLegend.h"
#include "AcclaimClustering.h"

ClassImp(Acclaim::InterferometricMap);

std::vector<Double_t> coarseBinEdgesPhi; // has size NUM_BINS_PHI+1
std::vector<Double_t> fineBinEdgesPhi; // has size NUM_BINS_PHI_ZOOM_TOTAL+1

std::vector<Double_t> coarseBinEdgesTheta; // has size NUM_BINS_THETA+1
std::vector<Double_t> fineBinEdgesTheta; // has size NUM_BINS_THETA_ZOOM_TOTAL+1
Double_t bin0PhiDeg = -9999;

TDecompSVD theSVD;
bool doneInitSVD = false;


std::vector<Double_t> phiCenterCenterDegs;


/** 
 * Utility function for quadratic peak interpolation.
 * 
 * Get the grid definition for the SVD method once and only once.
 * 
 * @param dPhiBinLow is the number of bins below the peak bin in the phi direction
 * @param dPhiBinHigh is the number of bins above the peak bin in the phi direction
 * @param dThetaBinLow is the number of bins below the peak bin in the theta direction
 * @param dThetaBinHigh is the number of bins below the peak bin in the theta direction
 */

void Acclaim::InterferometricMap::getDeltaBinRangeSVD(Int_t& dPhiBinLow, Int_t& dPhiBinHigh, Int_t& dThetaBinLow, Int_t& dThetaBinHigh){
  dPhiBinLow = -NUM_BINS_QUAD_FIT_PHI/2;
  dPhiBinHigh = dPhiBinLow + NUM_BINS_QUAD_FIT_PHI;
  dThetaBinLow = -NUM_BINS_QUAD_FIT_THETA/2;
  dThetaBinHigh = dThetaBinLow + NUM_BINS_QUAD_FIT_THETA;  
}





/** 
 * @brief Lazily prepares the SVD solver for the quadratic peak interpolation.
 * 
 * We are using the SVD method to find the best bit of the peak to a 2D quadratic.
 * There are 6 coefficients for a 2D quadratic which I imaginitvely call a,b,c,d,e,f
 * F(x,y) = ax^{2} + by^{2} + cxy + dx + ey + f
 * We have a bunch of values for F(x, y) which is the zoomed map values.
 * Currently the size of the grid is hard coded in the header with the NUM_BINS_QUAD_FIT values.
 * 
 * @return the populated SVD object
 */
const TDecompSVD& Acclaim::InterferometricMap::getSVD(){

  if(!doneInitSVD){

    // 2D quadratic defined by ax^{2} + by^{2} + cxy + dx + ey + f, so 6 coefficients
    const int nCoeffs = 6;

    // get grid offsets
    int firstBinPhi, lastBinPhi, firstBinTheta, lastBinTheta;    
    getDeltaBinRangeSVD(firstBinPhi, lastBinPhi, firstBinTheta, lastBinTheta);

    // one row per data point, one column per coefficient    
    TMatrixD A(NUM_BINS_QUAD_FIT_PHI*NUM_BINS_QUAD_FIT_THETA, nCoeffs);

    int row=0;
    for(int i=firstBinPhi; i < lastBinPhi; i++){
      const double dPhi = ZOOM_BIN_SIZE_PHI * i;
      for(int j=firstBinTheta; j < lastBinTheta; j++){
        const double dTheta = ZOOM_BIN_SIZE_THETA * j;
        A(row, 0) = dPhi*dPhi;
        A(row, 1) = dTheta*dTheta;
        A(row, 2) = dPhi*dTheta;
        A(row, 3) = dPhi;
        A(row, 4) = dTheta;
        A(row, 5) = 1;
        row++;
      }
    }
    theSVD.SetMatrix(A);
    doneInitSVD = true;
  }
  return theSVD;
}





/** 
 * Draw the interferometric map only (use DrawGroup to also draw TGraphInteractives)
 * 
 * @param opt 
 */
void Acclaim::InterferometricMap::Draw(Option_t* opt){

  // trigger the creation of the contained graphs upon drawing
  getSunGraph();
  getTruthGraph();
  if(fIsZoomMap){
    getEdgeBoxGraph();
  }
  TH2D::Draw(opt);
}


/** 
 * Method required by ROOT interative behaviour.
 * Draws a copy of the map when double clicked.
 * 
 * @param event is the interactive event
 * @param x is the mouse position in x
 * @param y is the mouse position in y
 */
void Acclaim::InterferometricMap::ExecuteEvent(int event, int x, int y){
 
  if(event == kButton1Double){
    (void) x;
    (void) y;
    new TCanvas();
    DrawGroup("colz");
  }
  TH2D::ExecuteEvent(event, x, y);
}




/** 
 * Get the center of the phi-sectors relative to ADU5 aft-fore.
 * For the purposes of the interferometric map they are exactly PHI_RANGE apart
 * 
 * @param phi is the phi-sector 0-15
 * 
 * @return the angle in degrees
 */
Double_t Acclaim::InterferometricMap::getPhiSectorCenterPhiDeg(int phi){
  if(phiCenterCenterDegs.size()==0){
    AnitaGeomTool* geom = AnitaGeomTool::Instance();
    Double_t aftForeOffset = geom->aftForeOffsetAngleVertical*TMath::RadToDeg();
    
    Double_t phi0 = -aftForeOffset;
    if(phi0 < -DEGREES_IN_CIRCLE/2){
      phi0+=DEGREES_IN_CIRCLE;
    }
    else if(phi0 >= DEGREES_IN_CIRCLE/2){
      phi0-=DEGREES_IN_CIRCLE;
    }
    
    for(int phi=0; phi < NUM_PHI; phi++){
      phiCenterCenterDegs.push_back(phi0 + PHI_RANGE*phi);
    }
  }
  
  return phiCenterCenterDegs.at(phi);
}



/** 
 * Get the position of the first bin in phi
 * 
 * @return the position of the first bin in phi in degrees
 */
Double_t Acclaim::InterferometricMap::getBin0PhiDeg(){

  if(bin0PhiDeg == -9999){
    double phi0 = getPhiSectorCenterPhiDeg(0);
    bin0PhiDeg = phi0 - PHI_RANGE/2;
  }
  return bin0PhiDeg;
}




/** 
 * Lazily populate the bin edges in theta for the coarse map.
 * Currently the bin sizes are scales by sin(theta) where theta=0 is the horizontal
 * @return the coase map bin edges in theta.
 */
const std::vector<Double_t>& Acclaim::InterferometricMap::getCoarseBinEdgesTheta(){

  if(coarseBinEdgesTheta.size()==0) // then not initialized so do it here...
  {
    
    // funk up the theta bin spacing...  
    UInt_t nBinsTheta = NUM_BINS_THETA;
    Double_t minTheta = MIN_THETA;
    Double_t maxTheta = MAX_THETA;
  
    // calculate the bin spaces in sin(theta)
    Double_t sinThetaMin = sin(minTheta*TMath::DegToRad());
    Double_t sinThetaMax = sin(maxTheta*TMath::DegToRad());
    Double_t sinThetaRange = sinThetaMax - sinThetaMin;
    Double_t dSinTheta = sinThetaRange/nBinsTheta;

    // std::vector<Double_t> binEdges(nBinsTheta+1);
    coarseBinEdgesTheta.reserve(nBinsTheta+1);
    for(unsigned bt = 0; bt <= nBinsTheta; bt++){
      Double_t thisSinTheta = sinThetaMin + bt*dSinTheta;
      Double_t thisTheta = TMath::RadToDeg()*TMath::ASin(thisSinTheta);
      // coarseBinEdgesTheta.at(bt) = thisTheta;
      coarseBinEdgesTheta.push_back(thisTheta);
    }
  }
  return coarseBinEdgesTheta;
}


/** 
 * Lazily populate the bin edges in theta for the fine map.
 * Currently the bin sizes are just integer units of theta in degrees.
 * @return the fine map bin edges in theta.
 */
const std::vector<Double_t>& Acclaim::InterferometricMap::getFineBinEdgesTheta(){

  if(fineBinEdgesTheta.size()==0) // then not initialized so do it here...
  {

    UInt_t nBinsTheta = NUM_BINS_THETA_ZOOM_TOTAL;
    Double_t minTheta = MIN_THETA - THETA_RANGE_ZOOM/2;
    Double_t maxTheta = MAX_THETA + THETA_RANGE_ZOOM/2;
  
    Double_t dTheta = (maxTheta - minTheta) / nBinsTheta;
    fineBinEdgesTheta.reserve(nBinsTheta+1);
    for(unsigned bt = 0; bt <= nBinsTheta; bt++){
      double thisTheta = minTheta + bt*dTheta;
      fineBinEdgesTheta.push_back(thisTheta);
    }
  }
  return fineBinEdgesTheta;
}





/** 
 * Lazily populate the bin edges in phi for the coarse map.
 * 
 * @return the coarse map bin edges in phi.
 */
const std::vector<Double_t>& Acclaim::InterferometricMap::getCoarseBinEdgesPhi(){

  if(coarseBinEdgesPhi.size()==0) // then not initialized so do it here...
  {
    UInt_t nBinsPhi = NUM_BINS_PHI*NUM_PHI;
    Double_t minPhi = getBin0PhiDeg();
    Double_t dPhi = double(DEGREES_IN_CIRCLE)/nBinsPhi;
    // std::vector<Double_t> binEdges(nBinsTheta+1);
    coarseBinEdgesTheta.reserve(nBinsPhi+1);
    for(unsigned bp = 0; bp <= nBinsPhi; bp++){
      Double_t thisPhi = minPhi + dPhi*bp;
      coarseBinEdgesPhi.push_back(thisPhi);
    }
  }
  return coarseBinEdgesPhi;
}



/** 
 * Lazily populate the bin edges in phi for the fine map.
 * 
 * @return the fine map bin edges in phi.
 */
const std::vector<Double_t>& Acclaim::InterferometricMap::getFineBinEdgesPhi(){

  if(fineBinEdgesPhi.size()==0) // then not initialized so do it here...
  {

    // funk up the theta bin spacing...  
    UInt_t nBinsPhi = NUM_BINS_PHI_ZOOM_TOTAL;
    Double_t minPhi = getBin0PhiDeg() - PHI_RANGE_ZOOM/2;
    // need some extra degrees here to account for the fact that you might get a coarse peak
    // on the edge of a map, therefore the fine peak might extend off the left/right edge 
    Double_t dPhi = double(DEGREES_IN_CIRCLE + PHI_RANGE_ZOOM)/nBinsPhi;
    // std::vector<Double_t> binEdges(nBinsTheta+1);
    fineBinEdgesTheta.reserve(nBinsPhi+1);
    for(unsigned bp = 0; bp <= nBinsPhi; bp++){
      Double_t thisPhi = minPhi + dPhi*bp;
      fineBinEdgesPhi.push_back(thisPhi);
    }
  }
  return fineBinEdgesPhi;
}




/** 
 * Sets the name of the histogram from an incrementing counter.
 * This avoids warnings in ROOT about things histograms having the same name.
 */
void Acclaim::InterferometricMap::setDefaultName(){

  static unsigned defaultCounter = 0;
  fName = TString::Format("hDefault%u", defaultCounter);
  defaultCounter++;
}



/** 
 * Sets an informative name and title using the member variables.
 */
void Acclaim::InterferometricMap::setNameAndTitle(){


  if(fIsZoomMap){
    fName = "hFine";
    fName += pol == AnitaPol::kVertical ? "ImageV" : "ImageH";
    fName += TString::Format("%u_%u", peakIndex, eventNumber);

    fTitle = TString::Format("Event %u ", eventNumber);
    fTitle += (pol == AnitaPol::kVertical ? "VPOL" : "HPOL");
    fTitle += TString::Format(" Zoom Map - Peak %u", peakIndex);
  }
  else{
    fName = "hCoarse";
    fName += pol == AnitaPol::kVertical ? "V" : "H";
    fName += TString::Format("%u", eventNumber);

    fTitle = TString::Format("Event %u ", eventNumber);
    fTitle += (pol == AnitaPol::kVertical ? "VPOL" : "HPOL");
  }
}













/** 
 * Fine map constructor
 * 
 * @param peakInd is the index of the peak in the coarse map (0 is the largest, 1 the second largest, etc.)
 * @param phiSector is the closest phi-sector to the centered zoom peak
 * @param zoomCentrePhi is the position on which to centre the map in phi (degrees)
 * @param phiRange is the size of the zoom map in phi (degrees)
 * @param zoomCentreTheta is the position on which to centre the map in theta (degrees)
 * @param thetaRange is the size of the zoom map in theta (degrees)
 */
Acclaim::InterferometricMap::InterferometricMap(Int_t peakInd, Int_t phiSector, Double_t zoomCentrePhi, Double_t phiRange, Double_t zoomCentreTheta, Double_t thetaRange)
    : fUsefulPat(NULL), truthLat(-9999), truthLon(-9999), truthAlt(-9999)
{
  fIsZoomMap = true;
  fPeakPhiSector = phiSector;
  peakIndex = peakInd;
  Double_t minPhiDesired = zoomCentrePhi - phiRange/2;
  Double_t maxPhiDesired = zoomCentrePhi + phiRange/2;
  const std::vector<double> fineBinsPhi = Acclaim::InterferometricMap::getFineBinEdgesPhi();
  // int minPhiBin, maxPhiBin;
  minPhiBin = -1;
  int maxPhiBin = -1;  
  getIndicesOfEdgeBins(fineBinsPhi, minPhiDesired, maxPhiDesired, minPhiBin, maxPhiBin);

  // std::cerr << "ctr 0" << "\t" << minPhiBin << "\t" << maxPhiBin << "\t" << fineBinsPhi.size() << std::endl;
  
  Double_t minThetaDesired = zoomCentreTheta - thetaRange/2;
  Double_t maxThetaDesired = zoomCentreTheta + thetaRange/2;
  const std::vector<double> fineBinsTheta = Acclaim::InterferometricMap::getFineBinEdgesTheta();
  // int minThetaBin, maxThetaBin;
  minThetaBin = -1;
  int maxThetaBin = -1;
  getIndicesOfEdgeBins(fineBinsTheta, minThetaDesired, maxThetaDesired, minThetaBin, maxThetaBin);

  // std::cerr << "ctr 1" << "\t" << minThetaBin << "\t" << maxThetaBin << "\t" << fineBinsTheta.size() << std::endl;
  
  SetBins(maxPhiBin - minPhiBin, &fineBinsPhi[minPhiBin], maxThetaBin - minThetaBin, &fineBinsTheta[minThetaBin]); // changes also errors array (if any)  
  initializeInternals();

  // std::cerr << "ctr 2" << std::endl;  
}



/** 
 * Default constructor, which creates a coarse map
 */
Acclaim::InterferometricMap::InterferometricMap()
    : fUsefulPat(NULL), truthLat(-9999), truthLon(-9999), truthAlt(-9999)
{
  fIsZoomMap = false;
  fPeakPhiSector = -1;
  minPhiBin = -1;
  minThetaBin = -1;
  peakIndex = -1;
  fSigmaTheta = -1;
  fSigmaPhi = -1;
  const std::vector<double> coarseBinsPhi = Acclaim::InterferometricMap::getCoarseBinEdgesPhi();
  const std::vector<double> coarseBinsTheta = Acclaim::InterferometricMap::getCoarseBinEdgesTheta();
  SetBins(coarseBinsPhi.size()-1, &coarseBinsPhi[0], coarseBinsTheta.size()-1, &coarseBinsTheta[0]);

  initializeInternals();  
}



/** 
 * Default destructor
 */
Acclaim::InterferometricMap::~InterferometricMap(){
  deletePat();
}


void Acclaim::InterferometricMap::Reset(Option_t* opt){
  TH2D::Reset(opt);
  deleteChildren();
}


/** 
 * If we have a copy of the usefulAdu5Pat, then delete it nicely.
 */
void Acclaim::InterferometricMap::deletePat(){
  if(fUsefulPat) {
    // std::cerr << fUsefulPat->heading << "\t" << std::endl;
    delete fUsefulPat;
  }
  fUsefulPat = NULL;
}



/** 
 * Create an internal UsefulAdu5Pat from the raw gps data
 * 
 * @param pat is a pointer to the Adu5Pat
 */
void Acclaim::InterferometricMap::addGpsInfo(const Adu5Pat* pat){
  // std::cerr << __PRETTY_FUNCTION__ << "\t" << this << "\t" << pat << std::endl;
  deletePat();

  fUsefulPat = new UsefulAdu5Pat(pat);
}



/** 
 * Create an internal UsefulAdu5Pat from the raw gps data
 * 
 * @param pat is a pointer to the Adu5Pat
 */
void Acclaim::InterferometricMap::addTruthInfo(const TruthAnitaEvent* truth){

  if(truth){
    truthLon = truth->sourceLon;
    truthLat = truth->sourceLat;
    truthAlt = truth->sourceAlt;
  }
}







/** 
 * Workhorse function to turn cross-correlations into an interferometric map.
 * 
 * @param thePol is the polarisation of the map to make
 * @param cc points to the cross-correlator which contains the cross-correlations
 * @param dtCache is a pointer to the object containing the cached set of deltaTs.
 */
void Acclaim::InterferometricMap::Fill(AnitaPol::AnitaPol_t thePol, CrossCorrelator* cc, InterferometryCache* dtCache){

  pol = thePol;
  eventNumber = cc->eventNumber[pol];
  setNameAndTitle();
  
  if(!fIsZoomMap){
    
    std::vector<Int_t>* combosToUse = NULL;
    Int_t binsPerPhiSector = GetNbinsPhi()/NUM_PHI;
    for(Int_t phiSector = 0; phiSector < NUM_PHI; phiSector++){
      combosToUse = &cc->combosToUseGlobal[phiSector];

      Int_t startPhiBin = phiSector*binsPerPhiSector;
      Int_t endPhiBin = startPhiBin + binsPerPhiSector;

      for(UInt_t comboInd=0; comboInd<combosToUse->size(); comboInd++){
	Int_t combo = combosToUse->at(comboInd);
	if(cc->kOnlyThisCombo >= 0 && combo!=cc->kOnlyThisCombo){
	  continue;
	}
	for(Int_t phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){
	  for(Int_t thetaBin = 0; thetaBin < GetNbinsTheta(); thetaBin++){

	    Int_t bin = (thetaBin+1)*(GetNbinsPhi()+2) + phiBin+1;
	    Double_t cInterp = cc->getCrossCorrelation(pol, combo, dtCache->coarseDt(pol, combo, phiBin, thetaBin));	  
	    AddBinContent(bin,cInterp);
	    fEntries++; // otherwise drawing don't work at all
	  }
	}
      }
    }

    for(Int_t phiSector = 0; phiSector < NUM_PHI; phiSector++){
      combosToUse = &cc->combosToUseGlobal[phiSector];

      Double_t normFactor = cc->kOnlyThisCombo < 0 && combosToUse->size() > 0 ? combosToUse->size() : 1;
      // absorb the removed inverse FFT normalization
      normFactor*=(cc->numSamples*cc->numSamples);

      Int_t startPhiBin = phiSector*binsPerPhiSector;
      Int_t endPhiBin = (phiSector+1)*binsPerPhiSector;
      for(Int_t phiBin = startPhiBin; phiBin < endPhiBin; phiBin++){
	for(Int_t thetaBin = 0; thetaBin < GetNbinsTheta(); thetaBin++){

	  Double_t val = GetBinContent(phiBin+1, thetaBin+1);
	  SetBinContent(phiBin+1, thetaBin+1, val/normFactor);

	  if(val > fPeakBinValue){	
            fPeakBinValue = val;
            // peakPhiDeg = GetPhiAxis()->GetBinLowEdge(phiBin+1);
            // peakThetaDeg = GetThetaAxis()->GetBinLowEdge(thetaBin+1);	  
            fPeakBinPhi = phiBin;
            fPeakBinTheta = thetaBin;
	    fPeakPhiSector = phiSector;
	  }
	}
      }
    }
    setPeakInfoJustFromPeakBin(fPeakBinPhi, fPeakBinTheta);
  }

  else{
    // std::cerr << "zmps " << fPeakPhiSector << std::endl;
    cc->doUpsampledCrossCorrelations(pol, fPeakPhiSector);    
    
    std::vector<Int_t>& combosToUse = cc->combosToUseGlobal[fPeakPhiSector];

    const Int_t offset = cc->numSamplesUpsampled/2;
    for(UInt_t comboInd=0; comboInd<combosToUse.size(); comboInd++){
      Int_t combo = combosToUse.at(comboInd);
      if(cc->kOnlyThisCombo >= 0 && combo!=cc->kOnlyThisCombo){
	continue;
      }
      int ant1 = cc->comboToAnt1s[combo];
      int ant2 = cc->comboToAnt2s[combo];
      // for(Int_t thetaBin = 0; thetaBin < NUM_BINS_THETA_ZOOM; thetaBin++){
      for(Int_t thetaBin = 0; thetaBin < GetNbinsTheta(); thetaBin++){	
	Int_t zoomThetaInd = minThetaBin + thetaBin;
	// Double_t zoomThetaWave = zoomedThetaWaves[zoomThetaInd];
	// Double_t partBA = partBAsZoom[pol][combo][zoomThetaInd];
	Double_t partBA = dtCache->fPartBAsZoom[dtCache->partBAsIndex(pol, combo, zoomThetaInd)]; //)[pol][combo][zoomThetaInd];      
	// Double_t dtFactor = dtFactors[zoomThetaInd];
	Double_t dtFactor = dtCache->fDtFactors[zoomThetaInd];

	// for(Int_t phiBin = 0; phiBin < NUM_BINS_PHI_ZOOM; phiBin++){
	for(Int_t phiBin = 0; phiBin < GetNbinsPhi(); phiBin++){
	  Int_t zoomPhiInd = minPhiBin + phiBin;
	  // Double_t zoomPhiWave = zoomedPhiWaveLookup[zoomPhiInd];

	  int p21 = dtCache->part21sIndex(pol, combo, zoomPhiInd);
	  Double_t offsetLowDouble = dtFactor*(partBA - dtCache->fPart21sZoom[p21]);//[pol][combo][zoomPhiInd]);		

	  offsetLowDouble += dtCache->fUseOffAxisDelay > 0 ? dtCache->fOffAxisDelaysDivided[p21] : 0;	  
	
	  offsetLowDouble += cc->startTimes[pol][ant1]/cc->correlationDeltaT;
	  offsetLowDouble -= cc->startTimes[pol][ant2]/cc->correlationDeltaT;
	

	  // hack for floor()
	  Int_t offsetLow = (int) offsetLowDouble - (offsetLowDouble < (int) offsetLowDouble);

	  Double_t deltaT = (offsetLowDouble - offsetLow);
	  offsetLow += offset;
	  Double_t c1 = cc->crossCorrelationsUpsampled[pol][combo][offsetLow];
	  Double_t c2 = cc->crossCorrelationsUpsampled[pol][combo][offsetLow+1];
	  Double_t cInterp = deltaT*(c2 - c1) + c1;

          Int_t bin = (thetaBin+1)*(GetNbinsPhi()+2) + phiBin+1;
	  AddBinContent(bin, cInterp);
	  fEntries++;
	}
      }
    }

    Double_t normFactor = cc->kOnlyThisCombo < 0 && combosToUse.size() > 0 ? combosToUse.size() : 1;
    // absorb the removed inverse FFT normalization
    normFactor*=(cc->numSamples*cc->numSamples);

    // set peak finding variables
    fPeakBinValue = -DBL_MAX;
  
    for(Int_t thetaBin = 0; thetaBin < GetNbinsTheta(); thetaBin++){
      for(Int_t phiBin = 0; phiBin < GetNbinsPhi(); phiBin++){
	Double_t val = GetBinContent(phiBin+1, thetaBin+1);
	val /= normFactor;
	SetBinContent(phiBin+1, thetaBin+1, val);
	if(val > fPeakBinValue){	
	  fPeakBinValue = val;
	  // peakPhiDeg = GetPhiAxis()->GetBinLowEdge(phiBin+1);
	  // peakThetaDeg = GetThetaAxis()->GetBinLowEdge(thetaBin+1);	  
	  fPeakBinPhi = phiBin;
	  fPeakBinTheta = thetaBin;
	}
      }
    }
    fitPeakWithQuadratic(fPeakBinPhi, fPeakBinTheta);
    
  }
  
}



/** 
 * Useful if we're not fitting, or something goes wrong with the fit procedure.
 * 
 * @param peakPhiBin is the phi bin on which to centre the 2D quadratic
 * @param peakThetaBin is the theta bin on which to centre the 2D quadratic
 */
void Acclaim::InterferometricMap::setPeakInfoJustFromPeakBin(Int_t peakPhiBin, Int_t peakThetaBin){
  fPeakPhi = GetPhiAxis()->GetBinLowEdge(peakPhiBin+1);
  fPeakTheta = GetThetaAxis()->GetBinLowEdge(peakThetaBin+1);
  fPeakValue = GetBinContent(peakPhiBin+1, peakThetaBin+1);
}


/** 
 * @brief Fit a 2D quadratic function to the bins around the peak of the map.
 * Only implemented for zoomed maps.
 * 
 * @param peakPhiBin is the phi bin on which to centre the 2D quadratic
 * @param peakThetaBin is the theta bin on which to centre the 2D quadratic
 */

void Acclaim::InterferometricMap::fitPeakWithQuadratic(Int_t peakPhiBin, Int_t peakThetaBin){

  if(!fIsZoomMap){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", not implemented for coarse map." << std::endl;    
    return;
  }
  
  // Get range and check we're not too close to the edge...  
  int minDeltaPhiBin, maxDeltaPhiBin, minDeltaThetaBin, maxDeltaThetaBin;
  getDeltaBinRangeSVD(minDeltaPhiBin, maxDeltaPhiBin, minDeltaThetaBin, maxDeltaThetaBin);

  int firstPhiBin = minDeltaPhiBin + peakPhiBin;
  int lastPhiBin = maxDeltaPhiBin + peakPhiBin;
  int firstThetaBin = minDeltaThetaBin + peakThetaBin;
  int lastThetaBin = maxDeltaThetaBin + peakThetaBin;  

  if(firstPhiBin < 0 || lastPhiBin >= GetNbinsPhi() || firstThetaBin < 0 || lastThetaBin >= GetNbinsTheta()){
    // don't warn as this happens a lot and gets very verbose...
    // std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", too close to to the edge to do quadratic peak fit" << std::endl;
    setPeakInfoJustFromPeakBin(peakPhiBin, peakThetaBin);    
    return;    
  }

  // put bin values around peak into a vector for the matrix equation we're about to solve
  TVectorD peakData(NUM_BINS_QUAD_FIT_PHI*NUM_BINS_QUAD_FIT_THETA);
  int row=0;
  for(int phiBin=firstPhiBin; phiBin < lastPhiBin; phiBin++){
    for(int thetaBin=firstThetaBin; thetaBin < lastThetaBin; thetaBin++){
      peakData[row] = GetBinContent(phiBin+1, thetaBin+1);
      row++;
    }
  }
    
  // solve for the quadratic coefficients a,b,c,d,e,f
  TDecompSVD& svd = const_cast<TDecompSVD&>(getSVD()); // cast away const-ness
  Bool_t ok = false;
  TVectorD quadraticCoefficients = svd.Solve(peakData, ok);
  if(!ok){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__  << ": something went wrong with the peak interpolation. Will just use peak bin." << std::endl;
    setPeakInfoJustFromPeakBin(peakPhiBin, peakThetaBin);
    return;
  }

  // pick 2D quadratic coefficients out from the results vector
  double a = quadraticCoefficients(0);
  double b = quadraticCoefficients(1);
  double c = quadraticCoefficients(2);
  double d = quadraticCoefficients(3);
  double e = quadraticCoefficients(4);
  double f = quadraticCoefficients(5);


  // do some algebra to determine the centre position and peak value
  const double denom = (4*a*b - c*c);
  const double x_c = (-2*b*d + c*e)/denom;
  const double y_c = (-2*a*e + c*d)/denom;
  const double p_c = a*x_c*x_c + b*y_c*y_c + c*x_c*y_c + d*x_c + e*y_c + f;

  // Put fit results into member variables containing the results
  fPeakPhi = x_c + GetPhiAxis()->GetBinLowEdge(peakPhiBin+1);
  fPeakTheta = y_c + GetThetaAxis()->GetBinLowEdge(peakThetaBin+1);
  fPeakValue = p_c;

  // now also get uncertainties, correlation, chisquare
  TVectorD residual = svd.GetMatrix()*quadraticCoefficients - peakData;
  fPeakReducedChisquare = residual.Norm2Sqr()/(residual.GetNrows());

  // @todo, add uncertainties in phi/theta and covariance
  // would have been better to do a gaussian fit
  // by taking the logarithm beforehand...
  fPeakSigmaPhi = 0;
  fPeakSigmaTheta = 0;
  fPeakCovariance = 0;
}







void Acclaim::InterferometricMap::findPeakValues(Int_t numPeaks, std::vector<Double_t>& peakValues, std::vector<Double_t>& phiDegs, std::vector<Double_t>& thetaDegs) const{

  
  // In this function I want to find numPeak peaks and set an exclusion zone around each peak
  // so that the next peak isn't just a close neighbour of a true peak.

  // Set not crazy, but still debug visible values for peak values/location
  // -DBL_MAX was causing me some while loop issues in RootTools::getDeltaAngleDeg(...)
  // which has an unrestricted while loop inside.

  peakValues.clear();
  phiDegs.clear();
  thetaDegs.clear();
  
  peakValues.resize(numPeaks, -999);
  phiDegs.resize(numPeaks, -999);
  thetaDegs.resize(numPeaks, -999);

  
  // As we start, all regions are allowed.
  // Int_t allowedBins[NUM_BINS_PHI*NUM_PHI][NUM_BINS_THETA];
  // Int_t allowedBins[NUM_BINS_PHI*NUM_PHI][NUM_BINS_THETA];
  int nTheta = GetNbinsTheta();
  int nPhi = GetNbinsPhi();  
  std::vector<short> allowedBins(nPhi*nTheta, 1);

  
  for(Int_t peakInd=0; peakInd < numPeaks; peakInd++){
    Int_t gotHere = 0;
    for(Int_t phiBin=0; phiBin<nPhi; phiBin++){
      for(Int_t thetaBin = 0; thetaBin < nTheta; thetaBin++){

        // if(allowedBins[phiBin][thetaBin] > 0){
        if(allowedBins[phiBin*nTheta + thetaBin] > 0){
          // then we are in an allowed region

          // Double_t mapValue = coarseMap[pol][phiBin][thetaBin];
          Double_t mapValue = GetBinContent(phiBin+1, thetaBin+1);
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
            // if(coarseMap[pol][lastPhiBin][thetaBin] < mapValue &&
            //    coarseMap[pol][nextPhiBin][thetaBin] < mapValue){
            if(GetBinContent(lastPhiBin+1, thetaBin+1) < mapValue &&
               GetBinContent(nextPhiBin+1, thetaBin+1) < mapValue){

              gotHere = 2;
              // So far so good, now check theta neigbours
              // This doesn't wrap, so I will allow edge cases to pass as maxima
              Int_t nextThetaBin = thetaBin != (NUM_BINS_THETA) - 1 ? thetaBin+1 : thetaBin;
              Int_t lastThetaBin = thetaBin != 0 ? thetaBin-1 : thetaBin;

              // Check theta bins below
              gotHere = 3;
              if(thetaBin == 0 || (GetBinContent(lastPhiBin+1, lastThetaBin+1) < mapValue &&
                                   GetBinContent(phiBin+1, lastThetaBin+1) < mapValue &&
                                   GetBinContent(nextPhiBin+1, lastThetaBin+1) < mapValue)){

                // Check theta bins above
                gotHere = 4;
                if(thetaBin == NUM_BINS_THETA-1 || (GetBinContent(lastPhiBin+1, nextThetaBin+1) < mapValue &&
                                                    GetBinContent(phiBin+1, nextThetaBin+1) < mapValue &&
                                                    GetBinContent(nextPhiBin+1, nextThetaBin+1) < mapValue)){

                  // Okay okay okay... you're a local maxima, you can be my next peak.
                  peakValues[peakInd] = GetBinContent(phiBin+1, thetaBin+1);
                  // phiDegs[peakInd] = phiWaveLookup[phiBin]*TMath::RadToDeg();
                  // thetaDegs[peakInd] = thetaWaves[thetaBin]*TMath::RadToDeg();
                  phiDegs[peakInd] = GetPhiAxis()->GetBinLowEdge(phiBin+1);
                  thetaDegs[peakInd] = GetThetaAxis()->GetBinLowEdge(thetaBin+1);
                  // thetaDegs[peakInd] = Acclaim::InterferometricMap::getCoarseBinEdgesTheta()[thetaBin];
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

      for(Int_t phiBin=0; phiBin<nPhi; phiBin++){
        // Double_t phiDeg = phiWaveLookup[phiBin]*TMath::RadToDeg();
        Double_t phiDeg = GetPhiAxis()->GetBinLowEdge(phiBin+1);
        Double_t absDeltaPhi = TMath::Abs(RootTools::getDeltaAngleDeg(phiDegs[peakInd], phiDeg));
        if(absDeltaPhi < PEAK_PHI_DEG_RANGE){

          for(Int_t thetaBin = 0; thetaBin < nTheta; thetaBin++){

            // Double_t thetaDeg = thetaWaves[thetaBin]*TMath::RadToDeg();
            Double_t thetaDeg = GetThetaAxis()->GetBinLowEdge(thetaBin+1);
            Double_t absDeltaTheta = TMath::Abs(RootTools::getDeltaAngleDeg(thetaDegs[peakInd], thetaDeg));
            if(absDeltaTheta < PEAK_THETA_DEG_RANGE){

              // Now this region in phi/theta is disallowed on the next loop through...
              // allowedBins[phiBin][thetaBin] = 0;
              allowedBins[phiBin*nTheta + thetaBin] = 0;
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







void Acclaim::InterferometricMap::initializeInternals(){
  fThetaAxisInSinTheta = true;
  pol = AnitaPol::kNotAPol;
  eventNumber = 0;
  fPeakPhi = -9999;
  fPeakTheta = -9999;
  fPeakValue = -1;
  fSigmaTheta = -1;
  fSigmaPhi = -1;
  setDefaultName();
}



void Acclaim::InterferometricMap::getIndicesOfEdgeBins(const std::vector<double>& binEdges, Double_t lowVal, Double_t highVal, Int_t& lowIndex, Int_t& highIndex){


  // double minDiff = 1e9;
  
  // get last below
  lowIndex = 0;
  for(unsigned i=0; i < binEdges.size(); i++){
    // if(fabs(binEdges.at(i) - lowVal) < minDiff){
    //   minDiff = fabs(binEdges.at(i) - lowVal);
    //   lowIndex = i;
    // }
    
    if(binEdges.at(i) < lowVal){
      lowIndex = i;
    }
    else{
      break;
    }
  }
  
  // if(binEdges.at(0) == fineBinEdgesPhi.at(0)){
  //   highIndex = lowIndex + NUM_BINS_PHI_ZOOM;
  // }
  // else{
  //   highIndex = lowIndex + NUM_BINS_THETA_ZOOM;
  // }
  
  highIndex = binEdges.size() - 1;
  // get last phi above
  for(int i=binEdges.size()-1; i >= 0; i--){
    if(binEdges.at(i) > highVal){
      highIndex = i;
    }
    else{
      break;
    }
  }
}

void Acclaim::InterferometricMap::setResolutionEstimateFromWaveformSNR(double snr){
  Acclaim::Clustering::getAngularResolution(snr, fSigmaTheta, fSigmaPhi);
  // std::cout << fSigmaTheta << "\t"  << fSigmaPhi << std::endl;
}

TPad* Acclaim::InterferometricMap::makeProjectionCanvas(TPad* pad) {

  if(fUsefulPat){
    double phiWave = fPeakPhi*TMath::DegToRad();
    double thetaWave = -fPeakTheta*TMath::DegToRad();
    double sourceLon, sourceLat, sourceAlt, thetaAdj;
    
    int success = fUsefulPat->traceBackToContinent3(phiWave, thetaWave, &sourceLon, &sourceLat, &sourceAlt, &thetaAdj);
    (void) success;

    TGraphAntarctica* grSource = new TGraphAntarctica(1, &sourceLon, &sourceLat);
    grSource->SetName("Source");

    TPad* thisPad = pad ?  pad : new TCanvas();
    thisPad->cd();
    grSource->SetMarkerStyle(8);

    const TGraphInteractive* grPeak = getPeakPointGraph();
    grSource->SetMarkerColor(grPeak->GetMarkerColor());

    AntarcticaBackground* b = new AntarcticaBackground();
    b->Draw();
    
    grSource->Draw("psame");    
    grSource->SetBit(kCanDelete);

    TGraphAntarctica* grAnita = new TGraphAntarctica(1, &fUsefulPat->longitude,  &fUsefulPat->latitude);
    grAnita->SetName("ANITA");    

    grAnita->SetMarkerColor(kGreen);
    grAnita->SetMarkerStyle(8);
    grAnita->Draw("psame");
    grAnita->SetBit(kCanDelete);


    TArrowAntarctica* arr = new TArrowAntarctica(grAnita,  grSource);
    arr->SetBit(kCanDelete);
    arr->SetFillColor(grSource->GetMarkerColor());
    arr->SetLineColor(grSource->GetMarkerColor());
    arr->Draw();

    
    thisPad->Update();

    double meanEasting = 0.5*(grSource->GetEasting()[0] + grAnita->GetEasting()[0]);
    double meanNorthing = 0.5*(grSource->GetNorthing()[0] + grAnita->GetNorthing()[0]);    

    const double zoomRangeEastingNorthing = 500*1000;
    grSource->GetXaxis()->SetRangeUser(meanEasting - zoomRangeEastingNorthing, meanEasting + zoomRangeEastingNorthing);
    grSource->GetYaxis()->SetRangeUser(meanNorthing - zoomRangeEastingNorthing, meanNorthing + zoomRangeEastingNorthing);

    b->SetGrayScale(true);
    b->SetShowBases(true);
    
    TLegend* l = new TLegend(0.79, 0.69, 0.99, 0.99);
    l->AddEntry(grAnita, "ANITA", "p");
    double thetaHorizon = fUsefulPat->getThetaHorizon(fPeakPhi*TMath::DegToRad());
    double dThetaHorizon = -thetaHorizon*TMath::RadToDeg() - fPeakTheta;
    l->AddEntry(grSource, TString::Format("#delta#theta_{horizon}=%4.2lf", dThetaHorizon), "p");
    l->AddEntry(arr, "Projection", "l");
    l->Draw();

    if(fSigmaTheta > 0 && fSigmaPhi > 0){

      const int numSigma = 3;
      const int numPointsPerContour = 100;
      const double x0 = phiWave;
      const double y0 = thetaWave;
      
      for(int sigma = 1; sigma <= numSigma; sigma++){
	
	TGraphAntarctica* grSigma = new TGraphAntarctica();
	grSigma->SetName(TString::Format("%d_sigma_contour", sigma));
	grSigma->SetBit(kCanDelete);

	const double sigmaThetaRad = sigma*fSigmaTheta*TMath::DegToRad();
	const double sigmaPhiRad = sigma*fSigmaPhi*TMath::DegToRad();

	for(int point=0; point < numPointsPerContour; point++){

	  // trace out an ellipse...
	  // 1.5 rather than 0.5 as +ve  theta is down in the UsefulAdu5Pat
	  // so this will stop lines getting drawn across the top of the ellipse 
	  // if it goes over the horizon
	  double alpha = 1.5*TMath::Pi() + TMath::TwoPi()*point/(numPointsPerContour - 1); // -1 to make a complete loop
	  double x = sigmaThetaRad*cos(alpha);
	  double y = sigmaPhiRad*sin(alpha);
	  
	  double lon, lat, alt, dTheta = 0;
	  int success = fUsefulPat->traceBackToContinent3(x+x0, y+y0, &lon, &lat, &alt, &dTheta);
	  (void) success;
	  if(fabs(dTheta) <= 1e-8){
	    grSigma->SetPoint(grSigma->GetN(), lon, lat);
	  }
	}

	if(grSigma->GetN() > 0){
	  grSigma->SetLineStyle(1+sigma);
	  grSigma->Draw("lsame");
	}
	else{
	  delete grSigma;
	}
      }
    }
    
    return thisPad;
  }
  return pad;
}




const Acclaim::TGraphInteractive* Acclaim::InterferometricMap::getPeakPointGraph(){

  TString name = "grPeak";
  TGraphInteractive* grPeak = const_cast<TGraphInteractive*>(findChild(name));
  if(grPeak==NULL){
    grPeak = new TGraphInteractive(1, &fPeakPhi, &fPeakTheta, "p");
    grPeak->SetName(name);
    addGuiChild(grPeak);
  }
  grPeak->SetMarkerColor(GetLineColor());
  grPeak->SetMarkerStyle(8);
  grPeak->SetMarkerSize(1);
  return grPeak;
  
}


const Acclaim::TGraphInteractive* Acclaim::InterferometricMap::getSunGraph(){

  TString name = "grSun";
  TGraphInteractive* grSun = const_cast<TGraphInteractive*>(findChild(name));
  if(fUsefulPat && grSun==NULL){
    Double_t phiDeg, thetaDeg;
    fUsefulPat->getSunPosition(phiDeg, thetaDeg);

    const double phi0 = getBin0PhiDeg();
    phiDeg = phiDeg < phi0 ? phiDeg + DEGREES_IN_CIRCLE : phiDeg;
    phiDeg = phiDeg >= phi0 + DEGREES_IN_CIRCLE ? phiDeg - DEGREES_IN_CIRCLE : phiDeg;
    thetaDeg*=-1;

    grSun = new TGraphInteractive(1, &phiDeg, &thetaDeg, "p");
    grSun->SetName(name);
    addGuiChild(grSun);
  }
  grSun->SetMarkerStyle(kOpenCircle);
  grSun->SetMarkerSize(1);
  return grSun;
}


const Acclaim::TGraphInteractive* Acclaim::InterferometricMap::getTruthGraph(){

  const char* name = "grTruth";
  TGraphInteractive* grTruth = const_cast<TGraphInteractive*>(findChild(name));

  if(fUsefulPat && grTruth == NULL && truthLon > -999 && truthLat > -999 && truthAlt > -999){
    Double_t thetaDeg, phiDeg;
    fUsefulPat->getThetaAndPhiWave(truthLon, truthLat, truthAlt, thetaDeg, phiDeg);

    thetaDeg*=-TMath::RadToDeg();
    phiDeg*=TMath::RadToDeg();

    const double phi0 = getBin0PhiDeg();
    phiDeg = phiDeg < phi0 ? phiDeg + DEGREES_IN_CIRCLE : phiDeg;
    phiDeg = phiDeg >= phi0 + DEGREES_IN_CIRCLE ? phiDeg - DEGREES_IN_CIRCLE : phiDeg;

    grTruth = new TGraphInteractive(1, &phiDeg, &thetaDeg, "p");
    addGuiChild(grTruth);    
    // std::cerr << phiDeg << "\t" << thetaDeg << "\t" << std::endl;
  }
  if(grTruth){
    grTruth->SetMarkerStyle(kFullStar);
    grTruth->SetMarkerSize(1.5);
    grTruth->SetMarkerColor(kRed);
  }
  return grTruth;
  
}


const Acclaim::TGraphInteractive* Acclaim::InterferometricMap::getHorizonGraph(){

  const char* name = "grHorizon";
  TGraphInteractive* grHorizon = const_cast<TGraphInteractive*>(findChild(name));

  if(grHorizon==NULL && fUsefulPat != nullptr){
    grHorizon = new TGraphInteractive(0, NULL, NULL, "l");

    AntarcticaBackground b;
    Geoid::Position anita(*fUsefulPat);
    Geoid::Position pos;

    for(int by=1; by < b.GetNbinsY(); by++){
      double northing = b.GetYaxis()->GetBinCenter(by);
      for(int bx=1; bx < b.GetNbinsX(); bx++){
	double easting = b.GetXaxis()->GetBinCenter(bx);

	double alt = RampdemReader::SurfaceAboveGeoidEN(easting, northing);
	pos.SetEastingNorthingAlt(easting, northing, alt);

	double d = pos.Distance(anita);

	if(d < 900e3 && d > 500e3){
	  // std::cout << easting << "\t" << northing << "\t" << d << std::endl;

	  double thetaWave, phiWave;
	  fUsefulPat->getThetaAndPhiWaveCart(&pos, thetaWave, phiWave);

	  thetaWave*=-TMath::RadToDeg();
	  phiWave*=TMath::RadToDeg();

	  grHorizon->SetPoint(grHorizon->GetN(), phiWave, thetaWave);

	}
      }
    }

    addGuiChild(grHorizon);
  }
  grHorizon->SetLineWidth(0);
  grHorizon->SetLineColor(GetLineColor());
  grHorizon->SetMarkerColor(kMagenta);
  grHorizon->SetMarkerSize(1);

  return grHorizon;

}



const Acclaim::TGraphInteractive* Acclaim::InterferometricMap::getEdgeBoxGraph(){

  const char* name = "grEdgeBox";
  TGraphInteractive* grEdgeBox = const_cast<TGraphInteractive*>(findChild(name));

  if(grEdgeBox==NULL){
    grEdgeBox = new TGraphInteractive(0, NULL, NULL, "l");
    const TAxis* phiAxis = GetPhiAxis();
    const TAxis* thetaAxis = GetThetaAxis();  
    // left->right across bottom edge
    for(int phiBin=1; phiBin <= GetNbinsPhi()+1; phiBin++){
      grEdgeBox->SetPoint(grEdgeBox->GetN(), phiAxis->GetBinLowEdge(phiBin), thetaAxis->GetBinLowEdge(1));
    }
    // bottom->top on right edge
    for(int thetaBin=1; thetaBin <= GetNbinsTheta()+1; thetaBin++){
      grEdgeBox->SetPoint(grEdgeBox->GetN(), phiAxis->GetBinLowEdge(GetNbinsPhi()+1), thetaAxis->GetBinLowEdge(thetaBin));
    }
    // right->left on top edge
    for(int phiBin=GetNbinsPhi()+1; phiBin > 0; phiBin--){
      grEdgeBox->SetPoint(grEdgeBox->GetN(), phiAxis->GetBinLowEdge(phiBin), thetaAxis->GetBinLowEdge(GetNbinsTheta()+1));
    }
    // top->bottom on left edge
    for(int thetaBin=GetNbinsTheta()+1; thetaBin > 0; thetaBin--){
      grEdgeBox->SetPoint(grEdgeBox->GetN(), phiAxis->GetBinLowEdge(1), thetaAxis->GetBinLowEdge(thetaBin));
    }
    addGuiChild(grEdgeBox);
  }
  grEdgeBox->SetLineWidth(3);
  grEdgeBox->SetLineColor(GetLineColor());  
  
  return grEdgeBox;
  
}



Int_t Acclaim::InterferometricMap::getPhiSectorFromPhiRough(double phiRough){
  Double_t phi0 = getBin0PhiDeg();
  return static_cast<Int_t>((phiRough - phi0)/PHI_RANGE);
}
