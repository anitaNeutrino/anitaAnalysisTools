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
    Draw("colz");

    TGraph& grPeak = getPeakPointGraph();
    grPeak.Draw("psame");

    if(fUsefulPat){
      
      TGraph& grSun = getSunGraph();
      grSun.Draw("psame");
    }
    
    
  }  
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
    : fUsefulPat(NULL)
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
Acclaim::InterferometricMap::InterferometricMap() : fUsefulPat(NULL){
  fIsZoomMap = false;
  fPeakPhiSector = -1;
  minPhiBin = -1;
  minThetaBin = -1;
  peakIndex = -1;

  const std::vector<double> coarseBinsPhi = Acclaim::InterferometricMap::getCoarseBinEdgesPhi();
  const std::vector<double> coarseBinsTheta = Acclaim::InterferometricMap::getCoarseBinEdgesTheta();
  SetBins(coarseBinsPhi.size()-1, &coarseBinsPhi[0], coarseBinsTheta.size()-1, &coarseBinsTheta[0]);

  initializeInternals();  
}



/** 
 * Default destructor
 */
Acclaim::InterferometricMap::~InterferometricMap(){
  // std::cerr << __PRETTY_FUNCTION__ << "\t" << this << "\t" << fUsefulPat << std::endl;  
  deletePat();
}


/** 
 * If we have a copy of the usefulAdu5Pat, then delete it nicely.
 */
void Acclaim::InterferometricMap::deletePat(){
  // std::cerr << __PRETTY_FUNCTION__ << "\t" << this << "\t" << fUsefulPat << std::endl;
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
 * @brief For making a heat map.
 * Loops over the TProfile2D projection bins (easting/northing).
 * For each bin, gets theta/phi at the payload and does bilinear interpolation of the map value.
 * This is currently VERY slow, and so I'm not recommending its use.
 *
 * @param proj should be a histogram of Antarctica in Easting/Northing
 * @param horizonKilometers is a cut-off distance at which to stop in km.
 */
void Acclaim::InterferometricMap::project(TProfile2D* proj, double horizonKilometers){

  if(!fUsefulPat){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", GPS info required to project interferometric map." << std::endl;
    return;
  }
  if(!proj){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", require projection histogram to be non-NULL." << std::endl;
    return;
  }


  Int_t nx = proj->GetNbinsX();
  Int_t ny = proj->GetNbinsY();

  const double closeNorthing = proj->GetYaxis()->GetBinLowEdge(2) - proj->GetYaxis()->GetBinLowEdge(1);
  const double closeEasting = proj->GetXaxis()->GetBinLowEdge(2) - proj->GetXaxis()->GetBinLowEdge(1);

  const int nPhi = GetNbinsPhi();
  const double phiMin = fXaxis.GetBinLowEdge(1);
  const double phiMax = fXaxis.GetBinLowEdge(nPhi);

  const double thetaMin = fYaxis.GetBinLowEdge(1);
  const double thetaMax = fYaxis.GetBinLowEdge(GetNbinsTheta());
  const double floatPointErrorEpsilon = 0.0001;
  for(int by=1; by <= ny; by++){
    double northing = proj->GetYaxis()->GetBinLowEdge(by);

    for(int bx=1; bx <= nx; bx++){
      double easting = proj->GetXaxis()->GetBinLowEdge(bx);

      double lon, lat;
      RampdemReader::EastingNorthingToLonLat(easting, northing, lon, lat);
      double alt = RampdemReader::SurfaceAboveGeoid(lon, lat);

      double dist_km = 1e-3*fUsefulPat->getDistanceFromSource(lat, lon, alt);
      if(dist_km < horizonKilometers){

        // first get the theta/phi at the payload
        Double_t thetaWave, phiWave;
        fUsefulPat->getThetaAndPhiWave(lon, lat, alt, thetaWave, phiWave);

        // this is a naive collision detection...
        // I'm going to trace those theta/phi back to the continent
        // if it hits the same place (or very close) then we will go ahead
        double lat2, lon2,  alt2, theta_adj;
        fUsefulPat->traceBackToContinent(phiWave, thetaWave, &lat2, &lon2, &alt2, &theta_adj);

        double easting2, northing2;
        RampdemReader::LonLatToEastingNorthing(lon2, lat2, easting2, northing2);

        const double deltaEasting = TMath::Abs(easting2 - easting);
        const double deltaNorthing = TMath::Abs(northing2 - northing);

        if(deltaEasting < closeEasting && deltaNorthing < closeNorthing){

          Double_t thetaDeg = -1*thetaWave*TMath::RadToDeg();

          if(thetaDeg >= thetaMin && thetaDeg <= thetaMax){
            Double_t phiDeg = phiWave*TMath::RadToDeg();
            phiDeg = phiDeg <  phiMin ? phiDeg + DEGREES_IN_CIRCLE : phiDeg;
            phiDeg = phiDeg >= phiMax ? phiDeg - DEGREES_IN_CIRCLE : phiDeg;

            if(phiDeg >= phiMin && phiDeg <= phiMax){


              // my axis binning is that the low edge of the bin is the actual value
              // so can't use TH2::Interpolate, which (sensibly) uses the bin centres
              // so here's my implementation.

              // Get the grid points and values
              Int_t phiBinLow = fXaxis.FindBin(phiDeg);
              Int_t thetaBinLow = fYaxis.FindBin(thetaDeg);
              Int_t phiBinHigh = phiBinLow == nPhi ? 1 : phiBinLow + 1;

              double x1 = fXaxis.GetBinLowEdge(phiBinLow);
              double x2 = fXaxis.GetBinLowEdge(phiBinHigh);

              double y1 = fYaxis.GetBinLowEdge(thetaBinLow);
              double y2 = fYaxis.GetBinLowEdge(thetaBinLow+1);

              double v11 = GetBinContent(phiBinLow, thetaBinLow);
              double v12 = GetBinContent(phiBinLow, thetaBinLow+1);
              double v21 = GetBinContent(phiBinHigh, thetaBinLow);
              double v22 = GetBinContent(phiBinHigh, thetaBinLow+1);

              // this image is very helpful to visualize what's going on here
              // https://en.wikipedia.org/wiki/Bilinear_interpolation#/media/File:Comparison_of_1D_and_2D_interpolation.svg
              // need to do three linear interpolations...
              double dx = phiDeg - x1;
              double dy = thetaDeg - y1;
              double deltaX = x2 - x1 > 0 ? x2 - x1 : DEGREES_IN_CIRCLE + x2 - x1;
              double deltaY = y2 - y1;

              double v11_v21_interp = v11 + dx*(v21 - v11)/(deltaX);
              double v12_v22_interp = v12 + dx*(v22 - v12)/(deltaX);
              double val = v11_v21_interp + dy*(v12_v22_interp - v11_v21_interp)/(deltaY);

              // std::cout << x1  << "\t" << x2  << "\t" << y1  << "\t" <<  y2 << std::endl;
              // std::cout << v11 << "\t" << v12 << "\t" << v21 << "\t" << v22 << std::endl;
              // std::cout << v11_v21_interp << "\t" << v12_v22_interp << "\t" << val << std::endl;

              proj->Fill(easting+floatPointErrorEpsilon, northing+floatPointErrorEpsilon, val);
            }
          }
        }
      }
    }
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
    
    std::vector<Int_t>* combosToUse = &cc->combosToUseGlobal[fPeakPhiSector];

    const Int_t offset = cc->numSamplesUpsampled/2;
    for(UInt_t comboInd=0; comboInd<combosToUse->size(); comboInd++){
      Int_t combo = combosToUse->at(comboInd);
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
	Double_t partBA = dtCache->partBAsZoom[dtCache->partBAsIndex(pol, combo, zoomThetaInd)]; //)[pol][combo][zoomThetaInd];      
	// Double_t dtFactor = dtFactors[zoomThetaInd];
	Double_t dtFactor = dtCache->dtFactors[zoomThetaInd];

	// for(Int_t phiBin = 0; phiBin < NUM_BINS_PHI_ZOOM; phiBin++){
	for(Int_t phiBin = 0; phiBin < GetNbinsPhi(); phiBin++){
	  Int_t zoomPhiInd = minPhiBin + phiBin;
	  // Double_t zoomPhiWave = zoomedPhiWaveLookup[zoomPhiInd];

	  int p21 = dtCache->part21sIndex(pol, combo, zoomPhiInd);
	  Double_t offsetLowDouble = dtFactor*(partBA - dtCache->part21sZoom[p21]);//[pol][combo][zoomPhiInd]);		

	  offsetLowDouble += dtCache->kUseOffAxisDelay > 0 ? dtCache->offAxisDelaysDivided[p21] : 0;	  
	
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

    Double_t normFactor = cc->kOnlyThisCombo < 0 && combosToUse->size() > 0 ? combosToUse->size() : 1;
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
  TDecompSVD& svd = (TDecompSVD&) getSVD(); // cast away const-ness
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

  // std::cerr << a << "\t" << b << "\t" << c << "\t" << d << "\t" << e << "\t" << f << std::endl;
  // std::cerr << fPeakPhi << "\t" << fPeakTheta << "\t" << fPeakValue << "\t" << GetBinContent(peakPhiBin+1, peakThetaBin+1) << std::endl;
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



TGraph& Acclaim::InterferometricMap::findOrMakeGraph(TString name){

  std::map<TString, TGraph>::iterator it = grs.find(name);
  if(it!=grs.end()){
    return it->second;
  }
  else{
    TGraph gr;
    gr.SetName(name);
    grs[name] = gr;
    return grs[name];
  }  
}


TGraph& Acclaim::InterferometricMap::getPeakPointGraph(){

  TGraph& gr = findOrMakeGraph("grPeak");
  if(gr.GetN()==0){
    gr.SetPoint(0, fPeakPhi, fPeakTheta);
  }
  return grs["grPeak"];
  
}


TGraph& Acclaim::InterferometricMap::getSunGraph(){

  TGraph& grSun = findOrMakeGraph("grSun");

  if(fUsefulPat){
    Double_t phiDeg, thetaDeg;
    fUsefulPat->getSunPosition(phiDeg, thetaDeg);

    const double phi0 = getBin0PhiDeg();
    phiDeg = phiDeg < phi0 ? phiDeg + DEGREES_IN_CIRCLE : phiDeg;
    phiDeg = phiDeg >= phi0 + DEGREES_IN_CIRCLE ? phiDeg - DEGREES_IN_CIRCLE : phiDeg;
    thetaDeg*=-1;

    grSun.SetPoint(0, phiDeg, thetaDeg);
  }
  return grs["grSun"];
  
}


TGraph& Acclaim::InterferometricMap::getEdgeBoxGraph(){

  TGraph& gr = findOrMakeGraph("grEdgeBox");
  const TAxis* phiAxis = GetPhiAxis();
  const TAxis* thetaAxis = GetThetaAxis();  
  // left->right across bottom edge
  for(int phiBin=1; phiBin <= GetNbinsPhi()+1; phiBin++){
    gr.SetPoint(gr.GetN(), phiAxis->GetBinLowEdge(phiBin), thetaAxis->GetBinLowEdge(1));
  }
  // bottom->top on right edge
  for(int thetaBin=1; thetaBin <= GetNbinsTheta()+1; thetaBin++){
    gr.SetPoint(gr.GetN(), phiAxis->GetBinLowEdge(GetNbinsPhi()+1), thetaAxis->GetBinLowEdge(thetaBin));
  }
  // right->left on top edge
  for(int phiBin=GetNbinsPhi()+1; phiBin > 0; phiBin--){
    gr.SetPoint(gr.GetN(), phiAxis->GetBinLowEdge(phiBin), thetaAxis->GetBinLowEdge(GetNbinsTheta()+1));
  }
  // top->bottom on left edge
  for(int thetaBin=GetNbinsTheta()+1; thetaBin > 0; thetaBin--){
    gr.SetPoint(gr.GetN(), phiAxis->GetBinLowEdge(1), thetaAxis->GetBinLowEdge(thetaBin));
  }
  
  return grs["grEdgeBox"];
  
}
