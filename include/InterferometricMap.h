/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             A Cross Correlator to interact with the ROOTified ANITA-3 data and do some interferometry.
***********************************************************************************************************/

#ifndef INTERFEROMETRIC_MAP_H
#define INTERFEROMETRIC_MAP_H

#define NUM_BINS_THETA 70
#define NUM_BINS_PHI 9

#define MIN_THETA -55
#define MAX_THETA 35
#define PHI_RANGE 22.5


#define DEGREES_IN_CIRCLE 360

#define MAX_NUM_PEAKS 5
#define PEAK_PHI_DEG_RANGE 10
#define PEAK_THETA_DEG_RANGE 10

#define NUM_BINS_THETA_ZOOM 40
#define NUM_BINS_PHI_ZOOM 40
#define ZOOM_BINS_PER_DEGREE_PHI 4
#define ZOOM_BINS_PER_DEGREE_THETA 4

#define ZOOM_BIN_SIZE_PHI (1./ZOOM_BINS_PER_DEGREE_PHI)
#define ZOOM_BIN_SIZE_THETA (1./ZOOM_BINS_PER_DEGREE_THETA)
#define THETA_RANGE_ZOOM (NUM_BINS_THETA_ZOOM*ZOOM_BIN_SIZE_THETA)
#define PHI_RANGE_ZOOM (NUM_BINS_PHI_ZOOM*ZOOM_BIN_SIZE_PHI)

#define NUM_BINS_PHI_ZOOM_TOTAL (DEGREES_IN_CIRCLE*ZOOM_BINS_PER_DEGREE_PHI + NUM_BINS_PHI_ZOOM)
#define NUM_BINS_THETA_ZOOM_TOTAL ((MAX_THETA - MIN_THETA)*ZOOM_BINS_PER_DEGREE_THETA + NUM_BINS_THETA_ZOOM)

#define NUM_BINS_QUAD_FIT_PHI 5
#define NUM_BINS_QUAD_FIT_THETA 5

#include "AnitaConventions.h"

#include "TH2D.h"
#include "TGraph.h"
#include <map>
#include "TDecompSVD.h"
#include "TGraphInteractive.h"

class Adu5Pat;
class UsefulAdu5Pat;
class TProfile2D;
class TruthAnitaEvent;

namespace Acclaim
{
  
class InterferometryCache;
class CrossCorrelator;


class InterferometricMap : public TH2D, public GuiParent {

 public:
  InterferometricMap(); /// Coarse map constructor (also default constructor for ROOT)
  InterferometricMap(Int_t peakInd, Int_t phiSector, Double_t centrePhi, Double_t phiRange, Double_t centreTheta, Double_t thetaRange); //// Fine map constructor
  virtual ~InterferometricMap();

  void Fill(AnitaPol::AnitaPol_t pol, CrossCorrelator* cc, InterferometryCache* dtCache);
  void findPeakValues(Int_t numPeaks, std::vector<Double_t>& peakValues, std::vector<Double_t>& phiDegs, std::vector<Double_t>& thetaDegs) const;
  void getPeakInfo(Double_t& peak, Double_t& phiDeg, Double_t& thetaDeg) const{
    peak = fPeakBinValue;
    phiDeg = fPeakPhi;
    thetaDeg = fPeakTheta;
  }
  void getPeakInfo(Double_t& peak, Double_t& phiDeg, Double_t& thetaDeg, Double_t& reducedChisquare) const{
    getPeakInfo(peak, phiDeg, thetaDeg);
    reducedChisquare = fPeakReducedChisquare;
  }
  
  virtual void Draw(Option_t* opt);
  virtual void ExecuteEvent(int event, int x, int y);

  void addGpsInfo(const Adu5Pat* pat);
  void addTruthInfo(const TruthAnitaEvent* truth);
  // void addGpsInfo(const UsefulAdu5Pat* usefulPat);


  /** 
   * @brief Project this peak direction to the surface and draw a summary on the pad
   * 
   * Requires the GPS data to work, i.e. #fUsefulPat is non-NULL.
   * @see addGpsInfo(const Adu5Pat* pat)
   * 
   * @param pad is the pad to draw on, makes a new TCanvas if null
   * 
   * @return pointer to the new canvas if pad is NULL, pad otherwise
   */
  TPad* makeProjectionCanvas(TPad* pad);

  /** 
   * Version of makeProjectionCanvas(TPad* pad) for ROOT context menu
   * 
   * @return same as makeProjectionCanvas(TPad* pad)
   */
  TPad* makeProjectionCanvas(){return makeProjectionCanvas(NULL);} //*MENU*

  inline Int_t GetNbinsPhi() const {return GetNbinsX();}
  inline Int_t GetNbinsTheta() const {return GetNbinsY();}
  inline const TAxis* GetPhiAxis() const {return GetXaxis();}
  inline const TAxis* GetThetaAxis() const {return GetYaxis();}

  const TGraphInteractive* getPeakPointGraph(); // for plotting
  const TGraphInteractive* getEdgeBoxGraph(); // for plotting
  const TGraphInteractive* getSunGraph(); // for plotting
  const TGraphInteractive* getTruthGraph(); // for plotting

  bool isAZoomMap() const {return fIsZoomMap;}
  Int_t getPeakPhiSector() const {return fPeakPhiSector;}
  AnitaPol::AnitaPol_t getPol() const {return pol;}

  virtual void Reset(Option_t* = "");

  /** 
   * Sets values for #fSigmaTheta and #fSigmaPhi which allows contours to be drawn in makeProjectionCanvas()
   * 
   * @param snr is the estimate of the waveform SNR
   * 
   * @see Acclaim::Clustering::ResolutionModel
   * @see makeProjectionCanvas(TPad* pad)
   */
  void setResolutionEstimateFromWaveformSNR(double snr);
  
    
 protected:

  UsefulAdu5Pat* fUsefulPat;
  void deletePat();
  Double_t truthLat;
  Double_t truthLon;
  Double_t truthAlt;

  void fitPeakWithQuadratic(Int_t peakPhiBin, Int_t peakThetaBin);
  void setPeakInfoJustFromPeakBin(Int_t peakPhiBin, Int_t peakThetaBin);

  TGraph& findOrMakeGraph(TString name);
  
  // name and title
  AnitaPol::AnitaPol_t pol;
  UInt_t eventNumber;    
  void setNameAndTitle(); // once we have the eventNumber/polarization
  void setDefaultName(); // before then...

  bool fThetaAxisInSinTheta;
  void initializeInternals();



  // doing the zoomed in maps requires knowing a little more information
  // isZoomMap = false, and all other = -1 if doing a coarse map
  bool fIsZoomMap;
  int fPeakPhiSector;
  int minThetaBin;
  int minPhiBin;
  int peakIndex;  // for name and title of zoomed map only
  
  // Get the appropriate bin edges for the zoom map
  void getIndicesOfEdgeBins(const std::vector<double>& binEdges, Double_t lowVal, Double_t highVal, Int_t& lowIndex, Int_t& highIndex); 

  Double_t fPeakBinValue; /// The value of the highest bin in the map
  Int_t fPeakBinPhi; /// The bin in phi (counting from 0) containing the peak value
  Int_t fPeakBinTheta; /// The bin in theta (counting from 0) containing the peak value

  Double_t fPeakValue; /// Value of the peak of quadratic fit to the region around the peak bin
  Double_t fPeakPhi; /// Location in phi of the peak of the quadratic fit to the region around the peak bin
  Double_t fPeakTheta; /// Location in theta of the peak of the quadratic fit to the region around the peak bin
  Double_t fPeakReducedChisquare; /// Residual of the peak fit... I hope
  Double_t fPeakSigmaPhi; /// Width of peak in phi
  Double_t fPeakSigmaTheta; /// Width of peak theta
  Double_t fPeakCovariance; /// Width of peak theta

  Double_t fSigmaTheta; /// Estimate of the pointing resolution from the waveform SNR @see setResolutionEstimateFromWaveformSNR(double snr)
  Double_t fSigmaPhi; /// Estimate of of the pointing resultion from the waveform SNR @see setResolutionEstimateFromWaveformSNR(double snr)

  ClassDef(InterferometricMap, 4);

  // static members, may end up elsewhere at some point
 public:
  static const std::vector<Double_t>& getCoarseBinEdgesTheta();
  static const std::vector<Double_t>& getFineBinEdgesTheta();
  static const std::vector<Double_t>& getCoarseBinEdgesPhi();
  static const std::vector<Double_t>& getFineBinEdgesPhi();
  static Double_t getBin0PhiDeg();
  static Double_t getPhiSectorCenterPhiDeg(int phi);
  static Int_t getPhiSectorFromPhiRough(double phiRough);

 private:
  static const TDecompSVD& getSVD();
  static void getDeltaBinRangeSVD(Int_t& dPhiBinLow, Int_t& dPhiBinHigh, Int_t& dThetaBinLow, Int_t& dThetaBinHigh);
    
};
}



#endif // INTERFEROMETRIC_MAP_H
