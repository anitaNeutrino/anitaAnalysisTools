/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             A Cross Correlator to interact with the ROOTified ANITA-3 data and do some interferometry.
***********************************************************************************************************/

#ifndef INTERFEROMETRIC_MAP_H
#define INTERFEROMETRIC_MAP_H

// #define NUM_BINS_THETA 60
// #define NUM_BINS_THETA 180
// #define NUM_BINS_PHI 50
// #define NUM_BINS_THETA 35
// #define NUM_BINS_PHI 9
#define NUM_BINS_THETA 70
#define NUM_BINS_PHI 9

// #define THETA_RANGE 150
#define MIN_THETA -55
#define MAX_THETA 35
#define PHI_RANGE 22.5


#define DEGREES_IN_CIRCLE 360

#define MAX_NUM_PEAKS 5
#define PEAK_PHI_DEG_RANGE 10
#define PEAK_THETA_DEG_RANGE 10

// Image definitions
// #define NUM_BINS_THETA 100
// #define NUM_BINS_PHI 15
// #define THETA_RANGE 150
// #define PHI_RANGE 22.5


// #define NUM_BINS_THETA_ZOOM 200
// #define NUM_BINS_PHI_ZOOM 200
// #define ZOOM_BINS_PER_DEGREE_PHI 20
// #define ZOOM_BINS_PER_DEGREE_THETA 20

// #define NUM_BINS_THETA_ZOOM 100
// #define NUM_BINS_PHI_ZOOM 200
// #define ZOOM_BINS_PER_DEGREE_PHI 20
// #define ZOOM_BINS_PER_DEGREE_THETA 20

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

#include "AnitaConventions.h"



#include "TH2D.h"
#include "TGraph.h"
#include <map>


namespace Acclaim
{
  
  class InterferometryCache;
  class CrossCorrelator;


  class InterferometricMap : public TH2D {

  public:
    InterferometricMap(); //!< Coarse map constructor (also default constructor for ROOT)
    InterferometricMap(Int_t peakInd, Int_t phiSector, Double_t centrePhi, Double_t phiRange, Double_t centreTheta, Double_t thetaRange); ///!< Fine map constructor


    void Fill(AnitaPol::AnitaPol_t pol, CrossCorrelator* cc, InterferometryCache* dtCache);
    void findPeakValues(Int_t numPeaks, std::vector<Double_t>& peakValues, std::vector<Double_t>& phiDegs, std::vector<Double_t>& thetaDegs) const;
    void getPeakInfo(Double_t& peak, Double_t& phiDeg, Double_t& thetaDeg) const{peak = imagePeak, phiDeg = peakPhiDeg, thetaDeg = peakThetaDeg;}

  
    inline Int_t GetNbinsPhi() const {return GetNbinsX();}
    inline Int_t GetNbinsTheta() const {return GetNbinsY();}  
    inline const TAxis* GetPhiAxis() const {return GetXaxis();}
    inline const TAxis* GetThetaAxis() const {return GetYaxis();}  

    TGraph& getPeakPointGraph(); // for plotting
    TGraph& getEdgeBoxGraph(); // for plotting

    bool isAZoomMap() const {return isZoomMap;}
    Int_t getPeakPhiSector() const {return peakPhiSector;}
    AnitaPol::AnitaPol_t getPol() const {return pol;}  

  protected:

    TGraph& findOrMakeGraph(TString name);
  
    // name and title
    AnitaPol::AnitaPol_t pol;
    UInt_t eventNumber;    
    void setNameAndTitle(); // once we have the eventNumber/polarization
    void setDefaultName(); // before then...
  
  
    bool thetaAxisInSinTheta;
    void initializeInternals();

    // doing the zoomed in maps requires knowing a little more information
    // isZoomMap = false, and all other = -1 if doing a coarse map
    bool isZoomMap;
    int peakPhiSector;
    int minThetaBin;
    int minPhiBin;
    int peakIndex;  // for name and title of zoomed map only
  
    // Get the appropriate bin edges for the zoom map
    void getIndicesOfEdgeBins(const std::vector<double>& binEdges, Double_t lowVal, Double_t highVal, Int_t& lowIndex, Int_t& highIndex); 

    Double_t imagePeak;
    Double_t peakPhiDeg;
    Double_t peakThetaDeg;

    std::map<TString, TGraph> grs;

    ClassDef(InterferometricMap, 1);

    // static members, may end up elsewhere at some point
  public:
    static const std::vector<Double_t>& getCoarseBinEdgesTheta();
    static const std::vector<Double_t>& getFineBinEdgesTheta();
    static const std::vector<Double_t>& getCoarseBinEdgesPhi();
    static const std::vector<Double_t>& getFineBinEdgesPhi();  
  
    static Double_t getBin0PhiDeg();
  
     
  };
}



#endif // INTERFEROMETRIC_MAP_H
