/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             Strictly templated cut functions for the AnalysisPlot class
***********************************************************************************************************/

#ifndef ANALYSIS_CUTS_H
#define ANALYSIS_CUTS_H

#include "TString.h"
#include "AnitaConventions.h"
#include "AnitaEventSummary.h"

#include <iostream>


namespace Acclaim
{



  /**
   * @namespace AnalysisCuts
   * @brief Contains a set of classes which implement an cuts on AnitaEventSummary objects
   */
  namespace AnalysisCuts
  {



    /**
     * @brief A class to perform a simple cut using the results stored in an AnitaEventSummary
     *
     * The base AnalysisCut class is designed to be inherited from.
     * It contains a purely virtual apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakIndex) const function.
     * The apply function is then overloaded in the derived classes where the actual cuts are implemented.
     * To work properly with the AnalysisPlot class, apply must return an integer from 0 to fMaxRetVal-1 (inclusive).
     * In addition to the expected return values, the name/title members are also used in the AnalysisPlot for histogram names and titles.
     * 
     * @todo Make this inherit from TCut, and have the functionality match apply() somehow.
     */
    class AnalysisCut {
    public:

      AnalysisCut(const char* name, const char* title, int numApplyRetVals=2);
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const = 0; // to be overloaded with actual cut
      inline int getMaximumReturnValue() const {return fMaxRetVal;}
      inline const char* getName() const {return fName.Data();}
      inline const char* getTitle() const {return fTitle.Data();}
    protected:
      /** 
       * Wrapper function to get the "trainingPeak", which is defined by the member AnitaEventSummary member function 
       * 
       * @param sum is the AnitaEventSummary on which to call the function
       * @param pol is set to the polarisation of interest
       * @param peakInd is set to the peak index of interest
       */
      inline static void handleDefaults(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t& pol, Int_t& peakInd){
	if(peakInd==-1){
	  peakInd = sum->trainingPeakInd();
	  pol = sum->trainingPol();
	}
      }
      TString fName;
      TString fTitle;
      int fMaxRetVal;
    };



    class IsAboveHorizontal : public AnalysisCut
    {
    public:
      IsAboveHorizontal() : AnalysisCut("isAboveHorizontal", "Above Horizontal") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
    };

    class IsTaggedAsWaisPulser : public AnalysisCut
    {
    public:
      IsTaggedAsWaisPulser() : AnalysisCut("isTaggedAsWaisPulser", "Tagged As WAIS Pulser") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
    };

    class IsTaggedAsLDBPulser : public AnalysisCut
    {
    public:
      IsTaggedAsLDBPulser() : AnalysisCut("isTaggedAsLDBPulser", "Tagged As LDB Pulser") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
    };

    class HigherPol : public AnalysisCut
    {
    public:
      HigherPol() : AnalysisCut("higherPol", "Polarization", AnitaPol::kNotAPol) {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns 0 for HPol, 1 for VPol
    };

    class HasSourceLocation : public AnalysisCut
    {
    public:
      HasSourceLocation() : AnalysisCut("hasSourceLocation", "Has Source Location") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
    };

    class IsOnContinent : public AnalysisCut
    {
    public:
      IsOnContinent() : AnalysisCut("IsOnContinent", "Reconstructs to Land") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
    };

    class IsTaggedAsPayloadBlast : public AnalysisCut
    {
    public:
      IsTaggedAsPayloadBlast() : AnalysisCut("isTaggedAsPayloadBlast", "Tagged as Payload Blast") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
    };

    class IsWithin20DegreesOfSunInPhi : public AnalysisCut{
    public:
      IsWithin20DegreesOfSunInPhi() : AnalysisCut("isWithin20DegreesOfSunInPhi", "|#delta#phi_{sun}| < 20") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
    };

    class IsGood : public AnalysisCut{
    public:
      IsGood() : AnalysisCut("isGood", "Is Good") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
    };


    class GoodGPS : public AnalysisCut{
    public:
      GoodGPS() : AnalysisCut("goodGPS", "Good GPS") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
    };

    class NonZeroStokesI : public AnalysisCut{
    public:
      NonZeroStokesI() : AnalysisCut("nonZeroStokesI", "Non-zero Stokes I") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
    };

    class RealSNR : public AnalysisCut{
    public:
      RealSNR() : AnalysisCut("realSNR", "Non-NaN SNR") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
    };


    class Anita3QuietTime : public AnalysisCut{
    public:
      Anita3QuietTime() : AnalysisCut("quietTime", "Quiet Time") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
    };


    class CloseToMC : public AnalysisCut{
    public:
      CloseToMC() : AnalysisCut("closeToMC", "Near MC") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
    };

    class CloseToWais : public AnalysisCut{
    public:
      CloseToWais() : AnalysisCut("closeToWais", "Near Wais") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
    };


    class IsRfTrigger : public AnalysisCut{
    public:
      IsRfTrigger() : AnalysisCut("isRfTrigger", "RF trigger") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns the peakIndex + 1 if VPol 
    };


    class SmallDeltaRough : public AnalysisCut{
    public:
      SmallDeltaRough() : AnalysisCut("smallDeltaRough", "Small angle between coarse/fine maps") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns the peakIndex + 1 if VPol 
    };


    class IsNotTaggedAsPulser : public AnalysisCut{
    public:
      IsNotTaggedAsPulser() : AnalysisCut("isNotTaggedAsPulser", "Not pulser") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns the peakIndex + 1 if VPol 
    };

    class SignalLikeFirstStandardizedPeakMoments : public AnalysisCut{
    public:
      SignalLikeFirstStandardizedPeakMoments() : AnalysisCut("signalLikeFirstStandardizedPeakMoments", "Signal-like first standardized peak moments") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns the peakIndex + 1 if VPol 
    };

    class PassesThesisCuts : public AnalysisCut{
    public:
      PassesThesisCuts() : AnalysisCut("passesThesisCuts", "Passes thesis thermal cut") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// true/false
    };


    class IsNotNorth : public AnalysisCut {
    public:
      IsNotNorth() : AnalysisCut("isNotNorth", "Is Not North") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// true/false
    };


    class HigherPeakHilbertAfterDedispersion : public AnalysisCut {
    public:
      HigherPeakHilbertAfterDedispersion() : AnalysisCut("HigherPeakHilbertAfterDedispersion", "Hilbert peak is higher after dedispersion") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const;
    };

    
    class HigherImpulsivityMeasureAfterDedispersion : public AnalysisCut {
    public:
      HigherImpulsivityMeasureAfterDedispersion() : AnalysisCut("HigherPeakHilbertAfterDedispersion", "ImpulsivityMeasure is higher after dedispersion") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const;
    };

    class LowerFracPowerWindowGradientAfterDedispersion : public AnalysisCut {
    public:
      LowerFracPowerWindowGradientAfterDedispersion() : AnalysisCut("HigherPeakHilbertAfterDedispersion", "Lower fracPowerWindowGradient after dedispersion") {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const;
    };



    // const globals so you don't need to instantiate these yourself
    const IsAboveHorizontal isAboveHorizontal;
    const IsTaggedAsWaisPulser isTaggedAsWaisPulser;
    const HigherPol higherPol;
    const HasSourceLocation hasSourceLocation;
    const IsOnContinent isOnContinent;
    const IsTaggedAsPayloadBlast isTaggedAsPayloadBlast;
    const IsWithin20DegreesOfSunInPhi isWithin20DegreesOfSunInPhi;
    const IsGood isGood;
    const GoodGPS goodGPS;
    const NonZeroStokesI nonZeroStokesI;
    const RealSNR realSNR;
    const Anita3QuietTime anita3QuietTime;
    const CloseToMC closeToMC;
    const CloseToWais closeToWais;
    const IsRfTrigger isRfTrigger;
    const SmallDeltaRough smallDeltaRough;
    const IsNotTaggedAsPulser isNotTaggedAsPulser;
    const SignalLikeFirstStandardizedPeakMoments signalLikeFirstStandardizedPeakMoments;
    const PassesThesisCuts passesThesisThermalCut;
    const IsNotNorth isNotNorth;
    const HigherPeakHilbertAfterDedispersion higherPeakHilbertAfterDedispersion;
    const HigherImpulsivityMeasureAfterDedispersion higherImpulsivityMeasureAfterDedispersion;
    const LowerFracPowerWindowGradientAfterDedispersion lowerFracPowerWindowGradientAfterDedispersion;  
  }
}



#endif
