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
#include "AnalysisSettings.h"

namespace Acclaim
{



  /**
   * @namespace AnalysisCuts
   * @brief Contains a set of classes which implement an cuts on AnitaEventSummary objects
   */
  namespace AnalysisCuts
  {

    /** 
     * @enum Mode defines the behaviour of AnalysisCut::handleDefaults(...)
     */
    enum Mode {
      kTraining,	/// Makes handleDefaults call the AnitaEventSummary::training family of functions 
      kAcclaimAnalysis  /// Makes handleDefaults call the AnitaEventSummary::mostImpulsive(1) function
    };

    Mode getMode();
    void setMode(Mode mode);

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
    class AnalysisCut : public TObject {
    public:

      AnalysisCut(int numApplyRetVals=2);
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const = 0; // to be overloaded with actual cut
      inline int getMaximumReturnValue() const {return fMaxRetVal;}
    protected:
      void handleDefaults(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t& pol, Int_t& peakInd) const;
      int fMaxRetVal;

      ClassDef(AnalysisCut, 0);
    };

    /**
     * @class IsAboveHorizontal
     * @brief Checks whether a given peak is above the horizontal
     */

    class IsAboveHorizontal : public AnalysisCut
    {
    public:
      IsAboveHorizontal() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
      ClassDef(IsAboveHorizontal, 0);
    };


    /**
     * @class IsTaggedAsWaisPulser
     * @brief Is the triggerTime stamp consistent with WAIS divide?
     */

    class IsTaggedAsWaisPulser : public AnalysisCut
    {
    public:
      IsTaggedAsWaisPulser() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
      ClassDef(IsTaggedAsWaisPulser, 0);
    };



    /**
     * @class IsTaggedAsPulser
     * @brief Is the triggerTime stamp consistent with any pulser?
     */

    class IsTaggedAsPulser : public AnalysisCut
    {
    public:
      IsTaggedAsPulser() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
      ClassDef(IsTaggedAsPulser, 0);
    };


    /**
     * @class IsTaggedAsLDBPulser
     * @brief Is the triggerTime stamp consistent with the LDB pulser?
     */
    class IsTaggedAsLDBPulser : public AnalysisCut
    {
    public:
      IsTaggedAsLDBPulser() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
      ClassDef(IsTaggedAsLDBPulser, 0);
    };



    /**
     * @class HigherPol
     * @brief Which polarisation has a higher map peak?
     */

    class HigherPol : public AnalysisCut
    {
    public:
      HigherPol() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns 0 for HPol, 1 for VPol
      ClassDef(HigherPol, 0);
    };

    /**
     * @class HasSourceLocation
     * @brief Did the event interest with the surface model of the Earth?
     */
    class HasSourceLocation : public AnalysisCut
    {
    public:
      HasSourceLocation() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
      ClassDef(HasSourceLocation, 0);
    };


    /**
     * @class IsOnContinent
     * @brief Did the event reconstruct to the Antarctic landmass (including ice shelfs)?
     */
    class IsOnContinent : public AnalysisCut
    {
    public:
      IsOnContinent() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
      ClassDef(IsOnContinent, 0);
    };


    /**
     * @class IsTaggedAsPayloadBlast
     * @brief Was this event tagged as a payload blast in my initial reconstruction?
     */
    class IsTaggedAsPayloadBlast : public AnalysisCut
    {
    public:
      IsTaggedAsPayloadBlast() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
      ClassDef(IsTaggedAsPayloadBlast, 0);
    };


    /**
     * @class IsWithin20DegreesOfSunInPhi
     * @brief Is this event within 20 degrees of the sun in payload phi?
     */
    class IsWithin20DegreesOfSunInPhi : public AnalysisCut{
    public:
      IsWithin20DegreesOfSunInPhi() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
      ClassDef(IsWithin20DegreesOfSunInPhi, 0);
    };

    /**
     * @class IsGood
     * @brief Is the is good flag set for this event?
     */
    class IsGood : public AnalysisCut{
    public:
      IsGood() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
      ClassDef(IsGood, 0);
    };


    /**
     * @class GoodGPS
     * @brief Is the GPS heading a sensible number?
     */
    class GoodGPS : public AnalysisCut{
    public:
      GoodGPS() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
      ClassDef(GoodGPS, 0);
    };


    /**
     * @class NonZeroStokesI
     * @brief Is the stokes I parameter non-zero?
     */
    class NonZeroStokesI : public AnalysisCut{
    public:
      NonZeroStokesI() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
      ClassDef(NonZeroStokesI, 0);
    };

    /**
     * @class RealSNR
     * @brief Is the SNR a real number (i.e. not a NaN)?
     */
    class RealSNR : public AnalysisCut{
    public:
      RealSNR() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
      ClassDef(RealSNR, 0);
    };

    /**
     * @class Anita3QuietTime
     * @brief Does this event fall in my choice of the quiet portion of the ANITA-3 flight time?
     */
    class Anita3QuietTime : public AnalysisCut{
    public:
      Anita3QuietTime() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
      ClassDef(Anita3QuietTime, 0);
    };

    /**
     * @class CloseToMC
     * @brief Is this event close to the monte carlo?
     */
    class CloseToMC : public AnalysisCut{
    public:
      CloseToMC() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
      ClassDef(CloseToMC, 0);
    };

    class CloseToWais : public AnalysisCut{
    public:
      CloseToWais() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns false(0) or true(1)
      ClassDef(CloseToWais, 0);
    };

    /**
     * @class IsRfTrigger
     * @brief Is this an event with an RF trigger?
     */
    class IsRfTrigger : public AnalysisCut{
    public:
      IsRfTrigger() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns the peakIndex + 1 if VPol
      ClassDef(IsRfTrigger, 0);
    };

    /**
     * @class SmallDeltaRough
     * @brief Do the rough and fine maps broadly agree in their reconstructed direction? 
     */
    class SmallDeltaRough : public AnalysisCut{
    public:
      SmallDeltaRough() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns the peakIndex + 1 if VPol
      ClassDef(SmallDeltaRough, 0);
    };


    /**
     * @class IsNotTaggedAsPulser
     * @brief Is the trigger time consistent with a pulser?
     */
    class IsNotTaggedAsPulser : public AnalysisCut{
    public:
      IsNotTaggedAsPulser() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns the peakIndex + 1 if VPol
      ClassDef(IsNotTaggedAsPulser, 0);
    };

    /**
     * @class SignalLikeFirstStandardizedPeakMoments
     * @brief An attempt at a cut using the moments of the coherently summed waves
     */
    class SignalLikeFirstStandardizedPeakMoments : public AnalysisCut{
    public:
      SignalLikeFirstStandardizedPeakMoments() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// Returns the peakIndex + 1 if VPol
      ClassDef(SignalLikeFirstStandardizedPeakMoments, 0);
    };

    /**
     * @class PassesThesisCuts
     * @brief Incomplete thesis cut
     */
    class PassesThesisCuts : public AnalysisCut{
    public:
      PassesThesisCuts() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// true/false
      ClassDef(PassesThesisCuts, 0);
    };

    /**
     * @class IsNotNorth
     * @brief Is the peak event pointing away from North?
     */
    class IsNotNorth : public AnalysisCut {
    public:
      IsNotNorth() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, Int_t peakInd = -1) const; /// true/false
      ClassDef(IsNotNorth, 0);
    };


    /**
     * @class HigherPeakHilbertAfterDedispersion
     * @brief Does the hilbert peak get higher after dedispersing?
     */
    class HigherPeakHilbertAfterDedispersion : public AnalysisCut {
    public:
      HigherPeakHilbertAfterDedispersion() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const;
      ClassDef(HigherPeakHilbertAfterDedispersion, 0);
    };


    /**
     * @class HigherImpulsivityMeasureAfterDedispersion
     * @brief Does the impulsivity measure get higher after dedispersion?
     */
    class HigherImpulsivityMeasureAfterDedispersion : public AnalysisCut {
    public:
      HigherImpulsivityMeasureAfterDedispersion() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const;
      ClassDef(HigherImpulsivityMeasureAfterDedispersion, 0);
    };


    /**
     * @class LowerFracPowerWindowGradientAfterDedispersion
     * @brief Does the fracPowerWindowGradient get lower after dedispersion?
     */
    class LowerFracPowerWindowGradientAfterDedispersion : public AnalysisCut {
    public:
      LowerFracPowerWindowGradientAfterDedispersion() : AnalysisCut() {;}
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const;
      ClassDef(LowerFracPowerWindowGradientAfterDedispersion, 0);
    };


    // /**
    //  * @class DedispersedFracPowerWindowGradientBelowThreshold
    //  * @brief An absolute cut on the fracPowerWindowGradient, rejecting events above the threshold.
    //  */
    // class DedispersedFracPowerWindowGradientBelowThreshold : public AnalysisCut {
    // public:
    //   DedispersedFracPowerWindowGradientBelowThreshold()
    // 	: AnalysisCut(), fThreshold(20)
    //   {
    // 	// fTitle = TString::Format("Dedispersed frac power window gradient > %lf", fThreshold);
    //   }
    //   virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const;
    // private:
    //   const Double_t fThreshold;
    //   ClassDef(DedispersedFracPowerWindowGradientBelowThreshold, 0);      
    // };


    /**
     * @class FisherScoreAboveThreshold
     * @brief An absolute cut on the fracPowerWindowGradient, rejecting events above the threshold.
     */
    class FisherScoreAboveThreshold : public AnalysisCut {
    public:
      FisherScoreAboveThreshold()
	: AnalysisCut(), fThreshold(7.35783)
      {
	// fTitle = TString::Format("FisherScore > %lf", fThreshold);
      }
      virtual int apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const;
    private:
      const Double_t fThreshold;
      ClassDef(FisherScoreAboveThreshold, 0);
    };

    const IsAboveHorizontal isAboveHorizontal;
    const IsTaggedAsPulser isTaggedAsPulser;
    const IsTaggedAsWaisPulser isTaggedAsWaisPulser;
    const IsTaggedAsLDBPulser isTaggedAsLdbPulser;
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
    // const DedispersedFracPowerWindowGradientBelowThreshold dedispersedFracPowerWindowGradientBelowThreshold;
    const FisherScoreAboveThreshold fisherScoreAboveThreshold;
  }
}



#endif
