#include "AnalysisCuts.h"
#include "AnitaEventSummary.h"
#include "RampdemReader.h"
#include "TMath.h"

ClassImp(Acclaim::AnalysisCuts::AnalysisCut);
ClassImp(Acclaim::AnalysisCuts::IsAboveHorizontal);
ClassImp(Acclaim::AnalysisCuts::IsTaggedAsPulser);
ClassImp(Acclaim::AnalysisCuts::IsTaggedAsWaisPulser);
ClassImp(Acclaim::AnalysisCuts::IsTaggedAsLDBPulser);
ClassImp(Acclaim::AnalysisCuts::HigherPol);
ClassImp(Acclaim::AnalysisCuts::HasSourceLocation);
ClassImp(Acclaim::AnalysisCuts::IsOnContinent);
ClassImp(Acclaim::AnalysisCuts::IsTaggedAsPayloadBlast);
ClassImp(Acclaim::AnalysisCuts::IsWithin20DegreesOfSunInPhi);
ClassImp(Acclaim::AnalysisCuts::IsGood);
ClassImp(Acclaim::AnalysisCuts::GoodGPS);
ClassImp(Acclaim::AnalysisCuts::NonZeroStokesI);
ClassImp(Acclaim::AnalysisCuts::RealSNR);
ClassImp(Acclaim::AnalysisCuts::Anita3QuietTime);
ClassImp(Acclaim::AnalysisCuts::CloseToMC);
ClassImp(Acclaim::AnalysisCuts::CloseToWais);
ClassImp(Acclaim::AnalysisCuts::IsRfTrigger);
ClassImp(Acclaim::AnalysisCuts::SmallDeltaRough);
ClassImp(Acclaim::AnalysisCuts::IsNotTaggedAsPulser);
ClassImp(Acclaim::AnalysisCuts::SignalLikeFirstStandardizedPeakMoments);
ClassImp(Acclaim::AnalysisCuts::PassesThesisCuts);
ClassImp(Acclaim::AnalysisCuts::IsNotNorth);
ClassImp(Acclaim::AnalysisCuts::HigherPeakHilbertAfterDedispersion);
ClassImp(Acclaim::AnalysisCuts::HigherImpulsivityMeasureAfterDedispersion);
ClassImp(Acclaim::AnalysisCuts::LowerFracPowerWindowGradientAfterDedispersion);
// ClassImp(Acclaim::AnalysisCuts::DedispersedFracPowerWindowGradientBelowThreshold);
ClassImp(Acclaim::AnalysisCuts::FisherScoreAboveThreshold);


static Acclaim::AnalysisCuts::Mode mode = Acclaim::AnalysisCuts::kTraining; /// Determines the behaviour of AnalysisCut::handleDefaults

void Acclaim::AnalysisCuts::setMode(Acclaim::AnalysisCuts::Mode m){
  mode = m;
}

Acclaim::AnalysisCuts::Mode Acclaim::AnalysisCuts::getMode(){
  return mode;
}



/** 
 * Function which determines which reconstructed direction is of interest.
 * 
 * Queries the AnalysisCuts::mode variable:
 *  if kTraining, then sets pol/peakInd to AnitaEventSummary::training family of functions
 *  if kAcclaimAnalysis, then called AnitaEventSummary::mostImpulsive(1) family of funtions
 * 
 * @param sum is the AnitaEventSummary on which to call the function
 * @param pol is set to the polarisation of interest
 * @param peakInd is set to the peak index of interest
 */
void Acclaim::AnalysisCuts::AnalysisCut::handleDefaults(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t& pol, Int_t& peakInd) const {
  if(peakInd==-1){
    if(mode == kTraining){
      peakInd = sum->trainingPeakInd();
      pol = sum->trainingPol();
    }
    else{
      const int metric = 1;
      peakInd = sum->mostImpulsiveInd(metric);
      pol = sum->mostImpulsivePol(metric);
    }
  }
}




/** 
 * I didn't realise isinf wasn't portable
 * 
 * @param x, is it infinite or not?
 * 
 * @return true if infinite, false otherwise
 */
bool portable_isinf(double x){
  volatile double xPlus1 = x + 1;
  return (x == xPlus1);
}



/**
 * Generic constructor, assigns name title and maximum return value
 *
 * @param name is the cut name
 * @param title is the cut title
 * @param mrv is the maximum value it is possible for apply() to return
 */
Acclaim::AnalysisCuts::AnalysisCut::AnalysisCut(int mrv)
    : fMaxRetVal(mrv)
{
  // just assign maximum return value
}



/**
 * Checks if the event was above horizontal
 *
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (default = AnitaPol::kNotAPol, see handleDefaults to see how this is handled)
 * @param peakInd is the peak index (default = -1, see handleDefaults to see how this is handled)
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::IsAboveHorizontal::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  handleDefaults(sum, pol, peakInd);
  return sum->peak[pol][peakInd].theta > 0 ? 1 : 0;
}




/**
 * Was the event tagged as any pulser (from timing info)?
 *
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (unused)
 * @param peakInd is the peak index (unused)
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::IsTaggedAsPulser::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  (void) pol;
  (void) peakInd;
  return isTaggedAsLdbPulser.apply(sum) || isTaggedAsWaisPulser.apply(sum);
}

/**
 * Was the event tagged as a WAIS pulser (from timing info)?
 *
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (unused)
 * @param peakInd is the peak index (unused)
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::IsTaggedAsWaisPulser::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  (void) pol;
  (void) peakInd;
  return sum->flags.pulser == AnitaEventSummary::EventFlags::WAIS;
}


/**
 * Was the event tagged as a LDB pulser (from timing info)
 *
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (unused)
 * @param peakInd is the peak index (unused)
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::IsTaggedAsLDBPulser::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  (void) pol;
  (void) peakInd;
  if(AnitaVersion::get()==3){
    if(sum->run >= 200){return false;}
  }
  return sum->flags.pulser == AnitaEventSummary::EventFlags::LDB;
}


/**
 * Get the highest polarisation, wraps the highestPeakPol() function
 *
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (unused)
 * @param peakInd is the peak index (unused)
 *
 * @return 0 for HPol, 1 for VPol
 */
int Acclaim::AnalysisCuts::HigherPol::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  (void) pol;
  (void) peakInd;
  return sum->highestPolAsInt();
}

/**
 * Does the event reconstruct to somewhere on the Earth (land or sea)?
 *
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (default = AnitaPol::kNotAPol, see handleDefaults to see how this is handled)
 * @param peakInd is the peak index (default = -1, see handleDefaults to see how this is handled)
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::HasSourceLocation::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  handleDefaults(sum, pol, peakInd);  
  bool didReconstruct = (sum->peak[pol][peakInd].latitude < -900 || TMath::Abs(sum->peak[pol][peakInd].theta_adjustment_needed) > 0) ? false : true;
  // if(!didReconstruct){
  //   std::cerr << sum->highestPeak().latitude << "\t" << sum->highestPeak().theta_adjustment_needed << std::endl;
  // }
  return didReconstruct;
}

/**
 * Does the event reconstruct to somewhere on the Antarctic land mass?
 * Wraps RampdemReader::isOnContinent(lon, lat)
 *
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (default = AnitaPol::kNotAPol, see handleDefaults to see how this is handled)
 * @param peakInd is the peak index (default = -1, see handleDefaults to see how this is handled)
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::IsOnContinent::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  handleDefaults(sum, pol, peakInd);
  return RampdemReader::isOnContinent(sum->peak[pol][peakInd].longitude, sum->peak[pol][peakInd].latitude);
}



/**
 * Did the reconstruction tag this event as a payload blast?
 *
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (unused)
 * @param peakInd is the peak index (unused)
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::IsTaggedAsPayloadBlast::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  (void) pol;
  (void) peakInd;
  return sum->flags.isPayloadBlast != 0 ? true : false;
}


/**
 * Did the reconstruction tag this event as a payload blast?
 *
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (default = AnitaPol::kNotAPol, see handleDefaults to see how this is handled)
 * @param peakInd is the peak index (default = -1, see handleDefaults to see how this is handled)
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::IsWithin20DegreesOfSunInPhi::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  handleDefaults(sum, pol, peakInd);
  return TMath::Abs(sum->peak[pol][peakInd].dPhiSun()) < 20 ? true : false;
}


/**
 * Did this event pass all quality cuts?
 *
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (unused)
 * @param peakInd is the peak index (unused)
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::IsGood::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  (void) pol;
  (void) peakInd;
  // here we handle the fact that the QualityCut, which sets this flag during reconstruction,
  // also check the clock channel, which is empty for MC generated events.
  // The num points quality cut is encoded as isVarner2, so we skip that check for the MC case.
  bool isGood = sum->mc.weight > 0 ? (sum->flags.isVarner == 0 && sum->flags.isPayloadBlast == 0) : sum->flags.isGood == 1;
  return isGood;
}


/**
 * Bad GPS will give NaN...
 *
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (unused)
 * @param peakInd is the peak index (unused)
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::GoodGPS::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  (void) pol;
  (void) peakInd;
  return !(TMath::IsNaN(sum->anitaLocation.heading) || portable_isinf(sum->anitaLocation.heading));
}


/**
 * Zero Stokes-I gives NaN linear/circular pol frac
 * This is probably a reconstruction bug I need to handle
 *
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (default = AnitaPol::kNotAPol, see handleDefaults to see how this is handled)
 * @param peakInd is the peak index (default = -1, see handleDefaults to see how this is handled)
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::NonZeroStokesI::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  handleDefaults(sum, pol, peakInd);  
  return (sum->deconvolved[pol][peakInd].I > 0);
}



/**
 * The first few events in the decimated sample will have NaN SNR
 * Because there's no numbers in the noise tree...
 * 
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (default = AnitaPol::kNotAPol, see handleDefaults to see how this is handled)
 * @param peakInd is the peak index (default = -1, see handleDefaults to see how this is handled)=
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::RealSNR::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  handleDefaults(sum, pol, peakInd);  
  return !(TMath::IsNaN(sum->deconvolved[pol][peakInd].snr) || portable_isinf(sum->deconvolved[pol][peakInd].snr));
}


/**
 * The first few events in the decimated sample will have NaN SNR
 * Because there's no numbers in the noise tree...
 * 
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (unused)
 * @param peakInd is the peak index (unused)
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::Anita3QuietTime::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  (void) pol;
  (void) peakInd;
  return (sum->realTime >= 1419320000 && sum->realTime < 1420250000);
}


/**
 * If an event is data, not MC event (weight >= 0), returns false
 * 
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (default = AnitaPol::kNotAPol, see handleDefaults to see how this is handled)
 * @param peakInd is the peak index (default = -1, see handleDefaults to see how this is handled)
 *
 * @return 1+peakInd if true, 0 if false
 */
int Acclaim::AnalysisCuts::CloseToMC::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  handleDefaults(sum, pol, peakInd);
  bool close = false;
  if(sum->mc.weight > 0){ // for non-MC events, return false
    close = false;
    // kinda arbitrary here... but should get use close enough
    const double dPhiClose = 5;
    const double dThetaClose = 3;
    
    double dPhi = sum->peak[pol][peakInd].dPhiMC();
    if(TMath::Abs(dPhi) < dPhiClose){
      double dTheta = sum->peak[pol][peakInd].dThetaMC();
      if(TMath::Abs(dTheta) < dThetaClose){
        close = true;
      }
    }
  }
  return close;
}




/**
 * If an event is data, not MC event (weight >= 0), returns false
 * 
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (default = AnitaPol::kNotAPol, see handleDefaults to see how this is handled)
 * @param peakInd is the peak index (default = -1, see handleDefaults to see how this is handled)
 *
 * @return 1+peakInd if true, 0 if false
 */
int Acclaim::AnalysisCuts::CloseToWais::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  handleDefaults(sum, pol, peakInd);
  bool close = false;
  if(sum->mc.weight > 0){ // for non-MC events, return false
    close = false;
    // kinda arbitrary here... but should get use close enough
    const double dPhiClose = 5;
    const double dThetaClose = 3;
    
    double dPhi = sum->peak[pol][peakInd].dPhiWais();
    if(TMath::Abs(dPhi) < dPhiClose){
      double dTheta = sum->peak[pol][peakInd].dThetaWais();
      if(TMath::Abs(dTheta) < dThetaClose){
        close = true;
      }
    }
  }
  return close;
}


/**
 * Was this an RF trigger?
 * 
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (default = AnitaPol::kNotAPol, see handleDefaults to see how this is handled)
 * @param peakInd is the peak index (default = -1, see handleDefaults to see how this is handled)
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::IsRfTrigger::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  (void) pol;
  (void) peakInd;
  return (sum->flags.isRF != 0);
}


/**
 * Was this an RF trigger?
 * 
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (unused)
 * @param peakInd is the peak index (unused)
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::SmallDeltaRough::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  handleDefaults(sum, pol, peakInd);
  const double dAngleCut = 4;
  return TMath::Abs(sum->peak[pol][peakInd].dphi_rough) < dAngleCut && TMath::Abs(sum->peak[pol][peakInd].dtheta_rough) < dAngleCut;
}


/**
 * Was this not tagged as a pulser?
 * 
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (unused)
 * @param peakInd is the peak index (unused)
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::IsNotTaggedAsPulser::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  (void) pol;
  (void) peakInd;
  return sum->flags.pulser == 0;
}



/**
 * Do the standardized peak moments match MC expectations?
 * 
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (default = AnitaPol::kNotAPol, see handleDefaults to see how this is handled)
 * @param peakInd is the peak index (default = -1, see handleDefaults to see how this is handled)
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::SignalLikeFirstStandardizedPeakMoments::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  handleDefaults(sum, pol, peakInd);
  return (sum->coherent_filtered[pol][peakInd].standardizedPeakMoment(1) > 60 && sum->deconvolved_filtered[pol][peakInd].standardizedPeakMoment(1) > 60);
}






/** 
 * Applies my thesis cuts in sequence, warts and all.
 * Hopefully this analysis will do better...
 * 
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (unused)
 * @param peakInd is the peak index (unused)
 *
 * @return 1 if true, 0 if false
 * 
 * @return 
 */
int Acclaim::AnalysisCuts::PassesThesisCuts::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const {
  (void) pol;
  (void) peakInd;
  static const int numFisherWeights = 3;
  static const Double_t fisherWeights[numFisherWeights] = {-2.80993, 14.6149, 0.0107283};
  static const Double_t fisherCutVal = -0.0270745;

  Double_t fisherScore = fisherWeights[0] + fisherWeights[1]*sum->trainingPeak().value + fisherWeights[2]*sum->trainingPeak().value;
  return fisherScore > fisherCutVal;
  

}


/** 
 * Applies my thesis cuts in sequence, warts and all.
 * Hopefully this analysis will do better...
 * 
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (default = AnitaPol::kNotAPol, see handleDefaults to see how this is handled)
 * @param peakInd is the peak index (default = -1, see handleDefaults to see how this is handled)
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::IsNotNorth::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const {
  handleDefaults(sum, pol, peakInd);  
  return (TMath::Abs(sum->peak[pol][peakInd].dPhiNorth()) > 90);

}



/** 
 * Checks to see if the hilbert peak of the coherently summed waveform is higher than it's dedispersed counterpart
 * 
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (default = AnitaPol::kNotAPol, see handleDefaults to see how this is handled)
 * @param peakInd is the peak index (default = -1, see handleDefaults to see how this is handled)
 * 
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::HigherPeakHilbertAfterDedispersion::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const {
  handleDefaults(sum, pol, peakInd);
  // std::cerr << pol << "\t" << peakInd << "\t" << sum->deconvolved_filtered[pol][peakInd].peakHilbert << "\t" << sum->coherent_filtered[pol][peakInd].peakHilbert << std::endl;
  return (sum->deconvolved_filtered[pol][peakInd].peakHilbert > sum->coherent_filtered[pol][peakInd].peakHilbert);
}



/** 
 * Checks to see if the impulsivity measure of the coherently summed waveform is higher than it's dedispersed counterpart
 * 
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (default = AnitaPol::kNotAPol, see handleDefaults to see how this is handled)
 * @param peakInd is the peak index (default = -1, see handleDefaults to see how this is handled)
 * 
 * @return 1 if true, 0 if false
 */
int Acclaim::AnalysisCuts::HigherImpulsivityMeasureAfterDedispersion::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const {
  handleDefaults(sum, pol, peakInd);
  return (sum->deconvolved_filtered[pol][peakInd].impulsivityMeasure > sum->coherent_filtered[pol][peakInd].impulsivityMeasure);
}



/** 
 * Checks to see if the fracPowerWindowGradient of the coherently summed waveform is lower than it's dedispersed counterpart
 * 
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (default = AnitaPol::kNotAPol, see handleDefaults to see how this is handled)
 * @param peakInd is the peak index (default = -1, see handleDefaults to see how this is handled)
 * 
 * @return 1 if true, 0 if false
 */

int Acclaim::AnalysisCuts::LowerFracPowerWindowGradientAfterDedispersion::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const {  
  handleDefaults(sum, pol, peakInd);
  return (sum->deconvolved_filtered[pol][peakInd].fracPowerWindowGradient() < sum->coherent_filtered[pol][peakInd].fracPowerWindowGradient());
}




// /** 
//  * Checks to see if the dedispersed fracPowerWindowGradient is below threshold
//  * Threshold value defined in AnalysisSettings.conf
//  * 
//  * @param sum is the AnitaEventSummary
//  * @param pol is the polarisation (default = AnitaPol::kNotAPol, see handleDefaults to see how this is handled)
//  * @param peakInd is the peak index (default = -1, see handleDefaults to see how this is handled)
//  * 
//  * @return 1 if true, 0 if false
//  */

// int Acclaim::AnalysisCuts::DedispersedFracPowerWindowGradientBelowThreshold::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const {  
//   handleDefaults(sum, pol, peakInd);
//   // std::cout << sum->deconvolved_filtered[pol][peakInd].fracPowerWindowGradient() << "\t" << fThreshold << std::endl;
//   return (sum->deconvolved_filtered[pol][peakInd].fracPowerWindowGradient() < fThreshold);
// }






/** 
 * Checks to see if the dedispersed fracPowerWindowGradient is below threshold
 * Threshold value defined in AnalysisSettings.conf
 * 
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (default = AnitaPol::kNotAPol, see handleDefaults to see how this is handled)
 * @param peakInd is the peak index (default = -1, see handleDefaults to see how this is handled)
 * 
 * @return 1 if true, 0 if false
 */

int Acclaim::AnalysisCuts::FisherScoreAboveThreshold::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const {  
  handleDefaults(sum, pol, peakInd);

  // zeroth pass
  // 10.534119+(0.005302*Abs_trainingPeak_dPhiSun)+(-0.011572*Abs_trainingPeak_minAbsHwAngle)+(-0.224413*trainingDeconvolvedFiltered_fracPowerWindowGradient)
  double fs = 10.534119;
  fs += 0.005302*sum->peak[pol][peakInd].dPhiSun();
  fs += -0.011572*TMath::Abs(sum->peak[pol][peakInd].minAbsHwAngle());
  fs += -0.224413*sum->deconvolved_filtered[pol][peakInd].fracPowerWindowGradient();

  return (fs > fThreshold);
}



