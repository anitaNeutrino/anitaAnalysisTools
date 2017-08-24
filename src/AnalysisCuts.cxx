#include "AnalysisCuts.h"
#include "AnitaEventSummary.h"
#include "RampdemReader.h"
#include "TMath.h"

/**
 * Generic constructor, assigns name title and maximum return value
 *
 * @param name is the cut name
 * @param title is the cut title
 * @param mrv is the maximum value it is possible for apply() to return
 */
Acclaim::AnalysisCut::AnalysisCut(const char* name, const char* title, int mrv)
    : fName(name), fTitle(title), fMaxRetVal(mrv)
{
  // just assign name, title, and maximum return value
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
int Acclaim::IsAboveHorizontal::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  handleDefaults(sum, pol, peakInd);
  return sum->peak[pol][peakInd].theta > 0 ? 1 : 0;
}




/**
 * Was the event tagged as a WAIS pulser (from timing info)
 *
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (unused)
 * @param peakInd is the peak index (unused)
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::IsTaggedAsWaisPulser::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
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
int Acclaim::IsTaggedAsLDBPulser::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  (void) pol;
  (void) peakInd;  
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
int Acclaim::HigherPol::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
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
int Acclaim::HasSourceLocation::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
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
int Acclaim::IsOnContinent::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
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
int Acclaim::IsTaggedAsPayloadBlast::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
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
int Acclaim::IsWithin20DegreesOfSunInPhi::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
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
int Acclaim::IsGood::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  (void) pol;
  (void) peakInd;
  return sum->flags.isGood != 0 ? true : false;
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
int Acclaim::GoodGPS::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  (void) pol;
  (void) peakInd;  
  return !(TMath::IsNaN(sum->anitaLocation.heading) || isinf(sum->anitaLocation.heading));
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
int Acclaim::NonZeroStokesI::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
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
int Acclaim::RealSNR::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  handleDefaults(sum, pol, peakInd);  
  return !(TMath::IsNaN(sum->deconvolved[pol][peakInd].snr) || isinf(sum->deconvolved[pol][peakInd].snr));
}


/**
 * The first few events in the decimated sample will have NaN SNR
 * Because there's no numbers in the noise tree...
 * 
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (default = AnitaPol::kNotAPol, see handleDefaults to see how this is handled)
 * @param peakInd is the peak index (default = -1, see handleDefaults to see how this is handled)
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::Anita3QuietTime::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  (void) pol;
  (void) peakInd;
  return (sum->realTime >= 1419320000 && sum->realTime < 1420250000);
}


/**
 * If an event is an MC event (weight >= 0)
 * 
 * @param sum is the AnitaEventSummary
 * @param pol is the polarisation (default = AnitaPol::kNotAPol, see handleDefaults to see how this is handled)
 * @param peakInd is the peak index (default = -1, see handleDefaults to see how this is handled)
 *
 * @return 1+peakInd if true, 0 if false
 */
int Acclaim::CloseToMC::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const
{
  handleDefaults(sum, pol, peakInd);
  if(sum->mc.weight <= 0){
    return true;
  }
  else{
    // kinda arbitrary here... but should get use close enough
    // I suppose I should record a signal sample
    const double dPhiClose = 5;
    const double dThetaClose = 3;
    AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;
    for(int peakInd=0; peakInd < sum->nPeaks[pol]; peakInd++){
      double dPhi = sum->peak[pol][peakInd].dPhiMC();
      if(TMath::Abs(dPhi) < dPhiClose){
        double dTheta = sum->peak[pol][peakInd].dThetaMC();
        if(TMath::Abs(dTheta) < dThetaClose){
          return peakInd+1;
        }
      }
    }
  }
  return 0;
}
