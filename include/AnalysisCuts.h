/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             A global set of TCuts
***********************************************************************************************************/

#ifndef ANALYSIS_CUTS_H
#define ANALYSIS_CUTS_H

#include "AnitaConventions.h"
#include "AnitaEventSummary.h"
#include "TString.h"
#include "TCut.h"

namespace Acclaim
{

  /**
   * @namespace AnalysisCuts
   * @brief Contains a set of TCuts which implement cuts on AnitaEventSummary objects
   * 
   * For arrays like the peak[][] or coherent_filtered[][] things, doing TTree::Draw("peak.value") 
   * will fill once per array entry, so a selection variable must be a cut e.g. highestPeak TCut
   * which selects one element by requiring it have the highest peak value.
   * 
   */
  namespace AnalysisCuts
  {

    // Replacement for the member functions in AnitaEventSummary...
    const TCut highestPeak = "Max$(peak.value)==peak.value";
    
    
    const TCut isRfTrigger = "flags.isRF > 0";
    const TCut smallDeltaRough = "TMath::Abs(peak.dphi_rough) < 4 && TMath::Abs(peak.dtheta_rough) < 4";
    const TCut isNotTaggedAsPulser = "flags.pulser == 0";
    const TCut isNotNorth = "TMath::Abs(peak.dPhiNorth()) > 90";
    const TCut acceptableHardwareAngle = "TMath::Abs(sum.minAbsHwAngle()) < 65.0";
    const TCut higherHilbertPeakAfterDedispersion = "deconvolved_filtered.peakHilbert > coherent_filtered.peakHilbert";
    const TCut higherImpulsivityMeasureAfterDedispersion = "deconvolved_filtered.impulsivityMeasure > coherent_filtered.impulsivityMeasure";    
    const TCut lowerFracPowerWindowGradientAfterDedispersion = "deconvolved_filtered.fracPowerWindowGradient() < coherent_filtered.fracPowerWindowGradient()";
    const TCut closeToWais = "mc.weight == 0 && peak.dPhiWais() < 5.5 && peak.dThetaWais() < 3.5"; /// always false for MC
    const TCut closeToMC = "mc.weight > 0 && peak.dPhiMC() < 5.5 && peak.dThetaMC() < 3.5"; /// always false for MC
    const TCut anita3QuietTime = "realTime >= 1419320000 && sum->realTime < 1420250000";
    const TCut realSNR = "(!TMath::IsNaN(deconvolved.snr) && TMath::Finite(deconvolved.snr));";
    const TCut goodGPS = "(!TMath::IsNaN(anitaLocation.heading) && TMath::Finite(anitaLocation.heading))";
    const TCut isGood = "((mc.weight > 0 && sum->flags.isVarner == 0 && sum->flags.isPayloadBlast == 0) || (mc.weight==0 && flags.isGood == 1))";
    const TCut isWithin20DegreesOfSunInPhi = "sum->peak.dPhiSun()) < 20";
    const TCut isTaggedAsPayloadBlast = "flags.isPayloadBlast > 0;";
    const TCut isOnContinent = "RampdemReader::isOnContinent(peak.longitude, peak.latitude)"; ///  untested and probably won't work...
    const TCut hasSourceLocation = "(peak.latitude < -900 || TMath::Abs(peak.theta_adjustment_needed) > 0";
    // if(AnitaVersion::get()==3){
    //   if(sum->run >= 200){return false;}
    // }
    // return sum->flags.pulser == AnitaEventSummary::EventFlags::LDB;
    const TCut isTaggedAsLDBPulser = TString::Format("flags.pulser==%d", AnitaEventSummary::EventFlags::LDB).Data(); /// @todo AnitaVersion::get()==3 ?
    const TCut isTaggedAsWaisPulser = TString::Format("flags.pulser==%d", AnitaEventSummary::EventFlags::WAIS).Data(); /// @todo AnitaVersion::get()==3 ?
    const TCut isTaggedAsPulser = isTaggedAsLDBPulser || isTaggedAsWaisPulser;
    const TCut isAboveHorizontal = "peak.theta > 0";
  }

}

#endif
