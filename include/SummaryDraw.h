/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             A global set of TStrings/TCuts for drawing on the command line
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
   * @namespace Draw
   * @brief Commonly used analysis variables in SummarySet::Draw() or TTree::Draw().
   *
   * You'll need to include Cuts on the command line.
   */
  namespace Draw {
    const TString dPhiWais = "FFTtools::wrap(peak.phi - wais.phi, 360, 0)";
    const TString dThetaWais = "(peak.theta - wais.theta)";
    const TString dPhiMC = "FFTtools::wrap(peak.phi - mc.phi, 360, 0)";
    const TString dThetaMC = "(peak.theta - mc.theta)";
  }


  /**
   * @namespace Cuts
   * @brief Contains a set of TCuts which implement cuts on AnitaEventSummary objects
   * 
   * For arrays like the peak[][] or coherent_filtered[][] things, doing TTree::Draw("peak.value") 
   * will fill once per array entry, so a selection variable must be a cut e.g. highestPeak TCut
   * which selects one element by requiring it have the highest peak value.
   * 
   */
  namespace Cuts
  {

    const double maxDeltaPhi = 5.5;
    const double maxDeltaTheta = 3.5;

    // Replacement for the member functions in AnitaEventSummary...
    const TCut highestPeak = "Max$(peak.value)==peak.value";    
    const TCut isRfTrigger("isRfTrigger", "flags.isRF > 0");
    const TCut smallDeltaRough("smallDeltaRough", "TMath::Abs(peak.dphi_rough) < 4 && TMath::Abs(peak.dtheta_rough) < 4");
    const TCut isNotTaggedAsPulser("isNotTaggedAsPulser", "flags.pulser == 0");
    const TCut isNotNorth("isNotNorth", "TMath::Abs(peak.dPhiNorth()) > 90");
    const TCut acceptableHardwareAngle("acceptableHardwareAngle", "TMath::Abs(sum.minAbsHwAngle()) < 65.0");
    const TCut higherHilbertPeakAfterDedispersion("higherHilbertPeakAfterDedispersion", "deconvolved_filtered.peakHilbert > coherent_filtered.peakHilbert");
    const TCut higherImpulsivityMeasureAfterDedispersion("higherImpulsivityMeasureAfterDedispersion", "deconvolved_filtered.impulsivityMeasure > coherent_filtered.impulsivityMeasure");
    const TCut lowerFracPowerWindowGradientAfterDedispersion("lowerFracPowerWindowGradientAfterDedispersion", "deconvolved_filtered.fracPowerWindowGradient() < coherent_filtered.fracPowerWindowGradient()");
    const TCut closeToWais("closeToWais", TString::Format("(mc.weight == 0 && %s < %lf && %s > %lf && %s < %lf)", Draw::dPhiWais.Data(), maxDeltaPhi, Draw::dPhiWais.Data(), -maxDeltaPhi, Draw::dThetaWais.Data(), maxDeltaTheta)); /// false for MC
    const TCut closeToMC("closeToMC",     TString::Format("(mc.weight >  0 && %s < %lf && %s > %lf && %s < %lf)", Draw::dPhiMC.Data(),   maxDeltaPhi, Draw::dPhiWais.Data(), -maxDeltaPhi, Draw::dThetaMC.Data(),   maxDeltaTheta)); /// false for non-MC
    const TCut anita3QuietTime("anita3QuietTime", "realTime >= 1419320000 && sum->realTime < 1420250000");
    const TCut realSNR("realSNR", "(!TMath::IsNaN(deconvolved.snr) && TMath::Finite(deconvolved.snr))");
    const TCut goodGPS("goodGPS", "(!TMath::IsNaN(anitaLocation.heading) && TMath::Finite(anitaLocation.heading))");
    const TCut isGood("isGood", "((mc.weight > 0 && sum->flags.isVarner == 0 && sum->flags.isPayloadBlast == 0) || (mc.weight==0 && flags.isGood == 1))");
    const TCut isWithin20DegreesOfSunInPhi("isWithin20DegreesOfSunInPhi", "sum->peak.dPhiSun()) < 20");
    const TCut isTaggedAsPayloadBlast("isTaggedAsPayloadBlast", "flags.isPayloadBlast > 0");
    const TCut isOnContinent("isOnContinent", "RampdemReader::isOnContinent(peak.longitude, peak.latitude)");
    const TCut hasSourceLocation("hasSourceLocation", "(peak.latitude < -900 || TMath::Abs(peak.theta_adjustment_needed) > 0");
    const TCut isTaggedAsLDBPulser("isTaggedAsLDBPulser", TString::Format("!(AnitaVersion::get()==3 && run >=200) && flags.pulser==%d", AnitaEventSummary::EventFlags::LDB).Data()); /// @todo AnitaVersion::get()==3 ?
    const TCut isTaggedAsWaisPulser("isTaggedAsWaisPulser", TString::Format("flags.pulser==%d", AnitaEventSummary::EventFlags::WAIS).Data());
    const TCut isTaggedAsPulser("isTaggedAsPulser", TString::Format("(%s) || (%s)", isTaggedAsLDBPulser.GetTitle(), isTaggedAsWaisPulser.GetTitle()));
    const TCut isAboveHorizontal("isAboveHorizontal", "peak.theta > 0");

    const TCut npbc0A("npbc0A", TString::Format("flags.middleOrBottomPower[0] < %lf * flags.topPower[0] + %lf", 30./7, 0.06));
    const TCut npbc0B("npbc0B", TString::Format("flags.middleOrBottomPower[0] > %lf * (flags.topPower[0] - %lf)", 7./30, 0.06));
    const TCut npbc1("npbc1", TString::Format("flags.middleOrBottomPower[1] < %lf * flags.topPower[1]", 1./0.28));
  }



}

#endif
