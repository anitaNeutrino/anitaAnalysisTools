/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             A global set of TStrings/TCuts for drawing in root.
             You'll need to #include this file on the root prompt to use them.
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
   * @namespace SumTree
   * @brief Should work with Draw/Scan on a tree of AnitaEventSummaries
   */
  namespace SumTree {

    const TString weight = "((mc.weight > 0)*mc.weight + 1*(mc.weight==0))"; ///Note: This sets the weight for data = 1
    const TString pol = TString::Format("floor(Iteration$/%d)", AnitaEventSummary::maxDirectionsPerPol); /// Converts the implied Draw loop iteration to polarisation
    const TString peakInd = TString::Format("(Iteration$ %% %d)", AnitaEventSummary::maxDirectionsPerPol); /// Converts the implied Draw loop iteration to peak index
    const TString deconvolved_filtered_fracPowerWindowGradient = TString::Format("10.0*(%s %s %s %s)",
										 "-0.2*(deconvolved_filtered[][].fracPowerWindowEnds[0] - deconvolved_filtered[][].fracPowerWindowBegins[0])",
										 "-0.1*(deconvolved_filtered[][].fracPowerWindowEnds[1] - deconvolved_filtered[][].fracPowerWindowBegins[1])",
										 "+0.1*(deconvolved_filtered[][].fracPowerWindowEnds[3] - deconvolved_filtered[][].fracPowerWindowBegins[3])",
										 "+0.2*(deconvolved_filtered[][].fracPowerWindowEnds[4] - deconvolved_filtered[][].fracPowerWindowBegins[4])");
    const TString coherent_filtered_fracPowerWindowGradient = TString::Format("10.0*(%s %s %s %s)",
									      "-0.2*(coherent_filtered[][].fracPowerWindowEnds[0] - coherent_filtered[][].fracPowerWindowBegins[0])",
									      "-0.1*(coherent_filtered[][].fracPowerWindowEnds[1] - coherent_filtered[][].fracPowerWindowBegins[1])",
									      "+0.1*(coherent_filtered[][].fracPowerWindowEnds[3] - coherent_filtered[][].fracPowerWindowBegins[3])",
									      "+0.2*(coherent_filtered[][].fracPowerWindowEnds[4] - coherent_filtered[][].fracPowerWindowBegins[4])");
    const TString hilbertPeakTimeShift = "coherent_filtered[][].peakTime - deconvolved_filtered[][].peakTime";
    const TString minAbsHwAngle = TString::Format("(%s + %s)",
						  "(TMath::Abs(peak[][].hwAngle) < TMath::Abs(peak[][].hwAngleXPol))*TMath::Abs(peak[][].hwAngle)",
						  "(TMath::Abs(peak[][].hwAngle) >= TMath::Abs(peak[][].hwAngleXPol))*TMath::Abs(peak[][].hwAngleXPol)");
    const TCut isRfTrigger("isRfTrigger", "flags.isRF > 0"); /// (N=1)
    const TCut isNotTaggedAsPulser("isNotTaggedAsPulser", "flags.pulser == 0"); /// (N=1)
    const TCut anita3QuietTime("anita3QuietTime", "realTime >= 1419320000 && realTime < 1420250000"); /// (N=1)
    // const TCut realSNR("realSNR", "(!TMath::IsNaN(deconvolved.snr) && TMath::Finite(deconvolved.snr))");
    const TCut goodGPS("goodGPS", "(!TMath::IsNaN(anitaLocation.heading) && TMath::Finite(anitaLocation.heading))"); /// (N=1)
    const TCut isGood("isGood", "((mc.weight > 0 && flags.isVarner == 0 && flags.isPayloadBlast == 0) || (mc.weight==0 && flags.isGood == 1))"); /// (N=1)
    const TCut isTaggedAsPayloadBlast("isTaggedAsPayloadBlast", "flags.isPayloadBlast > 0"); /// (N=1)
    const TCut isTaggedAsLDBPulser("isTaggedAsLDBPulser", TString::Format("!(AnitaVersion::get()==3 && run >=200) && flags.pulser==%d", AnitaEventSummary::EventFlags::LDB).Data()); /// (N=1)
    const TCut isTaggedAsWaisPulser("isTaggedAsWaisPulser", TString::Format("flags.pulser==%d", AnitaEventSummary::EventFlags::WAIS).Data()); /// (N=1)
    const TCut isTaggedAsPulser("isTaggedAsPulser", TString::Format("(%s) || (%s)", isTaggedAsLDBPulser.GetTitle(), isTaggedAsWaisPulser.GetTitle())); /// (N=1)
    const TCut npbc0A("npbc0A", TString::Format("flags.middleOrBottomPower[0] < %lf * flags.topPower[0] + %lf", 30./7, 0.06)); /// (N=1)
    const TCut npbc0B("npbc0B", TString::Format("flags.middleOrBottomPower[0] > %lf * (flags.topPower[0] - %lf)", 7./30, 0.06)); /// (N=1)
    const TCut npbc1("npbc1", TString::Format("flags.middleOrBottomPower[1] < %lf * flags.topPower[1]", 1./0.28)); /// (N=1)
    const TCut npbc2("npbc2", TString::Format("Max$(flags.maxBottomToTopRatio) < %lf", 3.0)); /// (N=1)


    const TCut highestPeak("highestPeak", "Max$(peak[][].value)==peak[][].value");
    const TCut mostImpulsivePeak("mostImpulsivePeak", "Max$(deconvolved_filtered[][].impulsivityMeasure)==deconvolved_filtered[][].impulsivityMeasure");

    const double maxDeltaPhiRough = 4.0;    
    const double maxDeltaThetaRough = 4.0;
    const TCut smallDeltaRough("smallDeltaRough",
			       TString::Format("TMath::Abs(peak[][].dphi_rough) < %lf && TMath::Abs(peak[][].dtheta_rough) < %lf",
					       maxDeltaPhiRough, maxDeltaThetaRough));
    const TCut higherHilbertPeakAfterDedispersion("higherHilbertPeakAfterDedispersion",
						  "deconvolved_filtered[][].peakHilbert > coherent_filtered[][].peakHilbert"); 
    const TCut higherImpulsivityMeasureAfterDedispersion("higherImpulsivityMeasureAfterDedispersion",
							 "deconvolved_filtered[][].impulsivityMeasure > coherent_filtered[][].impulsivityMeasure"); 
    const TCut lowerFracPowerWindowGradientAfterDedispersion("lowerFracPowerWindowGradientAfterDedispersion", 
							     TString::Format("(%s) < (%s)",
									     deconvolved_filtered_fracPowerWindowGradient.Data(),
									     coherent_filtered_fracPowerWindowGradient.Data()));
    const TCut reasonableHilbertPeakTimeShiftAfterDedispersion("reasonableHilbertPeakTimeShiftAfterDedispersion",
							       TString::Format("%s > %lf && %s < %lf",
									       hilbertPeakTimeShift.Data(), 0.0,
									       hilbertPeakTimeShift.Data(), 20.0)); 
    const TCut isOnContinent("isOnContinent",
			     "RampdemReader::isOnContinent(peak[][].longitude, peak[][].latitude)"); /// As yet untested
    
    const TCut hasSourceLocation("hasSourceLocation",
				 "(peak[][].latitude < -900 || TMath::Abs(peak[][].theta_adjustment_needed) > 0"); 
    const TCut isAboveHorizontal("isAboveHorizontal",
				 "peak[][].theta > 0"); 
    const TCut isBelowHorizontal("isBelowHorizontal",
				 "peak[][].theta < 0"); 
    const TCut npbc3("npbc3", TString::Format("%lf*deconvolved_filtered[][].peakHilbert > (1+flags.maxBottomToTopRatio[Iteration$/5])*flags.minBottomToTopRatio[Iteration$/5] - %lf"
					      ,14.0, 1000.0));
    const TCut isGood2("isGood2", TString::Format("(%s && %s && %s && %s && %s)", npbc0A.GetTitle(), npbc0B.GetTitle(), npbc1.GetTitle(), npbc2.GetTitle(), npbc3.GetTitle()));

  }
}

#endif
