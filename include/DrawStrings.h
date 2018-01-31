/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             A global set of TStrings/TCuts for drawing in root.
             You'll need to #include this file on the root prompt to use them.
***********************************************************************************************************/

#ifndef ACCLAIM_DRAW_STRINGS_H
#define ACCLAIM_DRAW_STRINGS_H

#include "AnitaConventions.h"
#include "AnitaEventSummary.h"
#include "FFTtools.h"
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
    const TCut newPayloadBlastCutPart0A("newPayloadBlastCutPart0A", TString::Format("flags.middleOrBottomPower[0] < %lf * flags.topPower[0] + %lf", 30./7, 0.06)); /// (N=1)
    const TCut newPayloadBlastCutPart0B("newPayloadBlastCutPart0B", TString::Format("flags.middleOrBottomPower[0] > %lf * (flags.topPower[0] - %lf)", 7./30, 0.06)); /// (N=1)
    const TCut newPayloadBlastCutPart1("newPayloadBlastCutPart1", TString::Format("flags.middleOrBottomPower[1] < %lf * flags.topPower[1]", 1./0.28)); /// (N=1)
    const TCut newPayloadBlastCutPart2("newPayloadBlastCutPart2", TString::Format("Max$(flags.maxBottomToTopRatio) < %lf", 2.67)); /// (N=1)

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
    const TCut newPayloadBlastCutPart3("newPayloadBlastCutPart3", TString::Format("%lf*deconvolved_filtered[][].peakHilbert > (1+flags.maxBottomToTopRatio[Iteration$/5])*flags.minBottomToTopRatio[Iteration$/5] - %lf"
					      ,14.0, 1000.0));
    const TCut isGood2("isGood2", TString::Format("(%s && %s && %s && %s && %s)", newPayloadBlastCutPart0A.GetTitle(), newPayloadBlastCutPart0B.GetTitle(), newPayloadBlastCutPart1.GetTitle(), newPayloadBlastCutPart2.GetTitle(), newPayloadBlastCutPart3.GetTitle()));

  }

  /**
   * @namespace ThermalTree
   * @brief Strings/cuts for use in the thermal trees (a reduction of the summary trees with variables slightly renamed)
   */
  namespace ThermalTree {

    const TString dPhiMC = "FFTtools::wrap(peak_phi - mc_phi, 360, 0)";
    const TString dThetaMC = "peak_theta + mc_theta";/// MC theta has down is +ve, Acclaim reco has down is -ve.
    const TString dPhiWais = "FFTtools::wrap(peak_phi - wais_phi, 360, 0)";
    const TString dThetaWais = "peak_theta + wais_theta";/// EventCorrelator theta has down is +ve, Acclaim reco has down is -ve.
    const TString dPhiHiCal = "FFTtools::wrap(peak_phi - hiCalPhi, 360, 0)";
    const TString dThetaHiCal = "peak_theta + hiCalTheta";/// EventCorrelator theta has down is +ve, Acclaim reco has down is -ve.

    const TCut weight(const TCut cut = "", double multiplier=1);
    
    const double dPhiClose = 5.5;
    const double dThetaClose = 3.5;
    const TCut closeToMC("closeToMC", TString::Format("TMath::Abs(%s) < %lf && TMath::Abs(%s) < %lf",
						      dPhiMC.Data(),   dPhiClose,
						      dThetaMC.Data(), dThetaClose));    
    const TCut closeToWais("closeToWais", TString::Format("TMath::Abs(%s) < %lf && TMath::Abs(%s) < %lf",
							  dPhiWais.Data(),   dPhiClose,
							  dThetaWais.Data(), dThetaClose));

    const TCut closeToHiCal("closeToHiCal", TString::Format("(duringHiCal > 0 && TMath::Abs(%s) < 5)",
							    dPhiHiCal.Data()));

    const TCut smallDeltaRough("smallDeltaRough","TMath::Abs(peak_dphi_rough) < 4.0 && TMath::Abs(peak_dtheta_rough) < 4.0");

    const TCut isAboveHorizontal("isAboveHorizontal", "peak_theta > 0");
    const TCut anita3QuietTime = SumTree::anita3QuietTime; // should work for both
    const TCut isNotTaggedAsPulser("isNotTaggedAsPulser", "flags_pulser==0");
    const TCut isTaggedAsWaisPulser("isTaggedAsWaisPulser", TString::Format("flags_pulser==%d && run >= 331 && run <= 354",
									    AnitaEventSummary::EventFlags::WAIS).Data());
    const TCut realSnr("realSnr", TString::Format("%s && %s && %s && %s",
						  "!TMath::IsNaN(coherent_filtered_snr)",
						  "!TMath::IsNaN(deconvolved_filtered_snr)",
						  "TMath::Finite(coherent_filtered_snr)",
						  "TMath::Finite(deconvolved_filtered_snr)"));
    const TCut isRF("isRF", "flags_isRF>0");
    const TCut notShortWaveform("notShortWaveform", "flags_isVarner2==0");
    const TCut goodGps("goodGps", "TMath::Finite(anitaLocation_heading)");

    const TCut newPayloadBlastCutPart0A("newPayloadBlastCutPart0A", TString::Format("flags_middleOrBottomPower_0 < %lf * flags_topPower_0 + %lf", 30./7, 0.06));
    const TCut newPayloadBlastCutPart0B("newPayloadBlastCutPart0B", TString::Format("flags_middleOrBottomPower_0 > %lf * (flags_topPower_0 - %lf)", 7./30, 0.06));
    const TCut newPayloadBlastCutPart1("newPayloadBlastCutPart1", TString::Format("flags_middleOrBottomPower_1 < %lf * flags_topPower_1", 1./0.28)); /// (N=1)
    const TCut newPayloadBlastCutPart2("newPayloadBlastCutPart2", TString::Format("TMath::Max(flags_maxBottomToTopRatio_0, flags_maxBottomToTopRatio_1) < %lf", 2.67)); /// (N=1)
    const TCut newPayloadBlastCutPart3("newPayloadBlastCutPart3", TString::Format("%lf*deconvolved_filtered_peakHilbert > (pol==0)*((1+flags_maxBottomToTopRatio_0)*flags_minBottomToTopRatio_0 - %lf) + (pol==1)*((1+flags_maxBottomToTopRatio_1)*flags_minBottomToTopRatio_1 - %lf)"
					      ,14.0, 1000.0, 1000.0));
					      // ,11.0, 1000.0, 1000.0));    
    const TCut notInBlastRotatedCrossCorrelationRegion("notInBlastRotatedCrossCorrelationRegion",
						       "!((coherent_filtered_peakHilbert > 100) && (peak_value < 0.15))");

    const TCut notPayloadBlast("notPayloadBlast", TString::Format("(%s) && (%s) && (%s) && (%s) && (%s) && (%s)",
								  newPayloadBlastCutPart0A.GetTitle(),
								  newPayloadBlastCutPart0B.GetTitle(),
								  newPayloadBlastCutPart1.GetTitle(),
								  newPayloadBlastCutPart2.GetTitle(),
								  newPayloadBlastCutPart3.GetTitle(),
								  notInBlastRotatedCrossCorrelationRegion.GetTitle()));

    const TCut passAllQualityCuts("passAllQualityCuts", TString::Format("(%s) && (%s) && (%s) && (%s) && (%s)",
									notPayloadBlast.GetTitle(),
									notShortWaveform.GetTitle(),
									realSnr.GetTitle(),
									goodGps.GetTitle(),
									smallDeltaRough.GetTitle()));


    const TCut analysisSample("analysisSample",   TString::Format("(%s) && !(%s) && (%s)", isRF.GetTitle(), isAboveHorizontal.GetTitle(), isNotTaggedAsPulser.GetTitle()));
    const TCut thermalSideBand("thermalSideBand", TString::Format("(%s) &&  (%s) && (%s)", isRF.GetTitle(), isAboveHorizontal.GetTitle(), isNotTaggedAsPulser.GetTitle()));

    const TCut higherHilbertPeakAfterDedispersion("higherHilbertPeakAfterDedispersion",
						  "deconvolved_filtered_peakHilbert > coherent_filtered_peakHilbert");
    const TCut higherImpulsivityMeasureAfterDedispersion("higherImpulsivityMeasureAfterDedispersion",
							 "deconvolved_filtered_impulsivityMeasure > coherent_filtered_impulsivityMeasure");
    const TCut lowerFracPowerWindowGradientAfterDedispersion("lowerFracPowerWindowGradientAfterDedispersion",
							     "deconvolved_filtered_fracPowerWindowGradient < coherent_filtered_fracPowerWindowGradient");

    const TString fisherDiscriminant = "0.898497+(1.929594*coherent_filtered_fracPowerWindowGradient/deconvolved_filtered_fracPowerWindowGradient)+(-0.195909*deconvolved_filtered_fracPowerWindowGradient)+(5.943355*coherent_filtered_impulsivityMeasure)+(0.826114*deconvolved_filtered_impulsivityMeasure)+(0.021763*coherent_filtered_peakHilbert)+(-0.012670*deconvolved_filtered_peakHilbert)+(-0.394201*peak_value)";

    const TCut fisherCut("fisherCut", ThermalTree::fisherDiscriminant + " > 5.800012");

    const TCut continentNotIceShelf("continentNotIceShelf", "RampdemReader::isOnContinent(longitude,latitude)");
  }

}

#endif
