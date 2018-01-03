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
    const TString dPhiWais = "FFTtools::wrap(peak[][].phi - wais.phi, 360, 0)";
    const TString dThetaWais = "(peak[][].theta - wais.theta)";

    const TString dPhiSun = "FFTtools::wrap(peak[][].phi - sun.phi, 360, 0)";
    const TString dThetaSun = "(peak[][].theta - sun.theta)";
    
    const TString dPhiMC = "FFTtools::wrap(peak[][].phi - mc.phi, 360, 0)";
    const TString dThetaMC = "(peak[][].theta - mc.theta)";

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

    /**
     * These cuts are just flags, only Length$ = 1
     */    
    const TCut isRfTrigger("isRfTrigger", "flags.isRF > 0"); /// (N=1)
    const TCut isNotTaggedAsPulser("isNotTaggedAsPulser", "flags.pulser == 0"); /// (N=1)
    const TCut anita3QuietTime("anita3QuietTime", "realTime >= 1419320000 && sum->realTime < 1420250000"); /// (N=1)
    // const TCut realSNR("realSNR", "(!TMath::IsNaN(deconvolved.snr) && TMath::Finite(deconvolved.snr))");
    const TCut goodGPS("goodGPS", "(!TMath::IsNaN(anitaLocation.heading) && TMath::Finite(anitaLocation.heading))"); /// (N=1)
    const TCut isGood("isGood", "((mc.weight > 0 && sum->flags.isVarner == 0 && sum->flags.isPayloadBlast == 0) || (mc.weight==0 && flags.isGood == 1))"); /// (N=1)
    const TCut isTaggedAsPayloadBlast("isTaggedAsPayloadBlast", "flags.isPayloadBlast > 0"); /// (N=1)
    const TCut isTaggedAsLDBPulser("isTaggedAsLDBPulser", TString::Format("!(AnitaVersion::get()==3 && run >=200) && flags.pulser==%d", AnitaEventSummary::EventFlags::LDB).Data()); /// (N=1)
    const TCut isTaggedAsWaisPulser("isTaggedAsWaisPulser", TString::Format("flags.pulser==%d", AnitaEventSummary::EventFlags::WAIS).Data()); /// (N=1)
    const TCut isTaggedAsPulser("isTaggedAsPulser", TString::Format("(%s) || (%s)", isTaggedAsLDBPulser.GetTitle(), isTaggedAsWaisPulser.GetTitle())); /// (N=1)
    const TCut npbc0A("npbc0A", TString::Format("flags.middleOrBottomPower[0] < %lf * flags.topPower[0] + %lf", 30./7, 0.06)); /// (N=1)
    const TCut npbc0B("npbc0B", TString::Format("flags.middleOrBottomPower[0] > %lf * (flags.topPower[0] - %lf)", 7./30, 0.06)); /// (N=1)
    const TCut npbc1("npbc1", TString::Format("flags.middleOrBottomPower[1] < %lf * flags.topPower[1]", 1./0.28)); /// (N=1)
    const TCut npbc2("npbc2", TString::Format("Max$(flags.maxBottomToTopRatio) < %lf", 4.0)); /// (N=1)
    



    // some constants
    const double maxDeltaPhi = 5.5;
    const double maxDeltaTheta = 3.5;


    
    /**
     * These cuts are per-direction (i.e. Iteration$ runs from 0-9 traversing the 2D array [2][5])
     */    
    const TCut highestPeak("highestPeak", "Max$(peak[][].value)==peak[][].value"); /// (N=10) Basically, you should almost always use something like this to just pick the peak direction.
    const TCut smallDeltaRough("smallDeltaRough", "TMath::Abs(peak[][].dphi_rough) < 4 && TMath::Abs(peak[][].dtheta_rough) < 4"); /// (N=10)
    // const TCut isNotNorth("isNotNorth", "TMath::Abs(peak.dPhiNorth()) > 90");
    // const TCut acceptableHardwareAngle("acceptableHardwareAngle", "TMath::Abs(sum.minAbsHwAngle()) < 65.0");
    const TCut higherHilbertPeakAfterDedispersion("higherHilbertPeakAfterDedispersion", "deconvolved_filtered[][].peakHilbert > coherent_filtered[][].peakHilbert"); ///(N=10)
    const TCut higherImpulsivityMeasureAfterDedispersion("higherImpulsivityMeasureAfterDedispersion", "deconvolved_filtered[][].impulsivityMeasure > coherent_filtered[][].impulsivityMeasure"); //////(N=10)

    const TCut lowerFracPowerWindowGradientAfterDedispersion("lowerFracPowerWindowGradientAfterDedispersion", ///(N=10)
							     TString::Format("(%s) < (%s)",
									     Draw::deconvolved_filtered_fracPowerWindowGradient.Data(),
									     Draw::coherent_filtered_fracPowerWindowGradient.Data()));
    const TCut closeToWais("closeToWais", /// (N=10)
			   TString::Format("(mc.weight == 0 && %s < %lf && %s > %lf && %s < %lf)",
					   Draw::dPhiWais.Data(), maxDeltaPhi, Draw::dPhiWais.Data(), -maxDeltaPhi,
					   Draw::dThetaWais.Data(), maxDeltaTheta));
    const TCut closeToMC("closeToMC", /// (N=10)
			 TString::Format("(mc.weight > 0 && %s < %lf && %s > %lf && %s < %lf %s > %lf)",
					 Draw::dPhiMC.Data(),   maxDeltaPhi,   Draw::dPhiMC.Data(),   -maxDeltaPhi,
					 Draw::dThetaMC.Data(), maxDeltaTheta, Draw::dThetaMC.Data(), -maxDeltaTheta));

    // const TCut isWithin20DegreesOfSunInPhi("isWithin20DegreesOfSunInPhi", "sum->peak.dPhiSun()) < 20");
    const TCut isOnContinent("isOnContinent", "RampdemReader::isOnContinent(peak[][].longitude, peak[][].latitude)"); /// (N=10)
    
    const TCut hasSourceLocation("hasSourceLocation", "(peak[][].latitude < -900 || TMath::Abs(peak[][].theta_adjustment_needed) > 0"); /// (N=10)
    const TCut isAboveHorizontal("isAboveHorizontal", "peak[][].theta > 0"); /// (N=10)

    const TCut npbc3("npbc3", TString::Format("%lf*deconvolved_filtered[][].peakHilbert > (1+flags.maxBottomToTopRatio[Iteration$/5])*flags.minBottomToTopRatio[Iteration$/5] - %lf", 1000.0, (15000.0-1000.0)/1000.0)); /// (N=10)
    const TCut isGood2("isGood2", TString::Format("(%s && %s && %s && %s && %s)", npbc0A.GetTitle(), npbc0B.GetTitle(), npbc1.GetTitle(), npbc2.GetTitle(), npbc3.GetTitle())); /// N(10)
    
  }
}

#endif
