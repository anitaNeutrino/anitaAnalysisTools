#include "SummarySet.h"
#include "ProgressBar.h"
#include "OutputConvention.h"
#include "AnitaEventSummary.h"

#include "TH2D.h"
#include "AnalysisPlot.h"
#include "AnalysisCuts.h"

using namespace Acclaim;


/** 
 * A script initially developed for looking at the 10% data set.
 * Although, I suppose it can look at any set of AnitaEventSummary files
 * 
 * @param summaryTreeFileGlob is a globbing expression for trees containing AnitaEventSummaries
 * (make sure you add a single quotes at the start and end to prevent the shell expanding your wildcard!)
 * 
 * @return 0 on success, 1 if too many arguments
 */

int main(int argc, char* argv[]){

  if(argc!=2){
    std::cerr << argv[0] << " 'globWithWildcardsLikeThis*.root'" << std::endl;
    std::cerr << "(make sure you enclose the glob inside single quotes)" << std::endl;
    return 1;
  }

  OutputConvention oc(1, argv);
  TFile* fOut = oc.makeFile();

  const char* summaryTreeFileGlob = argv[1];
  SummarySet ss(summaryTreeFileGlob);  

  const Long64_t N = ss.N();
  std::cout << "Processing " << summaryTreeFileGlob << " with " << N << " entries." << std::endl;

  AnalysisPlot* hPeakVsTime = ss.bookTimeAnalysisPlot("hPeakVsTime", "Higher map peak vs time; realTime; Map peak", 1024, 128, 0, 1);
  hPeakVsTime->addCut(&AnalysisCuts::isAboveHorizontal);
  hPeakVsTime->addCut(&AnalysisCuts::isTaggedAsWaisPulser);
  hPeakVsTime->addCut(&AnalysisCuts::isTaggedAsPayloadBlast);
  hPeakVsTime->addCut(&AnalysisCuts::isGood);
  hPeakVsTime->addCut(&AnalysisCuts::higherPol);
  hPeakVsTime->addCut(&AnalysisCuts::hasSourceLocation);

  AnalysisProf* hDeltaAngleSunProf = new AnalysisProf("hDeltaAngleSunProf", "Angle between peak and sun;#delta#phi (Degrees);#delta#theta (Degrees)", 128, -20, 20, 128, -20, 20);
  hDeltaAngleSunProf->addCut(&AnalysisCuts::higherPol);

  AnalysisPlot* hDeltaAngleSun = new AnalysisPlot("hDeltaAngleSun", "Angle between peak and sun;#delta#phi (Degrees);#delta#theta (Degrees)", 128, -180, 180);
  hDeltaAngleSun->addCut(&AnalysisCuts::higherPol);

  
  AnalysisPlot* hPeakThetaVsTime = ss.bookTimeAnalysisPlot("hPeakThetaVsTime", "Peak elevation vs time; realTime; Peak #theta (Degrees)", 1024, 128, -90, 90);
  hPeakThetaVsTime->addCut(&AnalysisCuts::isTaggedAsWaisPulser);
  hPeakThetaVsTime->addCut(&AnalysisCuts::higherPol);
  
  AnalysisPlot* hPeakBearingVsTime = ss.bookTimeAnalysisPlot("hPeakBearingVsTime", "Bearing of map peak vs time; realTime; Peak #phi (Degrees)", 1024, 256, 0, 360);
  hPeakBearingVsTime->addCut(&AnalysisCuts::isTaggedAsWaisPulser);
  hPeakBearingVsTime->addCut(&AnalysisCuts::higherPol);

  AnalysisPlot* hImagePeakVsCoherentHilbertPeak = new AnalysisPlot("hImagePeakVsCoherentHilbertPeak", "Coherent Rotated Cross Correlation", 128, 0, 1, 1024, 0, 1000);
  AnalysisPlot* hImagePeakVsCoherentFilteredHilbertPeak = new AnalysisPlot("hImagePeakVsCoherentFilteredHilbertPeak", "Coherent Filtered Rotated Cross Correlation", 128, 0, 1, 1024, 0, 1000);
  AnalysisPlot* hImagePeakVsDeconvolvedHilbertPeak = new AnalysisPlot("hImagePeakVsDeconvolvedHilbertPeak", "Deconvolved Rotated Cross Correlation", 128, 0, 1, 1024, 0, 1000);
  AnalysisPlot* hImagePeakVsDeconvolvedFilteredHilbertPeak = new AnalysisPlot("hImagePeakVsDeconvolvedFilteredHilbertPeak", "Deconvolved Filtered Rotated Cross Correlation", 128, 0, 1, 1024, 0, 1000);
  AnalysisPlot* hDeconvolvedPeakHilbertVsDeconvolvedPeakHilbertTime = new AnalysisPlot("hDeconvolvedPeakHilbertVsDeconvolvedPeakHilbertTime", "Deconvolved Hilbert Peak Vs. Deconvolved Hilbert Peak Time", 128, 0, 100, 1024, 0, 1000);

  const int nWaves = 5;
  AnalysisPlot* hs[nWaves] = {hImagePeakVsCoherentHilbertPeak, hImagePeakVsCoherentFilteredHilbertPeak, hImagePeakVsDeconvolvedHilbertPeak, hImagePeakVsDeconvolvedFilteredHilbertPeak,
                              hDeconvolvedPeakHilbertVsDeconvolvedPeakHilbertTime};
  for(int i=0; i < nWaves; i++){
    hs[i]->addCut(&AnalysisCuts::isAboveHorizontal);
    hs[i]->addCut(&AnalysisCuts::isTaggedAsWaisPulser);
    hs[i]->addCut(&AnalysisCuts::isTaggedAsPayloadBlast);
    hs[i]->addCut(&AnalysisCuts::isGood);
    hs[i]->addCut(&AnalysisCuts::isWithin20DegreesOfSunInPhi);
    hs[i]->addCut(&AnalysisCuts::higherPol);
    hs[i]->addCut(&AnalysisCuts::hasSourceLocation);
  }
  
  ProgressBar p(N);
  for(Long64_t entry=0; entry < N; entry++){

    ss.getEntry(entry);
    AnitaEventSummary* sum = ss.summary();

    hPeakVsTime->Fill(sum, sum->realTime, sum->higherPeak().value);

    hDeltaAngleSunProf->Fill(sum, sum->dPhiSun(), sum->dThetaSun(), sum->higherPeak().value);
    hDeltaAngleSun->Fill(sum, sum->dPhiSun());

    hPeakThetaVsTime->Fill(sum, sum->realTime, sum->higherPeak().theta);
    hPeakBearingVsTime->Fill(sum, sum->realTime, sum->peakBearing());

    hImagePeakVsCoherentHilbertPeak->Fill(sum, sum->higherPeak().value, sum->higherCoherent().peakHilbert);
    hImagePeakVsCoherentFilteredHilbertPeak->Fill(sum, sum->higherPeak().value, sum->higherCoherentFiltered().peakHilbert);
    hImagePeakVsDeconvolvedHilbertPeak->Fill(sum, sum->higherPeak().value, sum->higherDeconvolved().peakHilbert);
    hImagePeakVsDeconvolvedFilteredHilbertPeak->Fill(sum, sum->higherPeak().value, sum->higherDeconvolvedFiltered().peakHilbert);

    hDeconvolvedPeakHilbertVsDeconvolvedPeakHilbertTime->Fill(sum, sum->higherDeconvolvedFiltered().peakTime, sum->higherDeconvolvedFiltered().peakHilbert);
    p.inc(entry, N);
  }
  
  fOut->Write();
  fOut->Close();
  return 0;
}














