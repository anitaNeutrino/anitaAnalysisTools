#include "AnalysisFlow.h"

using namespace Acclaim;


/** 
 * All the complicated objects persist inside AnalysisReco until the next event is processed
 * This macro just analyses one event and extracts all the gory details.
 * Use this when inspecting the AnitaEventSummary simply won't do.
 * 
 * @param run is the run containing the event of interest
 * @param eventNumber is the event of interest
 */
void interactiveAnalysisOneEvent(int run=352, UInt_t eventNumber=60832108){

  AnalysisFlow flow(NULL, run, AnalysisFlow::kAll);
  flow.doAnalysis(eventNumber);
  AnalysisReco* reco = flow.getReco(); // this object is what is in magic display

  
  TCanvas* cCoarseMaps = new TCanvas("cCoarseMaps", "The coarse maps", 1200, 600);
  InterferometricMap* coarseMaps[AnitaPol::kNotAPol] = {reco->getMap(AnitaPol::kHorizontal), reco->getMap(AnitaPol::kVertical)};
  cCoarseMaps->Divide(AnitaPol::kNotAPol);
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    cCoarseMaps->cd(polInd + 1);
    coarseMaps[polInd]->Draw("colz");
  }  
  

  UInt_t np =  reco->GetNumPeaks();
  std::vector<InterferometricMap*>zoomMaps[AnitaPol::kNotAPol];
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    AnitaPol::AnitaPol_t pol = AnitaPol::AnitaPol_t(polInd);
    for(UInt_t peakInd = 0 ; peakInd < np; peakInd++){
      zoomMaps[polInd].push_back(reco->getZoomMap(pol, peakInd));
    }
  }
  TCanvas* cZoomMaps = new TCanvas("cFineMaps", "The fine maps", 1200, 600);
  cZoomMaps->Divide(AnitaPol::kNotAPol);
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    cZoomMaps->cd(polInd + 1);
    zoomMaps[polInd].at(0)->Draw("colz");
  }  


  TCanvas* cWaves = new TCanvas("cWaves", "The coherent waveforms", 1200, 600);
  cWaves->Divide(AnitaPol::kNotAPol);  
  AnalysisWaveform* coherentWaves[AnitaPol::kNotAPol];
  AnalysisWaveform* deconvolvedWaves[AnitaPol::kNotAPol];  
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    cWaves->cd(polInd+1);
    AnitaPol::AnitaPol_t pol = AnitaPol::AnitaPol_t(polInd);
    coherentWaves[pol] = reco->getCoherentFiltered(pol);
    TGraphAligned* gr = coherentWaves[pol]->updateEven();
    gr->Draw();

    deconvolvedWaves[pol] = reco->getDeconvolved(pol);
    TGraphAligned* gr2 = deconvolvedWaves[pol]->updateEven();
    gr2->SetLineColor(kRed);
    gr2->Draw("lsame");
  }

  
  
  
}










