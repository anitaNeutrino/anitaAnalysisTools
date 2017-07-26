#include "AnalysisFlow.h"
#include "AcclaimFilters.h"
#include "UCFilters.h"

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

  FilterStrategy* strat = new FilterStrategy();  
  UCorrelator::fillStrategyWithKey(strat,Acclaim::Filters::getCosminsFavouriteSineSubName());
  
  AnalysisFlow flow(NULL, run, AnalysisFlow::kAll, strat);
  flow.doEvent(eventNumber);
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


  auto c = new TCanvas();
  const int nEv = 2;
  FilteredAnitaEvent* fEvs[nEv] = {reco->getEvMin(), reco->getEvDeco()};
  c->Divide(nEv);
  for(int e=0; e < nEv; e++){
    c->cd(e+1);
    FilteredAnitaEvent* fEv = fEvs[e];

    
    Int_t phiSect = zoomMaps[AnitaPol::kHorizontal][0]->getPeakPhiSector();
    std::vector<int> ants = reco->phiSectorToCoherentAnts(phiSect);
    double peakPhiDeg, peakThetaDeg, peakVal;
    zoomMaps[AnitaPol::kHorizontal][0]->getPeakInfo(peakVal, peakPhiDeg, peakThetaDeg);
    std::vector<double> dts;  
    reco->directionAndAntennasToDeltaTs(ants, AnitaPol::kHorizontal, peakPhiDeg, peakThetaDeg, dts);
    std::vector<const AnalysisWaveform *> waves;
    for(auto& ant : ants){
      waves.push_back(fEv->getFilteredGraph(ant, AnitaPol::kHorizontal));
    }
    std::vector<TGraphAligned*> grs;
    reco->wavesInCoherent(waves, dts, grs);

    for(UInt_t i=0; i < grs.size(); i++){
      TString opt = i ==0 ? "al" : "lsame";      
      grs[i]->Draw(opt);
      grs[i]->SetLineColor(i+1);
      if(i==0){
        grs[i]->SetMaximum(100);
        grs[i]->SetMinimum(-100);
        TString title = e == 0 ? "Not dedispersed" : "Dedispersed";
        title += "; Time (ns); Amplitude (mV)";
        grs[i]->SetTitle(title);        
      }
    }
  }  
  
  // reco->wavesInCoherent(, std::vector<Double_t> &dts, std::vector<TGraphAligned *> &grs)
  
  

  TCanvas* cWaves = new TCanvas("cWaves", "The coherent waveforms", 1200, 600);
  // cWaves->Divide(AnitaPol::kNotAPol);  
  AnalysisWaveform* coherentWaves[AnitaPol::kNotAPol];
  AnalysisWaveform* deconvolvedWaves[AnitaPol::kNotAPol];  
  for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
    // cWaves->cd(polInd+1);
    AnitaPol::AnitaPol_t pol = AnitaPol::AnitaPol_t(polInd);
    coherentWaves[pol] = reco->getCoherentFiltered(pol);
    TGraphAligned* gr = coherentWaves[pol]->updateEven();
    gr->Draw();

    deconvolvedWaves[pol] = reco->getDeconvolved(pol);
    TGraphAligned* gr2 = deconvolvedWaves[pol]->updateEven();
    gr2->SetLineColor(kRed);
    gr2->Draw("lsame");
    break;
  }

  
  
  
}










