#include "AnalysisFlow.h"
#include "AcclaimFilters.h"
#include "UCFilters.h"
#include "ProgressBar.h"

using namespace Acclaim;


void findDeconvolveOffset(){
  int run=352;

  TH2D* hs[NUM_SEAVEYS][NUM_SEAVEYS] = {{NULL}};
  
  FilterStrategy* strat = new FilterStrategy();  
  UCorrelator::fillStrategyWithKey(strat,Acclaim::Filters::getCosminsFavouriteSineSubName());
  
  AnalysisFlow flow(NULL, run, AnalysisFlow::kWaisPulser, strat);
  AnalysisReco* reco = flow.getReco();

  CrossCorrelator* ccDeco = new CrossCorrelator();
  const double norm = 1./(ccDeco->numSamples*ccDeco->numSamples);

  
  ProgressBar p(flow.lastEntry()); 
  for(Long64_t entry=flow.firstEntry(); entry < flow.lastEntry(); entry++){

    AnitaEventSummary* sum = flow.doEntry(entry);
    if(sum){

      if(TMath::Abs(sum->dPhiWais(0, AnitaPol::kHorizontal)) < 5){
      
        CrossCorrelator* ccReco = reco->getCrossCorrelator();
  
        InterferometricMap* mFine = reco->getZoomMap(AnitaPol::kHorizontal);
        Int_t peakPhi =  mFine->getPeakPhiSector();

        FilteredAnitaEvent* fEvDeco = reco->getEvDeco();

        // get the upsampled cross-correlations of the deconvolved graphs
        ccDeco->correlateEvent(fEvDeco, AnitaPol::kHorizontal);
        ccDeco->doUpsampledCrossCorrelations(AnitaPol::kHorizontal, peakPhi);

        std::vector<int> ants = reco->phiSectorToCoherentAnts(peakPhi);


        for(UInt_t i=0; i < ants.size() - 1; i++){
          int ant1 = ants[i];
          for(UInt_t j=i+1; j < ants.size();  j++){
            int ant2 = ants[j];
            TGraph* grReco = ccReco->getUpsampledCrossCorrelationGraph(AnitaPol::kHorizontal, ant1, ant2);
            TGraph* grDeco = ccDeco->getUpsampledCrossCorrelationGraph(AnitaPol::kHorizontal, ant1, ant2);

            if(grReco && grDeco){

              if(!hs[ant1][ant2]){
                TString n = TString::Format("h_%d_%d",  ant1, ant2);
                hs[ant1][ant2] = new TH2D(n, n, 128, -10, 10, 128, -10, 10);
              }
        
              // std::cout << ant1 << "\t" << ant2 << "\t" << grReco << "\t" << grDeco << std::endl;

              Double_t cMax1, tMax1, cMin1, tMin1;
              RootTools::getMaxMin(grReco, cMax1, tMax1, cMin1, tMin1);

              Double_t cMax2, tMax2, cMin2, tMin2;
              RootTools::getMaxMin(grDeco, cMax2, tMax2, cMin2, tMin2);

              // std::cout << cMax1 * norm << std::endl;
            
              hs[ant1][ant2]->Fill(tMax1, tMax2, cMax1*norm);
        
              // std::cout << ant1 << "\t" << ant2 << "\t" << tMax1 - tMax2 << "\t" << cMax1 << "\t" << cMax2 << std::endl;
            }

            if(grReco) delete grReco;
            if(grDeco) delete grDeco;
          }
        }
        delete mFine;
      }
      delete sum;      
    }
    p.inc(entry, flow.lastEntry());
  }
  
}
