#include "AnitaEventSummary.h"
#include "TTree.h"
#include "TFile.h"
#include "SummarySet.h"
#include "AnitaDataset.h"
#include "ProgressBar.h"
#include "AnalysisReco.h"
#include "FilterStrategy.h"
#include "InterferometricMap.h"

using namespace Acclaim;


/** 
 * Script to make a copy of summary trees, and hackily edit a couple of values if possible
 */
int main(int argc, char* argv[]){

  if(argc!=4){
    std::cerr << "Usage: " <<  argv[0] << " inputDir inputGlob run" << std::endl;
    return 1;
  }

  // std::vector<UInt_t> eventNumbers {26323378, 27569662, 27897245, 28254843, 28627679, 28776748, 28930237, 28934295, 29006913, 29008144, 32879955, 32953817, 34745811, 34757629, 35020101, 35116217, 35118197, 36710262, 37009416, 37009533, 37009810, 38273261, 38619569, 39464799, 40633613, 40646231, 40996918, 41492058, 41654702, 41799439, 42256877, 42699600, 44675867, 45985415, 47512906, 47831713, 48225288, 48242449, 49313992, 50827460, 51950538, 53275871, 55496464, 55504417, 55861013, 56463512, 56829451, 56894100, 57108465, 57471113, 57519585, 57523926, 57654756, 57839882, 57948519, 58002861, 58008286, 58110716, 58226572, 58310817, 58326495, 58332241, 58421469, 59967256, 60072506, 61104509, 61173194, 61215330, 61238684, 61497702, 61534104, 61784557, 61876418, 62235301, 62235933, 62420819, 62732869, 62849802, 63024583, 63031985, 63484004, 63714468, 63845258, 63845261};
 
  TString inDir = argv[1];
  TString inFileBaseName = argv[2];
  int run = atoi(argv[3]);

  TString inGlob = inDir + "/" + inFileBaseName;// + TString::Format("_%d_*.root", run);
  // TString inGlob = inDir + "/" + inFileBaseName + TString::Format("_%d_*.root", run);  
  SummarySet ss(inGlob);

  const Long64_t N = ss.N();
  if(N > 0){    
    // TString outFileName = inFileBaseName + TString::Format("_%d.root", run);
    TString outFileName = TString::Format("editSumTree_%d.root", run);
    TFile* fOut = new TFile(outFileName, "recreate");
    TTree* sumTree = new TTree("sumTree", "sumTree");
    AnitaEventSummary* sum = NULL;
    sumTree->Branch("sum", &sum);

    AnalysisSettings settings;
    Acclaim::AnalysisReco reco;
    FilterStrategy doNoFiltering;
    AnitaDataset* d = NULL;    
    
    settings.apply(dynamic_cast<TObject*>(&reco));
    settings.write(fOut);
    
    ProgressBar p(N);
    // ProgressBar p(eventNumbers.size());
    for(Long64_t entry=0; entry < N; entry++){
      ss.getEntry(entry);

    // for(UInt_t eventInd=0; eventInd < eventNumbers.size(); eventInd++){
      // ss.getEvent(eventNumbers.at(eventInd));
      AnitaEventSummary* inSum = ss.summary();
      
      // if(std::find(eventNumbers.begin(), eventNumbers.end(), inSum->eventNumber)!=eventNumbers.end()){
      if(inSum->flags.isPayloadBlast || inSum->flags.pulser > 0){


	if(d && d->getCurrRun() != inSum->run){
	  delete d;
	  d = NULL;
	}
      
	if(!d){
	  d = new AnitaDataset(inSum->run);	
	}

	d->getEvent(inSum->eventNumber);
	FilteredAnitaEvent fEv(d->useful(), &doNoFiltering, d->gps(), d->header());


	// d.getEvent(inSum->eventNumber);

	*sum = *inSum;


	// for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
	//   for(int peakInd=0; peakInd < AnitaEventSummary::maxDirectionsPerPol; peakInd++){
	//     double temp = sum->peak[polInd][peakInd].latitude;
	//     sum->peak[polInd][peakInd].latitude = sum->peak[polInd][peakInd].longitude;
	//     sum->peak[polInd][peakInd].longitude = temp;
	//   }
	// }

	// const double phi0 = InterferometricMap::getBin0PhiDeg();
	// for(int polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){
	//   AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
	//   for(int peakInd=0; peakInd < sum->nPeaks[pol]; peakInd++){

	//     // hack to recover peakPhiSector...
	//     double phi_rough = sum->peak[pol][peakInd].phi - sum->peak[pol][peakInd].dphi_rough;
	//     int peakPhiSector = floor((phi_rough - phi0)/22.5);
	//     // std::cout << peakPhiSector << "\t" << phi_rough << std::endl;
	//     AnalysisReco::setTriggerInfoFromPeakPhi(d.header(), pol, peakPhiSector, sum->peak[pol][peakInd]);
	//   }
	// }

	const double ratioCutHigh = 2.8;
	if(sum->flags.maxBottomToTopRatio[0] > ratioCutHigh || sum->flags.maxBottomToTopRatio[1] > ratioCutHigh){
	  sum->flags.isPayloadBlast = 1;
	}
	else{
	  sum->flags.isPayloadBlast = 0;
	}

	reco.fillPowerFlags(&fEv, sum);
	// std::cout << sum->flags.meanPower[0] << "\t" << sum->flags.meanPower[1] << "\t" << sum->flags.meanPower[2]  << std::endl;


	sumTree->Fill();

      }

      p.inc(entry, N);
      // p.inc(eventInd, N);
    }
    fOut->Write();
    fOut->Close();
  }
  
  return 0;
}
