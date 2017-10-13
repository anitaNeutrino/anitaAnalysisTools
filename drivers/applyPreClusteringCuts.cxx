#include "AnalysisCuts.h"
#include "OutputConvention.h"
#include "SummarySet.h"


using namespace Acclaim;


/** 
 * Get a single file of AnitaEventSummaries passing all pre-clustering cuts, as defined here
 */
int main(int argc, char* argv[]){

  if(!(argc == 2)){
    std::cerr << argv[0] << " 'glob' " << std::endl;
    return 1;
  }
  const char* glob = argv[1];

  SummarySet ss(glob);

  OutputConvention oc(1, argv);
  TFile* fOut = oc.makeFile();
  TTree* sumTree = new TTree("sumTree", "Tree of AnitaEventSummaries passing pre-clustering cuts");
  AnitaEventSummary* sum = NULL;
  sumTree->Branch("sum", &sum);

  const int nCut = 9;
  const AnalysisCuts::AnalysisCut* preClusteringCuts[nCut] = {&AnalysisCuts::isGood, // not payload blast, SURF saturation, 
							      &AnalysisCuts::smallDeltaRough, // agreement between coarse/fine peak
							      &AnalysisCuts::goodGPS, // do we have GPS data?
							      &AnalysisCuts::realSNR,
							      &AnalysisCuts::isRfTrigger,
							      &AnalysisCuts::higherPeakHilbertAfterDedispersion,
							      &AnalysisCuts::higherImpulsivityMeasureAfterDedispersion,
							      &AnalysisCuts::lowerFracPowerWindowGradientAfterDedispersion,
							      &AnalysisCuts::dedispersedFracPowerWindowGradientBelowThreshold};
  
  Long64_t n = ss.N();
  for(Long64_t entry=0; entry < n; entry++){
    ss.getEntry(entry);
    
    sum = ss.summary();

    // AnitaPol::AnitaPol_t pol = sum->acclaimPol();
    // Int_t peakIndex = sum->acclaimPeakIndex();    
    AnitaPol::AnitaPol_t pol = AnitaPol::kVertical; ///@todo set to some acclaim specific function 
    Int_t peakIndex = 0; ///@todo set to some acclaim specific function

    Bool_t passesCuts = true;
    for(int i=0; i < nCut; i++){
      int cutRetVal = preClusteringCuts[i]->apply(sum, pol,  peakIndex);

      passesCuts = passesCuts && (cutRetVal > 0);

      if(!passesCuts){
	break;
      }

      if(passesCuts){
	sumTree->Fill();
      }
    }

    fOut->Write();
    fOut->Close();
    
  }

  return 0;
}
