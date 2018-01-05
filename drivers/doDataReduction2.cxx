#include "SummaryDraw.h"
#include "SumTreeReductionSelector.h"
#include "SummarySet.h"
#include "OutputConvention.h"
#include <iostream>
#include "TSystem.h"
#include "ProgressBar.h"

#include "TProof.h"
#include "TEntryList.h"
#include "TTreeFormula.h"

using namespace Acclaim;

int main(int argc, char* argv[]){

  if(argc != 2){
    std::cerr << argv[0] << " 'glob' " << std::endl;
    return 1;
  }
  OutputConvention oc(1, argv);  
  const char* glob = argv[1];

  // Cuts::setMode(Cuts::kAcclaimAnalysis);

  SummarySet ss(glob);
  ss.SetUseProof(true);

  TCut allCuts = Cuts::highestPeak;
  allCuts += Cuts::isRfTrigger;
  // allCuts += Cuts::isNotTaggedAsPulser;
  allCuts += Cuts::npbc0A;
  allCuts += Cuts::npbc0B;
  allCuts += Cuts::npbc1;
  allCuts += Cuts::npbc2;
  allCuts += Cuts::npbc3;
  allCuts += Cuts::smallDeltaRough;
  allCuts += Cuts::goodGPS;
  allCuts += Cuts::isBelowHorizontal;
  allCuts += Cuts::reasonableHilbertPeakTimeShiftAfterDedispersion;
  allCuts += Cuts::higherHilbertPeakAfterDedispersion;
  allCuts += Cuts::higherImpulsivityMeasureAfterDedispersion;
  allCuts += Cuts::lowerFracPowerWindowGradientAfterDedispersion;
  allCuts += TCut("(%s) > 5.0", Draw::fisherDiscriminant.Data()); // probably much lower than we will set the limit
  ss.Draw(">>elist1", allCuts, "entrylist");

  TEntryList* elist = (TEntryList*) gProof->GetOutputList()->FindObject("elist1");


  ProgressBar p(elist->GetN());


  // TFile* fOut = oc.makeFile();
  // TTree* sumTree = new TTree("sumTree", "sumTree");
  AnitaEventSummary* sum = NULL;
  // sumTree->Branch("sum", &sum);
  double fisherDiscriminant;
  // sumTree->Branch("fisherDiscriminant", &fisherDiscriminant);

  TTreeFormula fisherFormula("fisherFormula", Draw::fisherDiscriminant, ss.getChain());
  ss.getChain()->SetNotify(&fisherFormula);
  
  for(int i=0; i < elist->GetN(); i++){
    Long64_t entry = elist->Next();
    Int_t treenum = 0;
    elist->GetEntryAndTree(i, treenum);
    ss.getEntry(entry);
    sum = ss.summary();


    Int_t polInd = sum->highestPolAsInt();
    Int_t peakInd = sum->highestPeakInd();
    Int_t iteration = polInd*AnitaEventSummary::maxDirectionsPerPol + peakInd;

    fisherDiscriminant = fisherFormula.EvalInstance(iteration);

    std::cout << entry << "\t" << treenum << "\t" << sum->eventNumber << "\t" << polInd <<  "\t" << peakInd << "\t" << fisherDiscriminant << std::endl;
    // for(int j=0; j < 10; j++){
    //   std::cout << fisherFormula.EvalInstance(j) << "\t";
    // }
    // std::cout << std::endl;

    // sumTree->Fill();
    p.inc(i);
  }
    
  // fOut->Write();
  // fOut->Close();  


  gSystem->Exit(0); /// Must be called after initializing PROOF!
    

  return 0;
}
