#include "SummaryDraw.h"
//#include "SumTreeReductionSelector.h"
//#include "SummarySet.h"
#include "OutputConvention.h"
// #include <iostream>
// #include "TSystem.h"
#include "ProgressBar.h"

// #include "TProof.h"
// #include "TEntryList.h"
// #include "TTreeFormula.h"

// using namespace Acclaim;

void doDataReduction(const char* glob = NULL){

  if(!glob){
    std::cerr <<  "doDataReduction('glob')" << std::endl;
    return 1;
  }

  const char* fakeArgv[1] = {"doDataReduction2"};
  
  Acclaim::OutputConvention oc(1, (char**)fakeArgv);

  // Acclaim::Cuts::setMode(Acclaim::Cuts::kAcclaimAnalysis);

  Acclaim::SummarySet ss(glob);
  ss.SetUseProof(true);

  TCut allCuts = Acclaim::Cuts::highestPeak;
  allCuts += Acclaim::Cuts::isRfTrigger;
  // allCuts += Acclaim::Cuts::isNotTaggedAsPulser;
  allCuts += Acclaim::Cuts::npbc0A;
  allCuts += Acclaim::Cuts::npbc0B;
  allCuts += Acclaim::Cuts::npbc1;
  allCuts += Acclaim::Cuts::npbc2;
  allCuts += Acclaim::Cuts::npbc3;
  allCuts += Acclaim::Cuts::smallDeltaRough;
  allCuts += Acclaim::Cuts::goodGPS;
  allCuts += Acclaim::Cuts::isBelowHorizontal;
  allCuts += Acclaim::Cuts::reasonableHilbertPeakTimeShiftAfterDedispersion;
  allCuts += Acclaim::Cuts::higherHilbertPeakAfterDedispersion;
  allCuts += Acclaim::Cuts::higherImpulsivityMeasureAfterDedispersion;
  allCuts += Acclaim::Cuts::lowerFracPowerWindowGradientAfterDedispersion;
  allCuts += TCut("(%s) > 5.0", Acclaim::Draw::fisherDiscriminant.Data()); // probably much lower than we will set the limit
  ss.Draw(">>elist1", allCuts, "entrylist");

  TEntryList* elist = (TEntryList*) gProof->GetOutputList()->FindObject("elist1");


  Acclaim::ProgressBar p(elist->GetN());


  // TFile* fOut = oc.makeFile();
  // TTree* sumTree = new TTree("sumTree", "sumTree");
  AnitaEventSummary* sum = NULL;
  // sumTree->Branch("sum", &sum);
  double fisherDiscriminant;
  // sumTree->Branch("fisherDiscriminant", &fisherDiscriminant);

  TTreeFormula fisherFormula("fisherFormula", Acclaim::Draw::fisherDiscriminant, ss.getChain());
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


  // gSystem->Exit(0); /// Must be called after initializing PROOF!
    

  return 0;
}
