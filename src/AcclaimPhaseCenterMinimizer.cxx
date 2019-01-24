
#include "AcclaimPhaseCenterMinimizer.h"
#include "AcclaimParameterManager.h"
#include "AcclaimPhaseCenterFakeGeom.h"
#include "OutputConvention.h"
#include "TChain.h"
#include "TFile.h"
#include "UsefulAdu5Pat.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

const char* Acclaim::PhaseCenter::Minimizer::toCString(AllowedPairs ap){
  switch(ap){
  case AllowedPairs::All:
    return "AllowedPairs::All";
  case AllowedPairs::SamePhiSector:
    return "AllowedPairs::SamePhiSector";
  case AllowedPairs::SamePhiSectorOrHorizontalNeigbour:
    return "AllowedPairs::SamePhiSectorOrHorizontalNeigbour";
  }
  return nullptr;
}




Acclaim::PhaseCenter::Minimizer::Minimizer(const char* corrTreeFiles){

  fChain = std::make_shared<TChain>("corrTree");
  if(corrTreeFiles){
    fChain->Add(corrTreeFiles);
  }
  readInSummaries();

  fMin = std::shared_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2"));
  fMin->SetMaxFunctionCalls(1e5);
  fMin->SetTolerance(0.01);
  fMin->SetPrintLevel(0);  

}


void Acclaim::PhaseCenter::Minimizer::setParameterSpace(ParameterSpace ps){

  fParamManager = ps;
  fParamManager.setInputs(fMin.get(), fInputs);
  fFunc = ROOT::Math::Functor(this, &Acclaim::PhaseCenter::Minimizer::eval, fParamManager.N());
  fMin->SetFunction(fFunc);
}


bool Acclaim::PhaseCenter::Minimizer::allowedPair(int ant1, int ant2) const {
  return allowedPair(fAllowedPairs, ant1, ant2);
}


bool Acclaim::PhaseCenter::Minimizer::allowedPair(AllowedPairs ap, int ant1, int ant2) {

  // ensure ant1 <= ant2
  if(ant1 > ant2){
    std::swap(ant1, ant2);
  }
  
  auto samePhiSector = [](int ant1, int ant2)->bool {return (ant1%NUM_PHI)==(ant2%NUM_PHI);};
  auto horizontalNeighbour = [](int ant1, int ant2)->bool {
			       bool oneDifferent = ant2 - ant1 == 1;
			       if((ant1 % NUM_PHI)==0){ //  if ant1 is first number in row
				 // then ant2 = ant1 + 15 is a horizontal neighbour
				 bool fifteenDifferent = ant2 - ant1 == NUM_PHI - 1;
				 return oneDifferent || fifteenDifferent;
			       }
			       return oneDifferent;
			     };
  
  switch(ap){
  case AllowedPairs::All:
    return true;
  case AllowedPairs::SamePhiSector:
    return samePhiSector(ant1, ant2);
  case AllowedPairs::SamePhiSectorOrHorizontalNeigbour:
    return samePhiSector(ant1, ant2) || horizontalNeighbour(ant1, ant2);
  }
  
  // shouldn't get here but suppress compiler warning...
  return false;
}




const std::vector<double>& Acclaim::PhaseCenter::Minimizer::fit(AnitaPol::AnitaPol_t pol){

  std::cout << "Starting fit!" << std::endl;  
  fFitPol = pol;

  fApplyParams = false;
  fInitial = eval(&fInputs[0]); // get the start value.
  fApplyParams = true;

  std::vector<const char*> fakeArgv;
  TString baseName = PhaseCenter::toCString(fParamManager.space());
  baseName.ReplaceAll("::", "_");
  fakeArgv.push_back(baseName.Data());

  TString polName;
  switch(fFitPol){
  case AnitaPol::kHorizontal : polName = "HPol";     break;
  case AnitaPol::kVertical   : polName = "VPol";     break;
  default                    : polName = "BothPols"; break;
  }  
  fakeArgv.push_back(polName.Data());

  TString pairName = toCString(fAllowedPairs);
  pairName.ReplaceAll("::", "_");
  fakeArgv.push_back(pairName.Data());
  
  Acclaim::OutputConvention oc(fakeArgv.size(), const_cast<char**>(&fakeArgv[0]));
  TString outputFileName = oc.getOutputFileName();;
  TFile f(outputFileName, "recreate");

  std::vector<AnitaPol::AnitaPol_t> pols;
  if(fFitPol==AnitaPol::kNotAPol){
    pols.push_back(AnitaPol::kHorizontal);
    pols.push_back(AnitaPol::kVertical);
  }
  else{
    pols.push_back(fFitPol);
  }
  
  fMin->Minimize();
  fResults.clear();
  const int nDim = fParamManager.N();
  fResults.reserve(nDim);
  
  for(int i=0; i < nDim; i++){
    fResults.push_back(fMin->X()[i]);
  }
  fMinimum = fMin->MinValue();

  fSaveResults = true;
  eval(&fResults[0]);
  fSaveResults = false;
  
  std::cout << std::endl;

  f.Write();
  f.Close();
  
  return fResults;
}

void Acclaim::PhaseCenter::Minimizer::readInSummaries(){

  fChain->SetBranchAddress("correlationSummary", &fSum);
  Long64_t n = fChain->GetEntries();
  fSummaries.reserve(n);
    
  for(Long64_t entry=0; entry < n; entry++){
    fChain->GetEntry(entry);

    // std::cout << fSummaries.size() << std::endl;

    auto& sum = *fSum;
    auto pat = sum.fPat;
    UsefulAdu5Pat usefulPat(&pat);

    double thetaWave;
    double phiWave;
    usefulPat.getThetaAndPhiWaveWaisDivide(thetaWave, phiWave);
    std::vector<int> keepMe;
    auto oldPairs = sum.fPairs;
    for(auto& pair : oldPairs){
      double dtExpected = usefulPat.getDeltaTExpected(pair.ant2, pair.ant1, phiWave, thetaWave);

      double dtMeasured = pair.dt;

      if(fabs(dtMeasured - dtExpected) > fDeltaDeltaTThreshold ||
	 pair.correlation < fCorrelationThreshold){
	keepMe.push_back(0);
      }
      else{
	keepMe.push_back(1);
      }
    }

    sum.fPairs.clear();
    for(int p=0; p < oldPairs.size(); p++){
      if(keepMe[p]==1){
	sum.fPairs.emplace_back(oldPairs[p]);
      }
    }

    if(sum.fPairs.size() > 0){
      fSummaries.emplace_back(sum);
    }    
  }

  if(fSummaries.size()==0){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << " No summaries read in!" << std::endl;
  }
}


void Acclaim::PhaseCenter::Minimizer::printResults() const {

  std::cout << "Results:\n";
  std::cout  << "Before the fit, the value was " << fInitial << " with:\n";
  const int np = fParamManager.N();;
  
  for(int i=0; i < np; i++){
    std::cout << fMin->VariableName(i) << " = " << fInputs.at(i) << "\n";
  }

  std::cout << "\nAfter the fitting, minimum value is " << fMinimum << " with:\n";
  for(int i=0; i < np; i++){
    std::cout << fMin->VariableName(i) << " = " << fResults.at(i) << "\n";
  }
  std::cout << "\n\n" << std::flush;
}


double Acclaim::PhaseCenter::Minimizer::eval(const double* params){
 
  auto geom = AnitaGeomTool::Instance();
  geom->usePhotogrammetryNumbers(true);

  FakeGeomTool fakeGeom(geom);
  fakeGeom.overwritePhotogrammetryPositionsInAnitaGeomTool(geom);
  
  fParamManager.update(params);
  // std::cout << "Eval! " << fApplyParams << std::endl;
  if(fApplyParams){
    fParamManager.applyGeom(geom);
  }

  TTree* fOutputTree = nullptr;
  CorrelationSummary* outputSummary = nullptr;
  if(fSaveResults){
    fOutputTree = new TTree("expectedTree", "expectedTree");
    outputSummary = new CorrelationSummary();
    fOutputTree->Branch("expected", &outputSummary);
  }

  std::array<std::array<Int_t, NUM_SEAVEYS>, AnitaPol::kNotAPol> normalization;
  std::array<std::array<Double_t, NUM_SEAVEYS>, AnitaPol::kNotAPol> ddts;

  // reset delta delta T counters
  for(auto& a : normalization){a.fill(0);}
  for(auto& a : ddts){a.fill(0);}

  // bool firstOne = true;
#pragma omp parallel if(!fSaveResults)
  for(auto& cs : fSummaries){    

    if(fFitPol==AnitaPol::kNotAPol || cs.fPol==fFitPol){

      Adu5Pat pat = cs.fPat;
      if(fApplyParams){
	fParamManager.applyPat(&pat);
      }
      
      double thetaWave, phiWave;
      UsefulAdu5Pat usefulPat(&pat, false);
      usefulPat.getThetaAndPhiWaveWaisDivide(thetaWave, phiWave);

      int pairIndex = -1;
      if(fSaveResults){
	*outputSummary = cs;
	// std::cout << outputSummary->fEventNumber  <<  ", " << cs.fEventNumber << std::endl;
      }
      
      for(const auto& corrPair : cs.fPairs){
	pairIndex++;

	if(allowedPair(corrPair)==false
	   // Channel 46V (counting from 1, i.e. 45V counting from 0) is just bad
	   || (cs.fPol==AnitaPol::kVertical && (corrPair.ant1==45||corrPair.ant2==45))){
	  if(fSaveResults){
	    outputSummary->fPairs.at(pairIndex).dt = -999;
	    outputSummary->fPairs.at(pairIndex).dt_expected = -999;
	  }
	  continue;
	}

	double dtExpected = usefulPat.getDeltaTExpected(corrPair.ant2, corrPair.ant1, phiWave, thetaWave);

	// if(firstOne){
	//   std::cout  << corrPair.ant1 << "\t" << corrPair.ant2 << "\t" << pat.pitch << "\t" << pat.roll << "\t" << dtExpected << std::endl;
	//   firstOne = false;
	// }

	if(fSaveResults){
	  outputSummary->fPairs.at(pairIndex).dt_expected = dtExpected;
	}
	
	double dtMeasured = corrPair.dt;
	if(fApplyParams){
	  fParamManager.applyDelay(dtMeasured, static_cast<AnitaPol::AnitaPol_t>(cs.fPol),
				   corrPair.ant1, corrPair.ant2);
	}


	double ddt = (dtMeasured - dtExpected);

	normalization.at(cs.fPol).at(corrPair.ant1)++;
	normalization.at(cs.fPol).at(corrPair.ant2)++;	
	
	ddts.at(cs.fPol).at(corrPair.ant1) += ddt*ddt;
	ddts.at(cs.fPol).at(corrPair.ant2) += ddt*ddt;
     }
      if(fSaveResults){
	// std::cout << cs.fPairs.size() << "\t" << outputSummary->fPairs.size() << std::endl;
	fOutputTree->Fill();
      }
    }
  }

  
  double retVal = 0;
  for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      retVal += ddts.at(pol).at(ant)/normalization.at(pol).at(ant);
    }
  }

  fakeGeom.restorePhotogrammetryPositionsInAnitaGeomTool(geom);
  

  if(fPrintOnEval){
    std::cout << "Step: " << retVal << " from params: ";
    const int np = fParamManager.N();
    for(int i=0; i < np; i++){
      std::cout << fMin->VariableName(i) << " = " << params[i] << "\n";
    }
    std::cout << std::flush;
  }
  else{
    std::cout << "." << std::flush;
  }

  if(fOutputTree){
    fOutputTree->Write();
    delete fOutputTree;
    fOutputTree = nullptr;

  }
  if(outputSummary){
    delete outputSummary;
    outputSummary = nullptr;    
  }
  
  return retVal;
}






