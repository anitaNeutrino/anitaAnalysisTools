#include "AcclaimPhaseCenterFitter.h"
#include "TChain.h"
#include "TFile.h"
#include "UsefulAdu5Pat.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"


Acclaim::FakeGeomTool::FakeGeomTool(const AnitaGeomTool* geom){

  if(geom){
    for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
      for(int ant=0; ant < NUM_SEAVEYS; ant++){
	auto key = std::make_pair(pol, ant);
	fDefaultR[key] = geom->rPhaseCentreFromVerticalHornPhotogrammetry[ant][pol];
	fDefaultZ[key] = geom->zPhaseCentreFromVerticalHornPhotogrammetry[ant][pol];
	fDefaultPhi[key] = geom->azPhaseCentreFromVerticalHornPhotogrammetry[ant][pol];
      }
    }
  }
}


void Acclaim::FakeGeomTool::restorePhotogrammetryPositionsInAnitaGeomTool(AnitaGeomTool* geom) const {
  for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      auto key = std::make_pair(pol, ant);
      geom->rPhaseCentreFromVerticalHornPhotogrammetry[ant][pol]  = fDefaultR.at(key);
      geom->zPhaseCentreFromVerticalHornPhotogrammetry[ant][pol]  = fDefaultZ.at(key);
      geom->azPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = fDefaultPhi.at(key);
    }
  }
}

void Acclaim::FakeGeomTool::overwritePhotogrammetryPositionsInAnitaGeomTool(AnitaGeomTool* geom) const {
  for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      geom->rPhaseCentreFromVerticalHornPhotogrammetry[ant][pol]  = getAntR(ant, pol);
      geom->zPhaseCentreFromVerticalHornPhotogrammetry[ant][pol]  = getAntZ(ant, pol);
      geom->azPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = getAntPhiPositionRelToAftFore(ant, pol);
    }
  }
}



Acclaim::PhaseCenterFitter::PhaseCenterFitter(const char* corrTreeFiles){

  fChain = new TChain("corrTree");
  if(corrTreeFiles){
    fChain->Add(corrTreeFiles);
  }
  readInSummaries();
  makeFunctors();

  for(auto& a : fNormalization){a.fill(0);}
  for(auto& a : fDdts){a.fill(0);}  
  
  fMin = ROOT::Math::Factory::CreateMinimizer("Minuit2");
  
  fMin->SetMaxFunctionCalls(1e5); // for Minuit/Minuit2
  fMin->SetTolerance(0.01);
  fMin->SetPrintLevel(0);  

}


void Acclaim::PhaseCenterFitter::setFitParameterSpace(ParameterSpace ps){

  fParamSpace = ps;
  fMin->SetFunction(fFuncs[ps]);
  switch(fParamSpace){
  case ParameterSpace::PitchRollHeading:
    {
      fInputs.resize(3, 0);
      fMin->SetLimitedVariable(0, "pitch", 0, 0.01, -1, 1);
      fMin->SetLimitedVariable(1, "roll", 0, 0.01, -1, 1);
      fMin->SetLimitedVariable(2, "heading_offset", 0, 0.01, -5, 5);
      break;
    }
      
  case ParameterSpace::RingR:
    {
      auto geom = AnitaGeomTool::Instance();
      fInputs.resize(4);
      fInputs.at(0) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[0][0];
      fInputs.at(1) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[1][0];
      fInputs.at(2) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[NUM_PHI][0];
      fInputs.at(3) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[2*NUM_PHI][0];
      fMin->SetVariable(0, "Top 1 Ring R", fInputs.at(0), 0.01);
      fMin->SetVariable(1, "Top 2 Ring R", fInputs.at(1), 0.01);
      fMin->SetVariable(2, "Middle ring R",fInputs.at(2), 0.01);
      fMin->SetVariable(3, "Bottom ring R",fInputs.at(3), 0.01);
    }
  case ParameterSpace::RingZ:
    {
      auto geom = AnitaGeomTool::Instance();
      fInputs.resize(4);
      fInputs.at(0) = geom->zPhaseCentreFromVerticalHornPhotogrammetry[0][0];
      fInputs.at(1) = geom->zPhaseCentreFromVerticalHornPhotogrammetry[1][0];
      fInputs.at(2) = geom->zPhaseCentreFromVerticalHornPhotogrammetry[NUM_PHI][0];
      fInputs.at(3) = geom->zPhaseCentreFromVerticalHornPhotogrammetry[2*NUM_PHI][0];
      fMin->SetVariable(0, "Top 1 Ring Z", fInputs.at(0), 0.01);
      fMin->SetVariable(1, "Top 2 Ring Z", fInputs.at(1), 0.01);
      fMin->SetVariable(2, "Middle ring Z",fInputs.at(2), 0.01);
      fMin->SetVariable(3, "Bottom ring Z",fInputs.at(3), 0.01);
    }    
  }
}

void Acclaim::PhaseCenterFitter::makeFunctors(){
  fFuncs[ParameterSpace::PitchRollHeading] = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFitter::eval, 3);
  fFuncs[ParameterSpace::RingR]     = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFitter::eval, 4);
  fFuncs[ParameterSpace::RingZ]     = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFitter::eval, 4);
}


TH2D* Acclaim::PhaseCenterFitter::makeDdtHist(AnitaPol::AnitaPol_t pol, const TString& name, const char* title){
  if(pol==AnitaPol::kNotAPol){
    return nullptr;
  }
  
  TString hName(name);
  TString hTitle(title);
  hTitle += pol == AnitaPol::kVertical ? "VPol" : "HPol";
  hTitle += ";#Phi-sector; Ring";
  TH2D* h = new TH2D(hName, hTitle,
		     NUM_PHI, 0.5, NUM_PHI+0.5,
		     AnitaRing::kNotARing, 0, AnitaRing::kNotARing);
  h->GetYaxis()->SetBinLabel(1, "B");
  h->GetYaxis()->SetBinLabel(2, "M");
  h->GetYaxis()->SetBinLabel(3, "T");
  for(int phi=1; phi <= NUM_PHI; phi++){
    auto l = TString::Format("%d", phi);
    h->GetXaxis()->SetBinLabel(phi, l.Data());
  }
  for(int ant = 0; ant < NUM_SEAVEYS; ant++){
    int ring = ant/NUM_PHI;
    int phi = ant%NUM_PHI;
    double val = fNormalization.at(pol).at(ant) > 0 ? fDdts.at(pol).at(ant)/fNormalization.at(pol).at(ant) : 0;
    h->Fill(phi+1, 2 - ring, val);
  }
  return h;
}


Acclaim::PhaseCenterFitter::~PhaseCenterFitter(){
  delete fChain;
  delete fMin;
}


const std::vector<double>& Acclaim::PhaseCenterFitter::fit(AnitaPol::AnitaPol_t pol, const char* outFileName){

  std::cout << "Starting fit!" << std::endl;  
  fFitPol = pol;
  fInitial = eval(&fInputs[0]); // get the start value.

  TString outputFileName = outFileName ? TString(outFileName) : TString::Format("FitterResults_%d", fFitPol);
  TFile f(outputFileName, "recreate");

  std::vector<AnitaPol::AnitaPol_t> pols;
  if(fFitPol==AnitaPol::kNotAPol){
    pols.push_back(AnitaPol::kHorizontal);
    pols.push_back(AnitaPol::kVertical);
  }
  else{
    pols.push_back(fFitPol);
  }

  for(auto pol : pols){
    TString name = TString::Format("hInitial_%d", pol);
    auto h = makeDdtHist(pol, name);
    h->Write();
    delete h;
  }
  
  fMin->Minimize();
  fResults.clear();
  const int nDim = fFuncs.at(fParamSpace).NDim();
  fResults.reserve(nDim);
  
  for(int i=0; i < nDim; i++){
    fResults.push_back(fMin->X()[i]);
  }
  fMinimum = fMin->MinValue();

  eval(&fResults[0]);
  for(auto pol : pols){
    TString name = TString::Format("hFinal_%d", pol);
    auto h = makeDdtHist(pol, name);
    h->Write();
    delete h;
  }  

  std::cout << std::endl;

  f.Write();
  f.Close();
  
  return fResults;
}

void Acclaim::PhaseCenterFitter::readInSummaries(){

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
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << "No summaries read in!" << std::endl;
  }
}


void Acclaim::PhaseCenterFitter::printResults() const {

  std::cout << "Results:\n";
  std::cout  << "Before the fit, the value was " << fInitial << " with:\n";
  const int np = fFuncs.at(fParamSpace).NDim();
  
  for(int i=0; i < np; i++){
    std::cout << fMin->VariableName(i) << " = " << fInputs.at(i) << "\n";
  }

  std::cout << "\nAfter the fitting, minimum value is " << fMinimum << " with:\n";
  for(int i=0; i < np; i++){
    std::cout << fMin->VariableName(i) << " = " << fResults.at(i) << "\n";
  }
  std::cout << "\n\n" << std::flush;
  
}


double Acclaim::PhaseCenterFitter::eval(const double* params){
 
  auto geom = AnitaGeomTool::Instance();
  geom->usePhotogrammetryNumbers(true);
  FakeGeomTool fakeGeom(geom);
  fakeGeom.overwritePhotogrammetryPositionsInAnitaGeomTool(geom);

  // reset delta delta T counters
  for(auto& a : fNormalization){a.fill(0);}
  for(auto& a : fDdts){a.fill(0);}
  for(auto& cs : fSummaries){

    if(fFitPol==AnitaPol::kNotAPol || cs.fPol==fFitPol){

      Adu5Pat pat = cs.fPat;
      if(fParamSpace==ParameterSpace::PitchRollHeading){
	pat.pitch = params[0];
	pat.roll = params[1];
	pat.heading += params[2];
      }
      else{
	// these are the results of fitting pitch/roll with photogrammetry numbers
	// pitch = -0.199065
	// roll = -0.113464
	// heading_offset = -0.552301
	
	pat.pitch    = -0.199065;
	pat.roll     = -0.113464;
	pat.heading += -0.552301;
      }
    
      double thetaWave, phiWave;
      UsefulAdu5Pat usefulPat(&pat, false);
      usefulPat.getThetaAndPhiWaveWaisDivide(thetaWave, phiWave);

      for(const auto& corrPair : cs.fPairs){

	// this channel is just bad...
	if(cs.fPol==AnitaPol::kVertical && (corrPair.ant1==45||corrPair.ant2==45)){
	  continue;
	}

	double dtExpected = usefulPat.getDeltaTExpected(corrPair.ant2, corrPair.ant1, phiWave, thetaWave);
	double dtMeasured = corrPair.dt;
	double ddt = (dtMeasured - dtExpected);

	fNormalization.at(cs.fPol).at(corrPair.ant1)++;
	fNormalization.at(cs.fPol).at(corrPair.ant2)++;	
	
	fDdts.at(cs.fPol).at(corrPair.ant1) += ddt*ddt;
	fDdts.at(cs.fPol).at(corrPair.ant2) += ddt*ddt;

	// std::cout << normalization.at(cs.fPol).at(corrPair.ant1) << "\t"  << ddts.at(cs.fPol).at(corrPair.ant1) << "\t" << ddt << "\n";
      
      }
    }
  }

  double retVal = 0;
  for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      if(fNormalization.at(pol).at(ant) > 0){
	retVal += fDdts.at(pol).at(ant)/fNormalization.at(pol).at(ant);
      }
    }
  }

  fakeGeom.restorePhotogrammetryPositionsInAnitaGeomTool(geom);
  

  if(fPrintOnEval){
    std::cout << "Step: " << retVal << " from params: ";
    const int np = fFuncs[fParamSpace].NDim();
    for(int i=0; i < np; i++){
      std::cout << fMin->VariableName(i) << " = " << params[i] << "\n";
    }
    std::cout << std::flush;
  }
  else{
    std::cout << "." << std::flush;
  }
  
  return retVal;
}


