#include "AcclaimPhaseCenterFitter.h"
#include "AcclaimParameterManager.h"
#include "TChain.h"
#include "TFile.h"
#include "UsefulAdu5Pat.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"





Acclaim::PhaseCenterFit::FakeGeomTool::FakeGeomTool(const AnitaGeomTool* geom){

  if(geom){
    for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
      for(int ant=0; ant < NUM_SEAVEYS; ant++){
	auto photoKey = std::make_pair(pol, ant);
	fPhotoR[photoKey]   = geom->rPhaseCentreFromVerticalHornPhotogrammetry[ant][pol];
	fPhotoZ[photoKey]   = geom->zPhaseCentreFromVerticalHornPhotogrammetry[ant][pol];
	fPhotoPhi[photoKey] = geom->azPhaseCentreFromVerticalHornPhotogrammetry[ant][pol];
      }
    }
  }

}



void Acclaim::PhaseCenterFit::FakeGeomTool::restorePhotogrammetryPositionsInAnitaGeomTool(AnitaGeomTool* geom) const {
  for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      auto key = std::make_pair(pol, ant);
      geom->rPhaseCentreFromVerticalHornPhotogrammetry[ant][pol]  = fPhotoR.at(key);
      geom->zPhaseCentreFromVerticalHornPhotogrammetry[ant][pol]  = fPhotoZ.at(key);
      geom->azPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = fPhotoPhi.at(key);
    }
  }
}

void Acclaim::PhaseCenterFit::FakeGeomTool::overwritePhotogrammetryPositionsInAnitaGeomTool(AnitaGeomTool* geom) const {
  for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      geom->rPhaseCentreFromVerticalHornPhotogrammetry[ant][pol]  = getAntR(ant, pol);
      geom->zPhaseCentreFromVerticalHornPhotogrammetry[ant][pol]  = getAntZ(ant, pol);
      geom->azPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = getAntPhi(ant, pol);
    }
  }
}





Acclaim::PhaseCenterFit::Minimizer::Minimizer(const char* corrTreeFiles){

  fChain = std::make_shared<TChain>("corrTree");
  if(corrTreeFiles){
    fChain->Add(corrTreeFiles);
  }
  readInSummaries();
  makeFunctors();

  for(auto& a : fNormalization){a.fill(0);}
  for(auto& a : fDdts){a.fill(0);}  
  
  fMin = std::shared_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2"));
  
  fMin->SetMaxFunctionCalls(1e5); // for Minuit/Minuit2
  fMin->SetTolerance(0.01);
  fMin->SetPrintLevel(0);  

}


void Acclaim::PhaseCenterFit::Minimizer::setParameterSpace(ParameterSpace ps){

  fParamSpace = ps;
  fMin->SetFunction(fFuncs[ps]);
  switch(fParamSpace){
  case ParameterSpace::None:
    {
      fInputs.resize(0, 0);
    }
    break;
    
  case ParameterSpace::PitchRollHeading:
    {
      fInputs.resize(3, 0);
      fMin->SetLimitedVariable(0, "pitch", 0, 0.01, -1, 1);
      fMin->SetLimitedVariable(1, "roll", 0, 0.01, -1, 1);
      fMin->SetLimitedVariable(2, "heading_offset", 0, 0.01, -5, 5);
    }
    break;
      
  case ParameterSpace::RingR:
    {
      auto geom = AnitaGeomTool::Instance();
      fInputs.resize(4); // 12);

      fInputs.at(0) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[1][0];
      fInputs.at(1) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[0][0];
      fInputs.at(2) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[NUM_PHI][0];
      fInputs.at(3) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[2*NUM_PHI][0];
      
      fMin->SetVariable(0, "TopHigh Ring R", fInputs.at(0), 0.01); 
      fMin->SetVariable(1, "TopLow Ring R", fInputs.at(1), 0.01);
      fMin->SetVariable(2, "Middle Ring R", fInputs.at(2), 0.01);
      fMin->SetVariable(3, "Bottom Ring R", fInputs.at(3), 0.01);
    }
    break;

  case ParameterSpace::RingPhi:
    {
      auto geom = AnitaGeomTool::Instance();
      fInputs.resize(4); // 12);

      fInputs.at(0) = geom->azPhaseCentreFromVerticalHornPhotogrammetry[1][0] - TMath::TwoPi()/NUM_PHI;
      fInputs.at(1) = geom->azPhaseCentreFromVerticalHornPhotogrammetry[0][0];
      fInputs.at(2) = geom->azPhaseCentreFromVerticalHornPhotogrammetry[NUM_PHI][0];
      fInputs.at(3) = geom->azPhaseCentreFromVerticalHornPhotogrammetry[2*NUM_PHI][0];

      const double lim = 5*TMath::DegToRad();
      fMin->SetLimitedVariable(0, "TopHigh Ring Phi", fInputs.at(0), 0.01, -lim, lim); 
      fMin->SetLimitedVariable(1, "TopLow Ring Phi", fInputs.at(1), 0.01, -lim, lim);
      fMin->SetLimitedVariable(2, "Middle Ring Phi", fInputs.at(2), 0.01, -lim, lim);
      fMin->SetLimitedVariable(3, "Bottom Ring Phi", fInputs.at(3), 0.01, -lim, lim);
    }
    break;
    
  case ParameterSpace::RingZ:
    {
      auto geom = AnitaGeomTool::Instance();
      fInputs.resize(4);
      fInputs.at(0) = geom->zPhaseCentreFromVerticalHornPhotogrammetry[1][0];
      fInputs.at(1) = geom->zPhaseCentreFromVerticalHornPhotogrammetry[0][0];
      fInputs.at(2) = geom->zPhaseCentreFromVerticalHornPhotogrammetry[NUM_PHI][0];
      fInputs.at(3) = geom->zPhaseCentreFromVerticalHornPhotogrammetry[2*NUM_PHI][0];
      fMin->SetVariable(0, "TopHigh Ring Z", fInputs.at(0), 0.01);
      fMin->SetVariable(1, "TopLow Ring Z", fInputs.at(1), 0.01);
      fMin->SetVariable(2, "Middle ring Z",fInputs.at(2), 0.01);
      fMin->SetVariable(3, "Bottom ring Z",fInputs.at(3), 0.01);
    }
    break;
  case ParameterSpace::RingPhiRZ:
    {
      auto geom = AnitaGeomTool::Instance();
      fInputs.resize(12);

      fInputs.at(0) = geom->azPhaseCentreFromVerticalHornPhotogrammetry[1][0] - TMath::TwoPi()/NUM_PHI;
      fInputs.at(1) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[1][0];
      fInputs.at(2) = geom->zPhaseCentreFromVerticalHornPhotogrammetry[1][0];
      
      fInputs.at(3) = geom->azPhaseCentreFromVerticalHornPhotogrammetry[0][0];
      fInputs.at(4) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[0][0];
      fInputs.at(5) = geom->zPhaseCentreFromVerticalHornPhotogrammetry[0][0];
      
      fInputs.at(6) = geom->azPhaseCentreFromVerticalHornPhotogrammetry[NUM_PHI][0];
      fInputs.at(7) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[NUM_PHI][0];
      fInputs.at(8) = geom->zPhaseCentreFromVerticalHornPhotogrammetry[NUM_PHI][0];
      
      fInputs.at(9)  = geom->azPhaseCentreFromVerticalHornPhotogrammetry[2*NUM_PHI][0];
      fInputs.at(10) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[2*NUM_PHI][0];
      fInputs.at(11) = geom->zPhaseCentreFromVerticalHornPhotogrammetry[2*NUM_PHI][0];

      const double lim = 5*TMath::DegToRad();
      fMin->SetLimitedVariable(0, "TopHigh Ring Phi", fInputs.at(0), 0.01, -lim, lim);
      fMin->SetVariable(1, "TopHigh Ring R", fInputs.at(1), 0.01);
      fMin->SetVariable(2, "TopHigh Ring Z", fInputs.at(2), 0.01);

      fMin->SetLimitedVariable(3, "TopLow Ring Phi", fInputs.at(3), 0.01, -lim, lim);
      fMin->SetVariable(4, "TopLow Ring R", fInputs.at(4), 0.01);
      fMin->SetVariable(5, "TopLow Ring Z", fInputs.at(5), 0.01);      

      fMin->SetLimitedVariable(6, "Middle Ring Phi", fInputs.at(6), 0.01, -lim, lim);
      fMin->SetVariable(7, "Middle Ring R", fInputs.at(7), 0.01);
      fMin->SetVariable(8, "Middle Ring Z", fInputs.at(8), 0.01);
      
      fMin->SetLimitedVariable(9, "Bottom Ring Phi", fInputs.at(9), 0.01, -lim, lim);
      fMin->SetVariable(10, "Bottom Ring R", fInputs.at(10), 0.01);
      fMin->SetVariable(11, "Bottom Ring Z", fInputs.at(11), 0.01);      
    }
    break;
  case ParameterSpace::RingEllipse:
    {
      auto geom = AnitaGeomTool::Instance();

      const int n = EllipseParams::N();
      const int numPhysicalRings = 4;
      fInputs.resize(n*numPhysicalRings);
      
      int ringInd = 0;
      for(int ant : {1, 0,  NUM_PHI, 2*NUM_PHI}){
	fInputs.at(ringInd*n + 0) = 0; // x0
	fInputs.at(ringInd*n + 1) = 0; // y0
	fInputs.at(ringInd*n + 2) = 0; // alpha
	fInputs.at(ringInd*n + 3) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[ant][0]; // Ra
	fInputs.at(ringInd*n + 4) = 0; // eccentricity
	fInputs.at(ringInd*n + 5) = geom->zPhaseCentreFromVerticalHornPhotogrammetry[ant][0]; // z
	ringInd++;
      }

      // const double lim = 5*TMath::DegToRad();

      int varInd=0;
      for(int ringInd=0; ringInd < numPhysicalRings; ringInd++){
	std::stringstream ss;
	ss << static_cast<PhysicalRing>(ringInd);
	for(int p=0; p < EllipseParams::N(); p++){	  
	  TString parName = ss.str();
	  parName += TString::Format("_%s", EllipseParams::name(p));
	  if(p==4){
	    fMin->SetLimitedVariable(varInd, parName.Data(), fInputs.at(varInd), 0.01,  0,  1);	  
	  }
	  else if(p==2){
	    fMin->SetLimitedVariable(varInd, parName.Data(), fInputs.at(varInd), 0.01,  -TMath::Pi(),  TMath::Pi());
	  }
	  else{
	    fMin->SetVariable(varInd, parName.Data(), fInputs.at(varInd), 0.01);	  
	  }
	  varInd++;
	}
      }
    }
    break;
  case ParameterSpace::ExtraDeltaT:
    fInputs.resize(NUM_SEAVEYS*AnitaPol::kNotAPol, 0);
    for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
      const char* polName = pol == AnitaPol::kHorizontal ? "H" : "V";
      for(int ant=0; ant< NUM_SEAVEYS; ant++){
	const int varInd = NUM_SEAVEYS*pol + ant;
	fMin->SetVariable(varInd, TString::Format("%d%s", ant, polName).Data(), fInputs.at(varInd), 0.01);
	if(ant==0 || (AnitaPol::kVertical && ant==45)){
	  fMin->FixVariable(varInd);
	}
	
	
      }
    }
    break;
  }
}

void Acclaim::PhaseCenterFit::Minimizer::makeFunctors(){
  fFuncs[ParameterSpace::None]             = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFit::Minimizer::eval, 0);
  fFuncs[ParameterSpace::PitchRollHeading] = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFit::Minimizer::eval, 3);
  fFuncs[ParameterSpace::RingR]            = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFit::Minimizer::eval, 4);
  fFuncs[ParameterSpace::RingPhi]          = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFit::Minimizer::eval, 4);
  fFuncs[ParameterSpace::RingZ]            = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFit::Minimizer::eval, 4);
  fFuncs[ParameterSpace::RingPhiRZ]        = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFit::Minimizer::eval, 12);
  fFuncs[ParameterSpace::RingEllipse]      = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFit::Minimizer::eval, EllipseParams::N()*4);
  fFuncs[ParameterSpace::ExtraDeltaT]      = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFit::Minimizer::eval, NUM_SEAVEYS*AnitaPol::kNotAPol);
}


TH2D* Acclaim::PhaseCenterFit::Minimizer::makeDdtHist(AnitaPol::AnitaPol_t pol, const TString& name, const char* title){
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


const std::vector<double>& Acclaim::PhaseCenterFit::Minimizer::fit(AnitaPol::AnitaPol_t pol, const char* outFileName){

  std::cout << "Starting fit!" << std::endl;  
  fFitPol = pol;

  fApplyParams = false;
  fInitial = eval(&fInputs[0]); // get the start value.
  fApplyParams = true;

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

  fSaveResults = true;
  eval(&fResults[0]);
  fSaveResults = false;
  
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

void Acclaim::PhaseCenterFit::Minimizer::readInSummaries(){

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


void Acclaim::PhaseCenterFit::Minimizer::printResults() const {

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


double Acclaim::PhaseCenterFit::Minimizer::eval(const double* params){
 
  auto geom = AnitaGeomTool::Instance();
  geom->usePhotogrammetryNumbers(true);

  FakeGeomTool fakeGeom(geom);
  fakeGeom.overwritePhotogrammetryPositionsInAnitaGeomTool(geom);
  
  fParamManager.update(params);
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

  // reset delta delta T counters
  for(auto& a : fNormalization){a.fill(0);}
  for(auto& a : fDdts){a.fill(0);}
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

	// this channel is just bad...
	if(cs.fPol==AnitaPol::kVertical && (corrPair.ant1==45||corrPair.ant2==45)){
	  if(fSaveResults){
	    outputSummary->fPairs.at(pairIndex).dt = -999;
	    outputSummary->fPairs.at(pairIndex).dt_expected = -999;
	  }
	  continue;
	}

	double dtExpected = usefulPat.getDeltaTExpected(corrPair.ant2, corrPair.ant1, phiWave, thetaWave);
	if(fSaveResults){
	  outputSummary->fPairs.at(pairIndex).dt_expected = dtExpected;
	}
	
	double dtMeasured = corrPair.dt;
	if(fApplyParams && fParamSpace==ParameterSpace::ExtraDeltaT){
	  const int polOffset = cs.fPol*NUM_SEAVEYS;
	  dtMeasured += params[polOffset + corrPair.ant1];
	  dtMeasured -= params[polOffset + corrPair.ant2];
	}


	double ddt = (dtMeasured - dtExpected);
	fNormalization.at(cs.fPol).at(corrPair.ant1)++;
	fNormalization.at(cs.fPol).at(corrPair.ant2)++;	
	
	fDdts.at(cs.fPol).at(corrPair.ant1) += ddt*ddt;
	fDdts.at(cs.fPol).at(corrPair.ant2) += ddt*ddt;

	// std::cout << normalization.at(cs.fPol).at(corrPair.ant1) << "\t"  << ddts.at(cs.fPol).at(corrPair.ant1) << "\t" << ddt << "\n";

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






