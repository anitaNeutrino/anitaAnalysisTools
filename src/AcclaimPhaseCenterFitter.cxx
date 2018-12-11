#include "AcclaimPhaseCenterFitter.h"
#include "TChain.h"
#include "TFile.h"
#include "UsefulAdu5Pat.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

std::ostream& operator<<(std::ostream& os, const Acclaim::PhaseCenterFitter::PhysicalRing& r){
  switch(r){
  case Acclaim::PhaseCenterFitter::PhysicalRing::TopHigh:
    os << "PhysicalRing::TopHigh";
    return os;
  case Acclaim::PhaseCenterFitter::PhysicalRing::TopLow:
    os << "PhysicalRing::TopLow";
    return os;
  case Acclaim::PhaseCenterFitter::PhysicalRing::Middle:
    os << "PhysicalRing::Middle";
    return os;
  case Acclaim::PhaseCenterFitter::PhysicalRing::Bottom:
    os << "PhysicalRing::Bottom";
    return os;    
  }
  return os;
}


Acclaim::FakeGeomTool::FakeGeomTool(const AnitaGeomTool* geom){

  const int numPhysRing = 4;
  std::array<double, numPhysRing> regR_ring {{1.15333, 1.31832, 2.38709, 2.38694}};
  std::array<double, numPhysRing> regZ_ring {{2.73599, 3.70029, 0.0343191, -1.06195}};
  std::array<double, numPhysRing> regPhi_ring {{-0.00558104, -0.00626989, -0.00114436, -0.000663389}};

  regPhi_ring[0] += geom->azPhaseCentreFromVerticalHornPhotogrammetry[0][0];
  regPhi_ring[1] += geom->azPhaseCentreFromVerticalHornPhotogrammetry[1][0];
  regPhi_ring[2] += geom->azPhaseCentreFromVerticalHornPhotogrammetry[NUM_PHI][0];
  regPhi_ring[3] += geom->azPhaseCentreFromVerticalHornPhotogrammetry[2*NUM_PHI][0];  
  
  if(geom){
    for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
      for(int ant=0; ant < NUM_SEAVEYS; ant++){
	auto photoKey = std::make_pair(pol, ant);
	fPhotoR[photoKey]   = geom->rPhaseCentreFromVerticalHornPhotogrammetry[ant][pol];
	fPhotoZ[photoKey]   = geom->zPhaseCentreFromVerticalHornPhotogrammetry[ant][pol];
	fPhotoPhi[photoKey] = geom->azPhaseCentreFromVerticalHornPhotogrammetry[ant][pol];

	auto ring = PhaseCenterFitter::antToPhysicalRing(ant);
	int ringInt = static_cast<int>(ring);	
	int phiSector = ant%NUM_PHI;

	double phi = regPhi_ring.at(ringInt) + (phiSector-2)*TMath::TwoPi()/NUM_PHI;
	auto regKey = std::make_pair(pol, ring);

	fRegularPhi[regKey] = phi;
	fRegularZ[regKey]   = regZ_ring.at(ringInt);
	fRegularR[regKey]   = regR_ring.at(ringInt);
      }
    }
  }

  // print();
}

void Acclaim::FakeGeomTool::print() const {
  std::cout << "Ring\tR\tPhi\tZ" << std::endl;
  for(auto ring : {PhaseCenterFitter::PhysicalRing::TopHigh,
		     PhaseCenterFitter::PhysicalRing::TopLow,
		     PhaseCenterFitter::PhysicalRing::Middle,
		     PhaseCenterFitter::PhysicalRing::Bottom}){
    auto pol = AnitaPol::kVertical;
    auto key = std::make_pair(pol,  ring);
    std::cout << ring << "\t" << fRegularR.at(key) << "\t" << fRegularPhi.at(key) << "\t" << fRegularZ.at(key) << "\n";


    int ant;
    switch(ring){
    case PhaseCenterFitter::PhysicalRing::TopHigh: ant = 1; break;
    case PhaseCenterFitter::PhysicalRing::TopLow:  ant = 0; break;
    case PhaseCenterFitter::PhysicalRing::Middle: ant = NUM_PHI; break;
    case PhaseCenterFitter::PhysicalRing::Bottom: ant = 2*NUM_PHI; break;
    }
    
    auto k2 = std::make_pair(pol, ant);
    std::cout << ring << "\t" << fPhotoR.at(k2) << "\t" << fPhotoPhi.at(k2) << "\t" << fPhotoZ.at(k2) << "\n";
  }


  std::cout << "\n" << std::endl;
}


void Acclaim::FakeGeomTool::restorePhotogrammetryPositionsInAnitaGeomTool(AnitaGeomTool* geom) const {
  for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      auto key = std::make_pair(pol, ant);
      geom->rPhaseCentreFromVerticalHornPhotogrammetry[ant][pol]  = fPhotoR.at(key);
      geom->zPhaseCentreFromVerticalHornPhotogrammetry[ant][pol]  = fPhotoZ.at(key);
      geom->azPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = fPhotoPhi.at(key);
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
    }
    break;
      
  case ParameterSpace::RingR:
    {
      auto geom = AnitaGeomTool::Instance();
      int var=0;
      fInputs.resize(4); // 12);
      for(auto ring : {PhysicalRing::TopHigh, PhysicalRing::TopLow, PhysicalRing::Middle, PhysicalRing::Bottom}){
	switch(ring){
	case PhysicalRing::TopHigh:
	  fInputs.at(var) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[1][0];
	  fMin->SetVariable(var, "TopHigh Ring R", fInputs.at(var), 0.01); 
	  break;
	case PhysicalRing::TopLow:	  
	  fInputs.at(var) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[0][0];
	  fMin->SetVariable(var, "TopLow Ring R", fInputs.at(var), 0.01);
	  break;
	case PhysicalRing::Middle:
	  fInputs.at(var) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[NUM_PHI][0];
	  fMin->SetVariable(var, "Middle Ring R", fInputs.at(var), 0.01);
	  break;
	case PhysicalRing::Bottom:
	  fInputs.at(var) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[2*NUM_PHI][0];
	  fMin->SetVariable(var, "Bottom Ring R", fInputs.at(var), 0.01);
	  break;
	}
 	var++;
	// fInputs.at(var) = 0;
	// var++;
	// fInputs.at(var) = 0;
	// var++;
      }
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
  }
}

void Acclaim::PhaseCenterFitter::makeFunctors(){
  fFuncs[ParameterSpace::PitchRollHeading] = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFitter::eval, 3);
  fFuncs[ParameterSpace::RingR]            = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFitter::eval, 4);
  fFuncs[ParameterSpace::RingZ]            = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFitter::eval, 4);
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

  if(fParamSpace==ParameterSpace::RingR){
    for(auto pol : {AnitaPol::kHorizontal,  AnitaPol::kVertical}){
      for(int ant=0; ant < NUM_SEAVEYS; ant++){
	auto ring = antToPhysicalRing(ant);
	int ringInt = static_cast<int>(ring);
	geom->rPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = params[ringInt];
      }
    }
  }

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


