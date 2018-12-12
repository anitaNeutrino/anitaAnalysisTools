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

  // print();
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
      geom->azPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = getAntPhi(ant, pol);
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
      fInputs.resize(28);

      fInputs.at(0) = geom->azPhaseCentreFromVerticalHornPhotogrammetry[1][0] - TMath::TwoPi()/NUM_PHI;
      fInputs.at(1) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[1][0]; // Ra
      fInputs.at(2) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[1][0]; // Rb
      fInputs.at(3) = geom->zPhaseCentreFromVerticalHornPhotogrammetry[1][0]; // z
      fInputs.at(4) = 0; // x0
      fInputs.at(5) = 0; // y0
      fInputs.at(6) = 0; // angle of major axis
      
      
      fInputs.at(7) = geom->azPhaseCentreFromVerticalHornPhotogrammetry[0][0];
      fInputs.at(8) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[0][0];
      fInputs.at(9) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[0][0];      
      fInputs.at(10) = geom->zPhaseCentreFromVerticalHornPhotogrammetry[0][0];
      fInputs.at(11) = 0; // x0
      fInputs.at(12) = 0; // y0
      fInputs.at(13) = 0; // angle of major axis      
      
      fInputs.at(14) = geom->azPhaseCentreFromVerticalHornPhotogrammetry[NUM_PHI][0];
      fInputs.at(15) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[NUM_PHI][0];
      fInputs.at(16) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[NUM_PHI][0];      
      fInputs.at(17) = geom->zPhaseCentreFromVerticalHornPhotogrammetry[NUM_PHI][0];
      fInputs.at(18) = 0; // x0
      fInputs.at(19) = 0; // y0
      fInputs.at(20) = 0; // angle of major axis
      
      fInputs.at(21)  = geom->azPhaseCentreFromVerticalHornPhotogrammetry[2*NUM_PHI][0];
      fInputs.at(22) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[2*NUM_PHI][0];
      fInputs.at(23) = geom->rPhaseCentreFromVerticalHornPhotogrammetry[2*NUM_PHI][0];
      fInputs.at(24) = geom->zPhaseCentreFromVerticalHornPhotogrammetry[2*NUM_PHI][0];
      fInputs.at(25) = 0; // x0
      fInputs.at(26) = 0; // y0
      fInputs.at(27) = 0; // angle of major axis

      const double lim = 5*TMath::DegToRad();
      fMin->SetLimitedVariable(0, "TopHigh Ring Phi", fInputs.at(0), 0.01, -lim, lim);
      fMin->SetVariable(1, "TopHigh Ring Ra",         fInputs.at(1), 0.01);
      fMin->SetVariable(2, "TopHigh Ring Rb",         fInputs.at(2), 0.01);
      fMin->SetVariable(3, "TopHigh Ring Z",          fInputs.at(3), 0.01);
      fMin->SetVariable(4, "TopHigh x0",              fInputs.at(4), 0.01);
      fMin->SetVariable(5, "TopHigh y0",              fInputs.at(5), 0.01);
      fMin->SetVariable(6, "TopHigh alpha",           fInputs.at(6), 0.01);

      fMin->SetLimitedVariable(7, "TopLow Ring Phi", fInputs.at(7), 0.01, -lim, lim);
      fMin->SetVariable(8, "TopLow Ring Ra",         fInputs.at(8), 0.01);
      fMin->SetVariable(9, "TopLow Ring Rb",         fInputs.at(9), 0.01);
      fMin->SetVariable(10, "TopLow Ring Z",         fInputs.at(10), 0.01);
      fMin->SetVariable(11, "TopLow x0",             fInputs.at(11), 0.01);
      fMin->SetVariable(12, "TopLow y0",             fInputs.at(12), 0.01);
      fMin->SetVariable(13, "TopLow alpha",          fInputs.at(13), 0.01);
      
      fMin->SetLimitedVariable(14, "Middle Ring Phi", fInputs.at(14), 0.01, -lim, lim);
      fMin->SetVariable(15, "Middle Ring Ra",         fInputs.at(15), 0.01);
      fMin->SetVariable(16, "Middle Ring Rb",         fInputs.at(16), 0.01);
      fMin->SetVariable(17, "Middle Ring Z",          fInputs.at(17), 0.01);
      fMin->SetVariable(18, "Middle x0",              fInputs.at(18), 0.01);
      fMin->SetVariable(19, "Middle y0",              fInputs.at(19), 0.01);
      fMin->SetVariable(20, "Middle alpha",           fInputs.at(20), 0.01);

      fMin->SetLimitedVariable(21, "Bottom Ring Phi", fInputs.at(21), 0.01, -lim, lim);
      fMin->SetVariable(22, "Bottom Ring Ra",         fInputs.at(22), 0.01);
      fMin->SetVariable(23, "Bottom Ring Rb",         fInputs.at(23), 0.01);
      fMin->SetVariable(24, "Bottom Ring Z",          fInputs.at(24), 0.01);
      fMin->SetVariable(25, "Bottom x0",              fInputs.at(25), 0.01);
      fMin->SetVariable(26, "Bottom y0",              fInputs.at(26), 0.01);
      fMin->SetVariable(27, "Bottom alpha",           fInputs.at(27), 0.01);

    }
    break;
  }  
}

void Acclaim::PhaseCenterFitter::makeFunctors(){
  fFuncs[ParameterSpace::PitchRollHeading] = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFitter::eval, 3);
  fFuncs[ParameterSpace::RingR]            = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFitter::eval, 4);
  fFuncs[ParameterSpace::RingPhi]          = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFitter::eval, 4);
  fFuncs[ParameterSpace::RingZ]            = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFitter::eval, 4);
  fFuncs[ParameterSpace::RingPhiRZ]        = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFitter::eval, 12);
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

  if(fApplyParams){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      int phiSector = ant%NUM_PHI;
      for(auto pol : {AnitaPol::kHorizontal,  AnitaPol::kVertical}){
	auto ring = antToPhysicalRing(ant);
	int ringInt = static_cast<int>(ring);

	if(fParamSpace==ParameterSpace::RingR){
	  geom->rPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = params[ringInt];
	}
	else if(fParamSpace==ParameterSpace::RingPhi){
	  double phi = params[ringInt] + phiSector*TMath::DegToRad()*22.5;
	  geom->azPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = phi;
	}

	else if(fParamSpace==ParameterSpace::RingPhiRZ){
	  const int paramsPerRing = 3;
	  double phi = params[ringInt*paramsPerRing] + phiSector*TMath::DegToRad()*22.5;
	  geom->azPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = phi;
	  geom->rPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = params[ringInt*paramsPerRing+1];	
	  geom->zPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = params[ringInt*paramsPerRing+2];
	}
	else if(fParamSpace==ParameterSpace::RingEllipse){
	  const int paramsPerRing = 7;
	  double Ra = params[ringInt*paramsPerRing + 1];
	  double Rb = params[ringInt*paramsPerRing + 2];
	  double x0 = params[ringInt*paramsPerRing + 4];
	  double y0 = params[ringInt*paramsPerRing + 5];
	  
	  double alpha = params[ringInt*paramsPerRing + 6];
	  double phi = params[ringInt*paramsPerRing] + phiSector*TMath::DegToRad()*22.5 - alpha;

	  double x = Ra*cos(phi);
	  double y = Rb*sin(phi);
	  TVector2 v(x, y);
	  v = v.Rotate(alpha);

	  
	  
	  geom->azPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = phi;
	  geom->rPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = params[ringInt*paramsPerRing+1];	
	  geom->zPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = params[ringInt*paramsPerRing+2];
	}
	
      }
    }
  }

  // reset delta delta T counters
  for(auto& a : fNormalization){a.fill(0);}
  for(auto& a : fDdts){a.fill(0);}
  for(auto& cs : fSummaries){

    if(fFitPol==AnitaPol::kNotAPol || cs.fPol==fFitPol){

      Adu5Pat pat = cs.fPat;
      if(fApplyParams){
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






void Acclaim::PhaseCenterFitter::EllipseParams::phiToEllipseXY(double phi, double& x, double& y){
  double t = tFromPhi(phi);
  x = Ra*TMath::Cos(t);
  y = Rb*TMath::Sin(t);
  TVector2 v(x, y);
  v = v.Rotate(alpha);
  x = v.X() + x0;
  y = v.Y() + y0;
}

double Acclaim::PhaseCenterFitter::EllipseParams::tFromPhi(double phi){
  // here we find the parametric angle t from the angle phi.
  double localPhi = phi - alpha; // account for rotation of ellipse
  double tan_t = (Ra/Rb)*TMath::Tan(localPhi);
  double t = TMath::ATan(tan_t);

  // get correct quadrant, I hope...
  if(localPhi > 0.5*TMath::Pi() && localPhi <= 1.5*TMath::Pi()){
    t += TMath::Pi();
  }
  // std::cout << phi*TMath::RadToDeg() << "\t" << localPhi*TMath::RadToDeg() << "\t" << t*TMath::RadToDeg() << std::endl;
  
  return t;
}
