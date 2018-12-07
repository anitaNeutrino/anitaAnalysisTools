#include "AcclaimPhaseCenterFitter.h"
#include "TChain.h"
#include "UsefulAdu5Pat.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

Acclaim::PhaseCenterFitter::PhaseCenterFitter(const char* corrTreeFiles){

  fChain = new TChain("corrTree");
  if(corrTreeFiles){
    fChain->Add(corrTreeFiles);
  }
  readInSummaries();
  makeFunctors();

  fMin = ROOT::Math::Factory::CreateMinimizer("Minuit2");
  
  fMin->SetMaxFunctionCalls(1e5); // for Minuit/Minuit2
  fMin->SetTolerance(0.01);
  fMin->SetPrintLevel(0);

}


void Acclaim::PhaseCenterFitter::SetFitParameterSpace(ParameterSpace ps){

  fParamSpace = ps;
  fMin->SetFunction(fFuncs[ps]);
  switch(fParamSpace){
  case ParameterSpace::PitchRoll:
    fMin->SetLimitedVariable(0, "pitch", 0, 0.01, -1, 1);
    fMin->SetLimitedVariable(1, "roll", 0, 0.01, -1, 1);
    break;
  }  
}

void Acclaim::PhaseCenterFitter::makeFunctors(){

  fFuncs[ParameterSpace::PitchRoll] = ROOT::Math::Functor(this, &Acclaim::PhaseCenterFitter::eval, 2);
  
}



Acclaim::PhaseCenterFitter::~PhaseCenterFitter(){
  delete fChain;
  delete fMin;
}


const std::vector<double>& Acclaim::PhaseCenterFitter::fit(){

  fMin->Minimize();
  fResults.clear();
  const int nDim = fFuncs[fParamSpace].NDim();
  fResults.reserve(nDim);
  
  for(int i=0; i < nDim; i++){
    fResults.push_back(fMin->X()[i]);
  }
  return fResults;
}

void Acclaim::PhaseCenterFitter::readInSummaries(){

  fChain->SetBranchAddress("correlationSummary", &fSum);
   Long64_t n = fChain->GetEntries();
  fSummaries.reserve(n);
  for(Long64_t entry=0; entry < n; entry++){
    fChain->GetEntry(entry);

    // std::cout << fSummaries.size() << std::endl;
    fSummaries.emplace_back(*fSum);
  }

  if(fSummaries.size()==0){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << "No summaries read in!" << std::endl;
  }
}


void Acclaim::PhaseCenterFitter::printResults() const {
  
  std::cout << "Results = ";
  const int np = fFuncs.at(fParamSpace).NDim();
  for(int i=0; i < np; i++){
    std::cout << fMin->VariableName(i) << " = " << fResults.at(i);
    if(i < np - 1){
      std::cout << ", ";
    }
  }    
  
}


double Acclaim::PhaseCenterFitter::eval(const double* params){    

  // const double pitch = params[0];
  // const double roll = params[1];
 
  auto geom = AnitaGeomTool::Instance();
  geom->usePhotogrammetryNumbers(true);
  
  double retVal = 0;
  for(auto& cs : fSummaries){

    if(fParamSpace==ParameterSpace::PitchRoll){    
      cs.fPat.pitch = params[0];
      cs.fPat.roll = params[1];
    }

    double thetaWave,  phiWave;
    
    UsefulAdu5Pat usefulPat(&cs.fPat, false);
    usefulPat.getThetaAndPhiWaveWaisDivide(thetaWave, phiWave);

    double eventResidual = 0;
    int goodPairs = 0;
    for(const auto& corrPair : cs.fPairs){

      // this channel is just bad...
      if(cs.fPol==AnitaPol::kVertical && (corrPair.ant1==45||corrPair.ant2==45)){
	continue;
      }

      if(corrPair.correlation > fCorrelationThreshold){
	double dtExpected = usefulPat.getDeltaTExpected(corrPair.ant2, corrPair.ant1, phiWave, thetaWave);
	double dtMeasured = corrPair.dt;
	double ddt = (dtMeasured - dtExpected);

	eventResidual += ddt*ddt;
	goodPairs++;
      }
    }
    eventResidual = goodPairs > 0 ? eventResidual/goodPairs : 0;

    retVal += eventResidual;
    // p.inc(entry,  n);
  }

  std::cout << "Result: " << retVal << " from params: ";
  const int np = fFuncs[fParamSpace].NDim();
  for(int i=0; i < np; i++){
    std::cout << fMin->VariableName(i) << " = " << params[i];
    if(i < np - 1){
      std::cout << ", ";
    }
  }
  std::cout << std::endl;

  return retVal;
}


