#include "AcclaimPhaseCenterFitter.h"
#include "TChain.h"
#include "UsefulAdu5Pat.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

Acclaim::PhaseCenterFitter::PhaseCenterFitter(const char* corrTreeFiles){

  fChain = new TChain("corr Tree");
  if(corrTreeFiles){
    fChain->Add(corrTreeFiles);
  }
  readInSummaries();

  fMin = ROOT::Math::Factory::CreateMinimizer("Minuit2");
  
  fFunc = new ROOT::Math::Functor(this, &Acclaim::PhaseCenterFitter::evalPitchRoll, 2);

  fMin->SetMaxFunctionCalls(1e5); // for Minuit/Minuit2
  fMin->SetTolerance(0.0001);
  fMin->SetPrintLevel(0);
  fMin->SetFunction(*fFunc);

  fMin->SetLimitedVariable(0, "cos(pitch)", 0, 0.001, -1, 1);
  fMin->SetLimitedVariable(1, "cos(roll)", 0, 0.001, -1, 1);
}


Acclaim::PhaseCenterFitter::~PhaseCenterFitter(){
  delete fChain;
  delete fMin;
  delete fFunc;
}


void Acclaim::PhaseCenterFitter::fit(){
  fMin->Minimize();
}

void Acclaim::PhaseCenterFitter::readInSummaries(){

  fChain->SetBranchAddress("correlationSummary", &fSum);
  Long64_t n = fChain->GetEntries();
  fSummaries.reserve(n);
  for(Long64_t entry=0; entry < n; entry++){
    fChain->GetEntry(entry);

    std::cout << fSummaries.size() << std::endl;
    fSummaries.emplace_back(*fSum);
  }
}


double Acclaim::PhaseCenterFitter::evalPitchRoll(const double* params){    

  const double pitch = TMath::RadToDeg()*TMath::ACos(params[0]);
  const double roll = TMath::RadToDeg()*TMath::ACos(params[1]);
  
  double retVal = 0;
  for(auto& cs : fSummaries){

    double thetaWave,  phiWave;
    cs.fPat.pitch = pitch;
    cs.fPat.roll = roll;
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

  std::cout << pitch << "\t" << roll << "\t" << retVal << std::endl;
  return retVal;
}


