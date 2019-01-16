#include "AcclaimPhaseCenterParameters.h"

#include "AnitaGeomTool.h"
#include "Adu5Pat.h"

#include "TMath.h"
#include "Math/Minimizer.h"

#include <sstream>


std::ostream& operator<<(std::ostream& os, const Acclaim::PhaseCenter::PhysicalRing& r){
  switch(r){
  case Acclaim::PhaseCenter::PhysicalRing::TopHigh:
    os << "PhysicalRing::TopHigh";
    return os;
  case Acclaim::PhaseCenter::PhysicalRing::TopLow:
    os << "PhysicalRing::TopLow";
    return os;
  case Acclaim::PhaseCenter::PhysicalRing::Middle:
    os << "PhysicalRing::Middle";
    return os;
  case Acclaim::PhaseCenter::PhysicalRing::Bottom:
    os << "PhysicalRing::Bottom";
    return os;
  }
  return os;
}


Acclaim::PhaseCenter::PhysicalRing Acclaim::PhaseCenter::antToPhysicalRing(int ant){
  int ring = 1 + ant/NUM_PHI;
  if(ant < NUM_PHI && (ant % 2)==0){
    ring -= 1;
  }
  return static_cast<Acclaim::PhaseCenter::PhysicalRing>(ring);
}





Acclaim::PhaseCenter::ParameterManager::ParameterManager(ParameterSpace ps)
  : fParamSpace(ps)
{
  switch(fParamSpace){
  case ParameterSpace::None:             fN = 0;                              break;
  case ParameterSpace::PitchRoll:        fN = 2;                              break;
  case ParameterSpace::PitchRollHeading: fN = 3;                              break;
  case ParameterSpace::RingR:            fN = 4;                              break;
  case ParameterSpace::RingPhi:          fN = 4;                              break;
  case ParameterSpace::RingZ:            fN = 4;                              break;
  case ParameterSpace::RingPhiRZ:        fN = 12;                             break;
  case ParameterSpace::RingEllipse:      fN = EllipseParams::N()*4;           break;
  case ParameterSpace::ExtraDeltaT:      fN = NUM_SEAVEYS*AnitaPol::kNotAPol; break;
  }
}

Acclaim::PhaseCenter::ParameterManager::~ParameterManager(){;}


void Acclaim::PhaseCenter::ParameterManager::setInputs(ROOT::Math::Minimizer* fMin, std::vector<double>& fInputs) const {

  fInputs.resize(N(), 0);
  
  switch(fParamSpace){
  case ParameterSpace::None:
    break;
    
  case ParameterSpace::PitchRoll:
    {
      fMin->SetLimitedVariable(0, "pitch", 0, 0.1, -5, 5);
      fMin->SetLimitedVariable(1, "roll", 0, 0.1, -5, 5);
    }
    break;

  case ParameterSpace::PitchRollHeading:
    {
      fMin->SetLimitedVariable(0, "pitch", 0, 0.01, -1, 1);
      fMin->SetLimitedVariable(1, "roll", 0, 0.01, -1, 1);
      fMin->SetLimitedVariable(2, "heading_offset", 0, 0.01, -5, 5);
    }
    break;
    
      
  case ParameterSpace::RingR:
    {
      auto geom = AnitaGeomTool::Instance();
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
      const int numPhysicalRings = PhysicalRings.size();
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


void Acclaim::PhaseCenter::ParameterManager::update(const double* params){
  fParams = params;
  std::cout << "Update! " << fParams[0]  << "\t" << fParams[1]  <<  std::endl;

  if(fParamSpace==ParameterSpace::RingEllipse){
    fEllipseParams.clear();
    for(int ringInd=0; ringInd < PhysicalRings.size(); ringInd++){
      fEllipseParams.emplace_back(&fParams[ringInd*EllipseParams::N()]);
    }
  }
}


void Acclaim::PhaseCenter::ParameterManager::applyPat(Adu5Pat* pat) const {

  if(fParamSpace==ParameterSpace::PitchRoll){
    pat->pitch    = fParams[0];
    pat->roll     = fParams[1];
    std::cout << "apply!" << pat->pitch << "\t" << pat->roll << "\t"  << std::endl;
  }
  else if(fParamSpace==ParameterSpace::PitchRollHeading){
    pat->pitch    = fParams[0];
    pat->roll     = fParams[1];
    pat->heading += fParams[2];
  }
  else if(fParamSpace != ParameterSpace::None){
    // these are the results of fitting pitch/roll with photogrammetry numbers
    // pitch = -0.199065
    // roll = -0.113464
    // heading_offset = -0.552301

    pat->pitch    = -0.199065;
    pat->roll     = -0.113464;
    pat->heading += -0.552301;
  }
}

void Acclaim::PhaseCenter::ParameterManager::applyGeom(AnitaGeomTool* fGeom) const {

  for(int ant=0; ant < NUM_SEAVEYS; ant++){
    int phiSector = ant%NUM_PHI;
    for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
      auto ring = antToPhysicalRing(ant);
      int ringInt = static_cast<int>(ring);

      if(fParamSpace==ParameterSpace::RingR){
	fGeom->rPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = fParams[ringInt];
      }
      else if(fParamSpace==ParameterSpace::RingPhi){
	double phi = fParams[ringInt] + phiSector*TMath::DegToRad()*22.5;
	fGeom->azPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = phi;
      }
      else if(fParamSpace==ParameterSpace::RingPhiRZ){
	const int paramsPerRing = 3;
	double phi = fParams[ringInt*paramsPerRing] + phiSector*TMath::DegToRad()*22.5;
	fGeom->azPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = phi;
	fGeom->rPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = fParams[ringInt*paramsPerRing+1];
	fGeom->zPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = fParams[ringInt*paramsPerRing+2];
      }
      else if(fParamSpace==ParameterSpace::RingEllipse){

	for(int ant=0; ant < NUM_SEAVEYS; ant++){
	  PhysicalRing ring = antToPhysicalRing(ant);
	  int ringInd = static_cast<int>(ring);
	  double r, z;
	  double phi = fGeom->azPhaseCentreFromVerticalHornPhotogrammetry[ant][pol];
	  fEllipseParams.at(ringInd).phiToEllipseRZ(phi, r, z);
	  fGeom->rPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = r;
	  fGeom->zPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = z;
	}
      }
    }
  }
}













Acclaim::PhaseCenter::EllipseParams::EllipseParams(const double* params){
  if(params){
    fill(params);
  }
}

double Acclaim::PhaseCenter::EllipseParams::Rb() const {
  return Ra* TMath::Sqrt(1 - eccentricity*eccentricity);
}

void Acclaim::PhaseCenter::EllipseParams::fill(const double* params){
  x0           = params[0];
  y0           = params[1];
  alpha        = params[2];
  Ra           = params[3];
  eccentricity = params[4];
  z            = params[5];
}

const char* Acclaim::PhaseCenter::EllipseParams::name(int p) {
  switch(p){
  case 0: return "x0";
  case 1: return "y0";
  case 2: return "alpha";
  case 3: return "Ra";
  case 4: return "eccentricity";
  case 5: return "z";
  default: return "Unknown!";
  }
}


void Acclaim::PhaseCenter::EllipseParams::phiToEllipseXY(double phi, double& x, double& y) const {
  double t = tFromPhi(phi);
  x = Ra  *TMath::Cos(t);
  y = Rb()*TMath::Sin(t);
  TVector2 v(x, y);
  v = v.Rotate(alpha);
  x = v.X() + x0;
  y = v.Y() + y0;
}

void Acclaim::PhaseCenter::EllipseParams::phiToEllipseRZ(double phi, double& r, double& z) const {
  double x, y;
  phiToEllipseXY(phi, x, y);
  r = TMath::Sqrt(x*x + y*y);
  z = this->z;
}

double Acclaim::PhaseCenter::EllipseParams::tFromPhi(double phi) const {
  // here we find the parametric angle t from the angle phi.
  double localPhi = phi - alpha; // account for rotation of ellipse
  double tan_t = (Ra/Rb())*TMath::Tan(localPhi);
  double t = TMath::ATan(tan_t);

  // get correct quadrant, I hope...
  if(localPhi > 0.5*TMath::Pi() && localPhi <= 1.5*TMath::Pi()){
    t += TMath::Pi();
  }
  // std::cout << phi*TMath::RadToDeg() << "\t" << localPhi*TMath::RadToDeg() << "\t" << t*TMath::RadToDeg() << std::endl;

  return t;
}
