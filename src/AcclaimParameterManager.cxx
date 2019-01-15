#include "AcclaimParameterManager.h"
#include "AnitaGeomTool.h"
#include "Adu5Pat.h"
#include "TMath.h"


std::ostream& operator<<(std::ostream& os, const Acclaim::PhaseCenterFit::PhysicalRing& r){
  switch(r){
  case Acclaim::PhaseCenterFit::PhysicalRing::TopHigh:
    os << "PhysicalRing::TopHigh";
    return os;
  case Acclaim::PhaseCenterFit::PhysicalRing::TopLow:
    os << "PhysicalRing::TopLow";
    return os;
  case Acclaim::PhaseCenterFit::PhysicalRing::Middle:
    os << "PhysicalRing::Middle";
    return os;
  case Acclaim::PhaseCenterFit::PhysicalRing::Bottom:
    os << "PhysicalRing::Bottom";
    return os;
  }
  return os;
}



Acclaim::PhaseCenterFit::ParameterManager::ParameterManager(ParameterSpace ps,  int nParams, const double* params)
  : fN(nParams), fParams(params), fParamSpace(ps)
{

}

Acclaim::PhaseCenterFit::ParameterManager::~ParameterManager(){;}


void Acclaim::PhaseCenterFit::ParameterManager::update(const double* params){
  fParams = params;

  if(fParamSpace==ParameterSpace::RingEllipse){
    fEllipseParams.clear();
    const int numPhysicalRings = 4;
    for(int ringInd=0; ringInd < numPhysicalRings; ringInd++){
      fEllipseParams.emplace_back(&fParams[ringInd*EllipseParams::N()]);
    }
  }
}


void Acclaim::PhaseCenterFit::ParameterManager::applyPat(Adu5Pat* pat) const {
  if(fParamSpace==ParameterSpace::PitchRollHeading){
    pat->pitch = fParams[0];
    pat->roll = fParams[1];
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

void Acclaim::PhaseCenterFit::ParameterManager::applyGeom(AnitaGeomTool* fGeom) const {

  for(int ant=0; ant < NUM_SEAVEYS; ant++){
    int phiSector = ant%NUM_PHI;
    for(auto pol : {AnitaPol::kHorizontal,  AnitaPol::kVertical}){
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













Acclaim::PhaseCenterFit::EllipseParams::EllipseParams(const double* params){
  if(params){fill(params);}
}

double Acclaim::PhaseCenterFit::EllipseParams::Rb() const {
  return Ra* TMath::Sqrt(1 - eccentricity*eccentricity);
}

void Acclaim::PhaseCenterFit::EllipseParams::fill(const double* params){
  x0           = params[0];
  y0           = params[1];
  alpha        = params[2];
  Ra           = params[3];
  eccentricity = params[4];
  z            = params[5];
}

const char* Acclaim::PhaseCenterFit::EllipseParams::name(int p) {
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


void Acclaim::PhaseCenterFit::EllipseParams::phiToEllipseXY(double phi, double& x, double& y) const {
  double t = tFromPhi(phi);
  x = Ra  *TMath::Cos(t);
  y = Rb()*TMath::Sin(t);
  TVector2 v(x, y);
  v = v.Rotate(alpha);
  x = v.X() + x0;
  y = v.Y() + y0;
}

void Acclaim::PhaseCenterFit::EllipseParams::phiToEllipseRZ(double phi, double& r, double& z) const {
  double x, y;
  phiToEllipseXY(phi, x, y);
  r = TMath::Sqrt(x*x + y*y);
  z = this->z;
}

double Acclaim::PhaseCenterFit::EllipseParams::tFromPhi(double phi) const {
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
