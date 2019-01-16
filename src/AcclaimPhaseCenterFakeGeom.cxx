#include "AcclaimPhaseCenterFakeGeom.h"
#include "AnitaGeomTool.h"

Acclaim::PhaseCenter::FakeGeomTool::FakeGeomTool(const AnitaGeomTool* geom){

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



void Acclaim::PhaseCenter::FakeGeomTool::restorePhotogrammetryPositionsInAnitaGeomTool(AnitaGeomTool* geom) const {
  for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      auto key = std::make_pair(pol, ant);
      geom->rPhaseCentreFromVerticalHornPhotogrammetry[ant][pol]  = fPhotoR.at(key);
      geom->zPhaseCentreFromVerticalHornPhotogrammetry[ant][pol]  = fPhotoZ.at(key);
      geom->azPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = fPhotoPhi.at(key);
    }
  }
}

void Acclaim::PhaseCenter::FakeGeomTool::overwritePhotogrammetryPositionsInAnitaGeomTool(AnitaGeomTool* geom) const {
  for(auto pol : {AnitaPol::kHorizontal, AnitaPol::kVertical}){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      geom->rPhaseCentreFromVerticalHornPhotogrammetry[ant][pol]  = getAntR(ant, pol);
      geom->zPhaseCentreFromVerticalHornPhotogrammetry[ant][pol]  = getAntZ(ant, pol);
      geom->azPhaseCentreFromVerticalHornPhotogrammetry[ant][pol] = getAntPhi(ant, pol);
    }
  }
}

