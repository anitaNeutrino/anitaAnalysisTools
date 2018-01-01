#include "AnalysisCuts.h"





// /** 
//  * Does not point to a known moving source?
//  * 
//  * @param sum is the AnitaEventSummary
//  * @param pol is the polarisation (default = AnitaPol::kNotAPol, see handleDefaults to see how this is handled)
//  * @param peakInd is the peak index (default = -1, see handleDefaults to see how this is handled)
//  * 
//  * @return 1 (true) if does not point to known moving source, 0 (false) points to known moving source
//  */
// int Acclaim::AnalysisCuts::DoesNotPointToKnownMovingSource::apply(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd) const {
//   handleDefaults(sum, pol, peakInd);

//   Adu5Pat pat = sum->anitaLocation.pat();
//   UsefulAdu5Pat usefulPat(&pat);

//   const int numTraverses = BaseList::getNumPaths();

//   for(int flightInd=0; flightInd < numTraverses; flightInd++){
//     const BaseList::path p = BaseList::getPath(flightInd);

//     if(p.isValid(sum->realTime)){

//       Double_t distanceLimit = p.isFlight ? fCloseEnoughToFlightMetres : fCloseEnoughToTraverseMetres;
//       AntarcticCoord pos = p.getPosition(sum->realTime).as(AntarcticCoord::WGS84);
//       // x, y, z = lat/lon/alt ... I hate the AntarcticCoord naming convention!

//       Double_t separation = usefulPat.getDistanceFromSource(pos.x, pos.y, pos.z);

//       if(separation > distanceLimit){
// 	Double_t thetaWave, phiWave;
// 	usefulPat.getThetaAndPhiWave(pos.y, pos.x, pos.z, thetaWave, phiWave);
// 	Double_t phiDeg = phiWave*TMath::RadToDeg();

// 	Double_t deltaPhiDeg = RootTools::getDeltaAngleDeg(phiDeg, sum->peak[pol][peakInd].phi);

// 	Double_t deltaPhiCut = p.isFlight ? fDeltaPhiDegFlights : fDeltaPhiDegTraverses;
// 	if(TMath::Abs(deltaPhiDeg) < deltaPhiCut){
// 	  return false;
// 	}
//       }
//     }
//   }
//   return true;
// }

