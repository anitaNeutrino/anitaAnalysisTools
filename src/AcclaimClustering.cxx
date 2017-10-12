#include "AcclaimClustering.h"
#include "AnitaGeomTool.h"
#include "SummarySet.h"
#include "AnitaDataset.h"
#include "OutputConvention.h"
#include "AnitaEventSummary.h"
#include "AntarcticaMapPlotter.h"
#include "RootTools.h"
#include "ProgressBar.h"
#include "TGraphAntarctica.h"
#include "TH2DAntarctica.h"

const int nDim = 3;


Acclaim::Clustering::LogLikelihoodMethod::Point::Point(Adu5Pat* pat, Double_t lat, Double_t lon, Double_t alt,
						       Double_t theta, Double_t phi, Double_t sigmaTheta, Double_t sigmaPhi,
						       Int_t polIn){

  (void) pat;
  pol = (AnitaPol::AnitaPol_t) polIn;
  latitude = lat;
  longitude = lon;
  altitude = alt;

  // std::cout << __PRETTY_FUNCTION__ << std::endl << "latitude = " << latitude << "\tlongitude = " << longitude << std::endl << std::endl;

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->getCartesianCoords(lat, lon, alt, centre);

  thetaDeg = theta;
  phiDeg = phi;

  sigmaThetaDeg = sigmaTheta;
  sigmaPhiDeg = sigmaPhi;

  ll = DBL_MAX;
  inCluster = -1;
  llSecondBest = DBL_MAX;
  secondClosestCluster = -1;

}


Acclaim::Clustering::LogLikelihoodMethod::Point::Point(){
  latitude = 0;
  longitude = 0;
  altitude = 0;
  thetaDeg = -9999;
  phiDeg = -9999;

  // convert to degrees
  thetaDeg = -9999;
  phiDeg = -9999;
  // dTheta = 0;
  // dPhi = 0;
  sigmaThetaDeg = -9999;
  sigmaPhiDeg = -9999;
  ll = DBL_MAX;
  inCluster = -1;
  llSecondBest = DBL_MAX;
  secondClosestCluster = -1;
  for(int dim=0; dim < nDim; dim++){
    centre[dim] = 0;
  }

}



Acclaim::Clustering::LogLikelihoodMethod::McPoint::McPoint()
  : Point(){
  weight = 0; energy=0;
}



Acclaim::Clustering::LogLikelihoodMethod::McPoint::McPoint(Adu5Pat* pat, Double_t lat, Double_t lon, Double_t alt,
							   Double_t thetaDeg, Double_t phiDeg, Double_t sigmaTheta, Double_t sigmaPhi,
							   Int_t polIn, Double_t theWeight, Double_t theEnergy)
  : Point(pat,	lat, lon, alt, sigmaTheta, sigmaPhi, polIn){
  weight = theWeight;
  energy = theEnergy;
}


Acclaim::Clustering::LogLikelihoodMethod::Cluster::Cluster() {
  for(int dim=0; dim < nDim; dim++){
    centre[dim] = 0;
  }
  numEvents = 0;
  totalError = 0;
  latitude = 0;
  longitude = 0;
  altitude = 0;
  maxDist = 0;
  numPointsWithinMinLL = 0;
}

Acclaim::Clustering::LogLikelihoodMethod::Cluster::Cluster(const Point& seedPoint) {
  latitude = seedPoint.latitude;
  longitude = seedPoint.longitude;
  altitude = seedPoint.altitude;
  for(int dim=0; dim < nDim; dim++){
    centre[dim] = seedPoint.centre[dim];
  }
  numEvents = 1; // since the seed point will be in this cluster
  totalError = 0;
  maxDist = 0;
  numPointsWithinMinLL = 0;
}

Acclaim::Clustering::LogLikelihoodMethod::Cluster::Cluster(const BaseList::base& base) {
  latitude = base.latitude;
  longitude = base.longitude;
  altitude = base.altitude;

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->getCartesianCoords(latitude, longitude, altitude, centre);
  numEvents = 0;
  totalError = 0;
  maxDist = 0;
}



Acclaim::Clustering::LogLikelihoodMethod::LogLikelihoodMethod(){

  // numIter = numIterations;
  // numClusters = nClusters;
  // clusters.reserve(numClusters);

  // points.reserve(approxNumPoints);
  // pats.reserve(approxNumPoints);
  // eventNumbers.reserve(approxNumPoints);
  // runs.reserve(approxNumPoints);

  maxRetestClusterSize = 10;
  numIsolatedSmallClusters.resize(maxRetestClusterSize, 0);
  numIsolatedSmallBaseClusters.resize(maxRetestClusterSize, 0);

  llCut = 250;
  maxDistCluster = 800e3; // try 800km
  numCallsToRecursive = 0;
  
  // can use this to move cluster around surface, 3D positions correctly becomes 2D problem
  doneDefaultAssignment = false;
}

Acclaim::Clustering::LogLikelihoodMethod::~LogLikelihoodMethod(){

  // delete non-NULL antarctic histograms
  while(hUnclustereds.size() > 0){
    TH2DAntarctica* h = hUnclustereds.back();
    if(h){
      delete h;
    }
    hUnclustereds.pop_back();
  }
}




// for debugging
inline void prettyPrint(const int n, const double* array){
  for(int i=0; i < n; i++){
    std::cerr << array[i];
    if(i < n - 1){
      std::cerr << ", ";
    }
  }
  std::cerr << std::endl;
}



inline void prettyPrintConvert(const int n, const double* array){
  Double_t lat, lon, alt;
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  std::vector<Double_t> temp(array, array+n*sizeof(double));
  geom->getLatLonAltFromCartesian(&temp[0], lat, lon, alt);
  std::cerr << lat << "\t" << lon << "\t" << alt << std::endl;
}



// utility function hopefully this one gets inlined
inline Double_t getDistSq(const Acclaim::Clustering::LogLikelihoodMethod::Point& point, const Acclaim::Clustering::LogLikelihoodMethod::Cluster& cluster){
  Double_t d2=0;
  for(int dim=0; dim < nDim; dim++){
    Double_t d = cluster.centre[dim] - point.centre[dim];
    d2 += d*d;
  }
  return d2;
}


// utility function hopefully this one gets inlined
inline Double_t getDistSq(Acclaim::Clustering::LogLikelihoodMethod::Cluster& cluster1, const Acclaim::Clustering::LogLikelihoodMethod::Cluster& cluster2){
  Double_t d2=0;
  for(int dim=0; dim < nDim; dim++){
    Double_t d = cluster1.centre[dim] - cluster2.centre[dim];
    d2 += d*d;
  }
  return d2;
}




inline void getDeltaThetaDegDeltaPhiDegCluster(const Acclaim::Clustering::LogLikelihoodMethod::Point& point, const Acclaim::Clustering::LogLikelihoodMethod::Cluster& cluster, UsefulAdu5Pat& usefulPat, Double_t& deltaThetaDeg, Double_t& deltaPhiDeg){

  Double_t thetaWave, phiWave;
  usefulPat.getThetaAndPhiWave(cluster.longitude, cluster.latitude, cluster.altitude, thetaWave, phiWave);
  Double_t thetaDeg = -1*TMath::RadToDeg()*thetaWave;
  Double_t phiDeg = TMath::RadToDeg()*phiWave;

  deltaThetaDeg = (thetaDeg - point.thetaDeg);
  deltaPhiDeg = Acclaim::RootTools::getDeltaAngleDeg(phiDeg, point.phiDeg);
}






inline void getDeltaThetaDegDeltaPhiDegCluster(const Acclaim::Clustering::LogLikelihoodMethod::Point& point, const Acclaim::Clustering::LogLikelihoodMethod::Cluster& cluster, const Adu5Pat* pat, Double_t& deltaThetaDeg, Double_t& deltaPhiDeg){

  UsefulAdu5Pat usefulPat(pat);

  getDeltaThetaDegDeltaPhiDegCluster(point, cluster, usefulPat, deltaThetaDeg, deltaPhiDeg);
}







inline Double_t getAngDistSq(const Acclaim::Clustering::LogLikelihoodMethod::Point& point, const Acclaim::Clustering::LogLikelihoodMethod::Cluster& cluster, UsefulAdu5Pat& usefulPat){

  Double_t deltaThetaDeg, deltaPhiDeg;
  getDeltaThetaDegDeltaPhiDegCluster(point, cluster, usefulPat, deltaThetaDeg, deltaPhiDeg);

  // normalized
  Double_t dThetaNorm = deltaThetaDeg/point.sigmaThetaDeg;
  Double_t dPhiNorm = deltaPhiDeg/point.sigmaPhiDeg;

  // std::cerr << point.thetaDeg << "\t" << point.phiDeg << "\t" << deltaThetaDeg << "\t" << deltaPhiDeg << "\t" << dThetaNorm << "\t" << dPhiNorm << std::endl;  

  Double_t angSq =  dThetaNorm*dThetaNorm + dPhiNorm*dPhiNorm;
  
  return angSq;
}


inline Double_t getAngDistSq(const Acclaim::Clustering::LogLikelihoodMethod::Point& point, const Acclaim::Clustering::LogLikelihoodMethod::Cluster& cluster, const Adu5Pat* pat){
  UsefulAdu5Pat usefulPat(pat);
  return getAngDistSq(point, cluster, usefulPat);
}





Double_t Acclaim::Clustering::LogLikelihoodMethod::getSumOfMcWeights(){
  Double_t sumOfWeights = 0;
  for(int pointInd=0; pointInd < (int)mcPoints.size(); pointInd++){
    sumOfWeights += mcPoints.at(pointInd).weight;
  }
  return sumOfWeights;
}



Int_t Acclaim::Clustering::LogLikelihoodMethod::histogramUnclusteredEvents(Int_t& globalMaxBin){

  Int_t addClusterIter =  hUnclustereds.size();
  TString name = TString::Format("hUnclustered%d", addClusterIter);
  // const int nBinsLat = 30;
  // const int nBinsLon = 360;
  // const int nBins = 1024;
  // amp->addHistogram(name, name, nBins, nBins);

  TH2DAntarctica* h = new TH2DAntarctica(name, name);

  hUnclustereds.push_back(h);

  // find unclustered points
  ampBinNumbers.clear(); // doesn't clear memory
  ampBinNumbers.resize((Int_t) points.size()); // so this shouldn't reallocate memory
  for(int pointInd=0; pointInd < (Int_t) points.size(); pointInd++){
    const Point& point = points.at(pointInd);
    if(point.inCluster < 0){
      ampBinNumbers.at(pointInd) = h->Fill(point.longitude, point.latitude);
    }
  }
  globalMaxBin = hUnclustereds.back()->GetMaximumBin();

  Int_t maxEventsInBin = hUnclustereds.back()->GetBinContent(globalMaxBin);

  return maxEventsInBin;

}



void Acclaim::Clustering::LogLikelihoodMethod::recursivelyAddClusters(Int_t minBinContent){

  static int lastIntegral = 0;

  if(doneDefaultAssignment==false){
    assignEventsToDefaultClusters();
  }

  numCallsToRecursive++;
  Int_t globalMaxBin;
  Int_t maxBinVal = histogramUnclusteredEvents(globalMaxBin);

  // std::cout << minLat << "\t" << maxLat << "\t" << minLon << "\t" << maxLon << std::endl;

  // std::cout << "Info in " << __PRETTY_FUNCTION__ << ",  " << numCallsToRecursive << "\t" << hUnclustereds.back()->Integral() << std::endl;

  Int_t counter=0;
  Cluster cluster;
  for(int pointInd=0; pointInd < (Int_t) points.size(); pointInd++){
    const Point& point = points.at(pointInd);
    if(point.inCluster < 0){
      if(ampBinNumbers.at(pointInd)==globalMaxBin){
	for(int dim=0; dim < nDim; dim++){
	  cluster.centre[dim] += point.centre[dim];
	}
	counter++;
      }
    }
  }

  if(counter!=maxBinVal){
    std::cerr << "Now what?" << std::endl;
  }
  if(counter > minBinContent){

    for(int dim=0; dim < nDim; dim++){
      cluster.centre[dim]/=counter;
    }
    AnitaGeomTool* geom = AnitaGeomTool::Instance();
    geom->getLatLonAltFromCartesian(cluster.centre, cluster.latitude,\
				    cluster.longitude, cluster.altitude);
    clusters.push_back(cluster);
    numClusters = (Int_t) clusters.size();

    const int isMC = 0;
    for(int pointInd=0; pointInd < (Int_t) points.size(); pointInd++){
      assignSinglePointToCloserCluster(pointInd, isMC, numClusters-1);
      // std::cout << pointInd << "\t" << points.at(pointInd).ll << "\t" << std::endl;
    }

    const int numEventsUnclustered = hUnclustereds.back()->Integral();
    if(numEventsUnclustered > 0 && lastIntegral != numEventsUnclustered){
      lastIntegral = hUnclustereds.back()->Integral();
      recursivelyAddClusters(minBinContent);
    }
    else{
      delete hUnclustereds.back();
      hUnclustereds.pop_back();
    }

  }
}




void Acclaim::Clustering::LogLikelihoodMethod::findClosestPointToClustersOfSizeOne(){

  for(int i=0; i < maxRetestClusterSize; i++){
    numIsolatedSmallClusters.at(i) = 0;
    numIsolatedSmallBaseClusters.at(i) = 0;
  }

  // loop through clusters and find the small ones
  for(int clusterInd=0; clusterInd < numClusters; clusterInd++){
    if(clusters.at(clusterInd).numEvents > 0 && clusters.at(clusterInd).numEvents <= maxRetestClusterSize){


      Cluster& cluster = clusters.at(clusterInd);

      Int_t numNearPoints = 0;

      // now loop through and see how close this isolated cluster is close to any other point...
      Double_t minLL = 9e99;
      Int_t minPointInd = -1;
      for(int pointInd=0; pointInd < (int)points.size(); pointInd++){

	// skip points in the cluster under test
	if(points.at(pointInd).inCluster!=clusterInd){

	  UsefulAdu5Pat usefulPat(&pats.at(pointInd));
	  Double_t dist = usefulPat.getDistanceFromSource(cluster.latitude,
							  cluster.longitude,
							  cluster.altitude);

	  // this is where proper collision detection would be useful...
	  if(dist < maxDistCluster){
	    Double_t ll = getAngDistSq(points.at(pointInd), cluster, usefulPat);

	    if(ll < llCut){

	      numNearPoints++;
	      cluster.numPointsWithinMinLL++;
	      if(ll < minLL){
		minLL = ll;
		minPointInd = pointInd;
	      }
	    }
	  }
	}
      }

      // std::cout << runs.at(isolatedPointInd) << "\t" << eventNumbers.at(isolatedPointInd) << "\t"
      // 		<< clusterInd << "\t" << numNearPoints << std::endl;

      if(numNearPoints > 0){
	std::cout << "reassign " << clusterInd << "\t" << cluster.numEvents << " to " << points.at(minPointInd).inCluster << "\t" << clusters.at(points.at(minPointInd).inCluster).numEvents << std::endl;
	// re-assign points...
	for(int pointInd=0; pointInd < (int)points.size(); pointInd++){
	  if(points.at(pointInd).inCluster==clusterInd){
	    clusters.at(points.at(minPointInd).inCluster).numEvents++;
	    cluster.numEvents--;
	    points.at(pointInd).inCluster = points.at(minPointInd).inCluster;
	    points.at(pointInd).secondClosestCluster = points.at(minPointInd).secondClosestCluster;
	  }
	}
      }
      // else{
      // 	std::cout << "don't reassign " << clusterInd << "\t" << cluster.numEvents << std::endl;
      // }
    }
  }

  const int numBases = BaseList::getNumBases();
  for(int clusterInd=0; clusterInd < numClusters; clusterInd++){
    if(clusters.at(clusterInd).numEvents > 0 &&
       clusters.at(clusterInd).numEvents <= maxRetestClusterSize){

      if(clusterInd < numBases){
	// if(clusters.at(clusterInd).numEvents==1){
	//   std::cout << clusterInd << "\t" << clusters.at(clusterInd).numEvents  << std::endl;
	// }
	numIsolatedSmallBaseClusters.at(clusters.at(clusterInd).numEvents-1)++;
      }
      else{
	numIsolatedSmallClusters.at(clusters.at(clusterInd).numEvents-1)++;
      }

    }
  }
}







void Acclaim::Clustering::LogLikelihoodMethod::assignSinglePointToCloserCluster(Int_t pointInd, Int_t isMC, Int_t clusterInd){
  Point& point = isMC==0 ? points.at(pointInd) : mcPoints.at(pointInd);
  const Adu5Pat& pat = isMC==0 ? pats.at(pointInd)   : mcPats.at(pointInd);
  
  Cluster& cluster = clusters.at(clusterInd);

  UsefulAdu5Pat usefulPat(&pat);

  Double_t distM = usefulPat.getDistanceFromSource(cluster.latitude, cluster.longitude, cluster.altitude);

  if(distM < maxDistCluster){ // are we even close?
    
    Double_t ll = getAngDistSq(point, cluster, usefulPat);


    if(ll < llCut){ // within the cut?

      if(!isMC){
	std::cout << ll << std::endl;
      }
      
      if(ll < point.ll){ // the best point?
 
	point.secondClosestCluster = point.inCluster;
	point.llSecondBest = point.ll;

	point.ll = ll;
	point.inCluster = clusterInd;

	if(isMC==0){
	  if(point.secondClosestCluster >= 0){
	    clusters.at(point.secondClosestCluster).numEvents--;
	  }
	  cluster.numEvents++;
	}
      }
    }
    
  }
}




void Acclaim::Clustering::LogLikelihoodMethod::assignMcPointsToClusters(){

  const int isMC = 1;
  Double_t numSinglets = 0;
  for(int pointInd=0; pointInd < (Int_t)mcPoints.size(); pointInd++){
    for(int clusterInd=0; clusterInd < (Int_t)clusters.size(); clusterInd++){
      assignSinglePointToCloserCluster(pointInd, isMC, clusterInd);
    }
    if(mcPoints.at(pointInd).inCluster < 0){
      numSinglets+=mcPoints.at(pointInd).weight;
    }
  }
  numMcIsolatedSinglets = numSinglets;
}










size_t Acclaim::Clustering::LogLikelihoodMethod::addMcPoint(Adu5Pat* pat, Double_t latitude, Double_t longitude, Double_t altitude, Int_t run, UInt_t eventNumber, Double_t thetaDeg, Double_t phiDeg, Double_t sigmaThetaDeg, Double_t sigmaPhiDeg, AnitaPol::AnitaPol_t pol, Double_t weight, Double_t energy){

  mcPoints.push_back(McPoint(pat, latitude, longitude, altitude, thetaDeg, phiDeg, sigmaThetaDeg, sigmaPhiDeg, pol, weight, energy));
  mcEventNumbers.push_back(eventNumber);
  mcRuns.push_back(run);
  mcPats.push_back(*pat);

  return points.size();
}








size_t Acclaim::Clustering::LogLikelihoodMethod::addPoint(Adu5Pat* pat, Double_t latitude, Double_t longitude, Double_t altitude, Int_t run, UInt_t eventNumber, Double_t thetaDeg, Double_t phiDeg, Double_t sigmaThetaDeg, Double_t sigmaPhiDeg, AnitaPol::AnitaPol_t pol){

  points.push_back(Point(pat, latitude, longitude, altitude, thetaDeg, phiDeg, sigmaThetaDeg, sigmaPhiDeg, pol));
  // std::cout << __PRETTY_FUNCTION__ << std::endl << "latitude = " << latitude << "\tlongitude = " << longitude << std::endl << std::endl;  
  eventNumbers.push_back(eventNumber);
  runs.push_back(run);
  pats.push_back(*pat);

  return points.size();
}







void Acclaim::Clustering::LogLikelihoodMethod::assignEventsToDefaultClusters(){
  std::cout << "Info in " << __PRETTY_FUNCTION__ << " assigning points to closest cluster!" << std::endl;
  ProgressBar p(numClusters);
  const int isMC = 0;
  for(Long64_t clusterInd=0; clusterInd < numClusters; clusterInd++){
    for(int pointInd=0; pointInd < (int) points.size(); pointInd++){
      assignSinglePointToCloserCluster(pointInd, isMC, clusterInd);
    }
    p.inc(clusterInd, numClusters);
  }
  doneDefaultAssignment = true;
}


void Acclaim::Clustering::LogLikelihoodMethod::initializeEmptyBaseList(){
  numClusters=0;
  BaseList::makeEmptyBaseList();
}


void Acclaim::Clustering::LogLikelihoodMethod::initializeBaseList(){

  numClusters = (int) BaseList::getNumBases();
  std::cout << "Info in " << __PRETTY_FUNCTION__ << ": Initializing base list..." << std::endl;
  for(int clusterInd=0; clusterInd < numClusters; clusterInd++){
    const BaseList::base& base = BaseList::getBase(clusterInd);
    clusters.push_back(Cluster(base));
  }
}



void Acclaim::Clustering::LogLikelihoodMethod::resetClusters(){

  for(int pointInd=0; pointInd < (int) points.size(); pointInd++){
    points.at(pointInd).ll = DBL_MAX;
    points.at(pointInd).llSecondBest = DBL_MAX;
    points.at(pointInd).inCluster = -1;
    points.at(pointInd).secondClosestCluster = -1;
  }

  for(int pointInd=0; pointInd < (int) mcPoints.size(); pointInd++){
    mcPoints.at(pointInd).ll = DBL_MAX;
    mcPoints.at(pointInd).llSecondBest = DBL_MAX;
    mcPoints.at(pointInd).inCluster = -1;
    mcPoints.at(pointInd).secondClosestCluster = -1;
  }

  const int numBases = BaseList::getNumBases();
  while((int) clusters.size() >= numBases && numBases>0){
    clusters.pop_back();
  }
  
  numClusters = clusters.size();
  for(int clusterInd=0; clusterInd < (int) numClusters; clusterInd++){
    clusters.at(clusterInd).numEvents = 0;
  }
  doneDefaultAssignment = false;
}


































TGraphAntarctica* Acclaim::Clustering::LogLikelihoodMethod::makeClusterSummaryTGraph(Int_t clusterInd){

  TGraphAntarctica* gr = NULL;
  if(clusterInd >= 0 && clusterInd < numClusters){

    TString name  = TString::Format("grCluster%d", clusterInd);
    TString title  = TString::Format("Cluster %d; Easting (m); Northing (m)", clusterInd);
    gr = new TGraphAntarctica();
    gr->SetName(name);
    gr->SetTitle(title);

    // AnitaGeomTool* geom = AnitaGeomTool::Instance();

    for(int pointInd=0; pointInd < (Int_t) points.size(); pointInd++){
      if(points.at(pointInd).inCluster==clusterInd){

	// Double_t lat, lon, alt;
	// geom->getLatLonAltFromCartesian(points.at(pointInd).centre, lat, lon, alt);
	gr->SetPoint(gr->GetN(), points.at(pointInd).latitude, points.at(pointInd).longitude);
      }
    }
  }
  return gr;
}





// TTree* Acclaim::Clustering::LogLikelihoodMethod::makeClusterSummaryTree(TFile* fOut, TFile* fSignalBox){

//   fOut->cd();

//   TTree* clusterTree = new TTree("clusterTree", "Tree of clustered ANITA events");
//   fSignalBox->cd();
//   TTree* clusterTreeHidden = new TTree("clusterTreeHidden", "Tree of clustered hidden ANITA events");

//   Event* clusteredEvent;
//   clusterTree->Branch("clusteredEvent", &clusteredEvent);

//   clusterTreeHidden->Branch("clusteredEvent", &clusteredEvent);

//   //  std::cout << "Num points " << (Int_t)points.size() << std::endl;
//   // AnitaGeomTool* geom = AnitaGeomTool::Instance();
//   const int numBases = BaseList::getNumBases();

//   for(Int_t pointInd=0; pointInd < (Int_t)points.size(); pointInd++){

//     const Point& point = points.at(pointInd);

//     clusteredEvent = new Event();

//     Int_t clusterInd = points.at(pointInd).inCluster;
//     clusteredEvent->inCluster = clusterInd;
//     clusteredEvent->isBase = 0;
//     if(clusterInd >= 0){
//       clusteredEvent->eventNumber = eventNumbers.at(pointInd);
//       clusteredEvent->run = runs.at(pointInd);
//       clusteredEvent->pol = points.at(pointInd).pol;

//       // convert from km back to m for conversion to lat/lon/alt
//       // for(int dim=0; dim < nDim; dim++){
//       //   // clusteredEvent->eventPosition[dim] = point.centre[dim];
//       //   clusteredEvent->eventPosition[dim] = point.centre[dim];
//       // }
//       clusteredEvent->eventLat = point.latitude;
//       clusteredEvent->eventLon = point.longitude;
//       clusteredEvent->eventAlt = point.altitude;


//       clusteredEvent->anitaLat = pats.at(pointInd)->latitude;
//       clusteredEvent->anitaLon = pats.at(pointInd)->longitude;
//       clusteredEvent->anitaAlt = pats.at(pointInd)->altitude;

//       clusteredEvent->thetaDeg = point.thetaDeg;
//       clusteredEvent->phiDeg = point.phiDeg;
//       clusteredEvent->isMC = 0;
//       clusteredEvent->weight = 1;

//       const Cluster& cluster = clusters.at(clusterInd);

//       clusteredEvent->isBase = clusterInd < numBases ? 1 : 0;

//       clusteredEvent->clusterLat = cluster.latitude;
//       clusteredEvent->clusterLon = cluster.longitude;
//       clusteredEvent->clusterAlt = cluster.altitude;

//       // Double_t deltaThetaDeg, deltaPhiDeg;
//       getDeltaThetaDegDeltaPhiDegCluster(point, cluster, pats.at(pointInd), \
// 					 clusteredEvent->deltaThetaDeg, \
// 					 clusteredEvent->deltaPhiDeg);

//       // clusteredEvent->distanceToClusterCentroid = TMath::Sqrt(getDistSq(point, cluster));
//       // clusteredEvent->distanceToClusterCentroid = get(point, cluster));
//       clusteredEvent->minusTwoLogLikelihood = getAngDistSq(point, cluster, pats.at(pointInd));

//       UsefulAdu5Pat usefulPat(pats.at(pointInd));
//       clusteredEvent->distanceToClusterCentroid = usefulPat.getDistanceFromSource(cluster.latitude, \
// 										  cluster.longitude, \
// 										  cluster.altitude);
//       clusteredEvent->numEventsInCluster = cluster.numEvents;

//       if ( clusteredEvent->secondClosestCluster>=(int)clusters.size() ||  clusteredEvent->secondClosestCluster>=(int)points.size()  ){
// 	//	std::cout <<  "Second closest cluster problem? " << clusteredEvent->secondClosestCluster<< " " << clusterInd << " " << pointInd << " " << (int)clusters.size()  << " " << points.size() << std::endl;
//       } else    if(clusteredEvent->secondClosestCluster >= 0){
	
// 	const Cluster& cluster2 = clusters.at(clusteredEvent->secondClosestCluster);
	
// 	clusteredEvent->distanceToClusterCentroidSecondBest = usefulPat.getDistanceFromSource(cluster2.latitude, \
// 											      cluster2.longitude, \
// 											      cluster2.altitude);
//       }
    
//       clusteredEvent->numIterations = numIter;

//       clusteredEvent->numClusters = numClusters;

//     }
//     else{
//       std::cerr << "clusterInd = " << clusterInd << ". ";
//       std::cerr << "pointInd = " << pointInd << ". ";
//       std::cerr << "This shouldn't be possible!" << std::endl;
//     }

//     // inviting out of range exception if my clustering doesn't work properly
//     if(clusteredEvent->isBase == 0 &&
//        clusters.at(clusteredEvent->inCluster).numEvents==1 &&
//        clusters.at(clusteredEvent->inCluster).numPointsWithinMinLL == 0){
//       // the box of goodies
//       clusterTreeHidden->Fill();
//     }
//     else{
//       clusterTree->Fill();
//     }


//     delete clusteredEvent;
//   }


//   for(Int_t pointInd=0; pointInd < (Int_t)mcPoints.size(); pointInd++){

//     const McPoint& point = mcPoints.at(pointInd);

//     clusteredEvent = new Event();
//     clusteredEvent->eventNumber = mceventNumbers.at(pointInd);
//     clusteredEvent->run = mcruns.at(pointInd);
//     clusteredEvent->pol = point.pol;
//     clusteredEvent->isMC = 1;
//     clusteredEvent->weight = point.weight;
//     clusteredEvent->energy = point.energy;

//     // convert from km back to m for conversion to lat/lon/alt
//     // for(int dim=0; dim < nDim; dim++){
//     //   // clusteredEvent->eventPosition[dim] = point.centre[dim];
//     //   clusteredEvent->eventPosition[dim] = point.centre[dim];
//     // }
//     clusteredEvent->eventLat = point.latitude;
//     clusteredEvent->eventLon = point.longitude;
//     clusteredEvent->eventAlt = point.altitude;

//     clusteredEvent->anitaLat = mcpats.at(pointInd)->latitude;
//     clusteredEvent->anitaLon = mcpats.at(pointInd)->longitude;
//     clusteredEvent->anitaAlt = mcpats.at(pointInd)->altitude;

//     clusteredEvent->thetaDeg = point.thetaDeg;
//     clusteredEvent->phiDeg = point.phiDeg;

//     Int_t clusterInd = mcPoints.at(pointInd).inCluster;
//     clusteredEvent->inCluster = clusterInd;
//     clusteredEvent->isBase = 0;

//     if(clusterInd >= 0  && clusterInd<(int)clusters.size() ){
//       const Cluster& cluster = clusters.at(clusterInd);

//       clusteredEvent->isBase = clusterInd < numBases ? 1 : 0;
//       clusteredEvent->numPointsWithinMinLL = cluster.numPointsWithinMinLL;

//       clusteredEvent->clusterLat = cluster.latitude;
//       clusteredEvent->clusterLon = cluster.longitude;
//       clusteredEvent->clusterAlt = cluster.altitude;

//       // Double_t deltaThetaDeg, deltaPhiDeg;
//       getDeltaThetaDegDeltaPhiDegCluster(point, cluster, mcpats.at(pointInd), \
// 					 clusteredEvent->deltaThetaDeg, \
// 					 clusteredEvent->deltaPhiDeg);

//       // clusteredEvent->distanceToClusterCentroid = TMath::Sqrt(getDistSq(point, cluster));
//       // clusteredEvent->distanceToClusterCentroid = get(point, cluster));
//       clusteredEvent->minusTwoLogLikelihood = getAngDistSq(point, cluster, mcpats.at(pointInd));

//       UsefulAdu5Pat usefulPat(mcpats.at(pointInd));

//       clusteredEvent->distanceToClusterCentroid = usefulPat.getDistanceFromSource(cluster.latitude, \
// 										  cluster.longitude, \
// 										  cluster.altitude);


//       if(clusteredEvent->secondClosestCluster >= 0){
// 	const Cluster& cluster2 = clusters.at(clusteredEvent->secondClosestCluster);


// 	clusteredEvent->distanceToClusterCentroidSecondBest = usefulPat.getDistanceFromSource(cluster2.latitude, \
// 											      cluster2.longitude, \
// 											      cluster2.altitude);
//       }

//       clusteredEvent->numEventsInCluster = cluster.numEvents;

//       clusteredEvent->numClusters = numClusters;
//       clusteredEvent->numIterations = numIter;
//     }
//     clusterTree->Fill();
    
//     delete clusteredEvent;
    
//   }

  
//   clusterTreeHidden->Write();
//   fSignalBox->Write();
//   fSignalBox->Close();

//   fOut->cd();
//   return clusterTree;

// }


Long64_t Acclaim::Clustering::LogLikelihoodMethod::readInSummaries(const char* summaryGlob){

  Long64_t n = 0;
  if(summaryGlob){
    SummarySet ss(summaryGlob);
    n = ss.N();

    std::cout << "Info in " << __PRETTY_FUNCTION__ << ": reading in summaries: " << summaryGlob << std::endl;
    ProgressBar p(n);
    for(Long64_t entry=0; entry < n; entry++){
    
      ss.getEntry(entry);
      AnitaEventSummary* sum = ss.summary();

      AnitaPol::AnitaPol_t pol = sum->trainingPol();
      Int_t peakIndex = sum->trainingPeakInd();
    
      Double_t sourceLat = sum->peak[pol][peakIndex].latitude;
      Double_t sourceLon = sum->peak[pol][peakIndex].longitude;
      Double_t sourceAlt = sum->peak[pol][peakIndex].altitude;
      if(sourceLat > -999 && sourceLon > -999){

	Adu5Pat pat = static_cast<Adu5Pat>(sum->anitaLocation);

	Double_t thetaDeg = sum->peak[pol][peakIndex].theta;
	Double_t phiDeg = sum->peak[pol][peakIndex].phi;
	
	/// @todo update this
	Double_t sigmaTheta = 0.25;	
	Double_t sigmaPhi = 0.5;

	if(sum->mc.weight > 0){
	  addMcPoint(&pat, sourceLat, sourceLon, sourceAlt,
		     sum->run, sum->eventNumber, thetaDeg, phiDeg,
		     sigmaTheta, sigmaPhi, pol, sum->mc.weight, sum->mc.energy);
	
	}
	else{
	  addPoint(&pat, sourceLat,sourceLon,sourceAlt, sum->run, sum->eventNumber,
	  	   thetaDeg, phiDeg, sigmaTheta, sigmaPhi, pol);
	  // std::cout << __PRETTY_FUNCTION__ <<  "latitude = " << sourceLat << "\tlongitude = " << sourceLon <<  std::endl;
	  // addPoint(&pat, sourceLon,sourceLat,sourceAlt, sum->run, sum->eventNumber,
	  // 	   thetaDeg, phiDeg, sigmaTheta, sigmaPhi, pol);      
	}
      }
      p.inc(entry, n);
    }
  }
  return n;
}



void Acclaim::Clustering::LogLikelihoodMethod::doClustering(const char* dataGlob, const char* mcGlob, const char* outFileName){  
  
  readInSummaries(dataGlob);
  initializeBaseList();
  readInSummaries(mcGlob);

  char* fakeArgv0 = const_cast<char*>(outFileName);
  OutputConvention oc(1, &fakeArgv0);
  TFile* fOut = oc.makeFile();

  recursivelyAddClusters(0);

  for(UInt_t i=0; i < hUnclustereds.size(); i++){
    hUnclustereds.at(i)->Write();
    delete hUnclustereds.at(i);
    hUnclustereds.at(i) = NULL;
  }

  // TGraph* grNumPoints = new TGraph();
  // TGraph* grNumMcSinglets = new TGraph();
  // std::vector<TGraph*> grNumIsolateds(maxRetestClusterSize, NULL);
  // std::vector<TGraph*> grNumBaseIsolateds(maxRetestClusterSize, NULL);
  // for(int i=0; i < maxRetestClusterSize; i++){
  //   grNumIsolateds.at(i) = new TGraph();
  //   grNumBaseIsolateds.at(i) = new TGraph();
  // }

  // const int numLLs = 1;
  // double llCuts[numLLs] = {100};

  // for(int i=0; i < numLLs; i++){
  //   llCut = llCuts[i];
  //   resetClusters();
  //   recursivelyAddClusters(0);

  //   std::cout << "********************************" << std::endl;
  //   std::cout << "llCut = " << llCuts[i] << std::endl;
  //   std::cout << "********************************" << std::endl;
  //   findClosestPointToClustersOfSizeOne();

  //   assignMcPointsToClusters();

  //   for(int j=0; j < maxRetestClusterSize; j++){
  //     grNumIsolateds.at(j)->SetPoint(grNumIsolateds.at(j)->GetN(),
  // 				     llCut, numIsolatedSmallClusters.at(j));
      
  //     grNumBaseIsolateds.at(j)->SetPoint(grNumBaseIsolateds.at(j)->GetN(),
  // 					 llCut, numIsolatedSmallBaseClusters.at(j));
  //   }
  // }

  // for(int i=0; i < maxRetestClusterSize; i++){
  //   TString grName = TString::Format("grNumBaseIsolateds_%d", i+1);
  //   grNumBaseIsolateds.at(i)->SetName(grName);
  //   grNumBaseIsolateds.at(i)->Write();

  //   if(i > 0){
  //     grName = TString::Format("grNumIsolateds_%d", i+1);
  //     grNumIsolateds.at(i)->SetName(grName);
  //     grNumIsolateds.at(i)->Write();
  //   }
  // }

  fOut->Write();
  fOut->Close();

}
