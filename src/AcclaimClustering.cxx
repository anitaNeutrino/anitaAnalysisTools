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
							   Double_t theta, Double_t phi, Double_t sigmaTheta, Double_t sigmaPhi,
							   Int_t polIn, Double_t theWeight, Double_t theEnergy)
  : Point(pat,	lat, lon, alt, theta, phi, sigmaTheta, sigmaPhi, polIn){
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
  doneBaseClusterAssignment = false;
}

Acclaim::Clustering::LogLikelihoodMethod::~LogLikelihoodMethod(){

  // delete non-NULL antarctic histograms, tgraphs

  for(UInt_t i=0; i < hBaseClusteredEvents.size(); i++){
    if(hBaseClusteredEvents.at(i)){
      delete hBaseClusteredEvents.at(i);
      hBaseClusteredEvents.at(i) = NULL;
    }
  }
  for(UInt_t i=0; i < grBaseClusterCenters.size(); i++){
    if(grBaseClusterCenters.at(i)){
      delete grBaseClusterCenters.at(i);
      grBaseClusterCenters.at(i) = NULL;
    }
  }
  
  for(UInt_t i=0; i < hNonBaseClusteredEvents.size(); i++){
    if(hNonBaseClusteredEvents.at(i)){
      delete hNonBaseClusteredEvents.at(i);
      hNonBaseClusteredEvents.at(i) = NULL;
    }
  }

  for(UInt_t i=0; i < grNonBaseClusterCenters.size(); i++){
    if(grNonBaseClusterCenters.at(i)){
      delete grNonBaseClusterCenters.at(i);
      grNonBaseClusterCenters.at(i)= NULL;
    }
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

  TString name = TString::Format("hNonBaseClusteredEvents_%lu", clusters.size());
  // const int nBinsLat = 30;
  // const int nBinsLon = 360;
  // const int nBins = 1024;
  // amp->addHistogram(name, name, nBins, nBins);

  TH2DAntarctica* h = new TH2DAntarctica(name, name);

  hNonBaseClusteredEvents.push_back(h);

  // find unclustered points
  ampBinNumbers.clear(); // doesn't clear memory
  ampBinNumbers.resize((Int_t) points.size()); // so this shouldn't reallocate memory
  for(int pointInd=0; pointInd < (Int_t) points.size(); pointInd++){
    const Point& point = points.at(pointInd);
    if(point.inCluster < 0){
      ampBinNumbers.at(pointInd) = h->Fill(point.longitude, point.latitude);
    }
  }
  globalMaxBin = hNonBaseClusteredEvents.back()->GetMaximumBin();

  Int_t maxEventsInBin = hNonBaseClusteredEvents.back()->GetBinContent(globalMaxBin);

  return maxEventsInBin;

}



void Acclaim::Clustering::LogLikelihoodMethod::recursivelyAddClustersFromData(Int_t minBinContent){

  static int lastIntegral = 0;

  // in case this hasn't been done already
  if(doneBaseClusterAssignment==false){
    assignEventsToBaseClusters();
  }
  
  numCallsToRecursive++; // for debugging
  

  Int_t globalMaxBin;
  Int_t maxBinVal = histogramUnclusteredEvents(globalMaxBin);

  
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

  if(counter > minBinContent){ // then we create a new cluster and do the whole thing again

    for(int dim=0; dim < nDim; dim++){
      cluster.centre[dim] /= counter;
    }
    AnitaGeomTool* geom = AnitaGeomTool::Instance();
    geom->getLatLonAltFromCartesian(cluster.centre, cluster.latitude,
				    cluster.longitude, cluster.altitude);
    clusters.push_back(cluster);
    
    numClusters = (Int_t) clusters.size();

    const int isMC = 0;
    for(int pointInd=0; pointInd < (Int_t) points.size(); pointInd++){
      assignSinglePointToCloserCluster(pointInd, isMC, numClusters-1);
    }

    const int numEventsUnclustered = hNonBaseClusteredEvents.back()->Integral();
    if(numEventsUnclustered > 0 && lastIntegral != numEventsUnclustered){
      lastIntegral = hNonBaseClusteredEvents.back()->Integral();

      TGraphAntarctica* gr = new TGraphAntarctica(1, &cluster.longitude, &cluster.latitude);
      TString name = TString::Format("grClusterCenters%lu", clusters.size()-1);
      TString title = TString::Format("Cluster %lu", clusters.size()-1);
      gr->SetNameTitle(name, title);
      grNonBaseClusterCenters.push_back(gr);

      recursivelyAddClustersFromData(minBinContent); // recursion
    }
    else{
      delete hNonBaseClusteredEvents.back();
      hNonBaseClusteredEvents.pop_back();
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

  Long64_t nMc = mcPoints.size();

  std::cerr  << "Info in " << __PRETTY_FUNCTION__ << ": starting!" << std::endl;
  ProgressBar p(nMc);  
  for(Long64_t j=0; j < nMc; j++){
    for(int clusterInd=0; clusterInd < (Int_t)clusters.size(); clusterInd++){
      assignSinglePointToCloserCluster(j, isMC, clusterInd);
    }
    if(mcPoints.at(j).inCluster < 0){
      numSinglets+=mcPoints.at(j).weight;
    }
    p.inc(j, nMc);
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







void Acclaim::Clustering::LogLikelihoodMethod::assignEventsToBaseClusters(){
  std::cout << "Info in " << __PRETTY_FUNCTION__ << " assigning points to closest cluster!" << std::endl;
  ProgressBar p(numClusters);
  const int isMC = 0;

  
  for(Long64_t clusterInd=0; clusterInd < numClusters; clusterInd++){
    for(int pointInd=0; pointInd < (int) points.size(); pointInd++){
      assignSinglePointToCloserCluster(pointInd, isMC, clusterInd);
    }
    p.inc(clusterInd, numClusters);
  }

  
  hBaseClusteredEvents.resize(numClusters, NULL);
  grBaseClusterCenters.resize(numClusters, NULL);

  for(Long64_t clusterInd=0; clusterInd < numClusters; clusterInd++){
    if(clusters.at(clusterInd).numEvents > 0){
      TString name = TString::Format("hBaseClusteredEvents%d", (int)clusterInd);
      TH2DAntarctica* h = new TH2DAntarctica(name, name);
      hBaseClusteredEvents[clusterInd] = h;      

      for(int pointInd=0; pointInd < (int)points.size(); pointInd++){
	if(points.at(pointInd).inCluster == clusterInd){
	  h->Fill(points.at(pointInd).longitude, points.at(pointInd).latitude);
	}
      }

      grBaseClusterCenters.at(clusterInd) = new TGraphAntarctica(BaseList::getBase(clusterInd));
      grBaseClusterCenters.at(clusterInd)->SetName(TString::Format("grClusterCenter%d", (int)clusterInd));
    }
  }
  
  doneBaseClusterAssignment = true;
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
  doneBaseClusterAssignment = false;
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
	gr->SetPoint(gr->GetN(), points.at(pointInd).longitude, points.at(pointInd).latitude);
      }
    }
  }
  return gr;
}





void Acclaim::Clustering::LogLikelihoodMethod::makeSummaryTrees(){

  TTree* clusteredDataTree = new TTree("clusteredDataTree", "Tree of clustered ANITA events");
  Point* point = NULL;
  clusteredDataTree->Branch("point", &point);

  for(Int_t i=0; i < (Int_t)points.size(); i++){
    point = &points.at(i);
    clusteredDataTree->Fill();
  }
  
  TTree* clusteredMcTree = new TTree("clusteredMcTree", "Tree of clustered Monte Carlo ANITA events");
  McPoint* mcPoint = NULL;
  clusteredMcTree->Branch("mcPoint", &mcPoint);
  for(Int_t j=0; j < (Int_t)mcPoints.size(); j++){
    mcPoint = &mcPoints.at(j);
    clusteredMcTree->Fill();
  }

  TTree* clusterTree = new TTree("clusterTree", "Tree of clusters");
  Cluster* cluster = NULL;
  clusterTree->Branch("cluster", &cluster);
  for(int k=0; k < (Int_t)clusters.size(); k++){
    cluster = &clusters.at(k);
    clusterTree->Fill();
  }
}


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
	Double_t sigmaTheta = default_sigma_theta;	
	Double_t sigmaPhi = default_sigma_phi;

	if(sum->mc.weight > 0){
	  addMcPoint(&pat, sourceLat, sourceLon, sourceAlt,
		     sum->run, sum->eventNumber, thetaDeg, phiDeg,
		     sigmaTheta, sigmaPhi, pol, sum->mc.weight, sum->mc.energy);
	
	}
	else{
	  addPoint(&pat, sourceLat,sourceLon,sourceAlt, sum->run, sum->eventNumber,
	  	   thetaDeg, phiDeg, sigmaTheta, sigmaPhi, pol);
	}
      }
      p.inc(entry, n);
    }
  }
  return n;
}


void Acclaim::Clustering::LogLikelihoodMethod::writeAllGraphsAndHists(){
  for(UInt_t i=0; i < hBaseClusteredEvents.size(); i++){
    if(hBaseClusteredEvents.at(i)){
      hBaseClusteredEvents.at(i)->Write();
    }
  }
  for(UInt_t i=0; i < grBaseClusterCenters.size(); i++){
    if(grBaseClusterCenters.at(i)){
      grBaseClusterCenters.at(i)->Write();
    }
  }  
  for(UInt_t i=0; i < hNonBaseClusteredEvents.size(); i++){
    if(hNonBaseClusteredEvents.at(i)){
      hNonBaseClusteredEvents.at(i)->Write();
    }
  }
  for(UInt_t i=0; i < grNonBaseClusterCenters.size(); i++){
    if(grNonBaseClusterCenters.at(i)){
      grNonBaseClusterCenters.at(i)->Write();
    }
  }
}




void Acclaim::Clustering::LogLikelihoodMethod::doClustering(const char* dataGlob, const char* mcGlob, const char* outFileName){  

  readInSummaries(dataGlob);
  initializeBaseList();
  readInSummaries(mcGlob);

  char* fakeArgv0 = const_cast<char*>(outFileName);
  OutputConvention oc(1, &fakeArgv0);
  TFile* fOut = oc.makeFile();

  assignEventsToBaseClusters();
  recursivelyAddClustersFromData(0);


  assignMcPointsToClusters();

  

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
  //   recursivelyAddClustersFromData(0);

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

  writeAllGraphsAndHists();
  makeSummaryTrees();

  fOut->Write();
  fOut->Close();

}
