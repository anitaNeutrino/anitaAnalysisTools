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
#include "TF1.h"
#include "TCanvas.h"

ClassImp(Acclaim::Clustering::Event);
ClassImp(Acclaim::Clustering::McEvent);
ClassImp(Acclaim::Clustering::Cluster);

const int nDim = 3;


/**
 * @namespace ResolutionModel
 * @brief Parameters defining the resolution model
 * 
 * See the macro plotCalPulserResolution.C for the derivation of these numbers
 */
namespace ResolutionModel{
  const int n = 3;
  const double phiParams[n]   = {-2.50414e-01,  3.02406e-01, 2.43376e-01};
  const double thetaParams[n] = {-3.83773e-01, -3.00964e-01, 1.64537e-01};
  TString formula = "exp([0]*x + [1]) + [2]";
}



/**
 * @brief Wrapper function to calculate the angular resolution for clustering
 *
 * @param sum is the AnitaEventSummary
 * @param pol the polarisation of interest
 * @param peakInd the peak of the map of interest
 * @param sigma_theta the calculated theta resolution (degrees)
 * @param sigma_phi the calculated phi resolution (degrees)
 */
void Acclaim::Clustering::getAngularResolution(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd, double& sigma_theta, double& sigma_phi){
  const double x = sum->coherent_filtered[pol][peakInd].snr;
  getAngularResolution(x, sigma_theta, sigma_phi);
}


/**
 * @brief Calculate the angular resolution for clustering
 * @todo Currently derived from WAIS pulses, but should probably be from MC
 * 
 * 
 * @param x the parameterization variable
 * @param sigma_theta the calculated theta resolution (degrees)
 * @param sigma_phi the calculated phi resolution (degrees)
 */
void Acclaim::Clustering::getAngularResolution(double x, double& sigma_theta, double& sigma_phi){

  sigma_phi = exp(ResolutionModel::phiParams[0]*x + ResolutionModel::phiParams[1]) + ResolutionModel::phiParams[2];

  sigma_theta = exp(ResolutionModel::thetaParams[0]*x + ResolutionModel::thetaParams[1]) + ResolutionModel::thetaParams[2];

}

TCanvas* Acclaim::Clustering::drawAngularResolutionModel(double maxSnr){
  TCanvas* c1 = new TCanvas();

  TF1* fTheta = new TF1("fThetaResolutionModel", ResolutionModel::formula, 0, maxSnr);
  TF1* fPhi = new TF1("fThetaResolutionModel", ResolutionModel::formula, 0, maxSnr);

  for(int i=0; i < ResolutionModel::n; i++){
    fTheta->SetParameter(i, ResolutionModel::thetaParams[i]);
    fPhi->SetParameter(i, ResolutionModel::phiParams[i]);
  }

  fPhi->Draw();
  fPhi->SetLineColor(kRed);
  fTheta->Draw("lsame");
  fTheta->SetLineColor(kBlue);
  fPhi->SetBit(kCanDelete);
  fTheta->SetBit(kCanDelete);

  fPhi->SetMinimum(0.01);
  c1->Modified();
  c1->Update();
  return c1;
}

Acclaim::Clustering::Event::Event(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd){

  this->eventNumber = sum->eventNumber;
  this->run = sum->run;

  this->pol = pol;
  peakIndex = peakInd;
  const AnitaEventSummary::PointingHypothesis& peak = sum->peak[pol][peakInd];
  latitude = peak.latitude;
  longitude = peak.longitude;
  altitude = peak.altitude;
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->getCartesianCoords(latitude, longitude, altitude, centre);
  theta = peak.theta;
  phi = peak.phi;
  anita = sum->anitaLocation;

  // sigmaTheta = default_sigma_theta;
  // sigmaPhi = default_sigma_phi;
  getAngularResolution(sum, pol, peakInd, sigmaTheta, sigmaPhi);

  logLikelihood = DBL_MAX;
  cluster = -1;
  logLikelihood2 = DBL_MAX;
  cluster2 = -1;

  antarcticaHistBin = -1;

  dThetaCluster = -999;
  dPhiCluster = -999;
}


Acclaim::Clustering::Event::Event(){
  latitude = 0;
  longitude = 0;
  altitude = 0;
  theta = -9999;
  phi = -9999;
  anita.reset();

  eventNumber = 0;
  run = 0;
  pol = AnitaPol::kNotAPol;
  peakIndex = -1;

  // convert to degrees
  sigmaTheta = default_sigma_theta;
  sigmaPhi = default_sigma_phi;
  logLikelihood = DBL_MAX;
  cluster = -1;
  logLikelihood2 = DBL_MAX;
  cluster2 = -1;

  for(int dim=0; dim < nDim; dim++){
    centre[dim] = 0;
  }
  antarcticaHistBin = -1;
  dThetaCluster = -999;
  dPhiCluster = -999;

}



Acclaim::Clustering::McEvent::McEvent()
  : Event(){
  weight = 0; energy=0;
}



Acclaim::Clustering::Cluster::Cluster() {
  // std::cerr << __PRETTY_FUNCTION__ << std::endl;
  for(int dim=0; dim < nDim; dim++){
    centre[dim] = 0;
  }
  latitude = 0;
  longitude = 0;
  altitude = 0;
  knownBase = 0;
  numDataEvents = 0;
  sumMcWeights = 0;
  antarcticaHistBin = -1;
  seedEvent = -1;
}

Acclaim::Clustering::Cluster::Cluster(const BaseList::base& base) {
  // std::cerr << __PRETTY_FUNCTION__ << std::endl;
  AntarcticCoord ac = base.position.as(AntarcticCoord::WGS84);
  latitude = ac.x;
  longitude = ac.y;
  altitude = ac.z;
  knownBase = 1;

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->getCartesianCoords(latitude, longitude, altitude, centre);
  numDataEvents = 0;
  sumMcWeights = 0;
  antarcticaHistBin = -1;
  seedEvent = -1;
}

Acclaim::Clustering::LogLikelihoodMethod::LogLikelihoodMethod(){

  llCut = 100;
  maxDistCluster = 800e3; // try 800km
  numCallsToRecursive = 0;
  doneBaseClusterAssignment = false;
  hClusters = new TH2DAntarctica("hClusters", "hClusters");
  fDebug = false;
  fUseBaseList = true;
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

  if(hClusters){
    delete hClusters;
    hClusters = NULL;
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
Double_t Acclaim::Clustering::LogLikelihoodMethod::getDistSqEventCluster(Int_t eventInd, const Acclaim::Clustering::Cluster& cluster){
  // const Cluster& cluster = clusters.at(clusterInd);
  const Event& event = events.at(eventInd);
  Double_t d2=0;
  for(int dim=0; dim < nDim; dim++){
    Double_t d = cluster.centre[dim] - event.centre[dim];
    d2 += d*d;
  }
  return d2;
}


// Double_t Acclaim::Clustering::LogLikelihoodMethod::getDistSqClusterCluster(Int_t clusterInd1, Int_t clusterInd2){
//   const Cluster& cluster1 = clusters.at(clusterInd1);
//   const Cluster& cluster2 = clusters.at(clusterInd2);
//   Double_t d2=0;
//   for(int dim=0; dim < nDim; dim++){
//     Double_t d = cluster1.centre[dim] - cluster2.centre[dim];
//     d2 += d*d;
//   }
//   return d2;
// }




void Acclaim::Clustering::LogLikelihoodMethod::getDeltaThetaDegDeltaPhiDegEventCluster(Int_t eventInd, Int_t clusterInd, UsefulAdu5Pat& usefulPat, Double_t& deltaThetaDeg, Double_t& deltaPhiDeg){

  const Cluster& cluster = clusters.at(clusterInd);
  const Event& event = events.at(eventInd);

  Double_t thetaWave, phiWave;
  usefulPat.getThetaAndPhiWave(cluster.longitude, cluster.latitude, cluster.altitude, thetaWave, phiWave);
  Double_t thetaDeg = -1*TMath::RadToDeg()*thetaWave;
  Double_t phiDeg = TMath::RadToDeg()*phiWave;

  deltaThetaDeg = (thetaDeg - event.theta);
  deltaPhiDeg = Acclaim::RootTools::getDeltaAngleDeg(phiDeg, event.phi);

  if(fDebug){
    // if(event.cluster < 0 && event.antarcticaHistBin == cluster.antarcticaHistBin){
    if(event.cluster < 0 && eventInd == cluster.seedEvent){
      std::cout << "event theta = " << event.theta  << ", cluster theta = " << thetaDeg << std::endl;
      std::cout << "event phi = " << event.phi  << ", cluster phi = " << phiDeg << std::endl;
    }
  }
}










Double_t Acclaim::Clustering::LogLikelihoodMethod::getAngDistSqEventCluster(Int_t eventInd, Int_t clusterInd, UsefulAdu5Pat& usefulPat){

  const Cluster& cluster = clusters.at(clusterInd);
  const Event& event = events.at(eventInd);

  Double_t deltaThetaDeg, deltaPhiDeg;
  getDeltaThetaDegDeltaPhiDegEventCluster(eventInd, clusterInd, usefulPat, deltaThetaDeg, deltaPhiDeg);

  Double_t dThetaNorm = deltaThetaDeg/event.sigmaTheta;
  Double_t dPhiNorm = deltaPhiDeg/event.sigmaPhi;
  Double_t angSq =  dThetaNorm*dThetaNorm + dPhiNorm*dPhiNorm;
  if(fDebug){
    // if(event.cluster < 0 && event.antarcticaHistBin == cluster.antarcticaHistBin){
    if(event.cluster < 0 && eventInd == cluster.seedEvent){
      std::cout << event.eventNumber << std::endl;
      std::cout << "dTheta = " << deltaThetaDeg  << ", dPhi = " << deltaPhiDeg
		<< ", sigmaTheta = " << event.sigmaTheta << ", sigmaPhi = " << event.sigmaPhi
		<< ", dThetaNorm = " << dThetaNorm << ", dPhiNorm = " << dPhiNorm
		<< ", ll = " << angSq << std::endl;
      std::cout << "cluster lon/lat/alt = " << cluster.longitude << ", " << cluster.latitude << ", " << cluster.altitude << std::endl;
      std::cout << "event lon/lat/alt = " << event.longitude << ", " << event.latitude << ", " << event.altitude << std::endl;
      std::cout << "anita lon/lat/alt = " << usefulPat.longitude << ", " << usefulPat.latitude << ", " << usefulPat.altitude << std::endl;
    }
  }


  return angSq;
}


// Double_t Acclaim::Clustering::LogLikelihoodMethod::getAngDistSqEventCluster(Int_t eventInd, Int_t clusterInd, const Adu5Pat* pat){
//   UsefulAdu5Pat usefulPat(pat);
//   return getAngDistSqEventCluster(event, cluster, usefulPat);
// }





Double_t Acclaim::Clustering::LogLikelihoodMethod::getSumOfMcWeights(){
  Double_t sumOfWeights = 0;
  for(int i=0; i < (int)mcEvents.size(); i++){
    sumOfWeights += mcEvents.at(i).weight;
  }
  return sumOfWeights;
}



Int_t Acclaim::Clustering::LogLikelihoodMethod::histogramUnclusteredEvents(Int_t& globalMaxBin){

  TString name = TString::Format("hNonBaseClusteredEvents_%lu", clusters.size());
  TH2DAntarctica* h = new TH2DAntarctica(name, name);
  hNonBaseClusteredEvents.push_back(h);

  // find unclustered events
  for(int i=0; i < (Int_t) events.size(); i++){
    Event& event = events.at(i);
    if(event.cluster < 0){
      event.antarcticaHistBin = h->Fill(event.longitude, event.latitude);
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

  if(numCallsToRecursive==0){
    std::cout << "Beginning " << __PRETTY_FUNCTION__ << std::endl;
  }

  numCallsToRecursive++; // for debugging

  Int_t globalMaxBin;
  Int_t maxBinVal = histogramUnclusteredEvents(globalMaxBin);

  if(fDebug){
    std::cout << "maxBinVal = " << maxBinVal << std::endl;
  }

  Int_t counter=0;
  Cluster tempCluster;
  for(UInt_t i=0; i < events.size(); i++){
    const Event& event = events.at(i);
    if(event.cluster < 0 && event.antarcticaHistBin==globalMaxBin){
      for(int dim=0; dim < nDim; dim++){
	tempCluster.centre[dim] += event.centre[dim];
      }
      if(fDebug){
	std::cout << "The " << counter << "th event is at "
		  << "lon = " << event.longitude
		  << ", lat = " << event.latitude
		  << ", alt = " << event.altitude
		  << std::endl;
      }
      counter++;
    }
  }
  if(fDebug){
    std::cout << "There are " << counter << " events in the max bin " << std::endl;
  }

  if(counter!=maxBinVal){
    std::cerr << "Now what?" << std::endl;
  }

  if(counter > minBinContent){ // then we create a new cluster and do the whole thing again


    for(int dim=0; dim < nDim; dim++){
      tempCluster.centre[dim] /= counter;
    }



    AnitaGeomTool* geom = AnitaGeomTool::Instance();
    geom->getLatLonAltFromCartesian(tempCluster.centre,
				    tempCluster.latitude,
				    tempCluster.longitude,
				    tempCluster.altitude);

    // now loop over events and pick point closest to mean position of events
    Int_t eventIndClosestToBinMean = -1;
    Double_t minDSq = DBL_MAX;
    for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
      const Event& event = events.at(eventInd);
      if(event.cluster < 0 && event.antarcticaHistBin==globalMaxBin){
	Double_t dSq = getDistSqEventCluster(eventInd, tempCluster);
	if(dSq < minDSq){
	  minDSq = dSq;
	  eventIndClosestToBinMean = eventInd;
	}
      }
    }

    Cluster seedCluster;
    const Event& seedEvent = events.at(eventIndClosestToBinMean);
    seedCluster.latitude = seedEvent.latitude;
    seedCluster.longitude = seedEvent.longitude;
    seedCluster.altitude = seedEvent.altitude;
    seedCluster.seedEvent = eventIndClosestToBinMean;

    seedCluster.antarcticaHistBin = hClusters->Fill(seedCluster.longitude, seedCluster.latitude);

    clusters.push_back(seedCluster);

    if(fDebug){
      std::cout << "The average cluster location is "
		<< "lon = " << tempCluster.longitude
		<< ", lat = " << tempCluster.latitude
		<< ", alt = " << tempCluster.altitude
		<< std::endl;
      std::cout << "The seed cluster location is "
		<< "lon = " << seedCluster.longitude
		<< ", lat = " << seedCluster.latitude
		<< ", alt = " << seedCluster.altitude
		<< std::endl;

    }

    const int isMC = 0;
    for(int eventInd=0; eventInd < (Int_t) events.size(); eventInd++){
      assignSingleEventToCloserCluster(eventInd, isMC, clusters.size()-1);
    }

    const int numDataEventsUnclustered = hNonBaseClusteredEvents.back()->Integral();
    std::cout << "After " << numCallsToRecursive << ", " << numDataEventsUnclustered << " events are unclustered" << std::endl;
    if(numDataEventsUnclustered > 0 && lastIntegral != numDataEventsUnclustered){
      lastIntegral = hNonBaseClusteredEvents.back()->Integral();

      TGraphAntarctica* gr = new TGraphAntarctica(1, &seedCluster.longitude, &seedCluster.latitude);
      TString name = TString::Format("grClusterCenters%lu", clusters.size()-1);
      TString title = TString::Format("Cluster %lu", clusters.size()-1);
      gr->SetNameTitle(name, title);
      grNonBaseClusterCenters.push_back(gr);

      recursivelyAddClustersFromData(minBinContent); // recursion
    }
    else{
      delete hNonBaseClusteredEvents.back();
      hNonBaseClusteredEvents.pop_back();

      if(numDataEventsUnclustered > minBinContent && !fDebug){ // one last time through with debug on..
	fDebug = true;
	recursivelyAddClustersFromData(minBinContent); // recursion
	delete hNonBaseClusteredEvents.back();
	hNonBaseClusteredEvents.pop_back();
	clusters.pop_back();
      }

      if(fDebug){
	fDebug = false;
      }

      if(numDataEventsUnclustered > 0){
	std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", the clustering algorithm didn't manage to cluster "
		  << numDataEventsUnclustered << " events" << std::endl;
      }
    }
  }
}



void Acclaim::Clustering::LogLikelihoodMethod::redoSmallClusters(){
  std::cout << "In " << __PRETTY_FUNCTION__ << std::endl;

  // here we uncluster the small clusters...
  Int_t numEventsReset = 0;
  Int_t numClustersReset = 0;
  for(UInt_t clusterInd=0; clusterInd < clusters.size(); clusterInd++){
    Cluster& cluster = clusters.at(clusterInd);
    if(cluster.numDataEvents < SmallClusterSizeThreshold){
      cluster.numDataEvents = 0; // reset cluster data event counter
      numClustersReset++;

      for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
	if(events.at(eventInd).cluster==static_cast<Int_t>(clusterInd)){
	  // some reset function would be good here...

	  events.at(eventInd).logLikelihood = DBL_MAX;
	  events.at(eventInd).cluster = -1;
	  events.at(eventInd).logLikelihood2 = DBL_MAX;
	  events.at(eventInd).cluster2 = -1;
	  numEventsReset++;
	}
      }
    }
  }

  std::cout << numEventsReset << " events were reset over " << numClustersReset << " clusters!" << std::endl;
  // go through and empty non-base clusters?

  // now do the nSquared thing
  for(int eventInd=0; eventInd < (int)events.size(); eventInd++){
    if(events.at(eventInd).cluster < 0){
      const Int_t isMC = 0;
      for(int clusterInd=0; clusterInd < (Int_t)clusters.size(); clusterInd++){
	assignSingleEventToCloserCluster(eventInd, isMC, clusterInd);
      }
    }
  }
}


void Acclaim::Clustering::LogLikelihoodMethod::findClosestEventToClustersOfSizeOne(){

  // loop through clusters and find the small ones
  for(Int_t clusterInd=0; clusterInd < (Int_t)clusters.size(); clusterInd++){
    const Int_t maxRetestClusterSize = 10; ///@todo review this function
    if(clusters.at(clusterInd).numDataEvents > 0 && clusters.at(clusterInd).numDataEvents <= maxRetestClusterSize){


      Cluster& cluster = clusters.at(clusterInd);

      Int_t numNearEvents = 0;

      // now loop through and see how close this isolated cluster is close to any other event...
      Double_t minLL = 9e99;
      Int_t minEventInd = -1;
      for(int eventInd=0; eventInd < (int)events.size(); eventInd++){

	// skip events in the cluster under test
	if(events.at(eventInd).cluster!=clusterInd){
	  Adu5Pat pat = events.at(eventInd).anita.pat();
	  UsefulAdu5Pat usefulPat(&pat);
	  Double_t dist = usefulPat.getDistanceFromSource(cluster.latitude,
							  cluster.longitude,
							  cluster.altitude);

	  // this is where proper collision detection would be useful...
	  if(dist < maxDistCluster){
	    Double_t ll = getAngDistSqEventCluster(eventInd, clusterInd, usefulPat);

	    if(ll < llCut){

	      numNearEvents++;
	      // cluster.numEventsWithinMinLL++;
	      if(ll < minLL){
		minLL = ll;
		minEventInd = eventInd;
	      }
	    }
	  }
	}
      }

      if(numNearEvents > 0){
	std::cout << "reassign " << clusterInd << "\t" << cluster.numDataEvents << " to " << events.at(minEventInd).cluster << "\t" << clusters.at(events.at(minEventInd).cluster).numDataEvents << std::endl;
	// re-assign events...
	for(int eventInd=0; eventInd < (int)events.size(); eventInd++){
	  if(events.at(eventInd).cluster==clusterInd){
	    clusters.at(events.at(minEventInd).cluster).numDataEvents++;
	    cluster.numDataEvents--;
	    events.at(eventInd).cluster = events.at(minEventInd).cluster;
	    events.at(eventInd).cluster2 = events.at(minEventInd).cluster2;
	  }
	}
      }
      // else{
      // 	std::cout << "don't reassign " << clusterInd << "\t" << cluster.numDataEvents << std::endl;
      // }
    }
  }
}







void Acclaim::Clustering::LogLikelihoodMethod::assignSingleEventToCloserCluster(Int_t i, Int_t isMC, Int_t clusterInd){
  Event& event = isMC==0 ? events.at(i) : mcEvents.at(i);
  const Adu5Pat pat = event.anita.pat();
  UsefulAdu5Pat usefulPat(&pat);

  Cluster& cluster = clusters.at(clusterInd);

  Double_t distM = usefulPat.getDistanceFromSource(cluster.latitude, cluster.longitude, cluster.altitude);

  if(distM < maxDistCluster){ // are we even close?

    Double_t ll = getAngDistSqEventCluster(i, clusterInd, usefulPat);

    if(ll < llCut){ // within the cut?

      if(fDebug){
	if(event.cluster < 0){
	  std::cout << "clusterInd = " << clusterInd << ", i = " << i << ", ll = " << ll << ", event.ll = " << event.logLikelihood << ", event.cluster = " << event.cluster << std::endl;
	}
      }

      // if(!isMC){
      // 	std::cout << ll << std::endl;
      // }

      if(ll < event.logLikelihood){ // the best event?

	event.cluster2 = event.cluster;
	event.logLikelihood2 = event.logLikelihood;

	event.logLikelihood = ll;
	event.cluster = clusterInd;

	if(isMC==0){
	  if(event.cluster2 >= 0){
	    clusters.at(event.cluster2).numDataEvents--;
	  }
	  cluster.numDataEvents++;
	}
      }
    }
  }
}




void Acclaim::Clustering::LogLikelihoodMethod::assignMcEventsToClusters(){

  const int isMC = 1;
  Double_t numSinglets = 0;

  Long64_t nMc = mcEvents.size();

  std::cerr  << "Info in " << __PRETTY_FUNCTION__ << ": starting!" << std::endl;
  ProgressBar p(nMc);
  for(Long64_t j=0; j < nMc; j++){
    for(int clusterInd=0; clusterInd < (Int_t)clusters.size(); clusterInd++){
      assignSingleEventToCloserCluster(j, isMC, clusterInd);
    }
    if(mcEvents.at(j).cluster < 0){
      numSinglets+=mcEvents.at(j).weight;
    }
    p.inc(j, nMc);
  }
}









size_t Acclaim::Clustering::LogLikelihoodMethod::addMcEvent(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd){

  mcEvents.push_back(McEvent(sum,  pol, peakInd));

  return mcEvents.size();
}






size_t Acclaim::Clustering::LogLikelihoodMethod::addEvent(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd){

  events.push_back(Event(sum, pol, peakInd));

  return events.size();
}







void Acclaim::Clustering::LogLikelihoodMethod::assignEventsToBaseClusters(){
  std::cout << "Info in " << __PRETTY_FUNCTION__ << " assigning events to base clusters!" << std::endl;
  ProgressBar p(clusters.size());
  const int isMC = 0;
  for(Long64_t clusterInd=0; clusterInd < (Long64_t)clusters.size(); clusterInd++){
    for(int eventInd=0; eventInd < (int) events.size(); eventInd++){
      assignSingleEventToCloserCluster(eventInd, isMC, clusterInd);
    }

    clusters.at(clusterInd).antarcticaHistBin = hClusters->Fill(clusters.at(clusterInd).longitude, clusters.at(clusterInd).latitude);
    p.inc(clusterInd, clusters.size());
  }


  hBaseClusteredEvents.resize(clusters.size(), NULL);
  grBaseClusterCenters.resize(clusters.size(), NULL);

  for(Int_t clusterInd=0; clusterInd < (Int_t)clusters.size(); clusterInd++){
    if(clusters.at(clusterInd).numDataEvents > 0){
      TString name = TString::Format("hBaseClusteredEvents%d", (int)clusterInd);
      TH2DAntarctica* h = new TH2DAntarctica(name, name);
      hBaseClusteredEvents[clusterInd] = h;

      for(int eventInd=0; eventInd < (int)events.size(); eventInd++){
	if(events.at(eventInd).cluster == clusterInd){
	  h->Fill(events.at(eventInd).longitude, events.at(eventInd).latitude);
	}
      }

      grBaseClusterCenters.at(clusterInd) = new TGraphAntarctica(BaseList::getBase(clusterInd));
      grBaseClusterCenters.at(clusterInd)->SetName(TString::Format("grClusterCenter%d", (int)clusterInd));
    }
  }

  doneBaseClusterAssignment = true;
}


void Acclaim::Clustering::LogLikelihoodMethod::initializeEmptyBaseList(){
  BaseList::makeEmptyBaseList();
}


void Acclaim::Clustering::LogLikelihoodMethod::initializeBaseList(){

  std::cout << "Info in " << __PRETTY_FUNCTION__ << ": Initializing base list..." << std::endl;
  for(UInt_t clusterInd=0; clusterInd < BaseList::getNumBases(); clusterInd++){
    const BaseList::base& base = BaseList::getBase(clusterInd);
    clusters.push_back(Cluster(base));
  }
}



void Acclaim::Clustering::LogLikelihoodMethod::resetClusters(){

  for(int eventInd=0; eventInd < (int) events.size(); eventInd++){
    events.at(eventInd).logLikelihood = DBL_MAX;
    events.at(eventInd).logLikelihood2 = DBL_MAX;
    events.at(eventInd).cluster = -1;
    events.at(eventInd).cluster2 = -1;
  }

  for(int eventInd=0; eventInd < (int) mcEvents.size(); eventInd++){
    mcEvents.at(eventInd).logLikelihood = DBL_MAX;
    mcEvents.at(eventInd).logLikelihood2 = DBL_MAX;
    mcEvents.at(eventInd).cluster = -1;
    mcEvents.at(eventInd).cluster2 = -1;
  }

  const int numBases = BaseList::getNumBases();
  while((int) clusters.size() >= numBases && numBases>0){
    clusters.pop_back();
  }

  for(int clusterInd=0; clusterInd < (int) clusters.size(); clusterInd++){
    clusters.at(clusterInd).numDataEvents = 0;
  }
  doneBaseClusterAssignment = false;
}


































TGraphAntarctica* Acclaim::Clustering::LogLikelihoodMethod::makeClusterSummaryTGraph(Int_t clusterInd){

  TGraphAntarctica* gr = NULL;
  if(clusterInd >= 0 && clusterInd < (Int_t)clusters.size()){

    TString name  = TString::Format("grCluster%d", clusterInd);
    TString title  = TString::Format("Cluster %d; Easting (m); Northing (m)", clusterInd);
    gr = new TGraphAntarctica();
    gr->SetName(name);
    gr->SetTitle(title);

    // AnitaGeomTool* geom = AnitaGeomTool::Instance();

    for(int eventInd=0; eventInd < (Int_t) events.size(); eventInd++){
      if(events.at(eventInd).cluster==clusterInd){
	gr->SetPoint(gr->GetN(), events.at(eventInd).longitude, events.at(eventInd).latitude);
      }
    }
  }
  return gr;
}





void Acclaim::Clustering::LogLikelihoodMethod::makeSummaryTrees(){

  TTree* clusteredDataTree = new TTree("clusteredDataTree", "Tree of clustered ANITA events (all multiples and base singlets)");
  TTree* nonBaseSingletTree = new TTree("nonBaseSingletTree", "Tree of unclustered ANITA events (non-base singlets only)");
  Event* event = NULL;
  clusteredDataTree->Branch("event", &event);
  nonBaseSingletTree->Branch("event", &event);

  for(Int_t eventInd=0; eventInd < (Int_t)events.size(); eventInd++){
    event = &events.at(eventInd);
    if(event->cluster >= 0){
      Adu5Pat pat= event->anita.pat();
      UsefulAdu5Pat usefulPat(&pat);
      getDeltaThetaDegDeltaPhiDegEventCluster(eventInd, event->cluster, usefulPat,
					      event->dThetaCluster, event->dPhiCluster);

      if(clusters.at(event->cluster).numDataEvents == 1 && clusters.at(event->cluster).knownBase==0) {
	nonBaseSingletTree->Fill();
      }
      else{
	clusteredDataTree->Fill();
      }
    }
    else{
      // std::cerr << "Event " << event->eventNumber << " was unclustered!\n";
      // std::cerr << "It had theta = " << event->theta << ", phi = " << event->phi << "\n";
      // std::cerr << "lon = " << event->longitude << ", lat = " << event->latitude << std::endl;
    }
  }

  TTree* clusteredMcTree = new TTree("clusteredMcTree", "Tree of clustered Monte Carlo ANITA events");
  McEvent* mcEvent = NULL;
  clusteredMcTree->Branch("mcEvent", &mcEvent);

  for(Int_t j=0; j < (Int_t) mcEvents.size(); j++){
    mcEvent = &mcEvents.at(j);
    clusteredMcTree->Fill();
  }

  TTree* clusterTree = new TTree("clusterTree", "Tree of clusters");
  Cluster* cluster = NULL;
  clusterTree->Branch("cluster", &cluster);
  for(Int_t k=0; k < (Int_t)clusters.size(); k++){
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

      // AnitaPol::AnitaPol_t pol = sum->trainingPol();
      // Int_t peakIndex = sum->trainingPeakInd();
      AnitaPol::AnitaPol_t pol = sum->mostImpulsivePol(1);
      Int_t peakIndex = sum->mostImpulsiveInd(1);

      // std::cout << sum->eventNumber << "\t" << pol << "\t" << peakIndex << std::endl;

      Double_t sourceLat = sum->peak[pol][peakIndex].latitude;
      Double_t sourceLon = sum->peak[pol][peakIndex].longitude;
      Double_t snrHack = sum->deconvolved_filtered[pol][peakIndex].snr;
      Double_t thetaAdjustment = sum->peak[pol][peakIndex].theta_adjustment_needed;

      /// @todo remove the snr < 100 hack!!!
      if(sourceLat > -999 && sourceLon > -999 && snrHack < 100 && TMath::Abs(thetaAdjustment) < 1e-10){

	if(sum->mc.weight > 0){
	  addMcEvent(sum,  pol, peakIndex);
	}
	else{
	  addEvent(sum, pol, peakIndex);
	}
      }
      p.inc(entry, n);
    }
  }
  return n;
}


void Acclaim::Clustering::LogLikelihoodMethod::writeAllGraphsAndHists(){
  for(UInt_t eventInd=0; eventInd < hBaseClusteredEvents.size(); eventInd++){
    if(hBaseClusteredEvents.at(eventInd)){
      hBaseClusteredEvents.at(eventInd)->Write();
    }
  }
  for(UInt_t eventInd=0; eventInd < grBaseClusterCenters.size(); eventInd++){
    if(grBaseClusterCenters.at(eventInd)){
      grBaseClusterCenters.at(eventInd)->Write();
    }
  }
  for(UInt_t eventInd=0; eventInd < hNonBaseClusteredEvents.size(); eventInd++){
    if(hNonBaseClusteredEvents.at(eventInd)){
      hNonBaseClusteredEvents.at(eventInd)->Write();
    }
  }
  for(UInt_t eventInd=0; eventInd < grNonBaseClusterCenters.size(); eventInd++){
    if(grNonBaseClusterCenters.at(eventInd)){
      grNonBaseClusterCenters.at(eventInd)->Write();
    }
  }
}




void Acclaim::Clustering::LogLikelihoodMethod::doClustering(const char* dataGlob, const char* mcGlob, const char* outFileName){

  readInSummaries(dataGlob);

  if(fUseBaseList){
    initializeBaseList();
  }
  else{
    std::cout << "Info in " << __PRETTY_FUNCTION__ << ": not using base list!" << std::endl;
  }
  readInSummaries(mcGlob);

  char* fakeArgv0 = const_cast<char*>(outFileName);
  OutputConvention oc(1, &fakeArgv0);
  TFile* fOut = oc.makeFile();

  assignEventsToBaseClusters();
  recursivelyAddClustersFromData(0);
  redoSmallClusters();

  assignMcEventsToClusters();

  writeAllGraphsAndHists();
  makeSummaryTrees();

  fOut->Write();
  fOut->Close();

}
