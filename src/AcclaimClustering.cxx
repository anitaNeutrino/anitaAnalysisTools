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

Acclaim::Clustering::Event::Event(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd){

  this->pol = pol;
  const AnitaEventSummary::PointingHypothesis& peak = sum->peak[pol][peakInd];
  latitude = peak.latitude;
  longitude = peak.longitude;
  altitude = peak.altitude;
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->getCartesianCoords(latitude, longitude, altitude, centre);
  theta = peak.theta;
  phi = peak.phi;
  anita = sum->anitaLocation;


  /// @todo make this a function of coherently summed waveform SNR
  sigmaTheta = default_sigma_theta;
  sigmaPhi = default_sigma_phi;

  ll = DBL_MAX;
  inCluster = -1;
  llSecondBest = DBL_MAX;
  secondClosestCluster = -1;

}


Acclaim::Clustering::Event::Event(){
  latitude = 0;
  longitude = 0;
  altitude = 0;
  theta = -9999;
  phi = -9999;

  // convert to degrees
  sigmaTheta = default_sigma_theta;
  sigmaPhi = default_sigma_phi;
  ll = DBL_MAX;
  inCluster = -1;
  llSecondBest = DBL_MAX;
  secondClosestCluster = -1;
  for(int dim=0; dim < nDim; dim++){
    centre[dim] = 0;
  }

}



Acclaim::Clustering::McEvent::McEvent()
  : Event(){
  weight = 0; energy=0;
}



Acclaim::Clustering::Cluster::Cluster() {
  for(int dim=0; dim < nDim; dim++){
    centre[dim] = 0;
  }
  numDataEvents = 0;
  latitude = 0;
  longitude = 0;
  altitude = 0;
}

Acclaim::Clustering::Cluster::Cluster(const BaseList::base& base) {
  latitude = base.latitude;
  longitude = base.longitude;
  altitude = base.altitude;

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->getCartesianCoords(latitude, longitude, altitude, centre);
  numDataEvents = 0;
}



Acclaim::Clustering::LogLikelihoodMethod::LogLikelihoodMethod(){

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
inline Double_t getDistSq(const Acclaim::Clustering::Event& event, const Acclaim::Clustering::Cluster& cluster){
  Double_t d2=0;
  for(int dim=0; dim < nDim; dim++){
    Double_t d = cluster.centre[dim] - event.centre[dim];
    d2 += d*d;
  }
  return d2;
}


// utility function hopefully this one gets inlined
inline Double_t getDistSq(Acclaim::Clustering::Cluster& cluster1, const Acclaim::Clustering::Cluster& cluster2){
  Double_t d2=0;
  for(int dim=0; dim < nDim; dim++){
    Double_t d = cluster1.centre[dim] - cluster2.centre[dim];
    d2 += d*d;
  }
  return d2;
}




inline void getDeltaThetaDegDeltaPhiDegCluster(const Acclaim::Clustering::Event& event, const Acclaim::Clustering::Cluster& cluster, UsefulAdu5Pat& usefulPat, Double_t& deltaThetaDeg, Double_t& deltaPhiDeg){

  Double_t thetaWave, phiWave;
  usefulPat.getThetaAndPhiWave(cluster.longitude, cluster.latitude, cluster.altitude, thetaWave, phiWave);
  Double_t thetaDeg = -1*TMath::RadToDeg()*thetaWave;
  Double_t phiDeg = TMath::RadToDeg()*phiWave;

  deltaThetaDeg = (thetaDeg - event.theta);
  deltaPhiDeg = Acclaim::RootTools::getDeltaAngleDeg(phiDeg, event.phi);
}






inline void getDeltaThetaDegDeltaPhiDegCluster(const Acclaim::Clustering::Event& event, const Acclaim::Clustering::Cluster& cluster, const Adu5Pat* pat, Double_t& deltaThetaDeg, Double_t& deltaPhiDeg){

  UsefulAdu5Pat usefulPat(pat);

  getDeltaThetaDegDeltaPhiDegCluster(event, cluster, usefulPat, deltaThetaDeg, deltaPhiDeg);
}







inline Double_t getAngDistSq(const Acclaim::Clustering::Event& event, const Acclaim::Clustering::Cluster& cluster, UsefulAdu5Pat& usefulPat){

  Double_t deltaThetaDeg, deltaPhiDeg;
  getDeltaThetaDegDeltaPhiDegCluster(event, cluster, usefulPat, deltaThetaDeg, deltaPhiDeg);

  Double_t dThetaNorm = deltaThetaDeg/event.sigmaTheta;
  Double_t dPhiNorm = deltaPhiDeg/event.sigmaPhi;
  Double_t angSq =  dThetaNorm*dThetaNorm + dPhiNorm*dPhiNorm;
  
  return angSq;
}


inline Double_t getAngDistSq(const Acclaim::Clustering::Event& event, const Acclaim::Clustering::Cluster& cluster, const Adu5Pat* pat){
  UsefulAdu5Pat usefulPat(pat);
  return getAngDistSq(event, cluster, usefulPat);
}





Double_t Acclaim::Clustering::LogLikelihoodMethod::getSumOfMcWeights(){
  Double_t sumOfWeights = 0;
  for(int eventInd=0; eventInd < (int)mcEvents.size(); eventInd++){
    sumOfWeights += mcEvents.at(eventInd).weight;
  }
  return sumOfWeights;
}



Int_t Acclaim::Clustering::LogLikelihoodMethod::histogramUnclusteredEvents(Int_t& globalMaxBin){

  TString name = TString::Format("hNonBaseClusteredEvents_%lu", clusters.size());
  TH2DAntarctica* h = new TH2DAntarctica(name, name);
  hNonBaseClusteredEvents.push_back(h);

  // find unclustered events
  ampBinNumbers.clear(); // doesn't clear memory
  ampBinNumbers.resize((Int_t) events.size()); // so this shouldn't reallocate memory
  for(int eventInd=0; eventInd < (Int_t) events.size(); eventInd++){
    const Event& event = events.at(eventInd);
    if(event.inCluster < 0){
      ampBinNumbers.at(eventInd) = h->Fill(event.longitude, event.latitude);
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
  for(int eventInd=0; eventInd < (Int_t) events.size(); eventInd++){
    const Event& event = events.at(eventInd);
    if(event.inCluster < 0){
      if(ampBinNumbers.at(eventInd)==globalMaxBin){
	for(int dim=0; dim < nDim; dim++){
	  cluster.centre[dim] += event.centre[dim];
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
    for(int eventInd=0; eventInd < (Int_t) events.size(); eventInd++){
      assignSingleEventToCloserCluster(eventInd, isMC, numClusters-1);
    }

    const int numDataEventsUnclustered = hNonBaseClusteredEvents.back()->Integral();
    if(numDataEventsUnclustered > 0 && lastIntegral != numDataEventsUnclustered){
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




void Acclaim::Clustering::LogLikelihoodMethod::findClosestEventToClustersOfSizeOne(){

  for(int i=0; i < maxRetestClusterSize; i++){
    numIsolatedSmallClusters.at(i) = 0;
    numIsolatedSmallBaseClusters.at(i) = 0;
  }

  // loop through clusters and find the small ones
  for(int clusterInd=0; clusterInd < numClusters; clusterInd++){
    if(clusters.at(clusterInd).numDataEvents > 0 && clusters.at(clusterInd).numDataEvents <= maxRetestClusterSize){


      Cluster& cluster = clusters.at(clusterInd);

      Int_t numNearEvents = 0;

      // now loop through and see how close this isolated cluster is close to any other event...
      Double_t minLL = 9e99;
      Int_t minEventInd = -1;
      for(int eventInd=0; eventInd < (int)events.size(); eventInd++){

	// skip events in the cluster under test
	if(events.at(eventInd).inCluster!=clusterInd){
	  Adu5Pat pat = events.at(eventInd).anita.pat();
	  UsefulAdu5Pat usefulPat(&pat);
	  Double_t dist = usefulPat.getDistanceFromSource(cluster.latitude,
							  cluster.longitude,
							  cluster.altitude);

	  // this is where proper collision detection would be useful...
	  if(dist < maxDistCluster){
	    Double_t ll = getAngDistSq(events.at(eventInd), cluster, usefulPat);

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
	std::cout << "reassign " << clusterInd << "\t" << cluster.numDataEvents << " to " << events.at(minEventInd).inCluster << "\t" << clusters.at(events.at(minEventInd).inCluster).numDataEvents << std::endl;
	// re-assign events...
	for(int eventInd=0; eventInd < (int)events.size(); eventInd++){
	  if(events.at(eventInd).inCluster==clusterInd){
	    clusters.at(events.at(minEventInd).inCluster).numDataEvents++;
	    cluster.numDataEvents--;
	    events.at(eventInd).inCluster = events.at(minEventInd).inCluster;
	    events.at(eventInd).secondClosestCluster = events.at(minEventInd).secondClosestCluster;
	  }
	}
      }
      // else{
      // 	std::cout << "don't reassign " << clusterInd << "\t" << cluster.numDataEvents << std::endl;
      // }
    }
  }

  const int numBases = BaseList::getNumBases();
  for(int clusterInd=0; clusterInd < numClusters; clusterInd++){
    if(clusters.at(clusterInd).numDataEvents > 0 &&
       clusters.at(clusterInd).numDataEvents <= maxRetestClusterSize){

      if(clusterInd < numBases){
	// if(clusters.at(clusterInd).numDataEvents==1){
	//   std::cout << clusterInd << "\t" << clusters.at(clusterInd).numDataEvents  << std::endl;
	// }
	numIsolatedSmallBaseClusters.at(clusters.at(clusterInd).numDataEvents-1)++;
      }
      else{
	numIsolatedSmallClusters.at(clusters.at(clusterInd).numDataEvents-1)++;
      }

    }
  }
}







void Acclaim::Clustering::LogLikelihoodMethod::assignSingleEventToCloserCluster(Int_t eventInd, Int_t isMC, Int_t clusterInd){
  Event& event = isMC==0 ? events.at(eventInd) : mcEvents.at(eventInd);
  const Adu5Pat pat = event.anita.pat();
  UsefulAdu5Pat usefulPat(&pat);

  Cluster& cluster = clusters.at(clusterInd);

  Double_t distM = usefulPat.getDistanceFromSource(cluster.latitude, cluster.longitude, cluster.altitude);

  if(distM < maxDistCluster){ // are we even close?

    Double_t ll = getAngDistSq(event, cluster, usefulPat);

    if(ll < llCut){ // within the cut?

      // if(!isMC){
      // 	std::cout << ll << std::endl;
      // }
      
      if(ll < event.ll){ // the best event?
 
	event.secondClosestCluster = event.inCluster;
	event.llSecondBest = event.ll;

	event.ll = ll;
	event.inCluster = clusterInd;

	if(isMC==0){
	  if(event.secondClosestCluster >= 0){
	    clusters.at(event.secondClosestCluster).numDataEvents--;
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
    if(mcEvents.at(j).inCluster < 0){
      numSinglets+=mcEvents.at(j).weight;
    }
    p.inc(j, nMc);
  }
  numMcIsolatedSinglets = numSinglets;
}









size_t Acclaim::Clustering::LogLikelihoodMethod::addMcEvent(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd){
// size_t Acclaim::Clustering::LogLikelihoodMethod::addMcEvent(Adu5Pat* pat, Double_t latitude, Double_t longitude, Double_t altitude, Int_t run, UInt_t eventNumber, Double_t thetaDeg, Double_t phiDeg, Double_t sigmaThetasDeg, Double_t sigmaPhiDeg, AnitaPol::AnitaPol_t pol, Double_t weight, Double_t energy){
  

  mcEvents.push_back(McEvent(sum,  pol, peakInd));

  return mcEvents.size();
}






size_t Acclaim::Clustering::LogLikelihoodMethod::addEvent(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd){

  events.push_back(Event(sum, pol, peakInd));

  return events.size();
}







void Acclaim::Clustering::LogLikelihoodMethod::assignEventsToBaseClusters(){
  std::cout << "Info in " << __PRETTY_FUNCTION__ << " assigning events to closest cluster!" << std::endl;
  ProgressBar p(numClusters);
  const int isMC = 0;

  
  for(Long64_t clusterInd=0; clusterInd < numClusters; clusterInd++){
    for(int eventInd=0; eventInd < (int) events.size(); eventInd++){
      assignSingleEventToCloserCluster(eventInd, isMC, clusterInd);
    }
    p.inc(clusterInd, numClusters);
  }

  
  hBaseClusteredEvents.resize(numClusters, NULL);
  grBaseClusterCenters.resize(numClusters, NULL);

  for(Long64_t clusterInd=0; clusterInd < numClusters; clusterInd++){
    if(clusters.at(clusterInd).numDataEvents > 0){
      TString name = TString::Format("hBaseClusteredEvents%d", (int)clusterInd);
      TH2DAntarctica* h = new TH2DAntarctica(name, name);
      hBaseClusteredEvents[clusterInd] = h;      

      for(int eventInd=0; eventInd < (int)events.size(); eventInd++){
	if(events.at(eventInd).inCluster == clusterInd){
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

  for(int eventInd=0; eventInd < (int) events.size(); eventInd++){
    events.at(eventInd).ll = DBL_MAX;
    events.at(eventInd).llSecondBest = DBL_MAX;
    events.at(eventInd).inCluster = -1;
    events.at(eventInd).secondClosestCluster = -1;
  }

  for(int eventInd=0; eventInd < (int) mcEvents.size(); eventInd++){
    mcEvents.at(eventInd).ll = DBL_MAX;
    mcEvents.at(eventInd).llSecondBest = DBL_MAX;
    mcEvents.at(eventInd).inCluster = -1;
    mcEvents.at(eventInd).secondClosestCluster = -1;
  }

  const int numBases = BaseList::getNumBases();
  while((int) clusters.size() >= numBases && numBases>0){
    clusters.pop_back();
  }
  
  numClusters = clusters.size();
  for(int clusterInd=0; clusterInd < (int) numClusters; clusterInd++){
    clusters.at(clusterInd).numDataEvents = 0;
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

    for(int eventInd=0; eventInd < (Int_t) events.size(); eventInd++){
      if(events.at(eventInd).inCluster==clusterInd){
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

  for(Int_t i=0; i < (Int_t)events.size(); i++){
    event = &events.at(i);

    if(clusters.at(event->inCluster).numDataEvents > 1 || event->inCluster < (Int_t)BaseList::getNumBases()){
      clusteredDataTree->Fill();
    }
    else{
      nonBaseSingletTree->Fill();
    }
  }
  
  TTree* clusteredMcTree = new TTree("clusteredMcTree", "Tree of clustered Monte Carlo ANITA events");
  McEvent* mcEvent = NULL;
  clusteredMcTree->Branch("mcEvent", &mcEvent);
  for(Int_t j=0; j < (Int_t)mcEvents.size(); j++){
    mcEvent = &mcEvents.at(j);
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
      if(sourceLat > -999 && sourceLon > -999){

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


  assignMcEventsToClusters();

  

  // TGraph* grNumEvents = new TGraph();
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
  //   findClosestEventToClustersOfSizeOne();

  //   assignMcEventsToClusters();

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
