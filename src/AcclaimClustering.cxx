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
#include "TArrowAntarctica.h"
#include "TF1.h"
#include "TCanvas.h"
#include <limits.h>
#include "RampdemReader.h"
#include "Math/Factory.h"

ClassImp(Acclaim::Clustering::Event);
ClassImp(Acclaim::Clustering::McEvent);
ClassImp(Acclaim::Clustering::Cluster);

const int nDim = 3;

const double FITTER_INPUT_SCALING = 1e-4;
const double FITTER_OUTPUT_SCALING = 1./FITTER_INPUT_SCALING;

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

/** 
 * Since we're not saving the useful pat, 
 * this might need to be explicitly called 
 * when reading from a tree in the future
 */
void Acclaim::Clustering::Event::setupUsefulPat(){
  Adu5Pat pat = anita.pat();
  usefulPat = UsefulAdu5Pat(&pat);
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
  RampdemReader::LonLatToEastingNorthing(longitude, latitude, easting, northing);
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

  nearestNeighbourEventNumber = 0;
  nearestNeighbourLogLikelihood = DBL_MAX;

  setupUsefulPat();
}


Acclaim::Clustering::Event::Event(){
  latitude = 0;
  longitude = 0;
  altitude = 0;
  easting = DBL_MAX;
  northing = DBL_MAX;
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



/** 
 * Make a TArrowAntarctica that points from ANITA to the event's location on the continent
 * 
 * @return the TArrowAntarctica
 */
TArrowAntarctica* Acclaim::Clustering::Event::makeArrowFromAnitaToEvent(){
  return new TArrowAntarctica(anita.longitude, anita.latitude, longitude, latitude);
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

Acclaim::Clustering::LogLikelihoodMethod::LogLikelihoodMethod() :
  fMinimizer(NULL),
  fFunctor(this, &Acclaim::Clustering::LogLikelihoodMethod::evalPairLogLikelihoodAtLonLat, 2)
{
  grTestMinimizerWalk = NULL;
  llCut = 40;
  // fTestEvent1 = 55850024;
  // fTestEvent2 = 60424781;
  fTestEvent1 = 61259774;
  fTestEvent2 = 55730926;

  fFitHorizonDistM = 700e3;
  numCallsToRecursive = 0;
  doneBaseClusterAssignment = false;
  fKDTree = NULL;
  hClusters = new TH2DAntarctica("hClusters", "hClusters");
  fDebug = false;
  fUseBaseList = true;

  fMaxFitterAttempts = 1;
  fMinimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2");
  fMinimizer->SetMaxFunctionCalls(1e5); // for Minuit/Minuit2
  fMinimizer->SetTolerance(0.0001);
  fMinimizer->SetPrintLevel(0);
  fMinimizer->SetFunction(fFunctor);
  
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




void Acclaim::Clustering::LogLikelihoodMethod::getDeltaThetaDegDeltaPhiDegEventCluster(Int_t eventInd, Int_t clusterInd, Double_t& deltaThetaDeg, Double_t& deltaPhiDeg){

  const Cluster& cluster = clusters.at(clusterInd);
  const Event& event = events.at(eventInd);

  Double_t thetaWave, phiWave;
  event.usefulPat.getThetaAndPhiWave(cluster.longitude, cluster.latitude, cluster.altitude, thetaWave, phiWave);
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


Double_t Acclaim::Clustering::LogLikelihoodMethod::evalPairLogLikelihoodAtLonLat(const Double_t* params){

  // Double_t sourceLon = params[0];
  // Double_t sourceLat = params[1];
  Double_t sourceEasting = FITTER_OUTPUT_SCALING*params[0];
  Double_t sourceNorthing = FITTER_OUTPUT_SCALING*params[1];

  if(events.at(fFitEventInd1).eventNumber==fTestEvent1){
    if(events.at(fFitEventInd2).eventNumber==fTestEvent2){
      if(grTestMinimizerWalk){
	grTestMinimizerWalk->SetPoint(grTestMinimizerWalk->GetN(), sourceEasting, sourceNorthing);
      }
    }
  }

  Double_t sourceAlt = RampdemReader::SurfaceAboveGeoidEN(sourceEasting, sourceNorthing);
  Double_t sourceLon, sourceLat;
  RampdemReader::EastingNorthingToLonLat(sourceEasting, sourceNorthing, sourceLon, sourceLat);

  Double_t ll = 0;
  ll += dPoint(fFitEventInd1, sourceLon, sourceLat, sourceAlt, true);
  ll += dPoint(fFitEventInd2, sourceLon, sourceLat, sourceAlt, true);

  if(fMinimizer->PrintLevel() > 0){
    std::cout << sourceEasting << "\t" << sourceNorthing << "\t" << sourceLon << "\t" << sourceLat << "\t" << ll << std::endl;
  }
  return ll;
}



Double_t Acclaim::Clustering::LogLikelihoodMethod::dFit(Int_t eventInd1, Int_t eventInd2){
  
  fFitEventInd1 = eventInd1;
  fFitEventInd2 = eventInd2;

  Double_t ll12 = dAsym(eventInd1, eventInd2);
  Double_t ll21 = dAsym(eventInd2, eventInd1);

  // we actually want to start at the WORSE of the two sources, apparently
  // which is very unintuitive
  Int_t which = ll21 < ll12 ? eventInd1 : eventInd2;

  if(fDebug && events.at(eventInd1).eventNumber==fTestEvent1 && events.at(eventInd2).eventNumber==fTestEvent2){
    grTestMinimizerWalk = new TGraph();
  }

  fMinimizer->SetVariable(0, "sourceEasting", FITTER_INPUT_SCALING*events.at(which).easting, 1);
  fMinimizer->SetVariable(1, "sourceNorthing", FITTER_INPUT_SCALING*events.at(which).northing, 1);

  int old_level = gErrorIgnoreLevel; 
  gErrorIgnoreLevel = 1001;
  bool validMinimum = false;
  std::vector<double> minima;
  for(int attempt=0; attempt < fMaxFitterAttempts && !validMinimum; attempt++){
    validMinimum = fMinimizer->Minimize();
    minima.push_back(fMinimizer->MinValue());
    if(!validMinimum){
      fMinimizer->SetVariable(0, "sourceEasting", fMinimizer->X()[0], 1);
      fMinimizer->SetVariable(1, "sourceNorthing", fMinimizer->X()[1], 1);
    }
  }
  gErrorIgnoreLevel = old_level;

  if(grTestMinimizerWalk){
    grTestMinimizerWalk->SetName("grTestMinimizerWalk");
    grTestMinimizerWalk->Write();
    delete grTestMinimizerWalk;
    grTestMinimizerWalk = NULL;
  }


  if((!validMinimum && fDebug)){ // do the same thing again, this time with error messages!
    int old_print_level = fMinimizer->PrintLevel();
    fMinimizer->SetPrintLevel(3);
    fMinimizer->SetVariable(0, "sourceEasting", FITTER_INPUT_SCALING*events.at(which).easting, 1);
    fMinimizer->SetVariable(1, "sourceNorthing", FITTER_INPUT_SCALING*events.at(which).northing, 1);
    for(int attempt=0; attempt < fMaxFitterAttempts && !validMinimum; attempt++){
      validMinimum = fMinimizer->Minimize();
      if(!validMinimum){
	fMinimizer->SetVariable(0, "sourceEasting", fMinimizer->X()[0], 1);
	fMinimizer->SetVariable(1, "sourceNorthing", fMinimizer->X()[1], 1);
      }
    }
    fMinimizer->SetPrintLevel(old_print_level);

    AnitaVersion::set(3);
    double waisLon = AnitaLocations::getWaisLongitude();
    double waisLat = AnitaLocations::getWaisLatitude();
    double waisModelAlt = RampdemReader::SurfaceAboveGeoid(waisLon, waisLat);

    Double_t llWais1 = dPoint(fFitEventInd1, waisLon, waisLat, waisModelAlt);
    Double_t llWais2 = dPoint(fFitEventInd2, waisLon, waisLat, waisModelAlt);
    Double_t llWais1b = dPoint(fFitEventInd1, waisLon, waisLat, waisModelAlt, true);
    Double_t llWais2b = dPoint(fFitEventInd2, waisLon, waisLat, waisModelAlt, true);
    std::cerr << "Debug in " << __PRETTY_FUNCTION__ << ": The fitter minimum is " << fMinimizer->MinValue() << std::endl;
    std::cerr << "event 1: run " << events.at(eventInd1).run << ", eventNumber " << events.at(eventInd1).eventNumber << std::endl;
    std::cerr << "event 2: run " << events.at(eventInd2).run << ", eventNumber " << events.at(eventInd2).eventNumber << std::endl;
    std::cerr << "Just for fun trying true WAIS location..." << std::endl;
    std::cerr << "llWais1 = " << llWais1 << "\t" << llWais1b << std::endl;
    std::cerr << "llWais2 = " << llWais2 << "\t" << llWais2b << std::endl;
    std::cerr << "Seed was at " << events.at(which).easting << "\t" << events.at(which).northing << "\t"
	      << events.at(which).longitude << "\t" << events.at(which).latitude << "\t" << TMath::Max(ll21, ll12) << std::endl;

    const int e1 = eventInd1;
    std::cerr << "event 1 was at " << events.at(e1).easting << "\t" << events.at(e1).northing << "\t" << events.at(e1).longitude << "\t" << events.at(e1).latitude << "\t" << ll12 << std::endl;
    const int e2 = eventInd2;
    std::cerr << "event 2 was at " << events.at(e2).easting << "\t" << events.at(e2).northing << "\t" << events.at(e2).longitude << "\t" << events.at(e2).latitude << "\t" << ll21 << std::endl;        
    std::cerr << std::endl;
  }
  Double_t ll = fMinimizer->MinValue();

  return ll;

}


Double_t Acclaim::Clustering::LogLikelihoodMethod::dMin(Int_t eventInd1, Int_t eventInd2){
  return TMath::Min(dAsym(eventInd1, eventInd2), dAsym(eventInd2, eventInd1));
}


Double_t Acclaim::Clustering::LogLikelihoodMethod::dSum(Int_t eventInd1, Int_t eventInd2){
  return dAsym(eventInd1, eventInd2) + dAsym(eventInd2, eventInd1);
}

Double_t Acclaim::Clustering::LogLikelihoodMethod::dAsym(Int_t eventInd1, Int_t eventInd2){

  const Event& event2 = events.at(eventInd2);
  return dPoint(eventInd1, event2.longitude, event2.latitude, event2.altitude);
}




/** 
 * Evaluate the log-likelihood distance from any event to an arbitrary point
 * 
 * @param eventInd1 is the index of the event of interest
 * @param sourceLon is the source longitude
 * @param sourceLat is the source latitude
 * @param sourceAlt is the source altitude
 * @param bool is an option to add a penalty term in LL for (approximately) being over the horizon
 * 
 * @return the log-likelihood
 */
Double_t Acclaim::Clustering::LogLikelihoodMethod::dPoint(Int_t eventInd, Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt, bool addOverHorizonPenalty){

  const Event& event = events.at(eventInd);

  Double_t theta, phi;
  event.usefulPat.getThetaAndPhiWave(sourceLon, sourceLat, sourceAlt, theta, phi);
  theta = -1*TMath::RadToDeg()*theta;
  phi = TMath::RadToDeg()*phi;

  Double_t dTheta = (theta - event.theta)/event.sigmaTheta;
  Double_t dPhi = Acclaim::RootTools::getDeltaAngleDeg(phi, event.phi)/event.sigmaPhi;
  
  Double_t ll = dTheta*dTheta + dPhi*dPhi;

  if(addOverHorizonPenalty){
    Double_t distM = event.usefulPat.getDistanceFromSource(sourceLat, sourceLon, sourceAlt);    
    if(distM >= fFitHorizonDistM){
      double distOverHorizonM = fabs(distM - fFitHorizonDistM);
      ll += distOverHorizonM;
    }
  }

  return ll;
}





void Acclaim::Clustering::LogLikelihoodMethod::sortEventIndices(Int_t eventInd, std::vector<Int_t>& unsortedIndices){

  std::vector<std::pair<Double_t, Int_t> > sortable;
  sortable.reserve(unsortedIndices.size());

  for(UInt_t i=0; i < unsortedIndices.size(); i++){
    Int_t thisEventInd = unsortedIndices[i];
    sortable.push_back(std::make_pair(dSum(eventInd, thisEventInd), thisEventInd));
  }
  std::sort(sortable.begin(), sortable.end());

  for(UInt_t i=0; i < unsortedIndices.size(); i++){
    unsortedIndices[i] = sortable[i].second;
  }
}



Double_t Acclaim::Clustering::LogLikelihoodMethod::getAngDistSqEventCluster(Int_t eventInd, Int_t clusterInd){

  const Cluster& cluster = clusters.at(clusterInd);
  const Event& event = events.at(eventInd);

  Double_t deltaThetaDeg, deltaPhiDeg;
  getDeltaThetaDegDeltaPhiDegEventCluster(eventInd, clusterInd, deltaThetaDeg, deltaPhiDeg);

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
      std::cout << "anita lon/lat/alt = " << event.usefulPat.longitude << ", " << event.usefulPat.latitude << ", " << event.usefulPat.altitude << std::endl;
    }
  }


  return angSq;
}


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
    for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
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
	  Double_t dist = events.at(eventInd).usefulPat.getDistanceFromSource(cluster.latitude,
									      cluster.longitude,
									      cluster.altitude);

	  // this is where proper collision detection would be useful...
	  if(dist < fFitHorizonDistM){
	    Double_t ll = getAngDistSqEventCluster(eventInd, clusterInd);

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

  Cluster& cluster = clusters.at(clusterInd);

  Double_t distM = event.usefulPat.getDistanceFromSource(cluster.latitude, cluster.longitude, cluster.altitude);

  if(distM < maxDistCluster){ // are we even close?

    Double_t ll = getAngDistSqEventCluster(i, clusterInd);

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

  UInt_t nMc = mcEvents.size();

  std::cerr  << "Info in " << __PRETTY_FUNCTION__ << ": starting!" << std::endl;
  ProgressBar p(nMc);
  for(UInt_t j=0; j < nMc; j++){
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
  for(UInt_t clusterInd=0; clusterInd < clusters.size(); clusterInd++){
    for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
      assignSingleEventToCloserCluster(eventInd, isMC, clusterInd);
    }

    clusters.at(clusterInd).antarcticaHistBin = hClusters->Fill(clusters.at(clusterInd).longitude, clusters.at(clusterInd).latitude);
    p.inc(clusterInd, clusters.size());
  }


  hBaseClusteredEvents.resize(clusters.size(), NULL);
  grBaseClusterCenters.resize(clusters.size(), NULL);

  for(UInt_t clusterInd=0; clusterInd < clusters.size(); clusterInd++){
    if(clusters.at(clusterInd).numDataEvents > 0){
      TString name = TString::Format("hBaseClusteredEvents%u", clusterInd);
      TH2DAntarctica* h = new TH2DAntarctica(name, name);
      hBaseClusteredEvents[clusterInd] = h;

      for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
	if(events.at(eventInd).cluster == static_cast<Int_t>(clusterInd)){
	  h->Fill(events.at(eventInd).longitude, events.at(eventInd).latitude);
	}
      }

      grBaseClusterCenters.at(clusterInd) = new TGraphAntarctica(BaseList::getBase(clusterInd));
      grBaseClusterCenters.at(clusterInd)->SetName(TString::Format("grClusterCenter%u", clusterInd));
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

  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
    events.at(eventInd).logLikelihood = DBL_MAX;
    events.at(eventInd).logLikelihood2 = DBL_MAX;
    events.at(eventInd).cluster = -1;
    events.at(eventInd).cluster2 = -1;
  }

  for(UInt_t eventInd=0; eventInd < mcEvents.size(); eventInd++){
    mcEvents.at(eventInd).logLikelihood = DBL_MAX;
    mcEvents.at(eventInd).logLikelihood2 = DBL_MAX;
    mcEvents.at(eventInd).cluster = -1;
    mcEvents.at(eventInd).cluster2 = -1;
  }

  const UInt_t numBases = BaseList::getNumBases();
  while(clusters.size() >= numBases && numBases>0){
    clusters.pop_back();
  }

  for(UInt_t clusterInd=0; clusterInd < clusters.size(); clusterInd++){
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

    for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
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
      getDeltaThetaDegDeltaPhiDegEventCluster(eventInd, event->cluster, event->dThetaCluster, event->dPhiCluster);

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


/** 
 * Don't do this for large data sets!
 */
void Acclaim::Clustering::LogLikelihoodMethod::testTriangleInequality(){

  const int numENNeighbours = 100;
  std::vector<Int_t> indices;
  std::vector<Int_t> enSeparation;

  ProgressBar p(events.size());
  for(UInt_t i=0; i < events.size(); i++){

    for(UInt_t j=0; j < numENNeighbours; j++){

      Double_t d_ij = dMin(i,j);

      for(UInt_t k=0; k < numENNeighbours; k++){

	Double_t d_ik = dMin(i, k);
	Double_t d_jk = dMin(j, k);

	if(d_ik > d_ij + d_jk){
	  std::cerr << i << "\t" << j << "\t" << k << " fails!" << std::endl;
	}
      }
    }
    p.inc(i, events.size());
  }
}




void Acclaim::Clustering::LogLikelihoodMethod::initKDTree(){
  if(fDebug){
    std::cout << "About to build KDTree" << std::endl;
  }
  std::vector<double> eastings;
  std::vector<double> northings;
  for(UInt_t eventInd = 0; eventInd < events.size(); eventInd++){
    const Event& event = events.at(eventInd);
    eastings.push_back(event.easting);
    northings.push_back(event.northing);
  }

  const int binSize = 100000; // meters... too small?
  fKDTree = new TKDTreeID(events.size(), 2, binSize);
  fKDTree->SetData(0, &eastings[0]);
  fKDTree->SetData(1, &northings[0]);
  fKDTree->Build();

  if(fDebug){
    std::cout << "Built!" << std::endl;
  
    const int nTest = 10;
    std::vector<int> nns(nTest);
    std::vector<double> dists(nTest);
    fKDTree->FindNearestNeighbors(&events.at(0).easting, nTest, &nns[0], &dists[0]);
    std::cout << "The ten nearest neighbours of events[0] at " << events[0].longitude << "," << events[0].latitude << " are:" << std::endl;
    for(int i=0; i < nTest; i++){
      std::cout << "events[" << nns[i] << "] at " << events[nns[i]].longitude << ", " << events[nns[i]].longitude << std::endl;
    }


    std::vector<Int_t> neighbours;
    const double rangeEN = 700e3;
    fKDTree->FindInRange(&events.at(0).easting, rangeEN, neighbours);
    std::cout << "There are " << neighbours.size() << " events within a sphere of radius " << rangeEN << std::endl;
  }
}



void Acclaim::Clustering::LogLikelihoodMethod::rangeQueryEastingNorthing(Int_t eventInd, Int_t numNearbyEN, std::vector<Int_t>& seed){
  Event& event = events.at(eventInd);
  std::vector<int> neighbourEventInds(numNearbyEN);
  std::vector<double> distsEN(numNearbyEN);
  fKDTree->FindNearestNeighbors(&events.at(0).easting, numNearbyEN, &neighbourEventInds[0], &distsEN[0]);
  for(int i=0; i < numNearbyEN; i++){
    // Double_t ll = d(eventInd, neighbourEventInds[i]);
    Double_t ll = dFit(eventInd, neighbourEventInds[i]);

    if(ll < event.nearestNeighbourLogLikelihood){
      event.nearestNeighbourEventNumber = events.at(neighbourEventInds[i]).eventNumber;
      event.nearestNeighbourLogLikelihood = ll;
    }

    if(ll < llCut){
      seed.push_back(neighbourEventInds[i]);
    }
  }
}


void Acclaim::Clustering::LogLikelihoodMethod::makeAndWriteNSquaredEventEventHistograms(){

  ProgressBar p(events.size());
  UInt_t firstEventNumber = events.at(0).eventNumber;
  UInt_t lastEventNumber = events.at(0).eventNumber;
  for(UInt_t eventInd=1; eventInd < events.size(); eventInd++){
    UInt_t eventNumber = events.at(eventInd).eventNumber;
    if(eventNumber > lastEventNumber){
      lastEventNumber = eventNumber;
    }
    if(eventNumber < firstEventNumber){
      firstEventNumber = eventNumber;
    }
  }


  TH1D* hWais = new TH1D("hWais", "Log-likelihood - WAIS", 2048, 0, 2048);
  AnitaVersion::set(3);
  double waisLon = AnitaLocations::getWaisLongitude();
  double waisLat = AnitaLocations::getWaisLatitude();
  double waisModelAlt = RampdemReader::SurfaceAboveGeoid(waisLon, waisLat);

  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
    Double_t llWais = dPoint(eventInd, waisLon, waisLat, waisModelAlt);
    hWais->Fill(llWais);
  }

  TRandom3 randy(123);
  const Int_t nBins = 1024; //(lastEventNumber + 1) - firstEventNumber;
  // const Int_t stride = 50;
  std::cout << nBins << "\t" << firstEventNumber << "\t" << lastEventNumber << std::endl;
  TH2D* hUnfit = new TH2D("hUnfit_d_vs_pointingAngle", "", 1024, 0, 180, 1024, 0, 1024);
  TH2D* hFit = new TH2D("hFit_d_vs_pointingAngle", "", 1024, 0, 180, 1024, 0, 1024);
  TH2D* hFitVsUnfitted = new TH2D("hFit_vs_unfit", "Fitted vs. unfitted WAIS event-event log-likelihood; Unfitted log-likelihood; Fitted log-likelihood", 1024, 0, 1024, 1024, 0, 1024);
  TH1D* hUnfitMinusFit = new TH1D("hFit_minus_unfit", "(Unfitted - fitted) WAIS event-event log-likelihood; #delta (log likelihood)", 1024, -512, 512);
  
  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
    const Event& event1 = events.at(eventInd);

    // if(event1.eventNumber!=fTestEvent1) continue;

    TVector3 event1Pos = AntarcticCoord(AntarcticCoord::WGS84, event1.latitude, event1.longitude, event1.altitude).v();
    TVector3 anita1Pos = AntarcticCoord(AntarcticCoord::WGS84, event1.anita.latitude, event1.anita.longitude, event1.anita.altitude).v();
    TVector3 anitaToEvent1 = event1Pos - anita1Pos;
    
    // for(UInt_t eventInd2=eventInd; eventInd2 < events.size(); eventInd2+=stride){
    // for(UInt_t eventInd2=eventInd; eventInd2 < events.size(); eventInd2+=stride){
    const int eventInd2 = randy.Uniform(0, events.size());
    const Event& event2 = events.at(eventInd2);

    // if(event2.eventNumber!=fTestEvent2) continue;

    if(fDebug && event1.eventNumber==fTestEvent1 && event2.eventNumber==fTestEvent2){
      std::cerr << "Info in " << __PRETTY_FUNCTION__ << " mapping parameter space around WAIS for a single event test!" << std::endl;
      AnitaVersion::set(3);
      double waisEasting, waisNorthing;
      RampdemReader::LonLatToEastingNorthing(AnitaLocations::getWaisLongitude(), AnitaLocations::getWaisLatitude(), waisEasting, waisNorthing);
      double delta = 100e3;
      const int nBins = 1024;
      TH2D* hParams = new TH2D("hSingleEventTest", "param vs. e/n",
			       // nBins, -1090000, -1065000,
			       // nBins, -484200, -483600);
			       nBins, waisEasting-delta, waisEasting+delta,
			       nBins, waisNorthing-delta, waisNorthing+delta);

      fFitEventInd1 = eventInd;
      fFitEventInd2 = eventInd2;
      double params[2] = {0,0};
      for(int by=1; by <= hParams->GetNbinsY(); by++){
	params[1] = FITTER_INPUT_SCALING*hParams->GetYaxis()->GetBinCenter(by);
	for(int bx=1; bx <= hParams->GetNbinsX(); bx++){
	  params[0] = FITTER_INPUT_SCALING*hParams->GetXaxis()->GetBinCenter(bx);
	  double ll = evalPairLogLikelihoodAtLonLat(params);
	  hParams->Fill(FITTER_OUTPUT_SCALING*params[0], FITTER_OUTPUT_SCALING*params[1], ll);
	}
      }
      std::cout << "done!" << std::endl;
      hParams->Write();
      delete hParams;
    }

    TVector3 event2Pos = AntarcticCoord(AntarcticCoord::WGS84, event2.latitude, event2.longitude, event2.altitude).v();
    TVector3 anita2Pos = AntarcticCoord(AntarcticCoord::WGS84, event2.anita.latitude, event2.anita.longitude, event2.anita.altitude).v();
    TVector3 anitaToEvent2 = event2Pos - anita2Pos;

    double dist = dSum(eventInd, eventInd2);
    Double_t angleBetweenEvents = anitaToEvent1.Angle(anitaToEvent2);
      
    hUnfit->Fill(angleBetweenEvents*TMath::RadToDeg(), dist);

    double distFitted = dFit(eventInd, eventInd2);
    hFit->Fill(angleBetweenEvents*TMath::RadToDeg(), distFitted);

    hFitVsUnfitted->Fill(dist, distFitted);
    hUnfitMinusFit->Fill(dist - distFitted);

    p.inc(eventInd, events.size());
  }

  hWais->Write();
  delete hWais;
  
  hUnfit->Write();
  delete hUnfit;

  hFit->Write();
  delete hFit;

  hFitVsUnfitted->Write();
  delete hFitVsUnfitted;

  hUnfitMinusFit->Write();
  delete hUnfitMinusFit;
  
}


void Acclaim::Clustering::LogLikelihoodMethod::OPTICS(){

  // const int kUndefined = -2;
  // const int kNoise = -1;

  // coreDist = undefined if numNeighbours < fMinPts
  //          = distance to the fMinPts-th closest point


  // reachabilityDist defined from point point a to point b
  // reachabilityDist = undefined if numNeighbours < fMinPts
  //                  = max(coreDist, )
}

void Acclaim::Clustering::LogLikelihoodMethod::DBSCAN(){


  const int kUndefined = -2;
  const int kNoise = -1;
  
  const int fMinPts = 4;
  const int numNearbyEN = 2*fMinPts;
  int clusterInd = -1;

  // mark all clusters as undefined
  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
    events.at(eventInd).cluster = kUndefined;
  }

  
  std::cout << "Info in " << __PRETTY_FUNCTION__ << ": starting!" << std::endl;
  ProgressBar p(events.size());
  UInt_t numEventsProcessed = 0;
  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
    Event& event = events.at(eventInd);

    if(event.cluster == kUndefined){ // then not yet processed
    
      std::vector<Int_t> seed;
      rangeQueryEastingNorthing(eventInd, numNearbyEN, seed);
      p.inc(numEventsProcessed, events.size());
      numEventsProcessed++;
      
      if(seed.size() < fMinPts){
	event.cluster = kNoise;
      }
      else{
	Cluster cluster;
	cluster.latitude = event.latitude;
	cluster.longitude = event.longitude;
	cluster.altitude = event.altitude;
	clusters.push_back(cluster);
	clusterInd++;
      
	for(UInt_t seedInd=0; seedInd < seed.size(); seedInd++){
	  Int_t eventInd2 = seed.at(seedInd);
	  Event& event2 = events.at(eventInd2);

	  if(event2.cluster == kNoise){
	    event2.cluster = clusterInd;
	    clusters.at(clusterInd).numDataEvents++;
	  }
	  else if(event2.cluster == kUndefined){
	    event2.cluster = clusterInd;
	    clusters.at(clusterInd).numDataEvents++;

	    rangeQueryEastingNorthing(eventInd2, numNearbyEN, seed);
	    p.inc(numEventsProcessed, events.size());
	    numEventsProcessed++;
	  }
	}
      }
    }
  }
}



void Acclaim::Clustering::LogLikelihoodMethod::doClusteringDBSCAN(const char* dataGlob, const char* mcGlob, const char* outFileName){

  readInSummaries(dataGlob);
  initKDTree();
  setDebug(false);
  // return;
  // DBSCAN();
  const char* fakeArgv[2] = {outFileName, "DBSCAN"};
  OutputConvention oc(1, const_cast<char**>(fakeArgv));
  TFile* fOut = oc.makeFile();

  makeAndWriteNSquaredEventEventHistograms();
  // makeSummaryTrees();
  fOut->Write();
  fOut->Close();
  
  return;  
  
  
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
