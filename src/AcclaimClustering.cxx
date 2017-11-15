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

#include "AcclaimOpenMP.h"

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

  antarcticaHistBin = -1;

  resetClusteringNumbers();

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

  for(int dim=0; dim < nDim; dim++){
    centre[dim] = 0;
  }
  antarcticaHistBin = -1;
  resetClusteringNumbers();
}


void Acclaim::Clustering::Event::resetClusteringNumbers(){

  for(int z=0; z<2; z++){
    logLikelihood[z] = DBL_MAX;
    cluster[z] = -1;
    logLikelihood2[z] = DBL_MAX;
    cluster2[z] = -1;
  }
  nearestNeighbourEventNumber = 0;
  nearestNeighbourLogLikelihood = DBL_MAX;
  dThetaCluster = -999;
  dPhiCluster = -999;
  eventEventClustering = true;
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



Acclaim::Clustering::Cluster::Cluster(Int_t i) {
  for(int dim=0; dim < nDim; dim++){
    centre[dim] = 0;
  }
  latitude = 0;
  longitude = 0;
  altitude = 0;
  knownBase = 0;
  resetClusteringNumbers();
  antarcticaHistBin = -1;
  seedEvent = -1;
  index = i;
}

Acclaim::Clustering::Cluster::Cluster(const BaseList::base& base, Int_t i) {
  // std::cerr << __PRETTY_FUNCTION__ << std::endl;
  AntarcticCoord ac = base.position.as(AntarcticCoord::WGS84);
  latitude = ac.x;
  longitude = ac.y;
  altitude = ac.z;
  knownBase = 1;

  if(altitude == -999){
    altitude = RampdemReader::SurfaceAboveGeoid(longitude, latitude);
  }

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->getCartesianCoords(latitude, longitude, altitude, centre);
  resetClusteringNumbers();
  antarcticaHistBin = -1;
  seedEvent = -1;
  index = i;
}


Acclaim::Clustering::Cluster::Cluster(const Event& event, Int_t i) {
  // std::cerr << __PRETTY_FUNCTION__ << std::endl;
  latitude = event.longitude;
  longitude = event.latitude;
  altitude = event.altitude;
  knownBase = 0;

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->getCartesianCoords(latitude, longitude, altitude, centre);
  resetClusteringNumbers();
  index = i;
}


void Acclaim::Clustering::Cluster::resetClusteringNumbers(){
  numDataEvents = 0;
  sumMcWeights = 0;
}



Acclaim::Clustering::LogLikelihoodMethod::LogLikelihoodMethod() :
  fFitEvent1(NULL),
  fFitEvent2(NULL),
  fMinimizer(NULL),
  fFunctor(this, &Acclaim::Clustering::LogLikelihoodMethod::evalPairLogLikelihoodAtLonLat, 2)
{
  grTestMinimizerWalk = NULL;
  llClusterCut = 100;

  llEventCut = 20;
  llEventCuts.push_back(20);
  llEventCuts.push_back(40);
  fTestEvent1 = 58900975;
  fTestEvent2 = 55699027;

  fFitHorizonDistM = 700e3;
  numCallsToRecursive = 0;
  doneBaseClusterAssignment = false;
  fKDTree = NULL;
  hClusters = new TH2DAntarctica("hClusters", "hClusters");
  hEvents = NULL;
  // hEvents = makeUnevenlyBinnedEventHistogram();
  hEvents = new TH2DAntarctica("hEvents", "hEvents");

  fDebug = false;
  fUseBaseList = true;

  fMaxFitterAttempts = 1;
  fMinimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2");
  fMinimizer->SetMaxFunctionCalls(1e5); // for Minuit/Minuit2
  fMinimizer->SetTolerance(0.0001);
  fMinimizer->SetPrintLevel(0);
  fMinimizer->SetFunction(fFunctor);
}




TH2DAntarctica* Acclaim::Clustering::LogLikelihoodMethod::makeUnevenlyBinnedEventHistogram(){

  if(!hClusters){
    hClusters = new TH2DAntarctica("hClusters", "hClusters");
  }

  std::vector<Double_t> x(hClusters->GetNbinsX()+1);
  for(int i=0; i <= hClusters->GetNbinsX(); i++){
    x.at(i) = hClusters->GetXaxis()->GetBinUpEdge(i);
  }
  std::vector<Double_t> y(hClusters->GetNbinsX()+1);
  for(int i=0; i <= hClusters->GetNbinsY(); i++){
    y.at(i) = hClusters->GetYaxis()->GetBinUpEdge(i);
  }

  // now uneven-ize the bins
  double easting, northing;
  RampdemReader::LonLatToEastingNorthing(167, -78, easting, northing);
  const double fineBinWidthEasting = 5e2;
  const double fineBinWidthNorthing = 5e2;
  const double fineRangeEasting = 700e3;
  const double fineRangeNorthing = 700e3;
  for(double binEasting = easting-fineRangeEasting; binEasting <= easting + fineRangeEasting; binEasting += fineBinWidthEasting){
    x.push_back(binEasting);
  }
  for(double binNorthing = northing-fineRangeNorthing; binNorthing <= northing+fineRangeNorthing; binNorthing += fineBinWidthNorthing){
    y.push_back(binNorthing);
  }

  std::sort(x.begin(), x.end());
  std::sort(y.begin(), y.end());

  if(hEvents){
    delete hEvents;
  }
  hEvents = new TH2DAntarctica("hEvents", "hEvents", x, y);
  return hEvents;
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



Double_t Acclaim::Clustering::LogLikelihoodMethod::getDistSqEventCluster(const Event& event, const Acclaim::Clustering::Cluster& cluster){
  Double_t d2=0;
  for(int dim=0; dim < nDim; dim++){
    Double_t d = cluster.centre[dim] - event.centre[dim];
    d2 += d*d;
  }
  return d2;
}



void Acclaim::Clustering::LogLikelihoodMethod::getDeltaThetaDegDeltaPhiDegEventCluster(const Event& event, const Cluster& cluster, Double_t& deltaThetaDeg, Double_t& deltaPhiDeg){

  // const Cluster& cluster = clusters.at(clusterInd);
  // const Event& event = events.at(eventInd);

  Double_t thetaWave, phiWave;
  event.usefulPat.getThetaAndPhiWave(cluster.longitude, cluster.latitude, cluster.altitude, thetaWave, phiWave);
  Double_t thetaDeg = -1*TMath::RadToDeg()*thetaWave;
  Double_t phiDeg = TMath::RadToDeg()*phiWave;

  deltaThetaDeg = (thetaDeg - event.theta);
  deltaPhiDeg = Acclaim::RootTools::getDeltaAngleDeg(phiDeg, event.phi);

  // if(fDebug){
  //   if(event.cluster < 0 && eventInd == cluster.seedEvent){
  //     std::cout << "event theta = " << event.theta  << ", cluster theta = " << thetaDeg << std::endl;
  //     std::cout << "event phi = " << event.phi  << ", cluster phi = " << phiDeg << std::endl;
  //   }
  // }
}


Double_t Acclaim::Clustering::LogLikelihoodMethod::evalPairLogLikelihoodAtLonLat(const Double_t* params){

  Double_t sourceEasting = FITTER_OUTPUT_SCALING*params[0];
  Double_t sourceNorthing = FITTER_OUTPUT_SCALING*params[1];
  if(fFitEvent1->eventNumber==fTestEvent1){
    if(fFitEvent2->eventNumber==fTestEvent2){
      if(grTestMinimizerWalk){
	grTestMinimizerWalk->SetPoint(grTestMinimizerWalk->GetN(), sourceEasting, sourceNorthing);
      }
    }
  }

  Double_t sourceAlt = RampdemReader::SurfaceAboveGeoidEN(sourceEasting, sourceNorthing);
  Double_t sourceLon, sourceLat;
  RampdemReader::EastingNorthingToLonLat(sourceEasting, sourceNorthing, sourceLon, sourceLat);

  Double_t ll = 0;
  ll += dPoint(*fFitEvent1, sourceLon, sourceLat, sourceAlt, true);
  ll += dPoint(*fFitEvent2, sourceLon, sourceLat, sourceAlt, true);

  if(fMinimizer->PrintLevel() > 0){
    std::cout << sourceEasting << "\t" << sourceNorthing << "\t" << sourceLon << "\t" << sourceLat << "\t" << ll << std::endl;
  }
  return ll;
}



Double_t Acclaim::Clustering::LogLikelihoodMethod::dFit(const Event& event1, const Event& event2){
  fFitEvent1 = &event1;
  fFitEvent2 = &event2;

  Double_t ll12 = dAsym(event1, event2);
  Double_t ll21 = dAsym(event2, event1);

  // we actually want to start at the WORSE of the two sources, apparently
  // which is very unintuitive
  const Event* which = ll21 < ll12 ? fFitEvent1 : fFitEvent2;

  if(fFitEvent1->eventNumber==fTestEvent1 && fFitEvent2->eventNumber==fTestEvent2){
    grTestMinimizerWalk = new TGraph();
  }

  fMinimizer->SetVariable(0, "sourceEasting", FITTER_INPUT_SCALING*which->easting, 1);
  fMinimizer->SetVariable(1, "sourceNorthing", FITTER_INPUT_SCALING*which->northing, 1);

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
    fMinimizer->SetVariable(0, "sourceEasting", FITTER_INPUT_SCALING*which->easting, 1);
    fMinimizer->SetVariable(1, "sourceNorthing", FITTER_INPUT_SCALING*which->northing, 1);
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
    const Event& event1 = *fFitEvent1;
    const Event& event2 = *fFitEvent2;
    Double_t llWais1 = dPoint(event1, waisLon, waisLat, waisModelAlt);
    Double_t llWais2 = dPoint(event2, waisLon, waisLat, waisModelAlt);
    Double_t llWais1b = dPoint(event1, waisLon, waisLat, waisModelAlt, true);
    Double_t llWais2b = dPoint(event2, waisLon, waisLat, waisModelAlt, true);
    std::cerr << "Debug in " << __PRETTY_FUNCTION__ << ": The fitter minimum is " << fMinimizer->MinValue() << std::endl;
    std::cerr << "event 1: run " << event1.run << ", eventNumber " << event1.eventNumber << std::endl;
    std::cerr << "event 2: run " << event2.run << ", eventNumber " << event2.eventNumber << std::endl;
    std::cerr << "Just for fun trying true WAIS location..." << std::endl;
    std::cerr << "llWais1 = " << llWais1 << "\t" << llWais1b << std::endl;
    std::cerr << "llWais2 = " << llWais2 << "\t" << llWais2b << std::endl;
    std::cerr << "Seed was at " << which->easting << "\t" << which->northing << "\t"
	      << which->longitude << "\t" << which->latitude << "\t" << TMath::Max(ll21, ll12) << std::endl;

    std::cerr << "event 1 was at " << event1.easting << "\t" << event1.northing << "\t" << event1.longitude << "\t" << event1.latitude << "\t" << ll12 << std::endl;
    std::cerr << "event 2 was at " << event2.easting << "\t" << event2.northing << "\t" << event2.longitude << "\t" << event2.latitude << "\t" << ll21 << std::endl;
    std::cerr << std::endl;
  }

  fFitEasting  = fMinimizer->X()[0]*FITTER_OUTPUT_SCALING;
  fFitNorthing = fMinimizer->X()[1]*FITTER_OUTPUT_SCALING;
  Double_t ll = fMinimizer->MinValue();

  return ll;
}


Double_t Acclaim::Clustering::LogLikelihoodMethod::dMin(const Event& event1, const Event& event2){
  return TMath::Min(dAsym(event1, event2), dAsym(event2, event1));
}


Double_t Acclaim::Clustering::LogLikelihoodMethod::dSum(const Event& event1, const Event& event2){
  return dAsym(event1, event2) + dAsym(event2, event1);
}

Double_t Acclaim::Clustering::LogLikelihoodMethod::dAsym(const Event& event1, const Event& event2){
  return dPoint(event1, event2.longitude, event2.latitude, event2.altitude);
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
Double_t Acclaim::Clustering::LogLikelihoodMethod::dPoint(const Event& event, Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt, bool addOverHorizonPenalty){

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



Double_t Acclaim::Clustering::LogLikelihoodMethod::getAngDistSqEventCluster(const Event& event, const Cluster& cluster){
// Double_t Acclaim::Clustering::LogLikelihoodMethod::getAngDistSqEventCluster(Int_t eventInd, Int_t clusterInd){  

  Double_t deltaThetaDeg, deltaPhiDeg;
  getDeltaThetaDegDeltaPhiDegEventCluster(event, cluster, deltaThetaDeg, deltaPhiDeg);

  Double_t dThetaNorm = deltaThetaDeg/event.sigmaTheta;
  Double_t dPhiNorm = deltaPhiDeg/event.sigmaPhi;
  Double_t angSq =  dThetaNorm*dThetaNorm + dPhiNorm*dPhiNorm;
  // if(fDebug){
  //   // if(event.cluster < 0 && event.antarcticaHistBin == cluster.antarcticaHistBin){
  //   if(event.cluster < 0 && eventInd == cluster.seedEvent){
  //     std::cout << event.eventNumber << std::endl;
  //     std::cout << "dTheta = " << deltaThetaDeg  << ", dPhi = " << deltaPhiDeg
  // 		<< ", sigmaTheta = " << event.sigmaTheta << ", sigmaPhi = " << event.sigmaPhi
  // 		<< ", dThetaNorm = " << dThetaNorm << ", dPhiNorm = " << dPhiNorm
  // 		<< ", ll = " << angSq << std::endl;
  //     std::cout << "cluster lon/lat/alt = " << cluster.longitude << ", " << cluster.latitude << ", " << cluster.altitude << std::endl;
  //     std::cout << "event lon/lat/alt = " << event.longitude << ", " << event.latitude << ", " << event.altitude << std::endl;
  //     std::cout << "anita lon/lat/alt = " << event.usefulPat.longitude << ", " << event.usefulPat.latitude << ", " << event.usefulPat.altitude << std::endl;
  //   }
  // }


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
    if(event.cluster[0] < 0){
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
    if(event.cluster[0] < 0 && event.antarcticaHistBin==globalMaxBin){
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
      if(event.cluster[0] < 0 && event.antarcticaHistBin==globalMaxBin){
	Double_t dSq = getDistSqEventCluster(event, tempCluster);
	if(dSq < minDSq){
	  minDSq = dSq;
	  eventIndClosestToBinMean = eventInd;
	}
      }
    }

    Cluster seedCluster(clusters.size());
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
      assignSingleEventToCloserCluster(eventInd, isMC, clusters.back());
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
	if(events.at(eventInd).cluster[0]==static_cast<Int_t>(clusterInd)){
	  // some reset function would be good here...
	  events.at(eventInd).resetClusteringNumbers();
	  numEventsReset++;
	}
      }
    }
  }

  std::cout << numEventsReset << " events were reset over " << numClustersReset << " clusters!" << std::endl;
  // go through and empty non-base clusters?

  // now do the nSquared thing
  for(int eventInd=0; eventInd < (int)events.size(); eventInd++){
    if(events.at(eventInd).cluster[0] < 0){
      const Int_t isMC = 0;
      for(int clusterInd=0; clusterInd < (Int_t)clusters.size(); clusterInd++){
	Cluster& cluster = clusters.at(clusterInd);
	assignSingleEventToCloserCluster(eventInd, isMC, cluster);
      }
    }
  }
}





void Acclaim::Clustering::LogLikelihoodMethod::assignSingleEventToCloserCluster(Int_t i, Int_t isMC, Cluster& cluster){
  Event& event = isMC==0 ? events.at(i) : mcEvents.at(i);

  Double_t distM = event.usefulPat.getDistanceFromSource(cluster.latitude, cluster.longitude, cluster.altitude);

  if(distM < fFitHorizonDistM){ // are we even close?

    Double_t ll = getAngDistSqEventCluster(event, cluster);

    if(ll < llClusterCut){ // within the cut?

      // if(fDebug){
      // 	if(event.cluster < 0){
      // 	  std::cout << "clusterInd = " << clusterInd << ", i = " << i << ", ll = " << ll << ", event.ll = " << event.logLikelihood << ", event.cluster = " << event.cluster << std::endl;
      // 	}
      // }

      // if(!isMC){
      // 	std::cout << ll << std::endl;
      // }

      if(ll < event.logLikelihood[0]){ // the best event?

	event.cluster2[0] = event.cluster[0];
	event.logLikelihood2[0] = event.logLikelihood[0];

	event.logLikelihood[0] = ll;
	event.cluster[0] = cluster.index;

	if(isMC==0){
	  if(event.cluster2[0] >= 0){
	    clusters.at(event.cluster2[0]).numDataEvents--;
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
    for(UInt_t clusterInd=0; clusterInd < clusters.size(); clusterInd++){
      Cluster& cluster = clusters.at(clusterInd);
      assignSingleEventToCloserCluster(j, isMC, cluster);
    }
    if(mcEvents.at(j).cluster[0] < 0){
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
    Cluster& cluster = clusters.at(clusterInd);
    for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
      assignSingleEventToCloserCluster(eventInd, isMC, cluster);
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
	if(events.at(eventInd).cluster[0] == static_cast<Int_t>(clusterInd)){
	  h->Fill(events.at(eventInd).longitude, events.at(eventInd).latitude);
	}
      }

      grBaseClusterCenters.at(clusterInd) = new TGraphAntarctica(BaseList::getBase(clusterInd));
      grBaseClusterCenters.at(clusterInd)->SetName(TString::Format("grClusterCenter%u", clusterInd));
    }
  }

  doneBaseClusterAssignment = true;
}



void Acclaim::Clustering::LogLikelihoodMethod::initializeBaseList(){

  std::cout << "Info in " << __PRETTY_FUNCTION__ << ": Initializing base list..." << std::endl;
  for(UInt_t clusterInd=0; clusterInd < BaseList::getNumBases(); clusterInd++){
    const BaseList::base& base = BaseList::getBase(clusterInd);
    clusters.push_back(Cluster(base, clusters.size()));
  }
}



void Acclaim::Clustering::LogLikelihoodMethod::resetClusters(){

  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
    for(int z=0; z < 2; z++){
      events.at(eventInd).logLikelihood[z] = DBL_MAX;
      events.at(eventInd).logLikelihood2[z] = DBL_MAX;
      events.at(eventInd).cluster[z] = -1;
      events.at(eventInd).cluster2[z] = -1;
    }
  }

  for(UInt_t eventInd=0; eventInd < mcEvents.size(); eventInd++){
    for(int z=0; z < 2; z++){    
      mcEvents.at(eventInd).logLikelihood[z] = DBL_MAX;
      mcEvents.at(eventInd).logLikelihood2[z] = DBL_MAX;
      mcEvents.at(eventInd).cluster[z] = -1;
      mcEvents.at(eventInd).cluster2[z] = -1;
    }
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
      if(events.at(eventInd).cluster[0]==clusterInd){
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
    if(event->cluster[0] >= 0){
      const Cluster& cluster = clusters.at(event->cluster[0]);
      getDeltaThetaDegDeltaPhiDegEventCluster(*event, cluster, event->dThetaCluster, event->dPhiCluster);
    }
    if(event->cluster[0] < 0 || (clusters.at(event->cluster[0]).numDataEvents == 1 && clusters.at(event->cluster[0]).knownBase==0)) {
      nonBaseSingletTree->Fill();
    }
    else{
      clusteredDataTree->Fill();
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
      /// @todo maybe remove the run hack!!!
      if(inSandbox(sum->peak[pol][peakIndex]) && sourceLat > -999 && sourceLon > -999 && snrHack < 100 && TMath::Abs(thetaAdjustment) < 1e-10){
	// if(sum->run > 160 && sourceLat > -999 && sourceLon > -999 && snrHack < 100 && TMath::Abs(thetaAdjustment) < 1e-10){

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
    const Event& event_i = events.at(i);
    for(UInt_t j=0; j < numENNeighbours; j++){
      const Event& event_j = events.at(j);

      Double_t d_ij = dMin(event_i,event_j);

      for(UInt_t k=0; k < numENNeighbours; k++){
	const Event& event_k = events.at(k);
	Double_t d_ik = dMin(event_i, event_k);
	Double_t d_jk = dMin(event_j, event_k);

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

  fEventEastings.clear();
  fEventNorthings.clear();
  fEventEastings.reserve(events.size());
  fEventNorthings.reserve(events.size());
  for(UInt_t eventInd = 0; eventInd < events.size(); eventInd++){
    const Event& event = events.at(eventInd);
    if(event.eventEventClustering){
      fEventEastings.push_back(event.easting);
      fEventNorthings.push_back(event.northing);
    }
  }

  const int binSize = 100000; // meters... too small?
  fKDTree = new TKDTreeID(events.size(), 2, binSize);
  fKDTree->SetData(0, &fEventEastings[0]);
  fKDTree->SetData(1, &fEventNorthings[0]);
  fKDTree->Build();

  if(fDebug){
    std::cout << "Built!" << std::endl;

    const int nTest = 10;
    std::vector<int> nns(nTest);
    std::vector<double> dists(nTest);
    fKDTree->FindNearestNeighbors(&events.at(0).easting, nTest, &nns[0], &dists[0]);
    std::cout << "The ten nearest neighbours of events[0] at " << events[0].longitude << "," << events[0].latitude << " are:" << std::endl;
    for(int i=0; i < nTest; i++){
      std::cout << "events[" << nns[i] << "] at " << events[nns[i]].longitude << ", " << events[nns[i]].latitude << std::endl;
    }

    std::vector<Int_t> neighbours;
    const double rangeEN = 1000e3;
    double lookup[2] = {events.at(0).easting, events.at(0).northing};
    fKDTree->FindInRange(lookup, rangeEN, neighbours);
    std::cout << "There are " << neighbours.size() << " events within a sphere of radius " << rangeEN << std::endl;
  }
}


void Acclaim::Clustering::LogLikelihoodMethod::testAngleFindingSpeed(){

  std::cout << __PRETTY_FUNCTION__ << std::endl;
  TRandom3 randy;

  const double llEventCutSqrt = TMath::Sqrt(llEventCut);

  int grandTotal = 0;
  ProgressBar p(events.size());

  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){

    // if(eventInd < 200){

    Event& event = events.at(eventInd);
    const double deltaThetaMax = llEventCutSqrt*event.sigmaTheta; /// height of the error ellipse at phi=0
    const double deltaPhiMax = llEventCutSqrt*event.sigmaPhi;     /// height of the error ellipse at theta=0

    const int nSweepPoints = 10;

    double minEasting = DBL_MAX;
    double maxEasting = -DBL_MAX;
    double minNorthing = DBL_MAX;
    double maxNorthing = -DBL_MAX;

    double deltaPhi = -deltaPhiMax;
    double deltaPhiStep = 2*deltaPhiMax/(nSweepPoints-1);
    std::vector<Int_t> histBinsToConsider;

    for(int i=0; i < nSweepPoints; i++){
      double deltaPhiNorm = deltaPhi/event.sigmaPhi;
      double underSqrtSign = llEventCut - deltaPhiNorm*deltaPhiNorm;
      underSqrtSign = underSqrtSign < 0 ? 0 : underSqrtSign;
      double deltaTheta = event.sigmaTheta*TMath::Sqrt(underSqrtSign);

      const double thetaWaveHigh = -1*TMath::DegToRad()*(event.theta + deltaTheta);
      const double thetaWaveLow  = -1*TMath::DegToRad()*(event.theta - deltaTheta);
      const double phiWave = TMath::DegToRad()*(event.phi + deltaPhi);

      const int nQuadraticSolutions = 2;
      const double thetaWaves[nQuadraticSolutions] = {thetaWaveLow, thetaWaveHigh};

      for(int j=0; j < nQuadraticSolutions; j++){
	double lon, lat, alt, thetaAdjust;
	event.usefulPat.traceBackToContinent(phiWave, thetaWaves[j], &lon, &lat, &alt, &thetaAdjust, 0);

	if(lat < -60){
	  // std::cout << i << "\t" << bin << "\t" << deltaPhi << "\t" << deltaTheta << "\t" << lon << "\t" << lat << std::endl;
	  Double_t easting, northing;
	  RampdemReader::LonLatToEastingNorthing(lon, lat, easting, northing);
	  int bin = hEvents->FindBin(easting, northing);
	  if(std::find(histBinsToConsider.begin(), histBinsToConsider.end(), bin)==histBinsToConsider.end()){
	    histBinsToConsider.push_back(bin);
	  }

	  if(easting < minEasting){
	    minEasting = easting;
	  }
	  if(easting > maxEasting){
	    maxEasting = easting;
	  }
	  if(northing < minNorthing){
	    minNorthing = northing;
	  }
	  if(northing > maxNorthing){
	    maxNorthing = northing;
	  }
	}
      }

      deltaPhi += deltaPhiStep;
    }

    TAxis* xAxis = hEvents->GetXaxis();
    TAxis* yAxis = hEvents->GetYaxis();
    int bxLow = xAxis->FindBin(minEasting);
    int bxHigh = xAxis->FindBin(maxEasting);
    int byLow = yAxis->FindBin(minNorthing);
    int byHigh = yAxis->FindBin(maxNorthing);

    for(int by = byLow; by <= byHigh; by++){
      double northing = yAxis->GetBinCenter(by);
      for(int bx = bxLow; bx <= bxHigh; bx++){
	double easting = xAxis->GetBinCenter(bx);
	double lon, lat, alt;

	RampdemReader::EastingNorthingToLonLat(easting, northing, lon, lat);
	alt = RampdemReader::SurfaceAboveGeoid(lon, lat);

	Double_t theta, phi;
	event.usefulPat.getThetaAndPhiWave(lon, lat, alt, theta, phi);
	theta = -1*TMath::RadToDeg()*theta;
	phi = TMath::RadToDeg()*phi;

	double deltaTheta = (theta - event.theta)/event.sigmaTheta;
	double deltaPhi = RootTools::getDeltaAngleDeg(phi, event.phi)/event.sigmaPhi;
	double binLL = deltaPhi*deltaPhi + deltaTheta*deltaTheta;

	// std::cout << bx << "\t" << by << "\t" << event.easting - easting << "\t" << event.northing - northing << std::endl;
	// std::cout << theta << "\t" << phi << std::endl;
	// std::cout << deltaTheta*event.sigmaTheta << "\t" << deltaPhi*event.sigmaPhi << std::endl;
	// std::cout << event.sigmaTheta << "\t" << event.sigmaPhi << std::endl;
	// std::cout << binLL << std::endl;

	if(binLL < llEventCut){
	  int bin = hEvents->FindBin(easting, northing);
	  if(std::find(histBinsToConsider.begin(), histBinsToConsider.end(), bin)==histBinsToConsider.end()){
	    histBinsToConsider.push_back(bin);
	  }
	}
      }
    }

    int numEvents = 0;
    int nn = 0;
    for(UInt_t eventInd2=0; eventInd2 < events.size(); eventInd2++){
      Event& event2 = events.at(eventInd2);
      if(std::find(histBinsToConsider.begin(), histBinsToConsider.end(), event2.antarcticaHistBin) != histBinsToConsider.end()) {

	double ll = dSum(event, event2);
	if(ll > llEventCut){
	  ll = dFit(event, event2);
	}
	if(ll < llEventCut){
	  nn++;
	}

	numEvents++;
      }
    }
    p.inc(eventInd, events.size());
    std::cout << "eventInd = " << eventInd << ", must consider " << numEvents << " of " << events.size() << " over " << histBinsToConsider.size() << " bins " << std::endl;
    std::cout << "I found " << nn << " neighbours " << std::endl;
    std::cout << "the angular error deltaPhiMax = " << deltaPhiMax << ", deltaThetaMax = " << deltaThetaMax << std::endl;
    grandTotal += numEvents;
  }

  std::cout << "This brings to grand total to " << grandTotal << std::endl;
  std::cout << "done!" << std::endl;
  return;
}



void Acclaim::Clustering::LogLikelihoodMethod::nearbyEvents(Int_t eventInd, std::vector<Int_t>& nearbyEvents, std::vector<double>& nearbyEventLLs, double llRange, double llFitThreshold, double rangeEastingNorthing){

  // handle default arguments
  llFitThreshold = llFitThreshold < 0 ? llRange : llFitThreshold;

  Event& event = events.at(eventInd);
  event.cluster2[0] = 1;

  nearbyEvents.clear();
  nearbyEventLLs.clear();

  // first look up events close by in Easting/Northing
  std::vector<Int_t> nearbyEventIndsEastingNorthing;
  double lookup[2] = {event.easting, event.northing};
  fKDTree->FindInRange(lookup, rangeEastingNorthing, nearbyEventIndsEastingNorthing);

  // Loop through those and calculate the log-likelihoods
  for(UInt_t i=0; i < nearbyEventIndsEastingNorthing.size(); i++){
    int eventInd2 = nearbyEventIndsEastingNorthing[i];

    // only consider if eventInd2 > eventInd
    // this ensures that each pair only gets considered once
    // when 2 is in 1's lookup and not vice-versa
    if(eventInd2 > eventInd){
      Event& event2 = events.at(eventInd2);
      if(event2.eventEventClustering){

	double ll = dMin(event, event2); ///@todo make this handle non-direct events properly

	if(fDebug){
	  std::cout << "pair " << eventInd << ", " << eventInd2 << ", initial ll = " << ll;
	}
	if(ll > llFitThreshold){
	  ll = dFit(event, event2);
	  if(fDebug){
	    std::cout << ", fitted ll = " << ll;
	  }
	}

	// is this event's nearest neighbour?
	if(ll < event.nearestNeighbourLogLikelihood){
	  event.nearestNeighbourLogLikelihood = ll;
	  event.nearestNeighbourEventNumber = event2.eventNumber;
	}

	// is this event2's nearest neighbour?
	if(ll < event2.nearestNeighbourLogLikelihood){
	  event2.nearestNeighbourLogLikelihood = ll;
	  event2.nearestNeighbourEventNumber = event.eventNumber;
	}

	if(ll <= llRange){
	  nearbyEvents.push_back(eventInd2);
	  nearbyEventLLs.push_back(ll);
	  if(fDebug){
	    std::cerr  << ", added eventInd2 = " << eventInd2 << " with LL " << ll << " to nearbyEvents"  << std::endl;
	  }
	}
	if(fDebug){
	  std::cout << "\n";
	}
      }
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
    const Event& event = events.at(eventInd);
    Double_t llWais = dPoint(event, waisLon, waisLat, waisModelAlt);
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

  TH2D* hUnfitSqrt = new TH2D("hUnfitSqrt_d_vs_pointingAngle", "", 1024, 0, 180, 1024, 0, 1024);
  TH2D* hFitSqrt = new TH2D("hFitSqrt_d_vs_pointingAngle", "", 1024, 0, 180, 1024, 0, 1024);
  TH2D* hFitSqrtVsUnfittedSqrt = new TH2D("hFitSqrt_vs_unfitSqrt", "Fitted vs. unfitted WAIS event-event log-likelihood; Unfitted log-likelihood; Fitted log-likelihood", 1024, 0, 1024, 1024, 0, 1024);
  TH1D* hUnfitSqrtMinusFitSqrt = new TH1D("hFitSqrt_minus_unfitSqrt", "(Unfitted - fitted) WAIS event-event log-likelihood; #delta (log likelihood)", 1024, -512, 512);

  TH2DAntarctica* hEventPos = new TH2DAntarctica("hEventPos", "Event Position");
  TH2DAntarctica* hFittedPos = new TH2DAntarctica("hFittedPos", "Fitted Pair position");

  TGraphAntarctica* grAnita = new TGraphAntarctica();
  // TGraphAntarctica* grEventPos = new TGraphAntarctica();
  // TGraphAntarctica* grFittedPos = new TGraphAntarctica();

  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
    const Event& event1 = events.at(eventInd);

    // if(event1.eventNumber!=fTestEvent1) continue;

    TVector3 event1Pos = AntarcticCoord(AntarcticCoord::WGS84, event1.latitude, event1.longitude, event1.altitude).v();
    TVector3 anita1Pos = AntarcticCoord(AntarcticCoord::WGS84, event1.anita.latitude, event1.anita.longitude, event1.anita.altitude).v();
    TVector3 anitaToEvent1 = event1Pos - anita1Pos;

    const int eventInd2 = randy.Uniform(0, events.size());
    const Event& event2 = events.at(eventInd2);
    grAnita->SetPoint(grAnita->GetN(), event1.anita.longitude, event1.anita.latitude);
    // if(event2.eventNumber!=fTestEvent2) continue;

    if(event1.eventNumber==fTestEvent1 && event2.eventNumber==fTestEvent2){
      std::cerr << "Info in " << __PRETTY_FUNCTION__ << " mapping parameter space around WAIS for a single event test!" << std::endl;
      AnitaVersion::set(3);
      double waisEasting, waisNorthing;
      RampdemReader::LonLatToEastingNorthing(AnitaLocations::getWaisLongitude(), AnitaLocations::getWaisLatitude(), waisEasting, waisNorthing);
      double delta = 100e3;
      const int nBins = 1024;
      TH2D* hParams = new TH2D("hSingleEventTest", "Event-event fitted log likelihood; Easting (km); Northing (km); L_{sum}",
			       // nBins, -1090000, -1065000,
			       // nBins, -484200, -483600);
			       nBins, waisEasting-delta, waisEasting+delta,
			       nBins, waisNorthing-delta, waisNorthing+delta);

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

      TGraphAntarctica* grTestEvent1 = new TGraphAntarctica();
      grTestEvent1->SetPointEastingNorthing(0, event1.easting, event1.northing);
      grTestEvent1->SetName("grTestEvent1");
      grTestEvent1->Write();
      delete grTestEvent1;

      TGraphAntarctica* grTestEvent2 = new TGraphAntarctica();
      grTestEvent2->SetPointEastingNorthing(0, event2.easting, event2.northing);
      grTestEvent2->SetName("grTestEvent2");
      grTestEvent2->Write();
      delete grTestEvent2;

      TGraphAntarctica* grWaisTrue = new TGraphAntarctica();
      grWaisTrue->SetPoint(0, AnitaLocations::LONGITUDE_WAIS_A3, AnitaLocations::LATITUDE_WAIS_A3);
      grWaisTrue->SetName("grWaisTrue");
      grWaisTrue->Write();
      delete grWaisTrue;

    }

    TVector3 event2Pos = AntarcticCoord(AntarcticCoord::WGS84, event2.latitude, event2.longitude, event2.altitude).v();
    TVector3 anita2Pos = AntarcticCoord(AntarcticCoord::WGS84, event2.anita.latitude, event2.anita.longitude, event2.anita.altitude).v();
    TVector3 anitaToEvent2 = event2Pos - anita2Pos;

    double dist = dSum(event1, event2);
    Double_t angleBetweenEvents = anitaToEvent1.Angle(anitaToEvent2);


    hUnfit->Fill(angleBetweenEvents*TMath::RadToDeg(), dist);
    hUnfitSqrt->Fill(angleBetweenEvents*TMath::RadToDeg(), TMath::Sqrt(dist));

    // if(dist > 1000){
    //   std::cout << event1.eventNumber << "\t" << event2.eventNumber << std::endl;
    // }

    double distFitted = dFit(event1, event2);
    hFit->Fill(angleBetweenEvents*TMath::RadToDeg(), distFitted);
    hFitSqrt->Fill(angleBetweenEvents*TMath::RadToDeg(), TMath::Sqrt(distFitted));

    double fitLon, fitLat;
    RampdemReader::EastingNorthingToLonLat(fFitEasting, fFitNorthing, fitLon, fitLat);
    hFittedPos->Fill(fitLon, fitLat);
    hEventPos->Fill(event1.longitude, event1.latitude);
    // grFittedPos->SetPointEastingNorthing(grFittedPos->GetN(), fFitEasting, fFitNorthing);
    // grEventPos->SetPointEastingNorthing(grEventPos->GetN(), event1.easting, event1.northing);

    hFitVsUnfitted->Fill(dist, distFitted);
    hUnfitMinusFit->Fill(dist - distFitted);
    hFitSqrtVsUnfittedSqrt->Fill(TMath::Sqrt(dist), TMath::Sqrt(distFitted));
    hUnfitSqrtMinusFitSqrt->Fill(TMath::Sqrt(dist) - TMath::Sqrt(distFitted));

    p.inc(eventInd, events.size());
  }

  grAnita->SetName("grAnita");
  grAnita->Write();
  delete grAnita;

  hFittedPos->Write();
  delete hFittedPos;

  hEventPos->Write();
  delete hEventPos;

  // grFittedPos->SetName("grFittedPos");
  // grFittedPos->Write();
  // delete grFittedPos;

  // grEventPos->SetName("grEventPos");
  // grEventPos->Write();
  // delete grEventPos;


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





/**
 * An implementation of the event-to-event clustering with worst case O(N^2) efficiency
 */
void Acclaim::Clustering::LogLikelihoodMethod::doEventEventClustering(){

  ProgressBar p(events.size());
  for(UInt_t event1Ind=0; event1Ind < events.size(); event1Ind++){

    Event& event1 = events.at(event1Ind);
    if(event1.eventEventClustering){

      // look up nearby events
      std::vector<Int_t> event2Inds;
      std::vector<Double_t> eventEventLLs;
      nearbyEvents(event1Ind, event2Inds, eventEventLLs, llEventCut, llEventCut);

      std::vector<Int_t> matchedClusters;

      // Get clusters of matched events
      for(UInt_t i=0; i < event2Inds.size(); i++){
	UInt_t event2Ind = event2Inds.at(i);
	Event& event2 = events.at(event2Ind);
	if(event2.cluster[0] >= 0){
	  if(std::find(matchedClusters.begin(), matchedClusters.end(), event2.cluster[0])==matchedClusters.end()){
	    matchedClusters.push_back(event2.cluster[0]);
	  }
	}
      }

      // If none of the events are in a cluster, then add a new culster at event1
      if(matchedClusters.size()==0){

	// make a new cluster
	Cluster nc(event1, clusters.size());
	clusters.push_back(nc);

	matchedClusters.push_back(nc.index);
      }

      Int_t minClusterInd = TMath::LocMin(matchedClusters.size(), &matchedClusters[0]);
      Int_t thisCluster = matchedClusters.at(minClusterInd); // we will mark all events as in this cluster

      // mark event1 as in the minCluster
      if(event1.cluster[0] != thisCluster){
	if(event1.cluster[0] > 0){
	  clusters.at(event1.cluster[0]).numDataEvents--;
	}
	event1.cluster[0] = thisCluster;
	clusters.at(thisCluster).numDataEvents++;
      }

      // loop through all matched events, and put them in the appropriate cluster
      for(UInt_t i=0; i < event2Inds.size(); i++){
	Event& event2 = events.at(event2Inds.at(i));

	if(event2.cluster[0] >= 0){
	  clusters.at(event2.cluster[0]).numDataEvents--;
	}
	event2.cluster[0] = thisCluster;
	clusters.at(thisCluster).numDataEvents++;
      }

      // now loop through all events and if they match any of the matched clusters,
      // put them in the min cluster
      for(UInt_t event2Ind=0; event2Ind < events.size(); event2Ind++){
	Event& event2 = events.at(event2Ind);

	// if the event is in any of the matched clusters, reassign it
	if(std::find(matchedClusters.begin(), matchedClusters.end(), event2.cluster[0])!=matchedClusters.end()){
	  clusters.at(event2.cluster[0]).numDataEvents--;
	  event2.cluster[0] = thisCluster;
	  clusters.at(thisCluster).numDataEvents++;
	}
      }

      // quick check that we remove all events from the other clusters
      for(UInt_t m=0; m < matchedClusters.size(); m++){
	if(matchedClusters.at(m)!=thisCluster && clusters.at(matchedClusters.at(m)).numDataEvents > 0){
	  std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", expecting no events in cluster " << matchedClusters.at(m)
		    << ", have " << clusters.at(matchedClusters.at(m)).numDataEvents << std::endl;
	}
      }
    }

    p.inc(event1Ind);
  }
}



Int_t Acclaim::Clustering::LogLikelihoodMethod::removeLargeBasesNearMcMurdo(){

  fUseBaseList = true;
  if(fUseBaseList){
    initializeBaseList();
  }
  assignEventsToBaseClusters();
  const int manyEvents = 1000;
  int numRemoved=0;
  for(UInt_t clusterInd=0; clusterInd < clusters.size(); clusterInd++){
    Cluster& cluster = clusters.at(clusterInd);

    // remove these guys!
    if(cluster.numDataEvents >= manyEvents && isVaguelyNearMcMurdo(cluster)){
      for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
	events.at(eventInd).eventEventClustering = false;
	numRemoved++;
      }
    }
    else{
      cluster.resetClusteringNumbers();
      for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
	events.at(eventInd).resetClusteringNumbers();
      }
    }
  }
  return numRemoved;
}


void Acclaim::Clustering::LogLikelihoodMethod::doClusteringDBSCAN(const char* dataGlob, const char* mcGlob, const char* outFileName){

  readInSummaries(dataGlob);
  setDebug(false);
  // return;

  bool fRemoveLargeBasesNearMcMurdo = true;

  if(fRemoveLargeBasesNearMcMurdo){
    int numRemoved = removeLargeBasesNearMcMurdo();
    std::cout << "Removed " << numRemoved << " events from near McMurdo!" << std::endl;
  }
  fDebug=true;
  initKDTree();
  fDebug=false;

  // testAngleFindingSpeed();
  // DBSCAN();
  doEventEventClustering();

  // makeAndWriteNSquaredEventEventHistograms();
  const char* fakeArgv[2] = {outFileName, "DBSCAN"};
  OutputConvention oc(1, const_cast<char**>(fakeArgv));
  TFile* fOut = oc.makeFile();

  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
    Event& event = events.at(eventInd);
    event.antarcticaHistBin = hEvents->Fill(event.longitude, event.latitude);
  }
  hEvents->Write();

  writeAllGraphsAndHists();
  makeSummaryTrees();
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
  // readInSummaries(mcGlob);


  char* fakeArgv0 = const_cast<char*>(outFileName);
  OutputConvention oc(1, &fakeArgv0);
  TFile* fOut = oc.makeFile();

  assignEventsToBaseClusters();
  recursivelyAddClustersFromData(0);
  // redoSmallClusters();

  // assignMcEventsToClusters();

  writeAllGraphsAndHists();
  makeSummaryTrees();

  fOut->Write();
  fOut->Close();

}
