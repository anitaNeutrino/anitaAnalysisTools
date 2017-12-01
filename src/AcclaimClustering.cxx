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
#include "AcclaimOpenMP.h"

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
  usefulPat.setInterpSurfaceAboveGeoid(true);
  usefulPat.setSurfaceCloseEnoughInter(1e-3);
  usefulPat.setMaxLoopIterations(500); // make this arbitrarily large since it only happens once
  const double maxThetaAdjust = 8*TMath::DegToRad();
  usefulPat.traceBackToContinent(phi*TMath::DegToRad(), -theta*TMath::DegToRad(), &longitude, &latitude, &altitude, &thetaAdjustmentRequired, maxThetaAdjust, 100);
  if(latitude < -90){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", for eventNumber " << eventNumber << "\n";
    std::cerr << "Doing traceBackToContinenet again in debug mode!\n";
    usefulPat.setDebug(true);
    usefulPat.traceBackToContinent(phi*TMath::DegToRad(), -theta*TMath::DegToRad(), &longitude, &latitude, &altitude, &thetaAdjustmentRequired, maxThetaAdjust, 100);
    usefulPat.setDebug(false);
  }
  // selfLogLikelihood = logLikelihoodFromPoint(longitude, latitude, altitude, false);
  selfLogLikelihood = logLikelihoodFromPoint(longitude, latitude, altitude, true);
}


/**
 * Evaluate the log-likelihood distance from any event to an arbitrary point
 *
 * @param sourceLon is the source longitude
 * @param sourceLat is the source latitude
 * @param sourceAlt is the source altitude
 * @param addOverHorizonPenalty is an option to add a penalty term in LL for (approximately) being over the horizon
 *
 * @return the log-likelihood
 */
Double_t Acclaim::Clustering::Event::logLikelihoodFromPoint(Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt, bool addOverHorizonPenalty) const {

  Double_t thetaSource, phiSource;
  usefulPat.getThetaAndPhiWave2(sourceLon, sourceLat, sourceAlt, thetaSource, phiSource);
  thetaSource = -1*TMath::RadToDeg()*thetaSource;
  phiSource = TMath::RadToDeg()*phiSource;

  Double_t dTheta = (thetaSource - theta)/sigmaTheta;
  Double_t dPhi = Acclaim::RootTools::getDeltaAngleDeg(phiSource, phi)/sigmaPhi;

  Double_t ll = dTheta*dTheta + dPhi*dPhi;

  if(fDebug){
    std::cerr << __PRETTY_FUNCTION__ << " for " << eventNumber << ", dTheta = " << dTheta << ", dPhi = " << dPhi << ", ll = " << ll << std::endl;
  }

  if(addOverHorizonPenalty){
    Double_t distM = usefulPat.getDistanceFromSource(sourceLat, sourceLon, sourceAlt);
    if(distM >= default_horizon_distance){
      double distOverHorizonM = fabs(distM - default_horizon_distance);
      ll += distOverHorizonM;
    }
    if(fDebug){
      std::cerr << __PRETTY_FUNCTION__ << " for " << eventNumber << ", we are " << distM/1000 << "km from the source, after horizon penalty, ll = " << ll << std::endl;
    }    
  }

  return ll;
}




Acclaim::Clustering::Event::Event(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd, Int_t nT)
  : nThresholds(0), cluster(NULL),
    dThetaCluster(NULL), dPhiCluster(NULL)
{

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
  thetaAdjustmentRequired = peak.theta_adjustment_needed;
  selfLogLikelihood = -9999;
  anita = sum->anitaLocation;
  getAngularResolution(sum, pol, peakInd, sigmaTheta, sigmaPhi);
  antarcticaHistBin = -1;
  fDebug = false;
  nearestKnownBaseLogLikelihood = DBL_MAX;
  nearestKnownBaseCluster = -1;

  setNThresholds(nT);
  resetClusteringNumbers();
  setupUsefulPat();
}



Acclaim::Clustering::Event::Event(Int_t nT)
  : nThresholds(0), cluster(NULL),
    dThetaCluster(NULL), dPhiCluster(NULL)
{

  latitude = 0;
  longitude = 0;
  altitude = 0;
  easting = DBL_MAX;
  northing = DBL_MAX;
  theta = -9999;
  phi = -9999;
  thetaAdjustmentRequired = -9999;
  selfLogLikelihood = -9999;
  anita.reset();
  eventNumber = 0;
  run = 0;
  pol = AnitaPol::kNotAPol;
  peakIndex = -1;
  sigmaTheta = default_sigma_theta;
  sigmaPhi = default_sigma_phi;
  for(int dim=0; dim < nDim; dim++){
    centre[dim] = 0;
  }
  antarcticaHistBin = -1;
  fDebug = false;
  nearestKnownBaseLogLikelihood = DBL_MAX;
  nearestKnownBaseCluster = -1;

  setNThresholds(nT);
  resetClusteringNumbers();  
}


Acclaim::Clustering::Event::~Event(){
  deleteArrays();
}

Acclaim::Clustering::Event::Event(const Event& event){

  eventNumber = event.eventNumber;
  run = event.run;
  pol = event.pol;
  peakIndex = event.peakIndex;
  for(int i=0; i < 3; i++){
    centre[i] = event.centre[i];
  }
  latitude = event.latitude;
  longitude = event.longitude;
  altitude = event.altitude;
  easting = event.easting;
  northing = event.northing;
  anita = event.anita;
  theta = event.theta;
  phi = event.phi;
  sigmaTheta = event.sigmaTheta;
  sigmaPhi = event.sigmaPhi;
  thetaAdjustmentRequired = event.thetaAdjustmentRequired;
  selfLogLikelihood = event.selfLogLikelihood;
  nThresholds = event.nThresholds;
  cluster = new Int_t[nThresholds];
  dThetaCluster = new Double_t[nThresholds];
  dPhiCluster = new Double_t[nThresholds];
  for(int z=0; z < nThresholds; z++){
    cluster[z] = event.cluster[z];
    dThetaCluster[z] = event.dThetaCluster[z];
    dPhiCluster[z] = event.dPhiCluster[z];
  }
  // nearestNeighbourEventNumber = event.nearestNeighbourEventNumber;
  // nearestNeighbourLogLikelihood = event.nearestNeighbourLogLikelihood;
  eventEventClustering = event.eventEventClustering;
  nearestKnownBaseLogLikelihood = event.nearestKnownBaseLogLikelihood;
  nearestKnownBaseCluster = event.nearestKnownBaseCluster;
  antarcticaHistBin = event.antarcticaHistBin;
  fDebug = event.fDebug;
  usefulPat = event.usefulPat;
}


void Acclaim::Clustering::Event::deleteArrays(){
  if(cluster){
    delete [] cluster;
    cluster = NULL;
  }
  if(dPhiCluster){
    delete [] dPhiCluster;
    cluster = NULL;
  }
  if(dThetaCluster){
    delete [] dThetaCluster;
    dThetaCluster = NULL;
  }
  nThresholds = 0;
}



/** 
 * Set the N threshold arrays copying what was in there
 *  before into the new array, if long enough
 * 
 * @param n The size to set the internal arrays to
 */
void Acclaim::Clustering::Event::setNThresholds(Int_t n){

  std::vector<Int_t> tempCluster;
  std::vector<Double_t> tempDPhiCluster;
  std::vector<Double_t> tempDThetaCluster;
  
  if(nThresholds){
    tempCluster.insert(tempCluster.begin(), cluster, cluster+nThresholds);
    tempDPhiCluster.insert(tempDPhiCluster.begin(), dPhiCluster, dPhiCluster+nThresholds);
    tempDThetaCluster.insert(tempDThetaCluster.begin(), dThetaCluster, dThetaCluster+nThresholds);
  }

  deleteArrays();
  nThresholds = n;
  cluster = new Int_t[n];
  dPhiCluster = new Double_t[n];
  dThetaCluster = new Double_t[n];

  int nz = TMath::Min((int)tempCluster.size(), nThresholds);
  for(Int_t z=0; z < nz; z++){
    cluster[z] = tempCluster[z];
    dPhiCluster[z] = tempDPhiCluster[z];
    dThetaCluster[z] = tempDThetaCluster[z];
  }
  for(Int_t z=nz; z < nThresholds; z++){
    cluster[z] = -1;
    dPhiCluster[z] = -999;
    dThetaCluster[z] = -999;
  }
}



void Acclaim::Clustering::Event::resetClusteringNumbers(){

  for(Int_t z=0; z<nThresholds; z++){
    cluster[z] = -1;
    dThetaCluster[z] = -999;
    dPhiCluster[z] = -999;    
  }
  // nearestNeighbourEventNumber = 0;
  // nearestNeighbourLogLikelihood = DBL_MAX;
  // nearestKnownBaseLogLikelihood = DBL_MAX;
  // nearestKnownBaseCluster = -1;
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










Acclaim::Clustering::McEvent::McEvent(Int_t nT)
  : Event(nT){
  weight = 0;
  energy=0;
}

Acclaim::Clustering::McEvent::McEvent(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd, Int_t nT)
  : Event(sum, pol,  peakInd, nT)
{
  weight = sum->mc.weight;
  energy = sum->mc.energy;
  // std::cout << longitude << "\t" << latitude << "\t" << altitude << std::endl;
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
  AntarcticCoord ac = base.position.as(AntarcticCoord::WGS84);
  latitude = ac.x;
  longitude = ac.y;
  altitude = ac.z;
  knownBase = 1;

  if(altitude == -999){
    altitude = RampdemReader::BilinearInterpolatedSurfaceAboveGeoid(longitude, latitude);
  }

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->getCartesianCoords(latitude, longitude, altitude, centre);
  resetClusteringNumbers();
  antarcticaHistBin = -1;
  seedEvent = -1;
  index = i;
  llEventCutInd = 0;
}


Acclaim::Clustering::Cluster::Cluster(const Event& event, Int_t i) {
  latitude = event.longitude;
  longitude = event.latitude;
  altitude = event.altitude;
  knownBase = 0;

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->getCartesianCoords(latitude, longitude, altitude, centre);
  resetClusteringNumbers();
  index = i;
  llEventCutInd = 0;
}


void Acclaim::Clustering::Cluster::resetClusteringNumbers(){
  numDataEvents = 0;
  sumMcWeights = 0;
}



Acclaim::Clustering::LogLikelihoodMethod::LogLikelihoodMethod()
  : fROOTgErrorIgnoreLevel(gErrorIgnoreLevel)
{
  grTestMinimizerWalk = NULL;
  grTestMinimizerValue = NULL;  
  llClusterCut = 100;

  // llEventCuts.push_back(10);
  // llEventCuts.push_back(20);
  llEventCuts.push_back(40);
  llEventCuts.push_back(70);

  llEventCuts.push_back(100);
  llEventCuts.push_back(150);
  llEventCuts.push_back(200);
  llEventCuts.push_back(250);

  llEventCuts.push_back(350);
  llEventCuts.push_back(400);
  llEventCuts.push_back(450);
  llEventCuts.push_back(500);

  llEventCuts.push_back(600);
  llEventCuts.push_back(700);
  llEventCuts.push_back(800);
  llEventCuts.push_back(900);

  llEventCuts.push_back(1000);
  llEventCuts.push_back(1200);
  llEventCuts.push_back(1400);
  llEventCuts.push_back(1600);

  llEventCuts.push_back(2000);
  llEventCuts.push_back(2500);
  llEventCuts.push_back(3000);
  llEventCuts.push_back(4000);

  for(UInt_t z=0; z < llEventCuts.size(); z++){
    clusters.push_back(std::vector<Cluster>());
  }

  // both above horizon; good improvement; actual test
  fTestEvent1 = 55510391; // 61156660; //61033430;
  fTestEvent2 = 61338514; //55789194; //61424151;

  numCallsToRecursive = 0;
  doneBaseClusterAssignment = false;
  fKDTree = NULL;
  hClusters = new TH2DAntarctica("hClusters", "hClusters");
  hEvents = new TH2DAntarctica("hEvents", "hEvents");
  hMcEvents = new TH2DAntarctica("hMcEvents", "hMcEvents");

  fDebug = false;
  fUseBaseList = true;

  fMaxFitterAttempts = 1;

  int nThreads = Acclaim::OpenMP::getMaxThreads();
  fFitEvent1s.resize(nThreads, NULL);
  fFitEvent2s.resize(nThreads, NULL);
  fFitEastings.resize(nThreads, 0);
  fFitNorthings.resize(nThreads, 0);

  fFunctors.reserve(nThreads);
  fMinimizers.reserve(nThreads);
  for(int t=0; t < nThreads; t++){
    fMinimizers.push_back(ROOT::Math::Factory::CreateMinimizer("Minuit2"));
    fFunctors.push_back(ROOT::Math::Functor(this, &Acclaim::Clustering::LogLikelihoodMethod::evalPairLogLikelihoodAtLonLat, 2));    

    fMinimizers.at(t)->SetMaxFunctionCalls(1e5); // for Minuit/Minuit2
    fMinimizers.at(t)->SetTolerance(0.0001);
    fMinimizers.at(t)->SetPrintLevel(0);
    fMinimizers.at(t)->SetFunction(fFunctors.at(t));
  }
}





Acclaim::Clustering::LogLikelihoodMethod::~LogLikelihoodMethod(){

  // delete non-NULL antarctic histograms, tgraphs

  for(UInt_t i=0; i < hUnclusteredEvents.size(); i++){
    if(hUnclusteredEvents.at(i)){
      delete hUnclusteredEvents.at(i);
      hUnclusteredEvents.at(i) = NULL;
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

  Double_t thetaWave, phiWave;
  event.usefulPat.getThetaAndPhiWave2(cluster.longitude, cluster.latitude, cluster.altitude, thetaWave, phiWave);
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
  int t = OpenMP::thread();

  // Double_t sourceAlt = RampdemReader::SurfaceAboveGeoidEN(sourceEasting, sourceNorthing);
  Double_t sourceAlt = RampdemReader::BilinearInterpolatedSurfaceAboveGeoidEastingNorthing(sourceEasting, sourceNorthing);  
  Double_t sourceLon, sourceLat;
  RampdemReader::EastingNorthingToLonLat(sourceEasting, sourceNorthing, sourceLon, sourceLat);

  Double_t ll = 0;
  ll += fFitEvent1s.at(t)->logLikelihoodFromPoint(sourceLon, sourceLat, sourceAlt, true);
  ll += fFitEvent2s.at(t)->logLikelihoodFromPoint(sourceLon, sourceLat, sourceAlt, true);
  // ll += fFitEvent1s.at(t)->logLikelihoodFromPoint(sourceLon, sourceLat, sourceAlt, true);
  // ll += fFitEvent2s.at(t)->logLikelihoodFromPoint(sourceLon, sourceLat, sourceAlt, true);

  if(fFitEvent1s.at(t)->eventNumber==fTestEvent1){
    if(fFitEvent2s.at(t)->eventNumber==fTestEvent2){
      if(grTestMinimizerWalk){
	grTestMinimizerWalk->SetPoint(grTestMinimizerWalk->GetN(), sourceEasting, sourceNorthing);
      }
      if(grTestMinimizerValue){
	int n = grTestMinimizerValue->GetN();
	grTestMinimizerValue->SetPoint(n, n, ll);
	// std::cout << n << "\t" << sourceEasting << "\t" << sourceNorthing << "\t"
	// 	  << sourceAlt << "\t" << RampdemReader::SurfaceAboveGeoidEN(sourceEasting, sourceNorthing)
	// 	  << "\t" << ll << std::endl;
      }
    }
  }
  

  if(fMinimizers.at(t)->PrintLevel() > 0){
    std::cout << sourceEasting << "\t" << sourceNorthing << "\t" << sourceLon << "\t" << sourceLat << "\t" << ll << std::endl;
  }
  
  return ll;
}



Double_t Acclaim::Clustering::LogLikelihoodMethod::dFit(const Event* event1, const Event* event2){

  int t = Acclaim::OpenMP::thread();
  ROOT::Math::Minimizer* minimizer = fMinimizers.at(t);
  
  fFitEvent1s.at(t) = event1;
  fFitEvent2s.at(t) = event2;

  Double_t ll12 = dAsym(event1, event2);
  Double_t ll21 = dAsym(event2, event1);

  // const Event* which = ll12 < ll21 ? fFitEvent1s.at(t) : fFitEvent2s.at(t);
  const Event* which = ll21 < ll12 ? fFitEvent1s.at(t) : fFitEvent2s.at(t);  

  if(fFitEvent1s.at(t)->eventNumber==fTestEvent1 && fFitEvent2s.at(t)->eventNumber==fTestEvent2){
    grTestMinimizerWalk = new TGraph();
    grTestMinimizerValue = new TGraph();

    if(fDebug){
      event1->fDebug = true;
      event2->fDebug = true;
    }
  }

  
  minimizer->SetVariable(0, "sourceEasting", FITTER_INPUT_SCALING*which->easting, 1);
  minimizer->SetVariable(1, "sourceNorthing", FITTER_INPUT_SCALING*which->northing, 1);

  bool validMinimum = false;
  std::vector<double> minima;
  for(int attempt=0; attempt < fMaxFitterAttempts && !validMinimum; attempt++){
    validMinimum = minimizer->Minimize();
    minima.push_back(minimizer->MinValue());
    if(!validMinimum){
      if(fDebug){
	std::cerr << "Info in " << __PRETTY_FUNCTION__ << ", eventNumbers " << event1->eventNumber << " and "
		  << event2->eventNumber << " did not converge, the result was " << minimizer->MinValue() << std::endl;
      }
      minimizer->SetVariable(0, "sourceEasting", minimizer->X()[0], 1);
      minimizer->SetVariable(1, "sourceNorthing", minimizer->X()[1], 1);
    }
  }
  if(!validMinimum && fDebug){
    std::cerr << "Info in " << __PRETTY_FUNCTION__ << ", eventNumbers " << event1->eventNumber << " and "
	      << event2->eventNumber << " did not converge, the result was " << minimizer->MinValue() << std::endl;
  }

  if(fFitEvent1s.at(t)->eventNumber==fTestEvent1 && fFitEvent2s.at(t)->eventNumber==fTestEvent2){
    if(fDebug){
      event1->fDebug = false;
      event2->fDebug = false;
    }
  }

  if(grTestMinimizerWalk){
    grTestMinimizerWalk->SetName("grTestMinimizerWalk");
    grTestMinimizerWalk->Write();
    delete grTestMinimizerWalk;
    grTestMinimizerWalk = NULL;
  }
  if(grTestMinimizerValue){
    grTestMinimizerValue->SetName("grTestMinimizerValue");
    grTestMinimizerValue->Write();
    delete grTestMinimizerValue;
    grTestMinimizerValue = NULL;
  }

  if((!validMinimum && fDebug)){ // do the same thing again, this time with error messages!
    int old_print_level = minimizer->PrintLevel();
    gErrorIgnoreLevel = fROOTgErrorIgnoreLevel;
    minimizer->SetPrintLevel(3);
    minimizer->SetVariable(0, "sourceEasting", FITTER_INPUT_SCALING*which->easting, 1);
    minimizer->SetVariable(1, "sourceNorthing", FITTER_INPUT_SCALING*which->northing, 1);
    for(int attempt=0; attempt < fMaxFitterAttempts && !validMinimum; attempt++){
      validMinimum = minimizer->Minimize();
      if(!validMinimum){
	minimizer->SetVariable(0, "sourceEasting", minimizer->X()[0], 1);
	minimizer->SetVariable(1, "sourceNorthing", minimizer->X()[1], 1);
      }
    }
    minimizer->SetPrintLevel(old_print_level);

    AnitaVersion::set(3);
    double waisLon = AnitaLocations::getWaisLongitude();
    double waisLat = AnitaLocations::getWaisLatitude();
    double waisModelAlt = RampdemReader::BilinearInterpolatedSurfaceAboveGeoid(waisLon, waisLat);
    const Event& event1 = *fFitEvent1s.at(t);
    const Event& event2 = *fFitEvent2s.at(t);
    Double_t llWais1 = event1.logLikelihoodFromPoint(waisLon, waisLat, waisModelAlt, false);
    Double_t llWais2 = event2.logLikelihoodFromPoint(waisLon, waisLat, waisModelAlt, false);
    Double_t llWais1b = event1.logLikelihoodFromPoint(waisLon, waisLat, waisModelAlt, true);
    Double_t llWais2b = event2.logLikelihoodFromPoint(waisLon, waisLat, waisModelAlt, true);
    std::cerr << "Debug in " << __PRETTY_FUNCTION__ << ": The fitter minimum is " << minimizer->MinValue() << std::endl;
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
    gErrorIgnoreLevel = 1001;
  }

  fFitEastings.at(t) = minimizer->X()[0]*FITTER_OUTPUT_SCALING;
  fFitNorthings.at(t) = minimizer->X()[1]*FITTER_OUTPUT_SCALING;

  if(fFitEvent1s.at(t)->eventNumber==fTestEvent1 && fFitEvent2s.at(t)->eventNumber==fTestEvent2){
    grTestMinimumPosition = new TGraph(1, &fFitEastings.at(t), &fFitNorthings.at(t));
    grTestMinimumPosition->SetName("grTestMinimumPosition");
    grTestMinimumPosition->Write();
    delete grTestMinimumPosition;
    grTestMinimumPosition = NULL;
  }
  Double_t ll = minimizer->MinValue();

  return ll;
}




/** 
 * Return the Minimum of the self and asymetric log-likelihoods for both events.
 * i.e. the smaller of:
 * ANITA at the payload position during event1, viewing event1 and event2
 * and
 * ANITA at the payload position during event2, viewing event1 and event2
 * 
 * (self log-likelihood should be (very close to) 0 for events which reconstruct to the continent
 * 
 * @param event1 is a clustering event (order of arguments does not matter)
 * @param event2 is a clustering event (order of arguments does not matter)
 * 
 * @return Sum of log-likelihoods
 */
Double_t Acclaim::Clustering::LogLikelihoodMethod::dMin(const Event* event1, const Event* event2){
  return TMath::Min(dAsym(event1, event2) + event1->selfLogLikelihood, dAsym(event2, event1) + event2->selfLogLikelihood);
}


/** 
 * Return the SUM of the self and asymetric log-likelihoods for both events.
 * i.e.
 * ANITA at the payload position during event1, viewing event1 and event2
 * +
 * ANITA at the payload position during event2, viewing event1 and event2
 * 
 * (self log-likelihood should be (very close to) 0 for events which reconstruct to the continent
 * 
 * @param event1 is a clustering event (order of arguments does not matter)
 * @param event2 is a clustering event (order of arguments does not matter)
 * 
 * @return Sum of log-likelihoods
 */
Double_t Acclaim::Clustering::LogLikelihoodMethod::dSum(const Event* event1, const Event* event2){
  return dAsym(event1, event2) + event1->selfLogLikelihood + dAsym(event2, event1) + event2->selfLogLikelihood;
}


/** 
 * This the log likelihood of the source position of event2 in event1's frame of reference
 * 
 * @param event1 ANITA's position from here is used
 * @param event2 The source position from here is used.
 * 
 * @return The log-likelihood for seeing an event at event2's position, when ANITA is at event1.
 */
Double_t Acclaim::Clustering::LogLikelihoodMethod::dAsym(const Event* event1, const Event* event2){
  return event1->logLikelihoodFromPoint(event2->longitude, event2->latitude,event2->altitude, true);
}




Double_t Acclaim::Clustering::LogLikelihoodMethod::getAngDistSqEventCluster(const Event& event, const Cluster& cluster){

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



/** 
 * Here we find the next event to cluster.
 * This is done by histogramming unclustered events 
 * (unclustered with the smallest log likelihood threshold).
 * 
 * @return pointer to the next event to process, NULL when there are no remaining events
 */
Acclaim::Clustering::Event* Acclaim::Clustering::LogLikelihoodMethod::nextEvent(){

  // to do this, first we histogram all the events which aren't clustered (using the smallest log-likelihood threshold).
  TString name = TString::Format("hUnclusteredEvents_%lu", hUnclusteredEvents.size());
  TH2DAntarctica* h = new TH2DAntarctica(name, name);
  hUnclusteredEvents.push_back(h);

  h->setAcceptStereographic(true); // already have easting/northing, don't do unnecessary work

  // find unclustered events (using the smallest log-likelihood threshold)
  UInt_t numUnclustered=0;
  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
    Event& event = events.at(eventInd);
    if(event.cluster[0] < 0){
      event.antarcticaHistBin = h->Fill(event.easting, event.northing);
      numUnclustered++;
    }
  }
  h->setAcceptStereographic(false);

  Event* nextEvent = NULL;

  // then we're done...
  if(numUnclustered==0){
    delete h;
    hUnclusteredEvents.pop_back();
  }
  else{

    // here we pick the next event...
    // which is the one closest to the mean position
    // of all the unclustered events in this bin...
    // hopefully this will be close to the centre of
    // all unclustered events.
    
    Int_t globalMaxBin = h->GetMaximumBin();
    Int_t numUnclusteredInBin = h->GetBinContent(globalMaxBin);

    Double_t meanEasting=0;
    Double_t meanNorthing=0;

    for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
      const Event& event = events.at(eventInd);
      if(event.antarcticaHistBin==globalMaxBin && event.cluster[0] < 0){
	meanEasting += event.easting;
	meanNorthing += event.northing;
      }
    }
    meanNorthing/=numUnclusteredInBin;
    meanEasting/=numUnclusteredInBin;

    Double_t bestSurfaceSeparationSquared = DBL_MAX;
    for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
      Event& event = events.at(eventInd);
      if(event.antarcticaHistBin==globalMaxBin && event.cluster[0] < 0){
	Double_t dE = event.easting - meanEasting;
	Double_t dN = event.northing - meanNorthing;
	Double_t surfaceSeparationSquared = dE*dE + dN*dN;

	if(surfaceSeparationSquared < bestSurfaceSeparationSquared){
	  bestSurfaceSeparationSquared = surfaceSeparationSquared;
	  nextEvent = &event;
	}
      }
    }
  }
  return nextEvent;
}







void Acclaim::Clustering::LogLikelihoodMethod::assignSingleEventToCloserCluster(Int_t i, Int_t isMC, Cluster& cluster, Int_t z, double llCut){
  Event& event = isMC==0 ? events.at(i) : mcEvents.at(i);

  // std::cout << llCut << "\t";
  llCut = llCut < 0 ? cluster.llEventCut : llCut; // default llCut < 0 means use cluster stored value
  // if(cluster.index==19){
  //   std::cout << cluster.llEventCut << "\t" << llCut << std::endl;
  // }

  Double_t distM = event.usefulPat.getDistanceFromSource(cluster.latitude, cluster.longitude, cluster.altitude);

  if(distM < default_horizon_distance){ // are we even close?

    Double_t ll = getAngDistSqEventCluster(event, cluster);

    if(ll <= event.nearestKnownBaseLogLikelihood){
      event.nearestKnownBaseLogLikelihood = ll;
      event.nearestKnownBaseCluster = cluster.index;

      if(ll < llCut){
	if(!isMC){
	  int previousBestClusterIndex = event.cluster[z];
	  if(previousBestClusterIndex >= 0){
	    clusters.at(z).at(previousBestClusterIndex).numDataEvents--;
	    // if(clusters.at(z).at(previousBestClusterIndex).numDataEvents < 0){
	    //   std::cout << z << "\t" << clusters.at(z).at(previousBestClusterIndex).index << "\t" << clusters.at(z).at(previousBestClusterIndex).numDataEvents << std::endl;
	    // }
	  }
	  cluster.numDataEvents++;
	}
	event.cluster[z] = cluster.index;
      }
    }
  }
}




void Acclaim::Clustering::LogLikelihoodMethod::assignMcEventsToClusters(){

  const int isMC = 1;

  UInt_t nMc = mcEvents.size();

  std::cerr  << "Info in " << __PRETTY_FUNCTION__ << ": starting!" << std::endl;
  ProgressBar p(nMc);
  for(UInt_t j=0; j < nMc; j++){
    for(UInt_t z=0; z < llEventCuts.size(); z++){
      for(UInt_t clusterInd=0; clusterInd < clusters.at(z).size(); clusterInd++){
	Cluster& cluster = clusters.at(z).at(clusterInd);
	assignSingleEventToCloserCluster(j, isMC, cluster, z);
      }
    }
    p.inc(j, nMc);
  }
}









size_t Acclaim::Clustering::LogLikelihoodMethod::addMcEvent(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd){

  mcEvents.push_back(McEvent(sum,  pol, peakInd, llEventCuts.size()));

  return mcEvents.size();
}






size_t Acclaim::Clustering::LogLikelihoodMethod::addEvent(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd){

  events.push_back(Event(sum, pol, peakInd, llEventCuts.size()));

  return events.size();
}







void Acclaim::Clustering::LogLikelihoodMethod::forEachEventFindClosestKnownBase(int z){

  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
    Event& event = events.at(eventInd);
    for(UInt_t clusterInd=0; clusterInd < clusters.at(z).size(); clusterInd++){
      Cluster& cluster = clusters.at(z).at(clusterInd);

      Double_t distM = event.usefulPat.getDistanceFromSource(cluster.latitude, cluster.longitude, cluster.altitude);
      if(distM < default_horizon_distance){ // are we even close?

	Double_t ll = event.logLikelihoodFromPoint(cluster, true);

	if(ll <= event.nearestKnownBaseLogLikelihood){
	  event.nearestKnownBaseLogLikelihood = ll;
	  event.nearestKnownBaseCluster = cluster.index;
	}
      }
    }
  }
}




/** 
 * Puts an entry in each of the cluster[z] vectors for each of the known bases
 */
void Acclaim::Clustering::LogLikelihoodMethod::initializeBaseList(){

  std::cout << "Info in " << __PRETTY_FUNCTION__ << ": Initializing base list..." << std::endl;

  // make a copy for each llCut, just to ease the book keeping later
  for(UInt_t z=0; z < llEventCuts.size(); z++){
    for(UInt_t clusterInd=0; clusterInd < BaseList::getNumBases(); clusterInd++){
      const BaseList::base& base = BaseList::getBase(clusterInd);
      clusters.at(z).push_back(Cluster(base, clusters.at(z).size()));
      clusters.at(z).back().llEventCutInd = z;
      clusters.at(z).back().llEventCut = llEventCuts.at(z);
    }
  }
}



void Acclaim::Clustering::LogLikelihoodMethod::resetClusters(){

  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
    for(UInt_t z=0; z < llEventCuts.size(); z++){
      events.at(eventInd).cluster[z] = -1;
    }
  }

  for(UInt_t eventInd=0; eventInd < mcEvents.size(); eventInd++){
    for(UInt_t z=0; z < llEventCuts.size(); z++){    
      mcEvents.at(eventInd).cluster[z] = -1;
    }
  }

  clusters.resize(0);
  doneBaseClusterAssignment = false;
}


































// TGraphAntarctica* Acclaim::Clustering::LogLikelihoodMethod::makeClusterSummaryTGraph(Int_t clusterInd){

//   TGraphAntarctica* gr = NULL;
//   if(clusterInd >= 0 && clusterInd < (Int_t)clusters.size()){

//     TString name  = TString::Format("grCluster%d", clusterInd);
//     TString title  = TString::Format("Cluster %d; Easting (m); Northing (m)", clusterInd);
//     gr = new TGraphAntarctica();
//     gr->SetName(name);
//     gr->SetTitle(title);

//     // AnitaGeomTool* geom = AnitaGeomTool::Instance();

//     for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
//       if(events.at(eventInd).cluster[0]==clusterInd){
// 	gr->SetPoint(gr->GetN(), events.at(eventInd).longitude, events.at(eventInd).latitude);
//       }
//     }
//   }
//   return gr;
// }





void Acclaim::Clustering::LogLikelihoodMethod::makeSummaryTrees(){

  // first set the unknown cluster locations as the mean of the 
  // cartesian positions
  for(UInt_t z=0; z < clusters.size(); z++){
    for(UInt_t clusterInd=0; clusterInd < clusters.at(z).size(); clusterInd++){
      Cluster& cluster = clusters.at(z).at(clusterInd);

      if(cluster.knownBase==0){
    
	cluster.centre[0] = 0; // reset the
	cluster.centre[1] = 0; // cartesian
	cluster.centre[2] = 0; // coordinates
	int eventCounter = 0;
    
	for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
	  const Event& event = events.at(eventInd);
	  if(event.cluster[z]==(Int_t)clusterInd){
	    cluster.centre[0] += event.centre[0];
	    cluster.centre[1] += event.centre[1];
	    cluster.centre[2] += event.centre[2];
	    eventCounter++;
	  }
	}
	if(eventCounter > 0){
	  cluster.centre[0]/=eventCounter;
	  cluster.centre[1]/=eventCounter;
	  cluster.centre[2]/=eventCounter;
	}

	AnitaGeomTool* geom = AnitaGeomTool::Instance();
	geom->getLatLonAltFromCartesian(cluster.centre, cluster.latitude,
					cluster.longitude, cluster.altitude);

	if(eventCounter != cluster.numDataEvents){
	  std::cerr << "Error in " << __PRETTY_FUNCTION__
		    << ": was expecting " << cluster.numDataEvents
		    << " in cluster " << clusterInd << ", but counted "
		    << eventCounter << std::endl;
	}
      }
      // std::cout << "in summary making " << z << "\t" << clusterInd << "\t" << cluster.numDataEvents << std::endl;
    }
  }

  
  TTree* eventTree = new TTree("eventTree", "Tree of clustered ANITA events");
  Event* event = NULL;
  eventTree->Branch("event", &event);
  
  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
    event = &events.at(eventInd);

    for(UInt_t z = 0; z < llEventCuts.size(); z++){
      if(event->cluster[z] >= 0){
	const Cluster& cluster = clusters.at(z).at(event->cluster[z]);
	getDeltaThetaDegDeltaPhiDegEventCluster(*event, cluster, event->dThetaCluster[z], event->dPhiCluster[z]);
      }
      else {
	// std::cerr << "Error in " << __PRETTY_FUNCTION__
	// 	  << ": eventNumber " << event->eventNumber
	// 	  << " has a cluster[" << z << "] number "
	// 	  << event->cluster[z] << std::endl;
      }
    }
    eventTree->Fill();
  }

  TTree* mcEventTree = new TTree("mcEventTree", "Tree of clustered Monte Carlo ANITA events");
  McEvent* mcEvent = NULL;
  mcEventTree->Branch("mcEvent", &mcEvent);

  for(UInt_t j=0; j < mcEvents.size(); j++){
    mcEvent = &mcEvents.at(j);
    mcEventTree->Fill();
  }

  for(UInt_t z = 0; z < llEventCuts.size(); z++){
    TString treeName = TString::Format("clusterTree%u", z);
    TString treeTitle = TString::Format("Tree of clusters with llEventCut = %lf", llEventCuts[z]);
    TTree* clusterTree = new TTree(treeName, "Tree of clusters");
    Cluster* cluster = NULL;
    clusterTree->Branch("cluster", &cluster);
    std::vector<Cluster> &cs = clusters.at(z);
    for(Int_t k=0; k < (Int_t)clusters.at(z).size(); k++){
      cluster = &cs.at(k);
      clusterTree->Fill();
    }
  }
}


Long64_t Acclaim::Clustering::LogLikelihoodMethod::readInSummaries(const char* summaryGlob){

  Long64_t n = 0;
  if(summaryGlob){
    SummarySet ss(summaryGlob);
    n = ss.N();
    Int_t numReadIn = 0;
    std::cout << "Info in " << __PRETTY_FUNCTION__ << ": reading in summaries: " << summaryGlob << std::endl;

    bool useSandbox = false;
    // bool notUsingSandbox = TString(summaryGlob).Contains("wais");    

    ProgressBar p(n);

    events.reserve(n);    
    
    for(Long64_t entry=0; entry < n; entry++){

      ss.getEntry(entry);
      AnitaEventSummary* sum = ss.summary();

      // AnitaPol::AnitaPol_t pol = sum->trainingPol();
      // Int_t peakIndex = sum->trainingPeakInd();
      AnitaPol::AnitaPol_t pol = sum->trainingPol();
      Int_t peakIndex = sum->trainingPeakInd();
      // AnitaPol::AnitaPol_t pol = sum->mostImpulsivePol(1);
      // Int_t peakIndex = sum->mostImpulsiveInd(1);

      // std::cout << sum->eventNumber << "\t" << pol << "\t" << peakIndex << std::endl;
      Double_t snrHack = sum->deconvolved_filtered[pol][peakIndex].snr;

      /// @todo remove the snr > 160 hack???
      /// @todo remove the run hack!!!
      if(sum->run > 160
	 && (!useSandbox || inSandbox(sum->peak[pol][peakIndex]))
	 && sum->peak[pol][peakIndex].theta < 0
	 && snrHack < 100){
	
	if(sum->mc.weight > 0){
	  if(entry==0){
	    mcEvents.reserve(mcEvents.size() + n);
	  }
	  // if(sum->eventNumber==914599568){
	  // if((entry%10)==0){	    
	  addMcEvent(sum,  pol, peakIndex);
	  numReadIn++;	    
	  // }
	}
	else{
	  if(entry==0){
	    events.reserve(events.size() + n);
	  }	  
	  addEvent(sum, pol, peakIndex);	
	  numReadIn++;
	}
      }
      p.inc(entry, n);
    }
    std::cout << "Read in " << numReadIn << " summaries" << std::endl;
  }
  return n;
}


void Acclaim::Clustering::LogLikelihoodMethod::writeAllGraphsAndHists(){
  for(UInt_t eventInd=0; eventInd < hUnclusteredEvents.size(); eventInd++){
    if(hUnclusteredEvents.at(eventInd)){
      hUnclusteredEvents.at(eventInd)->Write();
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

      Double_t d_ij = dMin(&event_i, &event_j);

      for(UInt_t k=0; k < numENNeighbours; k++){
	const Event& event_k = events.at(k);
	Double_t d_ik = dMin(&event_i, &event_k);
	Double_t d_jk = dMin(&event_j, &event_k);

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

  const double llEventCutSqrt = TMath::Sqrt(llEventCuts[0]);

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
      double underSqrtSign = llEventCuts[0] - deltaPhiNorm*deltaPhiNorm;
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
	alt = RampdemReader::BilinearInterpolatedSurfaceAboveGeoid(lon, lat);

	Double_t theta, phi;
	event.usefulPat.getThetaAndPhiWave2(lon, lat, alt, theta, phi);
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

	if(binLL < llEventCuts[0]){
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

	// double ll = dSum(event, event2);
	double ll = dMin(&event, &event2);
	if(ll > llEventCuts[0]){
	  ll = dFit(&event, &event2);
	}
	if(ll < llEventCuts[0]){
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









void Acclaim::Clustering::LogLikelihoodMethod::nearbyEvents(Int_t eventInd, std::vector<Int_t>& nearbyEvents, std::vector<double>& nearbyEventLLs, bool mc, double llRange, double llFitThreshold, double rangeEastingNorthing){

  // handle default arguments
  llFitThreshold = llFitThreshold < 0 ? llRange : llFitThreshold;

  Event* event = mc ? &mcEvents.at(eventInd) : &events.at(eventInd);

  nearbyEvents.clear();
  nearbyEventLLs.clear();

  // first look up events close by in Easting/Northing
  std::vector<Int_t> nearbyEventIndsEastingNorthing;
  double lookup[2] = {event->easting, event->northing};
  fKDTree->FindInRange(lookup, rangeEastingNorthing, nearbyEventIndsEastingNorthing);

  // std::cout << rangeEastingNorthing << "\t" << nearbyEventIndsEastingNorthing.size() << std::endl;
  // /**
  //  * since we only get event2Ind > event1Ind, we need to allow the closest
  //  * event from previously to ensure you can get in the correct cluster
  //  * if we didn't do this, the final event in a real cluster of events on
  //  * the continent would't be considered by the algorithm.
  //  */
  // if(event->nearestNeighbourLogLikelihood < llRange){
  //   for(UInt_t j=0; j < events.size(); j++){
  //     if(events.at(j).eventNumber==event->nearestNeighbourEventNumber){
  // 	nearbyEvents.push_back(j);
  // 	nearbyEventLLs.push_back(event->nearestNeighbourLogLikelihood);
  // 	break;
  //     }
  //   }
  // }


  // move swapping global error level stuff out of parallelized loop
  // this used to be around the minimizer->minimize() 
  fROOTgErrorIgnoreLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = 1001;


  // // here we remove the events where eventNumber2 > eventNumber1 in advance of the threading
  // // to try and balance the workload between the threads...
  // std::vector<Int_t> nearbyEventIndsEastingNorthingTrimmed;
  // nearbyEventIndsEastingNorthingTrimmed.reserve(nearbyEventIndsEastingNorthing.size());
  // for(UInt_t i=0; i < nearbyEventIndsEastingNorthing.size(); i++){
  //   int eventInd2 = nearbyEventIndsEastingNorthing.at(i);
  //   if(events.at(eventInd2).eventNumber > event->eventNumber){
  //     nearbyEventIndsEastingNorthingTrimmed.push_back(nearbyEventIndsEastingNorthing.at(i));
  //   }
  // }

  int mt = OpenMP::getMaxThreads();
  std::vector<std::vector<Int_t> > nearbyEventsParallel(mt);
  std::vector<std::vector<Double_t> > nearbyEventLLsParallel(mt);
  // std::vector<Double_t > event1MinLL(mt, DBL_MAX);
  // std::vector<UInt_t> event1MinLLEvent2EventNumber(mt, 0);
#pragma omp parallel for
  // for(UInt_t i=0; i < nearbyEventIndsEastingNorthingTrimmed.size(); i++){
  //   int eventInd2 = nearbyEventIndsEastingNorthingTrimmed[i];
  for(UInt_t i=0; i < nearbyEventIndsEastingNorthing.size(); i++){
    int eventInd2 = nearbyEventIndsEastingNorthing[i];

    int t = OpenMP::thread();
    std::vector<Int_t>& nearbyEventsLoop = OpenMP::isEnabled ? nearbyEventsParallel.at(t) : nearbyEvents;
    std::vector<Double_t>& nearbyEventLLsLoop = OpenMP::isEnabled ? nearbyEventLLsParallel.at(t) : nearbyEventLLs;
    Event& event2 = events.at(eventInd2);

    if(event->eventNumber != event2.eventNumber){
      if(event2.eventEventClustering){

	// are we in the same cluster for even the smallest threshold?
	if(event->cluster[0] != event2.cluster[0] && event2.cluster[0] >= 0){

	  double ll = dMin(event, &event2);

	  if(fDebug){
	    std::cout << "pair " << eventInd << ", " << eventInd2 << ", initial ll = " << ll;
	  }
	  if(ll > llFitThreshold){
	    ll = dFit(event, &event2);
	    if(fDebug){
	      std::cout << ", fitted ll = " << ll;
	    }
	  }

	  // if(event->eventNumber==914599568){
	  //   std::cout << t << "\t" << event2.eventNumber << "\t" << ll << std::endl;
	  // }
	

	  // // is this event's nearest neighbour?
	  // if(ll < event1MinLL.at(t)){
	  //   event1MinLL.at(t) = ll;
	  //   event1MinLLEvent2EventNumber.at(t) = event2.eventNumber;
	  //   // std::cout << "new min in thread " << t << ",  ll = " << ll << std::endl;
	  // }

	  // is this event2's nearest neighbour?
	  // if(!mc && ll < event2.nearestNeighbourLogLikelihood){
	  //   event2.nearestNeighbourLogLikelihood = ll;
	  //   event2.nearestNeighbourEventNumber = event->eventNumber;
	  // }

	  if(ll <= llRange){
	    // nearbyEvents.push_back(eventInd2);
	    // nearbyEventLLs.push_back(ll);
	    nearbyEventsLoop.push_back(eventInd2);
	    nearbyEventLLsLoop.push_back(ll);

	    // std::cout << t << "\t" << eventInd2 << "\t" << ll << nearbyEventsLoop.size() << std::endl;
	  
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

  // need to copy the parallel stuff into resultant vector
  if(OpenMP::isEnabled){
    UInt_t numEvents = nearbyEvents.size();
    for(UInt_t t=0; t < nearbyEventsParallel.size(); t++){
      numEvents += nearbyEventsParallel.at(t).size();
    }
    nearbyEvents.reserve(numEvents);
    nearbyEventLLs.reserve(numEvents);

    for(UInt_t t=0; t < nearbyEventsParallel.size(); t++){
      nearbyEvents.insert(nearbyEvents.end(), nearbyEventsParallel.at(t).begin(), nearbyEventsParallel.at(t).end());
      nearbyEventLLs.insert(nearbyEventLLs.end(), nearbyEventLLsParallel.at(t).begin(), nearbyEventLLsParallel.at(t).end());
    }
  }

  // int minI = OpenMP::isEnabled ? TMath::LocMin(event1MinLL.size(), &event1MinLL[0]) : 0;
  // std::cout << "OpenMP::isEnabled = "  << OpenMP::isEnabled << ", minI = " <<  minI << std::endl;
  // event->nearestNeighbourLogLikelihood = event1MinLL.at(minI);
  // event->nearestNeighbourEventNumber = event1MinLLEvent2EventNumber.at(minI);

  // std::cout << "nearestNeighbourLogLikelihood = "  << event->nearestNeighbourLogLikelihood << ", nearestNeighbourEventNumber = " << event->nearestNeighbourEventNumber << std::endl;
  // std::cout << "after loop, eventInd = " << eventInd << ", nearbyEvents.size() = " << nearbyEvents.size() << ", nnEN ="
  // 	    << event->nearestNeighbourEventNumber << ", nnENLL = "
  // 	    << event->nearestNeighbourLogLikelihood << std::endl;

  
  // std::cout << "\n" << nearbyEvents.size() << "\t" << nearbyEventLLs.size() << std::endl;
  // for(unsigned i=0; i < nearbyEvents.size(); i++){
  //   std::cout << nearbyEvents[i] << "\t" << nearbyEventLLs[i] << std::endl;
  // }

  gErrorIgnoreLevel = fROOTgErrorIgnoreLevel;
}




void Acclaim::Clustering::LogLikelihoodMethod::makeAndWriteNSquaredEventEventHistograms(){


  fROOTgErrorIgnoreLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = 1001;
  
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
  double waisModelAlt = RampdemReader::BilinearInterpolatedSurfaceAboveGeoid(waisLon, waisLat);

  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
    const Event& event = events.at(eventInd);
    Double_t llWais = event.logLikelihoodFromPoint(waisLon, waisLat, waisModelAlt, true);
    hWais->Fill(llWais);
  }

  TRandom3 randy(123);
  std::vector<Int_t> event2Inds(events.size());
  for(UInt_t i=0; i < event2Inds.size(); i++){
    event2Inds.at(i) = randy.Uniform(0, events.size());
  }

  
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

    const int eventInd2 = 29695; //event2Inds.at(eventInd);
    const Event& event2 = events.at(eventInd2);
    grAnita->SetPoint(grAnita->GetN(), event1.anita.longitude, event1.anita.latitude);
    // if(event2.eventNumber!=fTestEvent2) continue;

    if(event1.eventNumber==fTestEvent1 && event2.eventNumber==fTestEvent2){
      std::cerr << "Info in " << __PRETTY_FUNCTION__ << " mapping parameter space around WAIS for event pair !"
		<< fTestEvent1 << " and " << fTestEvent2 << std::endl;
      AnitaVersion::set(3);
      double waisEasting, waisNorthing;
      RampdemReader::LonLatToEastingNorthing(AnitaLocations::getWaisLongitude(), AnitaLocations::getWaisLatitude(), waisEasting, waisNorthing);
      double delta = 700e3;
      const int nBins = 256;
      TH2D* hParams = new TH2D("hSingleEventTest", "Event-event fitted log likelihood; Easting (km); Northing (km); L_{sum}",
			       // nBins, -1090000, -1065000,
			       // nBins, -484200, -483600);
			       nBins, waisEasting-delta, waisEasting+delta,
			       nBins, waisNorthing-delta, waisNorthing+delta);

      fFitEvent1s.at(OpenMP::thread()) = &event1;
      fFitEvent2s.at(OpenMP::thread()) = &event2;
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

    double dist = dMin(&event1, &event2);
    Double_t angleBetweenEvents = TMath::RadToDeg()*anitaToEvent1.Angle(anitaToEvent2);

    
    hUnfit->Fill(angleBetweenEvents, dist);
    hUnfitSqrt->Fill(angleBetweenEvents, TMath::Sqrt(dist));

    double distFitted = dFit(&event1, &event2);
    hFit->Fill(angleBetweenEvents, distFitted);
    hFitSqrt->Fill(angleBetweenEvents, TMath::Sqrt(distFitted));

    // if(dist > 1000 && distFitted < 10){
    //   std::cout << dist << "\t" << distFitted << "\t" << angleBetweenEvents << "\t" << event1.eventNumber << "\t" << event2.eventNumber << std::endl;
    // }

    if(event1.theta > -5.5 && event2.theta > -5.5){
      std::cout << event1.eventNumber << "\t" << event2.eventNumber << std::endl;
    }
    // if(fabs(event1.thetaAdjustmentRequired) > 1e-2 && fabs(event2.thetaAdjustmentRequired) > 1e-2){
    //   std::cout << dist << "\t" << distFitted << "\t" << angleBetweenEvents << "\t" << event1.eventNumber << "\t" << event2.eventNumber << std::endl;
    // }
       

    double fitLon, fitLat;
    int t = OpenMP::thread();
    RampdemReader::EastingNorthingToLonLat(fFitEastings.at(t), fFitNorthings.at(t), fitLon, fitLat);
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

  gErrorIgnoreLevel = fROOTgErrorIgnoreLevel;
}





/**
 * An implementation of the event-to-event clustering with worst case O(N^2) efficiency
 */
void Acclaim::Clustering::LogLikelihoodMethod::doEventEventClustering(){

  double llFitThreshold = TMath::MinElement(llEventCuts.size(), &llEventCuts.at(0)); // fit if greater than this
  double llNearbyThreshold = TMath::MaxElement(llEventCuts.size(), &llEventCuts.at(0)); // ignore if greater than this

  std::cout << __PRETTY_FUNCTION__ << std::endl;
  std::cout << "llNearbyThreshold = " << llNearbyThreshold << ", llFitThreshold = " << llFitThreshold << std::endl;

  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
    std::cout << eventInd << "\t" << events.at(eventInd).cluster[0] << std::endl;
  }
  
  ProgressBar p(events.size());

  UInt_t numEventsProcessed = 0;
  Event* event1 = nextEvent();
  while(event1 && numEventsProcessed < events.size()){

    std::cout << event1->eventNumber << "\t" << event1->cluster[0] << std::endl;
    
    if(event1->eventEventClustering){

      // first look up events close by in Easting/Northing
      std::vector<Int_t> nearbyEventIndsEastingNorthing;
      double lookup[2] = {event1->easting, event1->northing};
      fKDTree->FindInRange(lookup, default_range_easting_northing, nearbyEventIndsEastingNorthing);	



      // move swapping global error level stuff out of parallelized loop
      // this used to be around the minimizer->minimize() 
      fROOTgErrorIgnoreLevel = gErrorIgnoreLevel;
      gErrorIgnoreLevel = 1001;



      // For storing the indices/and log likelihoods of the second events
      std::vector<Int_t> event2Inds;
      std::vector<Double_t> eventEventLLs;


      std::cerr << event2Inds.size() << "t" << eventEventLLs.size() << "\t" << nearbyEventIndsEastingNorthing.size() << std::endl;
      
      // For parallization
      int mt = OpenMP::getMaxThreads();
      std::vector<std::vector<Int_t> > event2IndsParallel(mt);
      std::vector<std::vector<Double_t> > eventEventLLsParallel(mt);
#pragma omp parallel for



      // Loop over nearby events in easting/northing
      for(UInt_t i=0; i < nearbyEventIndsEastingNorthing.size(); i++){
	int event2Ind = nearbyEventIndsEastingNorthing[i];


	int t = OpenMP::thread();

	// If not threaded, then put directly into event2Inds, otherwise parallel store.
	std::vector<Int_t>& nearbyEventsLoop = OpenMP::isEnabled ? event2IndsParallel.at(t) : event2Inds;
	std::vector<Double_t>& eventEventLLsLoop = OpenMP::isEnabled ? eventEventLLsParallel.at(t) : eventEventLLs;
	Event& event2 = events.at(event2Ind);


	// Here we try and save time..
	// 
	if(event1->eventNumber != event2.eventNumber && event2.eventEventClustering){
	  // are we in the same cluster for even the smallest threshold?
	  if(event1->cluster[0] != event2.cluster[0] || event2.cluster[0] >= 0){
	    double ll = dMin(event1, &event2);
	    if(ll > llFitThreshold){
	      ll = dFit(event1, &event2);
	    }
	    if(ll <= llNearbyThreshold){
	      nearbyEventsLoop.push_back(event2Ind);
	      eventEventLLsLoop.push_back(ll);
	    }
	  }
	}
      }

      // This block just combines the results of the threads...
      // blahBlahParallel -> blahBlah
      if(OpenMP::isEnabled){
	UInt_t numEvents = event2Inds.size();
	for(UInt_t t=0; t < event2IndsParallel.size(); t++){
	  numEvents += event2IndsParallel.at(t).size();
	}
	event2Inds.reserve(numEvents);
	eventEventLLs.reserve(numEvents);

	for(UInt_t t=0; t < event2IndsParallel.size(); t++){
	  event2Inds.insert(event2Inds.end(), event2IndsParallel.at(t).begin(), event2IndsParallel.at(t).end());
	  eventEventLLs.insert(eventEventLLs.end(), eventEventLLsParallel.at(t).begin(), eventEventLLsParallel.at(t).end());
	}
      }


      // Set fitter error level back to default
      gErrorIgnoreLevel = fROOTgErrorIgnoreLevel;



      // Now look at the clusters of the matched events (for each threshold)
      for(UInt_t z=0; z < llEventCuts.size(); z++){

	std::vector<Int_t> matchedClusters;
	if(event1->cluster[z] >= 0){
	  if(z==0){	    
	    std::cerr << "I think this is false for all events and this statement should never get printed..." << std::endl;
	    std::cerr << event1->eventNumber << "\t" << event1->cluster[0] << std::endl;
	  }
	  matchedClusters.push_back(event1->cluster[z]);
	}

	for(UInt_t i=0; i < event2Inds.size(); i++){
	  UInt_t event2Ind = event2Inds.at(i);
	  Event& event2 = events.at(event2Ind);

	  // within the cut value? and is the other event already associated with a cluster?
	  if(eventEventLLs.at(i) <= llEventCuts.at(z) && event2.cluster[z] >= 0){
	    // and if that other cluster isn't already in the list, then add it
	    if(std::find(matchedClusters.begin(), matchedClusters.end(), event2.cluster[z])==matchedClusters.end()){
	      matchedClusters.push_back(event2.cluster[z]);
	    }
	  }
	}

	// If none of the events are in a cluster, then add a new cluster!
	if(matchedClusters.size()==0){

	  // make a new cluster, it's intial position is this event location
	  // although that will be changed...
	  Cluster nc(*event1, clusters.at(z).size());
	  nc.llEventCutInd = z;
	  nc.llEventCut = llEventCuts.at(z);
	  clusters.at(z).push_back(nc);
	  matchedClusters.push_back(nc.index);
	}

	// If the list of matched clusters is greater than one
	// then we're going to merge clusters!
	// But what to call the new cluster?
	// We'll label it with the smallest index...
	Int_t thisCluster = TMath::MinElement(matchedClusters.size(), &matchedClusters[0]);


	// Mark event1 as in the minCluster
	Int_t oldCluster1Ind = event1->cluster[z];
	if(oldCluster1Ind >= 0){
	  if(z==0){
	    std::cerr << "Do we ever get here?" << std::endl;
	  }
	  clusters.at(z).at(oldCluster1Ind).numDataEvents--;
	}
	event1->cluster[z] = thisCluster;
	clusters.at(z).at(thisCluster).numDataEvents++;

	std::cerr << event2Inds.size() << std::endl;
	
	// Now reassign all events in the matched clusters to thisCluster
	for(UInt_t i=0; i < event2Inds.size(); i++){

	  // The only condition is that you're within the cut, close to this event
	  if(eventEventLLs.at(i) <= llEventCuts.at(z)){

	    Int_t event2Ind = event2Inds.at(i);
	    Event& event2 = events.at(event2Ind);
	    Int_t oldCluster2Ind = event2.cluster[z];
	    if(oldCluster2Ind!=thisCluster){
	      if(oldCluster2Ind >= 0){
		clusters.at(z).at(oldCluster2Ind).numDataEvents--;
	      }
	      event2.cluster[z] = thisCluster;
	      clusters.at(z).at(thisCluster).numDataEvents++;
	    }
	  }
	}
      }
    }

    numEventsProcessed = 0;
    for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
      if(events.at(eventInd).cluster[0] > -1){
	numEventsProcessed++;
      }
    }
    
    p.inc(numEventsProcessed);
    numEventsProcessed++;
    event1 = nextEvent();
  }
}






/**
 * An implementation of the event-to-event clustering with worst case O(N^2) efficiency
 */
void Acclaim::Clustering::LogLikelihoodMethod::doMcEventClustering(){

  double llFitThreshold = TMath::MinElement(llEventCuts.size(), &llEventCuts.at(0)); // fit if greater than this
  double llNearbyThreshold = TMath::MaxElement(llEventCuts.size(), &llEventCuts.at(0)); // ignore if greater than this

  std::cout << __PRETTY_FUNCTION__ << std::endl;
  std::cout << "llNearbyThreshold = " << llNearbyThreshold << ", llFitThreshold = " << llFitThreshold << std::endl;

  ProgressBar p(mcEvents.size());
  for(UInt_t event1Ind=0; event1Ind < mcEvents.size(); event1Ind++){
    McEvent& event1 = mcEvents.at(event1Ind);
    if(event1.eventEventClustering){

      // look up nearby events
      std::vector<Int_t> event2Inds;
      std::vector<Double_t> eventEventLLs;
      nearbyEvents(event1Ind, event2Inds, eventEventLLs, true, llNearbyThreshold, llFitThreshold);

// #pragma omp parallel for
    //   for(UInt_t z=0; z < llEventCuts.size(); z++){
    // 	event1.cluster[z] = -1;
    // 	for(UInt_t i=0; i < event2Inds.size(); i++){
    // 	  if(eventEventLLs.at(i) <= event1.nearestNeighbourLogLikelihood &&
    // 	     eventEventLLs.at(i) <= llEventCuts.at(z)){
    // 	    Int_t event2Ind = event2Inds.at(i);
    // 	    const Event& event2 = events.at(event2Ind);
    // 	    event1.cluster[z] = event2.cluster[z]; // put in same cluster
    // 	    clusters.at(z).at(event1.cluster[z]).sumMcWeights += event1.weight;
    // 	  }
    // 	}
    //   }
    }
    p.inc(event1Ind);
  }
}





Int_t Acclaim::Clustering::LogLikelihoodMethod::removeLargeBasesNearMcMurdo(){

  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
    Event& event = events.at(eventInd);
    if(event.nearestKnownBaseLogLikelihood < llClusterCut){
      int clusterInd = event.nearestKnownBaseCluster;
      for(UInt_t z=0; z < clusters.size(); z++){
	event.cluster[z] = clusterInd;
	clusters.at(z).at(clusterInd).numDataEvents++;
      }
    }
  }

  const int manyEvents = 1000;
  int numRemoved=0;
  for(UInt_t clusterInd=0; clusterInd < clusters.at(0).size(); clusterInd++){
    Cluster& cluster = clusters.at(0).at(clusterInd);

    // remove these guys!
    if(cluster.numDataEvents >= manyEvents && isVaguelyNearMcMurdo(cluster)){
      for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
	Event& event = events.at(eventInd);
	if(event.cluster[0]==(int)clusterInd){
	  events.at(eventInd).eventEventClustering = false;
	  numRemoved++;
	}
      }
    }
    else{
      for(UInt_t z=0; z < llEventCuts.size(); z++){
	clusters.at(z).at(clusterInd).resetClusteringNumbers();
      }
      for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
	Event& event = events.at(eventInd);
	if(event.cluster[0]==(int)clusterInd){
	  events.at(eventInd).resetClusteringNumbers();
	}
      }
    }
  }
  return numRemoved;
}




void Acclaim::Clustering::LogLikelihoodMethod::setInitialBaseClusters(){

  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
    Event& event = events.at(eventInd);
    Int_t clusterInd = event.nearestKnownBaseCluster;

    if(clusterInd >= 0){
      for(int z=0; z < event.nThresholds; z++){
	Cluster& cluster = clusters.at(z).at(clusterInd);
	if(event.nearestKnownBaseLogLikelihood < cluster.llEventCut){
	  cluster.numDataEvents++;
	  event.cluster[z] = cluster.index;
	  // std::cout << eventInd << "\t" << z << "\t" << cluster.index << "\t" << event.nearestKnownBaseLogLikelihood << "\t" << cluster.llEventCut << "\t" << cluster.numDataEvents << std::endl;
	}
      }
    }
  }
}


void Acclaim::Clustering::LogLikelihoodMethod::doMcBaseClustering(){

  std::cout << __PRETTY_FUNCTION__ << std::endl;

  ProgressBar p(mcEvents.size());
  for(UInt_t eventInd=0; eventInd < mcEvents.size(); eventInd++){
    McEvent* mcEvent = &mcEvents.at(eventInd);

    for(int clusterInd=0; clusterInd < (int)clusters.at(0).size(); clusterInd++){
      Cluster& cluster = clusters.at(0).at(clusterInd);
      if(cluster.knownBase){
	double distM = mcEvent->usefulPat.getDistanceFromSource(cluster.latitude, cluster.longitude, cluster.latitude);
	if(distM < fFitHorizonDistM){
	  double ll = mcEvent->logLikelihoodFromPoint(cluster.longitude, cluster.latitude, cluster.altitude, true);
	  for(int z=0; z < mcEvent->nThresholds; z++){
	    if(mcEvent->cluster[z] < 0){
	      if(ll < llEventCuts.at(z)){
		clusters[z][clusterInd].sumMcWeights += mcEvent->weight;
		mcEvent->cluster[z] += clusterInd;
	      }
	    }
	  }
	}
      }
    }
    p.inc(eventInd);
  }
}


void Acclaim::Clustering::LogLikelihoodMethod::doClustering(const char* dataGlob, const char* mcGlob, const char* outFileName){

  readInSummaries(dataGlob);
  readInSummaries(mcGlob);

  // initializeBaseList();
  // forEachEventFindClosestKnownBase();

  // bool fRemoveLargeBasesNearMcMurdo = true;
  // if(fRemoveLargeBasesNearMcMurdo){
  //   int numRemoved = removeLargeBasesNearMcMurdo();
  //   std::cout << "Removed " << numRemoved << " events from near McMurdo!" << std::endl;
  //   std::cout << (int)events.size() - numRemoved << " events remain for pairwise clustering!" << std::endl;
  // }

  // setInitialBaseClusters();
  initKDTree();
  doEventEventClustering();
  // doMcBaseClustering();
  // doMcEventClustering();

  const char* fakeArgv[1] = {outFileName};
  OutputConvention oc(1, const_cast<char**>(fakeArgv));
  TFile* fOut = oc.makeFile();

  // for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
  //   Event& event = events.at(eventInd);
  //   event.antarcticaHistBin = hEvents->Fill(event.longitude, event.latitude);
  // }
  // hEvents->Write();

  for(UInt_t eventInd=0; eventInd < mcEvents.size(); eventInd++){
    McEvent& mcEvent = mcEvents.at(eventInd);
    mcEvent.antarcticaHistBin = hMcEvents->Fill(mcEvent.longitude, mcEvent.latitude, mcEvent.weight);
  }
  hMcEvents->Write();
  
  // makeAndWriteNSquaredEventEventHistograms();

  writeAllGraphsAndHists();
  makeSummaryTrees();
  fOut->Write();
  fOut->Close();

  return;


}







