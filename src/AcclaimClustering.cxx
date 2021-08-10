#include "AcclaimClustering.h"
#include "AnitaGeomTool.h"
#include "SummarySet.h"
#include "AnitaDataset.h"
#include "OutputConvention.h"
#include "AnitaEventSummary.h"
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
#include "RootTools.h"
#include "DrawStrings.h"
#include "ThermalChain.h"
#include "Hical2.h"
#include "FFTtools.h"
#include "TBits.h"
#include <random>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"


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
  const int n = 6;
  const double phiParams[n]   = {-2.50414e-01,  3.02406e-01, 2.43376e-01, 5.09057,  8.01369e-01, 1.}; //A4 is the second set of 3 numbers, A3 is the first set
  const double thetaParams[n] = {-3.83773e-01, -3.00964e-01, 1.64537e-01, 1.34307, 7.09382e-01, 1.}; //A4 is the second set of 3 numbers, A3 is the first set
  TString formula = (AnitaVersion::get() == 3) ? "exp([0]*x + [1]) + [2]" : "[0]/(pow(x,[1]) + [2])";

}


/**
 * @namespace VarianceModel
 * @brief Parameters defining the variance model, where spherical curvature in the interferometric map is considered
 *
 * Derivation of these numbers is analogous to what is seen in the macro plotCalPulserResolution.C, except the formula has closer to do with the square of the argument
 */
namespace VarianceModel{
  const int n = 4;
  //  For now, the following numbers have to do with A4 only. First set of three are for coherent filtered SNR, second set are for deconvolved filtered SNR.
  const double phiParams[n]   = {8.34,  1.203, 13.54, 1.273};
  const double thetaParams[n] = {0.7114, 1.049, 1.472, 1.192};
  TString formula = "[0] / x^[1]";
//  const int n = 6;
//  const double phiParams[n]   = {25.19,  2.42, 2.876, 5.361, 3.674, 2.869};
//  const double thetaParams[n] = {9.14, 2.556, 0.1189, 2.685, 3.501, 0.1174};
//  TString formula = "[0] * TMath::Gaus(x, 0, [1]) + [2]";
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
  TString formula = "[0]/(pow(x,[1]) + [2])";
  sigma_phi = (AnitaVersion::get() == 3) ? exp(ResolutionModel::phiParams[0]*x + ResolutionModel::phiParams[1]) + ResolutionModel::phiParams[2] : ResolutionModel::phiParams[3]/(pow(x, ResolutionModel::phiParams[4]) + ResolutionModel::phiParams[5]);
  sigma_theta = (AnitaVersion::get() == 3) ? exp(ResolutionModel::thetaParams[0]*x + ResolutionModel::thetaParams[1]) + ResolutionModel::thetaParams[2] : ResolutionModel::thetaParams[3]/(pow(x, ResolutionModel::thetaParams[4]) + ResolutionModel::thetaParams[5]);
}


TCanvas* Acclaim::Clustering::drawAngularResolutionModel(double maxSnr){
  TCanvas* c1 = new TCanvas();

  TF1* fTheta = new TF1("fThetaResolutionModel", ResolutionModel::formula, 0, maxSnr);
  TF1* fPhi = new TF1("fPhiResolutionModel", ResolutionModel::formula, 0, maxSnr);
  int versionOffset = (AnitaVersion::get() == 3) ? 0 : 3;
  for(int i=0; i < ResolutionModel::n; i++){
    fTheta->SetParameter(i, ResolutionModel::thetaParams[i+versionOffset]);
    fPhi->SetParameter(i, ResolutionModel::phiParams[i+versionOffset]);
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
 * @brief Wrapper function to calculate the angular variance for clustering
 *
 * @param sum is the AnitaEventSummary
 * @param pol the polarisation of interest
 * @param peakInd the peak of the map of interest
 * @param var_theta the calculated theta variance (degrees)
 * @param var_phi the calculated phi variance (degrees)
 */
void Acclaim::Clustering::getAngularVariance(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd, double& var_theta, double& var_phi){
  const double x = sum -> coherent_filtered[pol][peakInd].snr;
  getAngularVariance(x, var_theta, var_phi);
}


/**
 * @brief Calculate the angular variance for clustering
 * @todo Currently derived from WAIS pulses, but should probably be from MC
 *
 * @param x the parameterization variable
 * @param var_theta the calculated theta variance (degrees)
 * @param var_phi the calculated phi variance (degrees)
 */
void Acclaim::Clustering::getAngularVariance(double x, double & var_theta, double & var_phi){
  var_phi = VarianceModel::phiParams[0] / pow(x, VarianceModel::phiParams[1]);
  var_theta = VarianceModel::thetaParams[0] / pow(x, VarianceModel::thetaParams[1]);
//  TString formula = "[0] * exp([1] * x) + [2]";
//  var_phi = VarianceModel::phiParams[0] * TMath::Gaus(x, 0, VarianceModel::phiParams[1]) + VarianceModel::phiParams[2];
//  var_theta = VarianceModel::thetaParams[0] * TMath::Gaus(x, 0, VarianceModel::thetaParams[1]) + VarianceModel::thetaParams[2];
}


TCanvas* Acclaim::Clustering::drawAngularVarianceModel(double maxSnr){
  TCanvas* c1 = new TCanvas();

  TF1* fTheta = new TF1("fThetaVarianceModel", VarianceModel::formula, 0, maxSnr);
  TF1* fPhi = new TF1("fPhiVarianceModel", VarianceModel::formula, 0, maxSnr);
  for(int i=0; i < VarianceModel::n; i++){
    fTheta->SetParameter(i, VarianceModel::thetaParams[i]);
    fPhi->SetParameter(i, VarianceModel::phiParams[i]);
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
 * 
 * @param recalculateNow a boolian for whether or not to fire traceBackToContinent when the usefulPat is initialized
 */
void Acclaim::Clustering::Event::setupUsefulPat(bool calculateSource){
  Adu5Pat pat = anita.pat();
  usefulPat = UsefulAdu5Pat(&pat);
  usefulPat.setInterpSurfaceAboveGeoid(true);
  // usefulPat.setSurfaceCloseEnoughInter(1e-3);
  usefulPat.setSurfaceCloseEnoughInter(1);
  usefulPat.setMaxLoopIterations(5000); // make this arbitrarily large since it only happens once
  // const double maxThetaAdjust = 8*TMath::DegToRad();

  if(calculateSource){
    //  if(eventNumber==10047816){ //9887706
    //   usefulPat.setDebug(true);
    // }

    usefulPat.traceBackToContinent3(phi*TMath::DegToRad(), -theta*TMath::DegToRad(), &longitude, &latitude, &altitude, &thetaAdjustmentRequired);//, maxThetaAdjust, 10);

    if(altitude < -999){
      usefulPat.setDebug(true);
      usefulPat.traceBackToContinent3(phi*TMath::DegToRad(), -theta*TMath::DegToRad(), &longitude, &latitude, &altitude, &thetaAdjustmentRequired);//, maxThetaAdjust, 10);
    }

    // std::cout << eventNumber << "\t" << longitude << "\t" << latitude << "\t" << altitude << "\t" << thetaAdjustmentRequired << std::endl;

    // if(usefulPat.getDebug()){
    //   exit(1);
    // }

    RampdemReader::LonLatToEastingNorthing(longitude, latitude, easting, northing);
    AnitaGeomTool* geom = AnitaGeomTool::Instance();
    geom->getCartesianCoords(latitude, longitude, altitude, centre);

    // std::cerr << eventNumber << "\t" << easting << "\t" << northing << std::endl;
    if(latitude < -90){
      std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", for eventNumber " << eventNumber << "\n";
      // std::cerr << "Doing traceBackToContinenet again in debug mode!\n";
      // usefulPat.setDebug(true);
      // usefulPat.traceBackToContinent(phi*TMath::DegToRad(), -theta*TMath::DegToRad(), &longitude, &latitude, &altitude, &thetaAdjustmentRequired, maxThetaAdjust, 100);
      // usefulPat.setDebug(false);
      // RampdemReader::LonLatToEastingNorthing(longitude, latitude, easting, northing);
      std::cerr << "Removing event " << eventNumber << " from event-event clustering!" << std::endl;
      eventEventClustering = false;
      easting = -99e20;
      northing = -99e20;
    }
    else{
      eventEventClustering = true;
    }

    // selfLogLikelihood = logLikelihoodFromPoint(longitude, latitude, altitude, false);
    // selfLogLikelihood = logLikelihoodFromPoint(longitude, latitude, altitude, true);
    selfLogLikelihood = logLikelihoodFromPoint(longitude, latitude, altitude, false);
  }
}


/**
 * Evaluate the log-likelihood distance from an arbitrary point relative to any event
 *
 * @param sourceLon is the source longitude
 * @param sourceLat is the source latitude
 * @param sourceAlt is the source altitude
 * @param addOverHorizonPenalty is an option to add a penalty term in LL for (approximately) being over the horizon
 *
 * @return the log-likelihood
 */
Double_t Acclaim::Clustering::Event::logLikelihoodFromPoint(Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt, bool addOverHorizonPenalty) const {

  Double_t phiSource, thetaSource, thetaMean;
  usefulPat.getThetaAndPhiWave2(sourceLon, sourceLat, sourceAlt, thetaSource, phiSource);
  phiSource = TMath::RadToDeg() * phiSource;
  thetaSource = -1 * TMath::RadToDeg() * thetaSource;
  thetaMean = (theta + thetaSource) / 2;

  Double_t dPhi = -1 * Acclaim::RootTools::getDeltaAngleDeg(phi, phiSource) * cos(TMath::DegToRad() * thetaMean) / sqrt(varPhi);
  //  Factor of -1 in front is ignorable for our purposes,
  //  but is there to drive home the ANITA angle convention in the geometric delays of our interferometic maps.
  //  dPhi originally weighted by cos(theta) as opposed to cos(thetaMean), but I think the article
  //  https://en.wikipedia.org/wiki/Geographical_distance#Spherical_Earth_projected_to_a_plane is on to something.
  Double_t dTheta = (theta - thetaSource) / sqrt(varTheta);

  Double_t ll = dPhi * dPhi + dTheta * dTheta;

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
      std::cerr << __PRETTY_FUNCTION__ << " for " << eventNumber << ", we are " << distM / 1000 << "km from the source, after horizon penalty, ll = " << ll << std::endl;
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
  // RampdemReader::LonLatToEastingNorthing(longitude, latitude, easting, northing);
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  geom->getCartesianCoords(latitude, longitude, altitude, centre);
  theta = peak.theta;
  phi = peak.phi;
  thetaAdjustmentRequired = peak.theta_adjustment_needed;
  selfLogLikelihood = -9999;
  anita = sum->anitaLocation;
  getAngularResolution(sum, pol, peakInd, sigmaTheta, sigmaPhi);
  getAngularVariance(sum, pol, peakInd, varTheta, varPhi);
  antarcticaHistBin = -1;
  fDebug = false;
  nearestKnownBaseLogLikelihood = DBL_MAX;
  nearestKnownBaseCluster = -1;

  setNThresholds(nT);
  resetClusteringNumbers();
  setupUsefulPat();

  // std::cout << eventNumber << "\t" << nThresholds << "\t" << cluster[0] << std::endl;
}









Acclaim::Clustering::Event::Event(int pol, int peakInd, double peak_phi, double peak_theta, int nT, UInt_t eventNumber, Int_t run,
    double anita_longitude, double anita_latitude, double anita_altitude, double anita_heading, double coherent_filtered_snr, double deconvolved_filtered_snr)// ,
// double longitude, double latitude, double altitude)
: nThresholds(0), cluster(NULL),
  dThetaCluster(NULL), dPhiCluster(NULL)
{

  this->eventNumber = eventNumber;
  this->run = run;
  this->pol = (AnitaPol::AnitaPol_t) pol;
  peakIndex = peakInd;
  // AnitaGeomTool* geom = AnitaGeomTool::Instance();
  // geom->getCartesianCoords(latitude, longitude, altitude, centre);
  theta = peak_theta;
  phi = peak_phi;
  // thetaAdjustmentRequired = peak.theta_adjustment_needed;
  selfLogLikelihood = -9999;
  anita.longitude = anita_longitude;
  anita.latitude = anita_latitude;
  anita.altitude = anita_altitude;
  anita.heading = anita_heading;

  getAngularResolution(coherent_filtered_snr, sigmaTheta, sigmaPhi);

  getAngularVariance(coherent_filtered_snr, varTheta, varPhi);

  // getAngularResolution(sum, pol, peakInd, sigmaTheta, sigmaPhi);
  antarcticaHistBin = -1;
  fDebug = false;
  nearestKnownBaseLogLikelihood = DBL_MAX;
  nearestKnownBaseCluster = -1;

  setNThresholds(nT);
  resetClusteringNumbers();
  setupUsefulPat();

  // getAngularResolution(sum, pol, peakInd, sigmaTheta, sigmaPhi);  

  // std::cout << eventNumber << "\t" << nThresholds << "\t" << cluster[0] << std::endl;
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
  varTheta = default_var_theta;
  varPhi = default_var_phi;
  for(int dim=0; dim < nDim; dim++){
    centre[dim] = 0;
  }
  antarcticaHistBin = -1;
  fDebug = false;


  nearestKnownBaseLogLikelihood = DBL_MAX;
  nearestKnownBaseCluster = -1;
  nearestKnownBaseSurfaceSeparationKm = DBL_MAX;
  nearestEventSurfaceDistanceKm = DBL_MAX;
  nearestEventSurfaceEventNumber = 0;
  nearestEventSurfaceLogLikelihood = DBL_MAX;



  setNThresholds(nT);
  resetClusteringNumbers();  
}


Acclaim::Clustering::Event::~Event(){
  deleteArrays();
}



Acclaim::Clustering::Event& Acclaim::Clustering::Event::operator=(const Event& event){

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

  thetaAdjustmentRequired = event.thetaAdjustmentRequired;

  sigmaTheta = event.sigmaTheta;
  sigmaPhi = event.sigmaPhi;

  varTheta = event.varTheta;
  varPhi = event.varPhi;

  if(nThresholds!=event.nThresholds){
    deleteArrays();
    nThresholds = event.nThresholds;
    cluster = new Int_t[nThresholds];
    dThetaCluster = new Double_t[nThresholds];
    dPhiCluster = new Double_t[nThresholds];
  }
  for(int z=0; z < nThresholds; z++){
    cluster[z] = event.cluster[z];
    dThetaCluster[z] = event.dThetaCluster[z];
    dPhiCluster[z] = event.dPhiCluster[z];
  }


  eventEventClustering = event.eventEventClustering;
  nearestKnownBaseLogLikelihood = event.nearestKnownBaseLogLikelihood;
  nearestKnownBaseSurfaceSeparationKm = event.nearestKnownBaseSurfaceSeparationKm;
  nearestKnownBaseCluster = event.nearestKnownBaseCluster;
  selfLogLikelihood = event.selfLogLikelihood;
  nearestEventSurfaceDistanceKm = event.nearestEventSurfaceDistanceKm;
  nearestEventSurfaceEventNumber = event.nearestEventSurfaceEventNumber;
  nearestEventSurfaceLogLikelihood = event.nearestEventSurfaceLogLikelihood;  

  antarcticaHistBin = event.antarcticaHistBin;
  usefulPat = event.usefulPat;
  fDebug = event.fDebug;

  return *this;
}


  Acclaim::Clustering::Event::Event(const Event& event)
: nThresholds(0), cluster(NULL),
  dThetaCluster(NULL), dPhiCluster(NULL)
{
  *this = event;
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
  nearestKnownBaseLogLikelihood = DBL_MAX;
  nearestKnownBaseSurfaceSeparationKm = DBL_MAX;
  nearestKnownBaseCluster = -1;
  selfLogLikelihood = 0;
  nearestEventSurfaceDistanceKm = DBL_MAX;
  nearestEventSurfaceEventNumber = 0;
  nearestEventSurfaceLogLikelihood = DBL_MAX;

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


Acclaim::Clustering::McEvent::McEvent(double weight, double energy, int pol, int peakInd, double peak_phi, double peak_theta, int nT, UInt_t eventNumber, Int_t run,
    double anita_longitude, double anita_latitude, double anita_altitude, double anita_heading, double coherent_filtered_snr, double deconvolved_filtered_snr)
:  Event(pol, peakInd, peak_phi, peak_theta, nT, eventNumber, run,
    anita_longitude, anita_latitude, anita_altitude, anita_heading, coherent_filtered_snr, deconvolved_filtered_snr)
{
  this->weight = weight;
  this->energy = energy;
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
: numMcDivisions(100), fEventsAlreadyClustered(false), fMyBackground(),
  fROOTgErrorIgnoreLevel(gErrorIgnoreLevel), fDrawNewNearbyEventsHistograms(true),
  fReadInBaseList(false), fStoreUnclusteredHistograms(true)

{
  const char* sgeTaskId = getenv("SGE_TASK_ID");
  if(!sgeTaskId){
    mcDivision = -1; // this means to all of it, ignoring the value of numMcDivisions 
  }
  else{
    mcDivision = atoi(sgeTaskId);
    std::cout << "Info in " << __PRETTY_FUNCTION__ << ", found SGE_TASK_ID=" << sgeTaskId
      << ", setting mcDivision = " << mcDivision << std::endl;
  }
  // fMyBackground.SetCoarseness(100);
  // fMyBackground.SetCoarseness(60);
  // fMyBackground.SetCoarseness(20);  // lower coarseness implies more bins...
  fMyBackground.SetCoarseness(120);  // lower coarseness implies more bins...
  std::cout << fMyBackground.GetXaxis()->GetBinWidth(1) << "\t" << fMyBackground.GetYaxis()->GetBinWidth(1) << std::endl;

  grTestMinimizerWalk = NULL;
  grTestMinimizerValue = NULL;
  
  for (int i = 0; i <= 60; ++i) llEventCuts.push_back(pow(10, 0.05 * i));

/*

  for(Int_t i = 1; i <= 20; ++i) {
  
    llEventCuts.push_back(i * i);
    llEventCuts.push_back(i * (i + 1));
    llEventCuts.push_back(i * (i + 2));
  }
  
  surfaceDistThresholdKm = 30;
  llEventCuts.push_back(1);
  llEventCuts.push_back(2);
  llEventCuts.push_back(4);
  llEventCuts.push_back(6);
  llEventCuts.push_back(8);

  llEventCuts.push_back(10);
  llEventCuts.push_back(12);
  llEventCuts.push_back(15);
  llEventCuts.push_back(20);
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
*/

  for(UInt_t z=0; z < llEventCuts.size(); z++){
    clusters.push_back(std::vector<Cluster>());
  }

  // both above horizon; good improvement; actual test
  fTestEvent1 = 61284883; //55510391; // 61156660; //61033430;
  fTestEvent2 = 55352006; //61338514; //55789194; //61424151;
  

  fKDTree = NULL;
  fDebug = false;
  fUseBaseList = true;
  fPermyriadOfMC = 0;
  fNumOfMC = 0;
  fCut = 0;
  fCutHical = 0;
  fSelfLLMax = -1;
  fEntryList = 0;

  fMaxFitterAttempts = 1;

  tr3 = new TRandom3(0);

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
  return TMath::Min(dAsym(event1, event2) + event2->selfLogLikelihood, dAsym(event2, event1) + event1->selfLogLikelihood);
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






bool Acclaim::Clustering::LogLikelihoodMethod::considerBin(const Event& event, Int_t bx, Int_t by, double& easting, double& northing){

  const double halfBinWidthNorthing = 0.5*fMyBackground.GetYaxis()->GetBinWidth(bx);
  const double halfBinWidthEasting = 0.5*fMyBackground.GetXaxis()->GetBinWidth(by);

  northing = fMyBackground.GetXaxis()->GetBinCenter(by);
  northing = northing > event.northing ? northing - halfBinWidthNorthing : northing + halfBinWidthNorthing;
  double dNorthing = northing - event.northing;
  
  easting = fMyBackground.GetXaxis()->GetBinCenter(bx);
  easting = easting > event.easting ? easting - halfBinWidthEasting : easting + halfBinWidthEasting;
  double dEasting = easting - event.easting;

  const double distSq = dNorthing * dNorthing + dEasting * dEasting;
  const double maxRangeSq = default_horizon_distance * default_horizon_distance;

  if(distSq < maxRangeSq){
    return true;
  }
  return false;
}



// void Acclaim::Clustering::LogLikelihoodMethod::nearbyEvents2(UInt_t eventInd, std::set<UInt_t>& nearbyEventInds){
void Acclaim::Clustering::LogLikelihoodMethod::nearbyEvents2(UInt_t eventInd, std::vector<UInt_t>& nearbyEventInds){  

  nearbyEventInds.clear();

  const int nx = fMyBackground.GetNbinsX();
  const int ny = fMyBackground.GetNbinsY();

  const Event& event = events.at(eventInd);
  std::vector<bool> found(events.size(), 0);


  for(Int_t by=1; by <= ny; by++){
    for(Int_t bx=1; bx <= nx; bx++){
      double easting, northing;
      if(considerBin(event, bx, by, easting, northing)){

        const std::vector<UInt_t>& eventsThisBin = fLookupEN[by-1][bx-1];

        for(UInt_t i=0; i < eventsThisBin.size(); i++){
          if(eventsThisBin[i]!=eventInd){ // (that's not the same event!)
            // if(RootTools::vectorContainsValue(eventsThisBin, eventInd)){
            // nearbyEventInds.insert(eventsThisBin[i]);
            found[eventsThisBin[i]] = 1;
            // nearbyEventInds.insert(eventsThisBin[i]);
            // }
          }
        }
      }
    }
  }

  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
    if(found[eventInd] > 0){
      nearbyEventInds.push_back(eventInd);
    }
  }


}






Double_t Acclaim::Clustering::LogLikelihoodMethod::getAngDistSqEventCluster(const Event& event, const Cluster& cluster){

  Double_t phiWave, thetaWave;
  event.usefulPat.getThetaAndPhiWave2(cluster.longitude, cluster.latitude, cluster.altitude, thetaWave, phiWave);
  Double_t phiDeg = TMath::RadToDeg() * phiWave;
  Double_t thetaDeg = -1 * TMath::RadToDeg() * thetaWave;
  Double_t thetaMean = (thetaDeg + event.theta) / 2;

  Double_t deltaPhiDeg, deltaThetaDeg;
  getDeltaThetaDegDeltaPhiDegEventCluster(event, cluster, deltaThetaDeg, deltaPhiDeg); 

  Double_t dPhiNorm = -1 * deltaPhiDeg * cos(TMath::DegToRad() * thetaMean) / sqrt(event.varPhi);
  //  Factor of -1 in front is ignorable for our purposes, but is there to drive home the ANITA angle convention
  //  in the geometric delays of our interferometic maps.
  //  dPhi originally weighted by cos(theta) as opposed to cos(thetaMean), but I think the article
  //  https://en.wikipedia.org/wiki/Geographical_distance#Spherical_Earth_projected_to_a_plane is on to something.
  Double_t dThetaNorm = deltaThetaDeg / sqrt(event.varTheta);

  Double_t angSq =  dThetaNorm * dThetaNorm + dPhiNorm * dPhiNorm;

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

  if(!fStoreUnclusteredHistograms){
    while(hUnclusteredEvents.size() > 0){
      delete hUnclusteredEvents.back();
      hUnclusteredEvents.pop_back();
    }
  }



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
      if(event.antarcticaHistBin < 0){
        std::cerr << event.eventNumber << "\t" << event.longitude << "\t" << event.latitude << "\t" << event.easting << "\t" << event.northing << std::endl;
      }
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
        Double_t dN = event.northing - meanNorthing;        
        Double_t dE = event.easting - meanEasting;

        Double_t surfaceSeparationSquared = dN * dN + dE * dE;

        if(surfaceSeparationSquared < bestSurfaceSeparationSquared){
          bestSurfaceSeparationSquared = surfaceSeparationSquared;
          nextEvent = &event;
        }
      }
    }
  }
  return nextEvent;
}




void Acclaim::Clustering::LogLikelihoodMethod::doBaseEventClustering(){

  std::cout << __PRETTY_FUNCTION__ << std::endl;

  const int nBases = BaseList::getNumBases();

  const Long64_t nEvents = events.size();
  ProgressBar p(nEvents);

  // I like to call this vector constructor bullshit
  std::vector< std::vector<std::vector<Int_t> > > matchedClusters(llEventCuts.size(), std::vector< std::vector<Int_t > >(nBases, std::vector<Int_t>()));

  for(Long64_t eventInd=0; eventInd < nEvents; eventInd++){
    Event* event = &events.at(eventInd);

    std::vector<std::vector<Int_t> > matchedClustersThisEvent(llEventCuts.size(), std::vector<Int_t>());
    for(int clusterInd=0; clusterInd < nBases; clusterInd++){
      Cluster& cluster = clusters.at(0).at(clusterInd);
      if(cluster.knownBase){
        double distM = event->usefulPat.getDistanceFromSource(cluster.latitude, cluster.longitude, cluster.latitude);
        if(distM < default_horizon_distance){
          double ll = event->logLikelihoodFromPoint(cluster.longitude, cluster.latitude, cluster.altitude, true);
          double surfaceSeparationKm = 1e-3*event->cartesianSeparation(cluster);

          if(surfaceSeparationKm < surfaceDistThresholdKm){ // then true for all cluster sizes
            for(int z=0; z < llEventCuts.size(); z++){
              matchedClustersThisEvent[z].push_back(clusterInd);
            }
          }
          else{
            for(int z=0; z < llEventCuts.size(); z++){
              if(ll < llEventCuts.at(z)){
                matchedClustersThisEvent[z].push_back(clusterInd);
              }
            }
          }

          if(ll < event->nearestKnownBaseLogLikelihood){
            event->nearestKnownBaseLogLikelihood = ll;
            event->nearestKnownBaseCluster = clusterInd;
          }
          if(surfaceSeparationKm < event->nearestKnownBaseSurfaceSeparationKm){
            event->nearestKnownBaseSurfaceSeparationKm = surfaceSeparationKm;
            event->nearestKnownBaseClusterSurface = clusterInd;
          }
        }
      }
      else{
        std::cerr << "You shouldn't get here!!!" << std::endl;
      }
    }
    for(int z=0; z < llEventCuts.size(); z++){
      if(matchedClustersThisEvent[z].size() > 0){

        // for all matched clusters
        for(int i = 0; i < matchedClustersThisEvent[z].size(); i++){
          Int_t matchedCluster = matchedClustersThisEvent[z][i];

          // add other matched clusters to their list...
          for(int j = 0; j < matchedClustersThisEvent[z].size(); j++){
            Int_t matchedCluster2 = matchedClustersThisEvent[z][j];

            if(!RootTools::vectorContainsValue(matchedClusters[z][matchedCluster], matchedCluster2)){
              matchedClusters[z][matchedCluster].push_back(matchedCluster2);
            }
          }
        }	
      }
    }

    p.inc(eventInd);
  }


  for(int z=0; z < llEventCuts.size(); z++){
    std::vector<Int_t> reassignedTo(0);
    reassignedTo.reserve(nBases);
    for(int i=0; i < nBases; i++){
      reassignedTo.push_back(i);
    }

    for(int b=0; b < nBases; b++){
      if(matchedClusters[z][b].size() > 0){

        // here we try and gather all the matched clusters...
        int lastNMatched = 0;
        int nMatches = 0;
        do{
          lastNMatched = nMatches;
          nMatches = matchedClusters[z][b].size();
          for(int i=lastNMatched; i < nMatches; i++){

            int b2 = matchedClusters[z][b][i];

            for(int j=0; j < matchedClusters[z][b2].size(); j++){
              int b3 = matchedClusters[z][b2][j];

              if(!RootTools::vectorContainsValue(matchedClusters[z][b], b3)){
                matchedClusters[z][b].push_back(b3);
              }
            }
          }
        }
        while(lastNMatched!=nMatches);

        if(fDebug){
          std::cout << "At llEventCut = " << llEventCuts[z] << ", base " << b << " was matched with: " << matchedClusters[z][b].size() << " bases" << std::endl;
          if(matchedClusters[z][b].size() < 20){
            std::cout << "They are: ";
            std::sort(matchedClusters[z][b].begin(), matchedClusters[z][b].end());
            for(UInt_t i=0; i < matchedClusters[z][b].size(); i++){
              std::cout << matchedClusters[z][b][i];
              if(i != matchedClusters[z][b].size() - 1){
                std::cout << ", ";
              }
            }
            std::cout<< std::endl;
          }
        }

        for(int i=0; i < matchedClusters[z][b].size(); i++){
          int cluster = matchedClusters[z][b][i];

          if(cluster < reassignedTo[b]){
            reassignedTo[b] = cluster;
          }
        }
      }
    }


    if(fDebug){
      bool printed = false;
      for(int b=0; b < nBases; b++){
        if(b!=reassignedTo[b]){
          if(printed == false){
            std::cout << "At llEventCut = " << llEventCuts[z] << ", the known base cluster reassignments are as follows:" << std::endl;
            printed = true;
          }
          std::cout << b << "->" << reassignedTo[b] << "\n";
        }
      }
    }

    for(int eventInd=0; eventInd < events.size(); eventInd++){
      Event& event = events.at(eventInd);
      int clusterInd = -1;
      if(event.nearestKnownBaseLogLikelihood < llEventCuts.at(z)){
        clusterInd = event.nearestKnownBaseCluster;
      }
      else if(event.nearestKnownBaseSurfaceSeparationKm < surfaceDistThresholdKm){
        clusterInd = event.nearestKnownBaseClusterSurface;
      }
      if(clusterInd >= 0){

        Int_t reassignedCluster = reassignedTo[clusterInd];
        event.cluster[z] = reassignedCluster;
        clusters.at(z).at(reassignedCluster).numDataEvents++;
      }
    }
  }
}


// void Acclaim::Clustering::LogLikelihoodMethod::doBaseEventClustering(){

//   std::cout << "Info in " << __PRETTY_FUNCTION__ << ", starting..." << std::endl;
//   ProgressBar p(clusters.at(0).size());

//   const Int_t nBases = BaseList::getNumBases();


//   for(Int_t clusterInd=0; clusterInd < nBases; clusterInd++){
//     const Cluster& cluster = clusters.at(0).at(clusterInd);
//     Double_t clusterEasting, clusterNorthing;
//     RampdemReader::LonLatToEastingNorthing(cluster.longitude, cluster.latitude, clusterEasting, clusterNorthing);

//     std::vector<Int_t> nearbyEventIndsEastingNorthing;
//     double lookup[2] = {clusterEasting, clusterNorthing};
//     fKDTree->FindInRange(lookup, default_horizon_distance, nearbyEventIndsEastingNorthing);

//     for(UInt_t i=0; i < nearbyEventIndsEastingNorthing.size(); i++){
//       UInt_t eventInd = nearbyEventIndsEastingNorthing.at(i);

//       Event& event = events.at(eventInd);

//       Double_t distToAnita = event.usefulPat.getDistanceFromSource(cluster.latitude, cluster.longitude, cluster.altitude);
//       if(distToAnita < default_horizon_distance){

// 	double ll = event.logLikelihoodFromPoint(clusters.at(0).at(clusterInd));
// 	double surfaceSeparationKm = 1e-3*event.cartesianSeparation(clusters.at(0).at(clusterInd));

// 	for(Int_t z=0; z < llEventCuts.size(); z++){
// 	  std::vector<Int_t> matchedClusters;
// 	  if(ll < llEventCuts.at(z)){
// 	    matchedClusters.push_back(clusterInd);
// 	  }
// 	  else if(surfaceSeparationKm < surfaceDistThresholdKm){
// 	    matchedClusters.push_back(clusterInd);
// 	  }


// 	  if(ll < event.nearestKnownBaseLogLikelihood){
// 	    event.nearestKnownBaseLogLikelihood = ll;
// 	    event.nearestKnownBaseCluster = clusterInd;
// 	  }

// 	  if(surfaceSeparationKm < event.nearestKnownBaseSurfaceSeparationKm){
// 	    event.nearestKnownBaseSurfaceSeparationKm = surfaceSeparationKm;
// 	    event.nearestKnownBaseClusterSurface = clusterInd;
// 	  }
// 	}
//       }
//     }
//     p.inc(clusterInd);    
//   }
//       // 	double ll = event.logLikelihoodFromPoint(clusters.at(0).at(clusterInd));
//       // 	double surfaceSeparationKm = 1e-3*event.cartesianSeparation(clusters.at(0).at(clusterInd));

//       // 	if(surfaceSeparationKm < event.nearestKnownBaseSurfaceSeparationKm){
//       // 	  event.nearestKnownBaseSurfaceSeparationKm = surfaceSeparationKm;

//       // 	  if(event.nearestKnownBaseSurfaceSeparationKm < surfaceDistThresholdKm){
//       // 	    event.nearestKnownBaseCluster = clusterInd;
//       // 	  }
//       // 	}

//       // 	if(ll < event.nearestKnownBaseLogLikelihood){
//       // 	  event.nearestKnownBaseLogLikelihood = ll;
//       // 	  if(event.nearestKnownBaseSurfaceSeparationKm >= surfaceDistThresholdKm){
//       // 	    event.nearestKnownBaseCluster = clusterInd;
//       // 	  }
//       // 	}

//   //   }
//   //   p.inc(clusterInd);
//   // }

//   for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
//     Event& event = events.at(eventInd);
//     // std::cout << event.eventNumber << "\t" << event.nearestKnownBaseSurfaceSeparationKm << std::endl;
//     for(int z=0; z < event.nThresholds; z++){
//       // if(event.nearestKnownBaseLogLikelihood < llEventCuts.at(z) || event.nearestKnownBaseSurfaceSeparationKm < surfaceDistThresholdKm){
//       // 	Int_t clusterInd = event.nearestKnownBaseCluster;
//       // 	clusters.at(z).at(clusterInd).numDataEvents++;
//       // 	event.cluster[z] = clusterInd;
//       // }

//       if(event.nearestKnownBaseLogLikelihood < llEventCuts.at(z)){
// 	Int_t clusterInd = event.nearestKnownBaseCluster;
// 	clusters.at(z).at(clusterInd).numDataEvents++;
// 	event.cluster[z] = clusterInd;
//       }
//       else if(event.nearestKnownBaseSurfaceSeparationKm < surfaceDistThresholdKm){
// 	Int_t clusterInd = event.nearestKnownBaseClusterSurface;
// 	clusters.at(z).at(clusterInd).numDataEvents++;
// 	event.cluster[z] = clusterInd;
//       }
//       // else{
//       // 	// double surfaceSeparationKm = 1e-3*event.cartesianSeparation(clusters.at(0).at(clusterInd));
//       // 	// if()
//       // }

//     }
//   }
// }












size_t Acclaim::Clustering::LogLikelihoodMethod::addMcEvent(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd){

  mcEvents.push_back(McEvent(sum,  pol, peakInd, llEventCuts.size()));

  return mcEvents.size();
}






size_t Acclaim::Clustering::LogLikelihoodMethod::addEvent(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd){

  events.push_back(Event(sum, pol, peakInd, llEventCuts.size()));

  return events.size();
}








/** 
 * Puts an entry in each of the cluster[z] vectors for each of the known bases
 */
void Acclaim::Clustering::LogLikelihoodMethod::readInBaseList(){

  if(!fReadInBaseList){
    std::cout << "Info in " << __FUNCTION__ << ": Initializing base list..." << std::endl;

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
  fReadInBaseList = true;
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


void Acclaim::Clustering::LogLikelihoodMethod::addToHistograms(TH2D* h, TH2D* h2){
  for(int i = 0; i < events.size(); i++)
  {
    const Event& event = events.at(i);
    h->Fill(event.nearestEventSurfaceDistanceKm, event.nearestEventSurfaceLogLikelihood);
  }
  for(int i = 0; i < clusters.size(); i++)
  {
    if(llEventCuts.at(i) > 40) continue;
    
    int n_singlets = 0;
    for(int j = 0; j < clusters.at(i).size(); j++)
    {
      Cluster& cluster = clusters.at(i).at(j);
      if(cluster.numDataEvents == 1) n_singlets++;
      h2->Fill(llEventCuts.at(i), double(n_singlets)/double(events.size()));
    }
  }
}




void Acclaim::Clustering::LogLikelihoodMethod::makeSummaryTrees(){


  if(!fEventsAlreadyClustered){

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

    eventTree->Write();
    delete eventTree;
    eventTree = NULL;
  }



  // monte carlo tree
  if(mcEvents.size() > 0){
    TTree* mcEventTree = new TTree("mcEventTree", "Tree of clustered Monte Carlo ANITA events");
    McEvent* mcEvent = NULL;
    mcEventTree->Branch("mcEvent", &mcEvent);

    for(UInt_t j=0; j < mcEvents.size(); j++){
      mcEvent = &mcEvents.at(j);
      mcEventTree->Fill();
    }
    mcEventTree->Write();
    delete mcEventTree;
  }



  for(UInt_t z = 0; z < llEventCuts.size(); z++){
    TString treeName = TString::Format("clusterTree%u", z);
    TString treeTitle = TString::Format("Tree of clusters with llEventCut = %lf", llEventCuts[z]);
    //TTree* clusterTree = new TTree(treeName, "Tree of clusters");
    TTree* clusterTree = new TTree(treeName, treeTitle);
    Cluster* cluster = NULL;
    clusterTree->Branch("cluster", &cluster);
    std::vector<Cluster> &cs = clusters.at(z);
    for(Int_t k=0; k < (Int_t)clusters.at(z).size(); k++){
      cluster = &cs.at(k);
      clusterTree->Fill();
    }
  }




  for(UInt_t eventInd=0; eventInd < hUnclusteredEvents.size(); eventInd++){
    if(hUnclusteredEvents.at(eventInd)){
      hUnclusteredEvents.at(eventInd)->Write();
    }
  }


  if(!fEventsAlreadyClustered){
    TH2DAntarctica* hEvents = new TH2DAntarctica("hEvents", "hEvents");
    for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
      const Event& event = events.at(eventInd);
      hEvents->Fill(event.longitude, event.latitude);
    }
    hEvents->Write();
    delete hEvents;
  }



  if(mcEvents.size() > 0){
    TH2DAntarctica* hMcEvents = new TH2DAntarctica("hMcEvents", "hMcEvents");

    for(UInt_t eventInd=0; eventInd < mcEvents.size(); eventInd++){
      const McEvent& mcEvent = mcEvents.at(eventInd);
      hMcEvents->Fill(mcEvent.longitude, mcEvent.latitude, mcEvent.weight);
    }
    hMcEvents->Write();
    delete hMcEvents;
  }

}


Long64_t Acclaim::Clustering::LogLikelihoodMethod::readInSummaries(const char* summaryGlob){

  Long64_t n = 0;
  if(summaryGlob){
    /**
     * First, let's try and see if we're reading in the output of a clustering!
     * This will be the case in the new MC clustering paradigm
     */

    TFile* f = TFile::Open(summaryGlob);
    TTree* eventTree = NULL;
    if(f){
      eventTree = (TTree*) f->Get("eventTree");
    }

    /**
     * Here we enter the viper's nest of resetting the cluster/event info.
     * This must overwrite the this->llEventCuts vector to match the clusters in this file.
     */
    if(eventTree){
      std::cout << "Info in " << __PRETTY_FUNCTION__ << ": reading in already clustered events: " << summaryGlob << std::endl;
      Event* event = NULL;
      eventTree->SetBranchAddress("event", &event);
      n = eventTree->GetEntries();
      events.reserve(n);
      ProgressBar p(n);
      for(Long64_t entry=0; entry < n; entry++){
        eventTree->GetEntry(entry);
        event->setupUsefulPat(false);
        events.push_back(*event);
        p.inc(entry);	
      }
      std::cout << "Read in " << n << " events." << std::endl;

      clusters.clear();
      llEventCuts.clear();

      Acclaim::Clustering::Cluster* cluster = NULL;
      TTree* clusterTree = NULL;
      Int_t treeInd=0;
      do {
        TString treeName = TString::Format("clusterTree%d", treeInd);
        clusterTree = (TTree*)f->Get(treeName);
        if(clusterTree){
          clusterTree->SetBranchAddress("cluster", &cluster);
          const Long64_t nC = clusterTree->GetEntries();

          clusters.push_back(std::vector<Cluster>());
          clusters.back().reserve(nC);

          for(Long64_t entry=0; entry < nC; entry++){
            clusterTree->GetEntry(entry);
            if(entry==0){
              llEventCuts.push_back(cluster->llEventCut);
            }
            clusters.back().push_back(*cluster);
          }
          treeInd++;
        }
      } while(clusterTree!=NULL);
      f->Close();



      std::cout << "Info in " << __PRETTY_FUNCTION__ << ", overwrote llEventCuts to: ";
      for(UInt_t z=0; z < llEventCuts.size(); z++){
        std::cout << llEventCuts[z];
        if(z < llEventCuts.size() - 1){
          std::cout << ", ";
        }
      }
      std::cout << std::endl;


      // make sure to not do event clustering again, or read in base list,
      // that hard work was already done...
      fEventsAlreadyClustered = true;
      fReadInBaseList = true;
    }

    else{
      ThermalChain tc(summaryGlob);
      ProgressBar pElist(1);

      // TCut hack("eventNumber==14545077||eventNumber==15202247");
      // TCut goodPosition = "(onContinent > 0 && onIceShelf==0)";
      // tc.setCut(goodPosition + !ThermalTree::isAboveHorizontal + ThermalTree::passAllQualityCuts + ThermalTree::isNotTaggedAsPulser + ThermalTree::fisherCut + !ThermalTree::closeToHiCal);
      // tc.setCut(!ThermalTree::isAboveHorizontal + ThermalTree::passAllQualityCuts + ThermalTree::isNotTaggedAsPulser + ThermalTree::fisherCut + !ThermalTree::closeToHiCal);
      // tc.setCut(!ThermalTree::isAboveHorizontal + ThermalTree::passAllQualityCuts + ThermalTree::isNotTaggedAsPulser + ThermalTree::fisherCut + !ThermalTree::closeToHiCal);      
      // tc.setCut(!ThermalTree::isAboveHorizontal + ThermalTree::passAllQualityCuts + ThermalTree::isNotTaggedAsPulser + ThermalTree::fisherCut + !ThermalTree::closeToHiCal + ThermalTree::closeToMC);
      // tc.setCut(!ThermalTree::isAboveHorizontal + ThermalTree::passAllQualityCuts + ThermalTree::isNotTaggedAsPulser + ThermalTree::fisherCut + ThermalTree::closeToMC);      
      tc.setCut(ThermalTree::isTaggedAsWaisPulser + ThermalTree::closeToWais);
      // tc.setCut(hack);

      n = tc.N();
      std::cout << "There are " << n << " entries matching the selection" << std::endl;

      ProgressBar p(n);
      for(Long64_t entry=0; entry < n; entry++){
        tc.getEntry(entry);

        if(entry==0){
          std::cout << tc.weight << std::endl;
          if(tc.weight == 1)
            events.reserve(events.size()+n);
          else{	    
            mcEvents.reserve(mcEvents.size()+n);
          }
        }

        if(tc.weight==1){
          events.push_back(Event(static_cast<int>(tc.pol), static_cast<int>(tc.peakInd),
                tc.peak_phi, tc.peak_theta,
                (int)llEventCuts.size(), tc.eventNumber, tc.run,
                tc.anita_longitude, tc.anita_latitude, tc.anita_altitude, tc.anita_heading,
                tc.coherent_filtered_snr, tc.deconvolved_filtered_snr));
//                tc.coherent_filtered_snr));
        }
        else{
          mcEvents.push_back(McEvent(tc.weight, tc.mc_energy, static_cast<int>(tc.pol), static_cast<int>(tc.peakInd),
                tc.peak_phi, tc.peak_theta,
                (int)llEventCuts.size(), tc.eventNumber, tc.run,
                tc.anita_longitude, tc.anita_latitude, tc.anita_altitude, tc.anita_heading,
                tc.coherent_filtered_snr, tc.deconvolved_filtered_snr));
//                tc.coherent_filtered_snr));
        }
        p.inc(entry);
      }
    }
  }
  return n;
}


Long64_t Acclaim::Clustering::LogLikelihoodMethod::readInTMVATreeSummaries(const char* summaryGlob, bool isMC){

  Long64_t n = 0;
  if(summaryGlob){
    /**
     * First, let's try and see if we're reading in the output of a clustering!
     * This will be the case in the new MC clustering paradigm
     */

    TFile* f = TFile::Open(summaryGlob);
    TTree* eventTree = NULL;
    if(f){
      eventTree = (TTree*) f->Get("eventTree");
    }

    /**
     * Here we enter the viper's nest of resetting the cluster/event info.
     * This must overwrite the this->llEventCuts vector to match the clusters in this file.
     */
    if(eventTree){
      std::cout << "Info in " << __PRETTY_FUNCTION__ << ": reading in already clustered events: " << summaryGlob << std::endl;
      Event* event = NULL;
      eventTree->SetBranchAddress("event", &event);
      n = eventTree->GetEntries();
      events.reserve(n);
      ProgressBar p(n);
      for(Long64_t entry=0; entry < n; entry++){
        eventTree->GetEntry(entry);
        event->setupUsefulPat(false);
        events.push_back(*event);
        p.inc(entry);	
      }
      std::cout << "Read in " << n << " events." << std::endl;

      clusters.clear();
      llEventCuts.clear();

      Acclaim::Clustering::Cluster* cluster = NULL;
      TTree* clusterTree = NULL;
      Int_t treeInd=0;
      do {
        TString treeName = TString::Format("clusterTree%d", treeInd);
        clusterTree = (TTree*)f->Get(treeName);
        if(clusterTree){
          clusterTree->SetBranchAddress("cluster", &cluster);
          const Long64_t nC = clusterTree->GetEntries();

          clusters.push_back(std::vector<Cluster>());
          clusters.back().reserve(nC);

          for(Long64_t entry=0; entry < nC; entry++){
            clusterTree->GetEntry(entry);
            if(entry==0){
              llEventCuts.push_back(cluster->llEventCut);
            }
            clusters.back().push_back(*cluster);
          }
          treeInd++;
        }
      } while(clusterTree!=NULL);
      f->Close();

      std::cout << "Info in " << __PRETTY_FUNCTION__ << ", overwrote llEventCuts to: ";
      for(UInt_t z=0; z < llEventCuts.size(); z++){
        std::cout << llEventCuts[z];
        if(z < llEventCuts.size() - 1){
          std::cout << ", ";
        }
      }
      std::cout << std::endl;

      // make sure to not do event clustering again, or read in base list,
      // that hard work was already done...
      fEventsAlreadyClustered = true;
      fReadInBaseList = true;
    } else {
      TChain* t = new TChain("sumTree");
      t->Add(summaryGlob);

      float decoImpulsivity, pol, peakInd, run, anita_longitude, anita_latitude, anita_altitude, anita_heading, peak_phi, peak_theta, coherent_filtered_snr, deconvolved_filtered_snr, F, lastFew, weight, mc_energy, isWais;
      UInt_t eventNumber;
      Int_t evNum;

      t->SetBranchAddress("pol", &pol);
      t->SetBranchAddress("ind", &peakInd);
      t->SetBranchAddress("weight", &weight);
      t->SetBranchAddress("energy", &mc_energy);
      t->SetBranchAddress("phi", &peak_phi);
      t->SetBranchAddress("theta", &peak_theta);
      t->SetBranchAddress("run", &run);
      t->SetBranchAddress("anita_latitude", &anita_latitude);
      t->SetBranchAddress("anita_longitude", &anita_longitude);
      t->SetBranchAddress("anita_altitude", &anita_altitude);
      t->SetBranchAddress("anita_heading", &anita_heading);
      t->SetBranchAddress("snr", &coherent_filtered_snr);
      t->SetBranchAddress("eventNumber", &evNum);
      t->SetBranchAddress("lastFewDigits", &lastFew);
      t->SetBranchAddress("F", &F);
      t->SetBranchAddress("isWais", &isWais);
      t->SetBranchAddress("decoImpulsivity", &decoImpulsivity);

      t->Draw(">>fEntryList", fCut, "entrylist");
      fEntryList = (TEntryList*) gDirectory->Get("fEntryList");
      t->SetEntryList(fEntryList);
      printf("%d entries loaded\n", fEntryList->GetN());

      for(Long64_t entry=0; entry < fEntryList->GetN(); entry++){
        n++;
        t->GetEntry(t->GetEntryNumber(entry));
        eventNumber = UInt_t( int(evNum / 10000) * 10000 + int(lastFew));
//        if(fCutHical && Hical2::isHical(eventNumber, FFTtools::wrap(anita_heading - peak_phi, 360, 0), coherent_filtered_snr)) continue;
        if(peak_theta > 0) {

          // switches theta convention (i used the UCorrelator convention for theta)
          peak_theta = -1* peak_theta;
          events.push_back(Event(static_cast<int>(pol), static_cast<int>(peakInd),
                (double)peak_phi, (double)peak_theta,
                (int)llEventCuts.size(), eventNumber, (int)run,
                (double)anita_longitude, (double)anita_latitude, (double)anita_altitude, (double)anita_heading,
		(double)coherent_filtered_snr, (double)deconvolved_filtered_snr));
//                (double)coherent_filtered_snr));
          if(fSelfLLMax > 0 && events.back().selfLogLikelihood > fSelfLLMax) events.pop_back();
        }
      }
      delete t;
    }
  }
  return n;
}


Long64_t Acclaim::Clustering::LogLikelihoodMethod::readInSampleSummaries(const char* summaryGlob, bool isMC){

  Long64_t n = 0;
  if(summaryGlob){
    /**
     * First, let's try and see if we're reading in the output of a clustering!
     * This will be the case in the new MC clustering paradigm
     */

    TFile* f = TFile::Open(summaryGlob);
    TTree* eventTree = NULL;
    if(f){
      eventTree = (TTree*) f->Get("eventTree");
    }

    /**
     * Here we enter the viper's nest of resetting the cluster/event info.
     * This must overwrite the this->llEventCuts vector to match the clusters in this file.
     */
    if(eventTree){
      std::cout << "Info in " << __PRETTY_FUNCTION__ << ": reading in already clustered events: " << summaryGlob << std::endl;
      Event* event = NULL;
      eventTree->SetBranchAddress("event", &event);
      n = eventTree->GetEntries();
      events.reserve(n);
      ProgressBar p(n);
      for(Long64_t entry=0; entry < n; entry++){
        eventTree->GetEntry(entry);
        event->setupUsefulPat(false);
        events.push_back(*event);
        p.inc(entry);	
      }
      std::cout << "Read in " << n << " events." << std::endl;

      clusters.clear();
      llEventCuts.clear();

      Acclaim::Clustering::Cluster* cluster = NULL;
      TTree* clusterTree = NULL;
      Int_t treeInd=0;
      do {
        TString treeName = TString::Format("clusterTree%d", treeInd);
        clusterTree = (TTree*)f->Get(treeName);
        if(clusterTree){
          clusterTree->SetBranchAddress("cluster", &cluster);
          const Long64_t nC = clusterTree->GetEntries();

          clusters.push_back(std::vector<Cluster>());
          clusters.back().reserve(nC);

          for(Long64_t entry=0; entry < nC; entry++){
            clusterTree->GetEntry(entry);
            if(entry==0){
              llEventCuts.push_back(cluster->llEventCut);
            }
            clusters.back().push_back(*cluster);
          }
          treeInd++;
        }
      } while(clusterTree!=NULL);
      f->Close();



      std::cout << "Info in " << __PRETTY_FUNCTION__ << ", overwrote llEventCuts to: ";
      for(UInt_t z=0; z < llEventCuts.size(); z++){
        std::cout << llEventCuts[z];
        if(z < llEventCuts.size() - 1){
          std::cout << ", ";
        }
      }
      std::cout << std::endl;


      // make sure to not do event clustering again, or read in base list,
      // that hard work was already done...
      fEventsAlreadyClustered = true;
      fReadInBaseList = true;
    } else {

      TString summaryGlobStr(summaryGlob);
//      TString sampleStr;
//
//      if (summaryGlobStr.Contains("minBias")) sampleStr = "minBias";
//      else if (summaryGlobStr.Contains("other")) sampleStr = "other";
//      else if (summaryGlobStr.Contains("payloadBlast")) sampleStr = "payloadBlast";
//      else if (summaryGlobStr.Contains("refinedHiCal2A")) sampleStr = "refinedHiCal2A";
//      else if (summaryGlobStr.Contains("refinedHiCal2B")) sampleStr = "refinedHiCal2B";
//      else if (summaryGlobStr.Contains("refinedWAISHPol")) sampleStr = "refinedWAISHPol";
//      else if (summaryGlobStr.Contains("refinedWAISVPol")) sampleStr = "refinedWAISVPol";
//      else if (summaryGlobStr.Contains("signal")) sampleStr = "signal";
//      else if (summaryGlobStr.Contains("strongCW")) sampleStr = "strongCW";
//      else if (summaryGlobStr.Contains("thermal")) sampleStr = "thermal";

      TChain * t = new TChain("sampleA4");
//      TChain * t = new TChain(sampleStr);
//      TChain* t = new TChain("sumTree");
      t -> Add(summaryGlob);

      AnitaEventSummary * sampleSum = 0;
      t -> SetBranchAddress("summary", & sampleSum);
      
      t -> GetEntry(0);

      int & run = sampleSum -> run;
      unsigned int & eventNumber = sampleSum -> eventNumber;
      float & anita_longitude = sampleSum -> anitaLocation.longitude;
      float & anita_latitude = sampleSum -> anitaLocation.latitude;
      float & anita_altitude = sampleSum -> anitaLocation.altitude;
      float & anita_heading = sampleSum -> anitaLocation.heading;
      double & mc_weight = sampleSum -> mc.weight;
      double & mc_energy = sampleSum -> mc.energy;

      //  Variables from AnitaEventSummary which can't be referenced like those above.
      int pol, peakInd;
      double peak_theta, peak_phi;
      double coherent_filtered_snr, deconvolved_filtered_snr;
      
//      float pol, peakInd, run, anita_longitude, anita_latitude, anita_altitude, anita_heading, peak_phi, peak_theta, coherent_filtered_snr, deconvolved_filtered_snr, lastFew;
//      float decoImpulsivity, pol, peakInd, run, anita_longitude, anita_latitude, anita_altitude, anita_heading, peak_phi, peak_theta, coherent_filtered_snr, F, lastFew, weight, mc_energy, isWais;
//      UInt_t eventNumber;
//      Int_t evNum;
//
//      t->SetBranchAddress("pol", &pol);
//      t->SetBranchAddress("ind", &peakInd);
//      t->SetBranchAddress("weight", &weight);
//      t->SetBranchAddress("energy", &mc_energy);
//      t->SetBranchAddress("phi", &peak_phi);
//      t->SetBranchAddress("theta", &peak_theta);
//      t->SetBranchAddress("run", &run);
//      t->SetBranchAddress("anita_latitude", &anita_latitude);
//      t->SetBranchAddress("anita_longitude", &anita_longitude);
//      t->SetBranchAddress("anita_altitude", &anita_altitude);
//      t->SetBranchAddress("anita_heading", &anita_heading);
//      t->SetBranchAddress("coherent_filtered_snr", &coherent_filtered_snr);
//      t->SetBranchAddress("deconvovled_filtered_snr", &deconvolved_filtered_snr);
//      t->SetBranchAddress("eventNumber", &evNum);
//      t->SetBranchAddress("lastFewDigits", &lastFew);
//      t->SetBranchAddress("F", &F);
//      t->SetBranchAddress("isWais", &isWais);
//      t->SetBranchAddress("decoImpulsivity", &decoImpulsivity);

      t->Draw(">>fEntryList", fCut, "entrylist");
      fEntryList = (TEntryList*) gDirectory->Get("fEntryList");
      t->SetEntryList(fEntryList);
      printf("%d entries loaded\n", fEntryList->GetN());
      
      //  Create vector of length fEntryList -> GetN() with randomly shuffled indices.
      std::vector<int> entryListIdx(fEntryList -> GetN());
      std::iota(std::begin(entryListIdx), std::end(entryListIdx), 0);
      std::shuffle(std::begin(entryListIdx), std::end(entryListIdx), std::mt19937_64(0));

      for(Long64_t entry=0; entry < fEntryList->GetN(); entry++){
        n++;
        t -> GetEntry(t->GetEntryNumber(entry));
//        eventNumber = UInt_t(int(evNum/10000) * 10000 + int(lastFew));
//        if(fCutHical && Hical2::isHical(eventNumber, FFTtools::wrap(anita_heading - peak_phi, 360, 0), coherent_filtered_snr)) continue;
//        if(peak_theta > 0) {

        pol = sampleSum -> mostImpulsivePolAsInt(2);
        peakInd = sampleSum -> mostImpulsiveInd(2);
        peak_theta = sampleSum -> mostImpulsivePeak(2).theta;
        peak_theta *= -1;  //  Switches theta convention from UCorrelator convention.
        peak_phi = sampleSum -> mostImpulsivePeak(2).phi;
        coherent_filtered_snr = sampleSum -> mostImpulsiveCoherent(2).snr;
        deconvolved_filtered_snr = sampleSum -> mostImpulsiveDeconvolved(2).snr;
//        // Switches theta convention (using the UCorrelator convention for theta)
//        peak_theta = -1* peak_theta;

	if (!isMC && peak_theta < 0) {

          events.push_back(Event(static_cast<int>(pol), static_cast<int>(peakInd),
                 peak_phi, peak_theta,
                 (int) llEventCuts.size(), eventNumber, run,
                 (double) anita_longitude, (double) anita_latitude, (double) anita_altitude, (double) anita_heading,
                 coherent_filtered_snr, deconvolved_filtered_snr));
//          events.push_back(Event(static_cast<int>(pol), static_cast<int>(peakInd),
//                (double)peak_phi, (double)peak_theta,
//                (int)llEventCuts.size(), eventNumber, (int)run,
//                (double)anita_longitude, (double)anita_latitude, (double)anita_altitude, (double)anita_heading,
//                (double)coherent_filtered_snr, (double)deconvolved_filtered_snr));
//                (double)coherent_filtered_snr));

          if(fSelfLLMax > 0 && events.back().selfLogLikelihood > fSelfLLMax) events.pop_back();
        }

	if (isMC) {

          if (fPermyriadOfMC && fNumOfMC) continue;
          if (tr3 -> Integer(10001) >= fPermyriadOfMC && !fNumOfMC) continue; 
          if (!fPermyriadOfMC && entryListIdx[entry] >= fNumOfMC) continue;
//          if (!fPermyriadOfMC && tr3 -> Integer(fEntryList -> GetN()) >= fNumOfMC) continue;
//          // switches theta convention
//          peak_theta = -1* peak_theta;
          mcEvents.push_back(McEvent((double) mc_weight, (double) mc_energy, static_cast<int>(pol), static_cast<int>(peakInd),
                peak_phi, peak_theta,
                (int)llEventCuts.size(), eventNumber, (int)run,
                (double)anita_longitude, (double)anita_latitude, (double)anita_altitude, (double)anita_heading,
                coherent_filtered_snr, deconvolved_filtered_snr));
//          mcEvents.push_back(McEvent((double)weight, (double)mc_energy, static_cast<int>(pol), static_cast<int>(peakInd),
//                (double)peak_phi, (double)peak_theta,
//                (int)llEventCuts.size(), eventNumber, (int)run,
//                (double)anita_longitude, (double)anita_latitude, (double)anita_altitude, (double)anita_heading,
//                (double)coherent_filtered_snr));
        }

      }
      delete t;
    }
  }
  return n;
}


void Acclaim::Clustering::LogLikelihoodMethod::readInSummariesForTesting(const char* summaryGlob){

  Long64_t n = 0;
  if(summaryGlob){
    /**
     * First, let's try and see if we're reading in the output of a clustering!
     * This will be the case in the new MC clustering paradigm
     */

      fChain = new TChain("sumTree");
      fChain->Add(summaryGlob);

      float decoImpulsivity, pol, peakInd, run, anita_longitude, anita_latitude, anita_altitude, anita_heading, peak_phi, peak_theta, coherent_filtered_snr, F, lastFew, weight, mc_energy, isWais;
      UInt_t eventNumber;
      Int_t evNum;

      fChain->SetBranchAddress("pol", &pol);
      fChain->SetBranchAddress("ind", &peakInd);
      fChain->SetBranchAddress("weight", &weight);
      fChain->SetBranchAddress("energy", &mc_energy);
      fChain->SetBranchAddress("phi", &peak_phi);
      fChain->SetBranchAddress("theta", &peak_theta);
      fChain->SetBranchAddress("run", &run);
      fChain->SetBranchAddress("anita_latitude", &anita_latitude);
      fChain->SetBranchAddress("anita_longitude", &anita_longitude);
      fChain->SetBranchAddress("anita_altitude", &anita_altitude);
      fChain->SetBranchAddress("anita_heading", &anita_heading);
      fChain->SetBranchAddress("snr", &coherent_filtered_snr);
      fChain->SetBranchAddress("eventNumber", &evNum);
      fChain->SetBranchAddress("lastFewDigits", &lastFew);
      fChain->SetBranchAddress("F", &F);
      fChain->SetBranchAddress("isWais", &isWais);
      fChain->SetBranchAddress("decoImpulsivity", &decoImpulsivity);

      fChain->Draw(">>fEntryList", fCut, "entrylist");
      fEntryList = (TEntryList*) gDirectory->Get("fEntryList");
      fChain->SetEntryList(fEntryList);
      printf("%d entries loaded\n", fEntryList->GetN());
  }
  return;
}


void Acclaim::Clustering::LogLikelihoodMethod::readInSampleSummariesForTesting(const char* summaryGlob){

  Long64_t n = 0;
  if(summaryGlob){
    /**
     * First, let's try and see if we're reading in the output of a clustering!
     * This will be the case in the new MC clustering paradigm
     */

      TString summaryGlobStr(summaryGlob);
//      TString sampleStr;
//
//      if (summaryGlobStr.Contains("iceMC")) sampleStr = "iceMC";
//      else if (summaryGlobStr.Contains("minBias")) sampleStr = "minBias";
//      else if (summaryGlobStr.Contains("other")) sampleStr = "other";
//      else if (summaryGlobStr.Contains("payloadBlast")) sampleStr = "payloadBlast";
//      else if (summaryGlobStr.Contains("HiCal2A")) sampleStr = "HiCal2A";
//      else if (summaryGlobStr.Contains("HiCal2B")) sampleStr = "HiCal2B";
//      else if (summaryGlobStr.Contains("WAISHPol")) sampleStr = "WAISHPol";
//      else if (summaryGlobStr.Contains("WAISVPol")) sampleStr = "WAISVPol";
//      else if (summaryGlobStr.Contains("signal")) sampleStr = "signal";
//      else if (summaryGlobStr.Contains("strongCW")) sampleStr = "strongCW";
//      else if (summaryGlobStr.Contains("thermal")) sampleStr = "thermal";

      TChain * fChain = new TChain("sampleA4");
//      TChain * fChain = new TChain(sampleStr);
      fChain -> Add(summaryGlob);

      AnitaEventSummary * sampleSum = 0;
      fChain -> SetBranchAddress("summary", & sampleSum);

      fChain -> GetEntry(0);

      fChain -> Draw(">>fEntryList", fCut, "entrylist");
      fEntryList = (TEntryList*) gDirectory->Get("fEntryList");
      fChain -> SetEntryList(fEntryList);
      printf("%d entries loaded\n", fEntryList->GetN());
  }
  return;
}


void Acclaim::Clustering::LogLikelihoodMethod::pickEventsFromList(int n_in_cluster)
{
  float decoImpulsivity, pol, peakInd, run, anita_longitude, anita_latitude, anita_altitude, anita_heading, peak_phi, peak_theta, coherent_filtered_snr, deconvolved_filtered_snr, F, lastFew, weight, mc_energy, isWais;
  UInt_t eventNumber;
  Int_t evNum;

  fChain->SetBranchAddress("pol", &pol);
  fChain->SetBranchAddress("ind", &peakInd);
  fChain->SetBranchAddress("weight", &weight);
  fChain->SetBranchAddress("energy", &mc_energy);
  fChain->SetBranchAddress("phi", &peak_phi);
  fChain->SetBranchAddress("theta", &peak_theta);
  fChain->SetBranchAddress("run", &run);
  fChain->SetBranchAddress("anita_latitude", &anita_latitude);
  fChain->SetBranchAddress("anita_longitude", &anita_longitude);
  fChain->SetBranchAddress("anita_altitude", &anita_altitude);
  fChain->SetBranchAddress("anita_heading", &anita_heading);
  fChain->SetBranchAddress("snr", &coherent_filtered_snr);
  fChain->SetBranchAddress("eventNumber", &evNum);
  fChain->SetBranchAddress("lastFewDigits", &lastFew);
  fChain->SetBranchAddress("F", &F);
  fChain->SetBranchAddress("isWais", &isWais);
  fChain->SetBranchAddress("decoImpulsivity", &decoImpulsivity);

  TBits * bits = new TBits(fEntryList->GetN());

  int i = 0;
  while(i < n_in_cluster){
    Int_t j = tr3->Uniform(0, fEntryList->GetN());
    if(bits->TestBitNumber(j)) continue;
    bits->SetBitNumber(j);
    fChain->GetEntry(fChain->GetEntryNumber(j));
    eventNumber = UInt_t(int(evNum/10000)*10000 + int(lastFew));
    if(peak_theta > 0)
    {
      i++;
      // switches theta convention (i used the UCorrelator convention for theta)
      peak_theta *= -1;
      events.push_back(Event(static_cast<int>(pol), static_cast<int>(peakInd),
            (double)peak_phi, (double)peak_theta,
            (int)llEventCuts.size(), eventNumber, (int)run,
            (double)anita_longitude, (double)anita_latitude, (double)anita_altitude, (double)anita_heading,
            (double)coherent_filtered_snr, (double)deconvolved_filtered_snr));
//            (double)coherent_filtered_snr));
    }
  }
  
//  std::cout << "Is this function called in standard processing?" << std::endl;
  
  delete bits;
}


void Acclaim::Clustering::LogLikelihoodMethod::pickSampleEventsFromList(int n_in_cluster)
{
  int pol, peakInd, run;
  unsigned int eventNumber;
  float anita_longitude, anita_latitude, anita_altitude, anita_heading;
  double peak_phi, peak_theta;
  double coherent_filtered_snr, deconvolved_filtered_snr;
  float weight, mc_energy;
//  float pol, peakInd, run, anita_longitude, anita_latitude, anita_altitude, anita_heading, peak_phi, peak_theta, coherent_filtered_snr, deconvolved_filtered_snr;
//  float decoImpulsivity, pol, peakInd, run, anita_longitude, anita_latitude, anita_altitude, anita_heading, peak_phi, peak_theta, coherent_filtered_snr, F, lastFew, weight, mc_energy, isWais;
//  UInt_t eventNumber;
//  Int_t evNum;

  fChain->SetBranchAddress("pol", &pol);
  fChain->SetBranchAddress("ind", &peakInd);
  fChain->SetBranchAddress("weight", &weight);
  fChain->SetBranchAddress("energy", &mc_energy);
  fChain->SetBranchAddress("phi", &peak_phi);
  fChain->SetBranchAddress("theta", &peak_theta);
  fChain->SetBranchAddress("run", &run);
  fChain->SetBranchAddress("anita_latitude", &anita_latitude);
  fChain->SetBranchAddress("anita_longitude", &anita_longitude);
  fChain->SetBranchAddress("anita_altitude", &anita_altitude);
  fChain->SetBranchAddress("anita_heading", &anita_heading);
  fChain->SetBranchAddress("coherent_filtered_snr", &coherent_filtered_snr);
  fChain->SetBranchAddress("deconvolved_filtered_snr", &deconvolved_filtered_snr);
  fChain->SetBranchAddress("eventNumber", & eventNumber);
//  fChain->SetBranchAddress("eventNumber", &evNum);
//  fChain->SetBranchAddress("lastFewDigits", &lastFew);
//  fChain->SetBranchAddress("F", &F);
//  fChain->SetBranchAddress("isWais", &isWais);
//  fChain->SetBranchAddress("decoImpulsivity", &decoImpulsivity);

  TBits * bits = new TBits(fEntryList->GetN());

  int i = 0;
  while(i < n_in_cluster){
    Int_t j = tr3->Uniform(0, fEntryList->GetN());
    if(bits->TestBitNumber(j)) continue;
    bits->SetBitNumber(j);
    fChain->GetEntry(fChain->GetEntryNumber(j));
//    eventNumber = UInt_t(int(evNum/10000)*10000 + int(lastFew));
    if(peak_theta > 0)
    {
      i++;
      // switches theta convention (i used the UCorrelator convention for theta)
      peak_theta *= -1;
      events.push_back(Event(static_cast<int>(pol), static_cast<int>(peakInd),
            peak_phi, peak_theta,
            (int)llEventCuts.size(), eventNumber, run,
            (double)anita_longitude, (double)anita_latitude, (double)anita_altitude, (double)anita_heading,
            coherent_filtered_snr, deconvolved_filtered_snr));
//            (double)coherent_filtered_snr));
    }
  }
  
//  std::cout << "Is this function called in standard processing?" << std::endl;
  
  delete bits;
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
  if(fKDTree){
    delete fKDTree;
  }  
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




  // const int n = 2;
  // std::vector<int> nns(n);
  // std::vector<double> dists(n);

  // Acclaim::ProgressBar p(events.size());
  // for(Long64_t eventInd=0; eventInd < events.size(); eventInd++){
  //   Event& event = events.at(eventInd);
  //   fKDTree->FindNearestNeighbors(&event.easting, n, &nns[0], &dists[0]);
  //   event.nearestEventSurfaceDistanceKm = 1e-3*dists.at(1);
  //   const Event& event2 = events.at(nns.at(1));
  //   event.nearestEventSurfaceEventNumber = event2.eventNumber;
  //   events.at(eventInd).nearestEventSurfaceLogLikelihood = dFit(&event, &event2);

  //   p.inc(eventInd);
  // }  
}






void Acclaim::Clustering::LogLikelihoodMethod::testSmallClustersFromPointSource(){

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
  TRandom3 randy(123);

  const int fakeClusterSize = 2;
  const int nTrials = 1000;

  TH2D* hClusterTrials = new TH2D("hClusterTrials", "clustering trials",
      llEventCuts.size(), 0, llEventCuts.size(),
      fakeClusterSize, 1, fakeClusterSize+1);


  for(int threshInd=0; threshInd < llEventCuts.size(); threshInd++){
    TString bl = TString::Format("%4.0lf", llEventCuts.at(threshInd));
    hClusterTrials->GetXaxis()->SetBinLabel(threshInd+1, bl);
  }


  std::vector<Event> tempEvents;
  tempEvents.reserve(events.size());
  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
    tempEvents.push_back(events.at(eventInd));
  }

  // pick random numbers
  ProgressBar p(nTrials);
  for(int trialInd=0; trialInd < nTrials; trialInd++){


    events.clear();

    for(UInt_t z=0; z < clusters.size(); z++){
      clusters.at(z).clear();
    }

    for(UInt_t i=0; i < fakeClusterSize; i++){
      UInt_t eventInd = randy.Uniform(tempEvents.size());
      events.push_back(tempEvents.at(eventInd));
    }    
    initKDTree();

    doEventEventClustering();

    for(UInt_t z=0; z < clusters.size(); z++){
      hClusterTrials->Fill(z, clusters.at(z).size());
    }

    p.inc(trialInd);
  }
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


  const Int_t nBins = 1024; //(lastEventNumber + 1) - firstEventNumber;
  // const Int_t stride = 50;
  std::cout << nBins << "\t" << firstEventNumber << "\t" << lastEventNumber << std::endl;
  TH2D* hUnfit = new TH2D("hUnfit_d_vs_pointingAngle", "", 1024, 0, 180, 1024, 0, 1024);
  TH2D* hFit = new TH2D("hFit_d_vs_pointingAngle", "", 1024, 0, 180, 1024, 0, 1024);

  TH2D* hFitVsSumOfSelfLLs = new TH2D("hFit_vs_sumOfSelfLLs", "", 1024, 0, 50, 1024, 0, 1024);

  TH2D* hFitVsUnfitted = new TH2D("hFit_vs_unfit", "Fitted vs. unfitted WAIS event-event log-likelihood; Unfitted log-likelihood; Fitted log-likelihood", 1024, 0, 1024, 1024, 0, 1024);
  TH1D* hUnfitMinusFit = new TH1D("hFit_minus_unfit", "(Unfitted - fitted) WAIS event-event log-likelihood; #delta (log likelihood)", 1024, -512, 512);

  TH2DAntarctica* hEventPos = new TH2DAntarctica("hEventPos", "Event Position");
  TH2DAntarctica* hFittedPos = new TH2DAntarctica("hFittedPos", "Fitted Pair position");

  TGraphAntarctica* grAnita = new TGraphAntarctica();
  // TGraphAntarctica* grEventPos = new TGraphAntarctica();
  // TGraphAntarctica* grFittedPos = new TGraphAntarctica();

  const int fakeClusterSize = 2;

  for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
    const Event& event1 = events.at(eventInd);

    // if(event1.eventNumber!=fTestEvent1) continue;

    TVector3 event1Pos = AntarcticCoord(AntarcticCoord::WGS84, event1.latitude, event1.longitude, event1.altitude).v();
    TVector3 anita1Pos = AntarcticCoord(AntarcticCoord::WGS84, event1.anita.latitude, event1.anita.longitude, event1.anita.altitude).v();
    TVector3 anitaToEvent1 = event1Pos - anita1Pos;

    std::vector<Int_t> event2Inds(fakeClusterSize-1);
    for(UInt_t eventInd2=0; eventInd2 < event2Inds.size(); eventInd2++){
      event2Inds.at(eventInd2) = randy.Uniform(0, events.size());
    }

    for(UInt_t eventInd2=0; eventInd2 < event2Inds.size(); eventInd2++){
      const Event& event2 = events.at(eventInd2);

      if(event1.eventNumber==fTestEvent1 && event2.eventNumber==fTestEvent2){
	std::cerr << "Info in " << __PRETTY_FUNCTION__ << " mapping parameter space around WAIS for event pair !"
		  << fTestEvent1 << " and " << fTestEvent2 << std::endl;
	AnitaVersion::set(3);
	double waisEasting, waisNorthing;
	RampdemReader::LonLatToEastingNorthing(AnitaLocations::getWaisLongitude(), AnitaLocations::getWaisLatitude(), waisEasting, waisNorthing);
	double delta = 700e3;
	const int nBins = 256;
	TH2D* hParams = new TH2D("hSingleEventTest", "Event-event fitted log likelihood; Easting (m); Northing (m); L_{sum}",
				 nBins, waisEasting-delta, waisEasting+delta,
				 nBins, waisNorthing-delta, waisNorthing+delta);
	grAnita->SetPoint(grAnita->GetN(), event1.usefulPat.longitude, event1.usefulPat.latitude);	
	grAnita->SetPoint(grAnita->GetN(), event2.usefulPat.longitude, event2.usefulPat.latitude);

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

      double distFitted = dFit(&event1, &event2);
      hFit->Fill(angleBetweenEvents, distFitted);

      hFitVsSumOfSelfLLs->Fill(event1.selfLogLikelihood+event2.selfLogLikelihood, distFitted);

      if(event1.theta > -5.5 && event2.theta > -5.5){
        std::cout << event1.eventNumber << "\t" << event2.eventNumber << std::endl;
      }

      double fitLon, fitLat;
      int t = OpenMP::thread();
      RampdemReader::EastingNorthingToLonLat(fFitEastings.at(t), fFitNorthings.at(t), fitLon, fitLat);
      hFittedPos->Fill(fitLon, fitLat);
      hEventPos->Fill(event1.longitude, event1.latitude);

      hFitVsUnfitted->Fill(dist, distFitted);
      hUnfitMinusFit->Fill(dist - distFitted);
    }
    p.inc(eventInd, events.size());
  }

  grAnita->SetName("grAnita");
  grAnita->Write();
  delete grAnita;

  hFittedPos->Write();
  delete hFittedPos;

  hEventPos->Write();
  delete hEventPos;

  hFitVsSumOfSelfLLs->Write();
  delete hFitVsSumOfSelfLLs;

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
  // move swapping global error level stuff out of parallelized loop
  // this used to be around the minimizer->minimize()
  fROOTgErrorIgnoreLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = 1001;

  TStopwatch timer;
  timer.Start(kTRUE);

  Event* event1 = nextEvent();
  Int_t loopCount = 0;

  while(event1){
    // update number done...
    UInt_t numEventsProcessed = 0;
    for(UInt_t eventInd=0; eventInd < events.size(); eventInd++){
      if(events.at(eventInd).cluster[0] > -1){
        numEventsProcessed++;
      }
    }
    std::cerr << "Have found a cluster for " << numEventsProcessed << " of " << events.size()
      << " events, " <<  (events.size() - numEventsProcessed) << " remaining.\n";
    std::cout << "Starting loop = " << loopCount << std::endl;

    if(event1->eventEventClustering){

      // first look up events close by in Easting/Northing
      std::vector<Int_t> nearbyEventsInds;
      double lookup[2] = {event1->easting, event1->northing};
      fKDTree->FindInRange(lookup, default_range_easting_northing, nearbyEventsInds);



      // prepare parallelized storage things
      const int mt = OpenMP::getMaxThreads();
      std::vector<std::vector<UInt_t> > event2Inds[mt]; // [t][z].size()==numMatchedClusters
      std::vector<std::vector<Int_t> > matchedClusters[mt]; // [t][z].size()==numMatchedClusters
      for(int t=0; t < mt; t++){
        matchedClusters[t] = std::vector<std::vector<Int_t> >(event1->nThresholds, std::vector<Int_t>());
        event2Inds[t] = std::vector<std::vector<UInt_t> >(event1->nThresholds, std::vector<UInt_t>());
      }



      // Add event1 event to the list of "matched clusters"
      for(int z=0; z < event1->nThresholds; z++){
        if(event1->cluster[z] > -1){
          if(z==0){
            std::cerr << "I think this is false for all events and this statement should never get printed..." << std::endl;
            std::cerr << event1->eventNumber << "\t" << event1->cluster[0] << std::endl;
          }
          for(int t=0; t < mt; t++){
            matchedClusters[t][z].push_back(event1->cluster[z]);
          }
        }
      }


      std::map<std::pair<Int_t, Int_t>, Double_t> bestUnfittedLogLikelihoodThisClusterThisHistBin[mt];

      // Loop over nearby events in easting/northing
      UInt_t n2 = OpenMP::isEnabled ? ceil(double(nearbyEventsInds.size())/mt) : nearbyEventsInds.size();
      Int_t innerLoopCounter = 0;
      std::cerr <<  "There are " << nearbyEventsInds.size() << " nearby events to process." << std::endl;
      ProgressBar p2(n2);

#pragma omp parallel for schedule(dynamic, 1)
      for(UInt_t i=0; i < nearbyEventsInds.size(); i++){
        int event2Ind = nearbyEventsInds[i];
        Event& event2 = events.at(event2Ind);

        int t = OpenMP::thread();

        if(event1->eventNumber != event2.eventNumber && event2.eventEventClustering){

          // you must check against this if event2 isn't already in a cluster
          // or event2 isn't in a cluster that you've already matched with at the lowest threshold

          // i.e. ONLY don't do an event if you've already matched against
          // another event from the same cluster at the lowest threshold
          // if(event2.cluster[0] < 0 || std::find(matchedClusters[t][0].begin(), matchedClusters[t][0].end(), event2.cluster[0])==matchedClusters[t][0].end())

          // If the clustering is working like it should

          if(event2.cluster[0] < 0 || !RootTools::vectorContainsValue(matchedClusters[t][0], event2.cluster[0])){

            Double_t ll = dMin(event1, &event2);
            Double_t surfaceDist = 1e-3*event1->cartesianSeparation(event2);

            // Double_t llAsym = dAsym(&event2, event1) + event1->selfLogLikelihood;

            // add extra condition here?
            // if that's smaller than the minLL for event2's cluster on that hist-bin...
            // then do the fit, otherwise, no	    
            // std::pair<Int_t, Int_t> key = std::make_pair(event2.cluster[0], event2.antarcticaHistBin);
            // std::map<std::pair<Int_t ,Int_t>, Double_t>::iterator it = bestUnfittedLogLikelihoodThisClusterThisHistBin[t].find(key);
            // bool tryTheFit = false;
            // if(it==bestUnfittedLogLikelihoodThisClusterThisHistBin[t].end()){
            //   bestUnfittedLogLikelihoodThisClusterThisHistBin[t][key] = llAsym;
            //   tryTheFit = true;
            // }
            // else if(it!=bestUnfittedLogLikelihoodThisClusterThisHistBin[t].end() && llAsym < it->second){
            //   it->second = llAsym;
            //   tryTheFit = true;
            // }	    
            // if(tryTheFit && llAsym > llFitThreshold)
            if(ll > llFitThreshold && surfaceDist > surfaceDistThresholdKm){
              ll = dFit(event1, &event2);
            }

            if(ll < event1->nearestEventSurfaceLogLikelihood){
              event1->nearestEventSurfaceLogLikelihood = ll;
              event1->nearestEventSurfaceLLEventNumber = event2.eventNumber;
            }

            if(surfaceDist < event1->nearestEventSurfaceDistanceKm){
              event1->nearestEventSurfaceDistanceKm = surfaceDist;
              event1->nearestEventSurfaceEventNumber = event2.eventNumber;
            }
            
            if(ll < event2.nearestEventSurfaceLogLikelihood){
              event2.nearestEventSurfaceLogLikelihood = ll;
              event2.nearestEventSurfaceLLEventNumber = event1->eventNumber;
            }

            if(surfaceDist < event2.nearestEventSurfaceDistanceKm){
              event2.nearestEventSurfaceDistanceKm = surfaceDist;
              event2.nearestEventSurfaceEventNumber = event1->eventNumber;
            }

            for(int z=0; z < event1->nThresholds; z++){
              if(surfaceDist < surfaceDistThresholdKm || ll <= llEventCuts.at(z)){
                event2Inds[t][z].push_back(event2Ind); // the events to merge
                if(event2.cluster[z] > -1){ // the clusters to merge
                  matchedClusters[t][z].push_back(event2.cluster[z]);
                }
              }
            }
          }
        }
        if(t==0){
          p2.inc(innerLoopCounter);
          innerLoopCounter++;
        }
      } // omp loop
      p2.inc(n2, n2); // finish the inner loop progress bar by hand, just in case I got the guestimate of the number of threads wrong

      if(OpenMP::isEnabled){
        for(int t=1; t < mt; t++){
          for(int z=0; z < event1->nThresholds; z++){
            for(UInt_t i=0; i < matchedClusters[t][z].size(); i++){
              if(!RootTools::vectorContainsValue(matchedClusters[0][z], matchedClusters[t][z][i])){
                matchedClusters[0][z].push_back(matchedClusters[t][z][i]);
              }
            }
            for(UInt_t i=0; i < event2Inds[t][z].size(); i++){
              event2Inds[0][z].push_back(event2Inds[t][z].at(i));
            }
            matchedClusters[t][z].clear();
            event2Inds[t][z].clear();
          }
        }
      }

      // Now look at the clusters of the matched events (for each threshold)
      for(UInt_t z=0; z < llEventCuts.size(); z++){

        // for speed, sort the matched event numbers
        std::sort(event2Inds[0][z].begin(), event2Inds[0][z].end());


        // Int_t canICount = 0;
        // for(UInt_t clusterInd=0; clusterInd < clusters.at(z).size(); clusterInd++){
        //   const Cluster& cluster = clusters.at(z).at(clusterInd);
        //   canICount += cluster.numDataEvents;
        // }

        // If none of the events are in a cluster, then add a new cluster!
        if(matchedClusters[0][z].size()==0){

          // make a new cluster, it's intial position is this event location
          // although that will be changed...
          Cluster nc(*event1, clusters.at(z).size());
          nc.llEventCutInd = z;
          nc.llEventCut = llEventCuts.at(z);
          clusters.at(z).push_back(nc);
          matchedClusters[0][z].push_back(nc.index);
        }

        /**
         * If the list of matched clusters is greater than one then we're going to merge clusters!
         * But what to call the new cluster? We'll label it with the smallest index.
         * this ensures that any merged base/non-base clusters will be considered part of a base,
         * since the bases populate the start of the list.
         */

        Int_t thisCluster = TMath::MinElement(matchedClusters[0][z].size(), &matchedClusters[0][z][0]);

        // Mark event1 as in the minCluster
        Int_t oldCluster1Ind = event1->cluster[z];
        if(oldCluster1Ind > -1){
          if(z==0){
            std::cerr << "Do we ever get here?" << std::endl; // no we don't
          }
          clusters.at(z).at(oldCluster1Ind).numDataEvents--;
        }
        event1->cluster[z] = thisCluster;
        clusters.at(z).at(thisCluster).numDataEvents++;

        // std::cerr << event2Inds.size() << std::endl;

        // Now reassign all matched events or events in matched clusters to the new cluster
        UInt_t matchedEventIndex=0;
        for(UInt_t event2Ind=0; event2Ind < events.size(); event2Ind++){
          Event& event2 = events.at(event2Ind);

          // check if this event is contained in the matched events or the matched clusters
          bool eventMatch = matchedEventIndex < event2Inds[0][z].size() && event2Ind==event2Inds[0][z][matchedEventIndex];
          if(eventMatch){
            // std::cout << matchedEventIndex << "\t" << event2Ind << "\t" << event2Inds[0][z][matchedEventIndex] << "\t" << event2Inds[0][z].size() << std::endl;
            matchedEventIndex++;
          }

          if(eventMatch || RootTools::vectorContainsValue(matchedClusters[0][z], event2.cluster[z])){

            Int_t oldCluster2Ind = event2.cluster[z];
            if(oldCluster2Ind > -1){
              clusters.at(z).at(oldCluster2Ind).numDataEvents--;
            }
            event2.cluster[z] = thisCluster;
            clusters.at(z).at(thisCluster).numDataEvents++;
          }
        }

        // Int_t canICount2 = 0;
        // for(UInt_t clusterInd=0; clusterInd < clusters.at(z).size(); clusterInd++){
        //   const Cluster& cluster = clusters.at(z).at(clusterInd);
        //   canICount2 += cluster.numDataEvents;
        // }

        // std::cout << z << "\t" << canICount << "\t" << canICount2 << std::endl;
      }
    }

    Int_t seconds = Int_t(timer.RealTime());
    Int_t hours = seconds / 3600;
    hours = hours < 0 ? 0 : hours;
    seconds = seconds - hours * 3600;
    Int_t mins = seconds / 60;
    mins = mins < 0 ? 0 : mins;
    seconds = seconds - mins * 60;

    std::cerr << "Finished loop = " << loopCount << ", in ";
    fprintf(stderr, ANSI_COLOR_MAGENTA); // because, why not?
    fprintf(stderr, "%02d:%02d:%02d\n", hours, mins, seconds);
    fprintf(stderr, ANSI_COLOR_RESET);
    timer.Start(kFALSE);    

    loopCount++;
    event1 = nextEvent();
  }

  // Set fitter error level back to default
  gErrorIgnoreLevel = fROOTgErrorIgnoreLevel;

}






/**
 * An implementation of the event-to-event clustering with worst case O(N^2) efficiency
 */
void Acclaim::Clustering::LogLikelihoodMethod::doMcEventClustering(){

  double llFitThreshold = TMath::MinElement(llEventCuts.size(), &llEventCuts.at(0)); // fit if greater than this
  double llNearbyThreshold = TMath::MaxElement(llEventCuts.size(), &llEventCuts.at(0)); // ignore if greater than this

  std::cout << __PRETTY_FUNCTION__ << std::endl;
  std::cout << "llNearbyThreshold = " << llNearbyThreshold << ", llFitThreshold = " << llFitThreshold << std::endl;

  fROOTgErrorIgnoreLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = 1001;

  ProgressBar p(mcEvents.size());
  
  for(UInt_t event1Ind=0; event1Ind < mcEvents.size(); event1Ind++) {
  
    McEvent& event1 = mcEvents.at(event1Ind);
    
    if (event1.eventEventClustering) {

      double lookup[2] = {event1.easting, event1.northing};
      // look up nearby DATA events, not MC, don't want to cluster MC to MC...
      std::vector<Int_t> event2Inds;
      std::vector<Double_t> event2EastingNorthingDistances;
      UInt_t lastNumNeighbours = 0;
      UInt_t numNeighbours = 2048;
      Double_t furthestConsidered = 0;
      Int_t numConsidered = 0;

      while (furthestConsidered < default_horizon_distance && event1.cluster[0] < 0) {
      
        event2Inds.resize(numNeighbours, -1);
        event2EastingNorthingDistances.resize(numNeighbours, -1);
        fKDTree->FindNearestNeighbors(lookup, numNeighbours, & event2Inds[0], & event2EastingNorthingDistances[0]);

        for (UInt_t i=lastNumNeighbours; i < event2Inds.size() && event1.cluster[0] < 0 && furthestConsidered < default_horizon_distance; i++) {
        
          UInt_t event2Ind = event2Inds.at(i);
          const Event& event2 = events.at(event2Ind);

          if (event2EastingNorthingDistances[i] > furthestConsidered) furthestConsidered = event2EastingNorthingDistances[i];

          double ll = dMin(& event1, & event2);
          Double_t surfaceDist = 1e-3 * event1.cartesianSeparation(event2);
          
          if (ll > llFitThreshold && surfaceDist > surfaceDistThresholdKm) ll = dFit(& event1, & event2);

          if (ll < event1.nearestEventSurfaceLogLikelihood) {
          
            event1.nearestEventSurfaceLogLikelihood = ll;
            event1.nearestEventSurfaceLLEventNumber = event2.eventNumber;
          }

          if (surfaceDist < event1.nearestEventSurfaceDistanceKm) {
          
            event1.nearestEventSurfaceDistanceKm = surfaceDist;
            event1.nearestEventSurfaceEventNumber = event2.eventNumber;
          }

          for (int z=0; z < event1.nThresholds; z++) {
          
            if (surfaceDist < surfaceDistThresholdKm || ll < llEventCuts.at(z)) {
            
              double eventWeight = event1.weight;
              event1.cluster[z] = event2.cluster[z];
              clusters[z][event1.cluster[z]].sumMcWeights += eventWeight;  //  Comparing to doMcBaseClustering, this should fill iceMC event weights where neccessary.             
            }
          }
          
          numConsidered++;
        }
        
        lastNumNeighbours = numNeighbours;
        numNeighbours *= 2;
        
        if (numNeighbours >= events.size()) continue;
      }
      
      const char* prefix = event1.cluster[0] < 0 ? "Did not find" : "Found";
      std::cerr << prefix << " a match after " << numConsidered << " events" << std::endl; 
    }
    
    p.inc(event1Ind);
  }

  gErrorIgnoreLevel = fROOTgErrorIgnoreLevel;
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
  const int nBases = BaseList::getNumBases();

  for(UInt_t eventInd=0; eventInd < mcEvents.size(); eventInd++){
    McEvent* mcEvent = &mcEvents.at(eventInd);
    for(int clusterInd=0; clusterInd < nBases; clusterInd++){
      Cluster& cluster = clusters.at(0).at(clusterInd);
      if(cluster.knownBase){
        double distM = mcEvent->usefulPat.getDistanceFromSource(cluster.latitude, cluster.longitude, cluster.latitude);
        if(distM < default_horizon_distance){
          double ll = mcEvent->logLikelihoodFromPoint(cluster.longitude, cluster.latitude, cluster.altitude, true);
          double surfaceSeparationKm = 1e-3*mcEvent->cartesianSeparation(cluster);

          if(surfaceSeparationKm < mcEvent->nearestKnownBaseSurfaceSeparationKm){
            mcEvent->nearestKnownBaseSurfaceSeparationKm = surfaceSeparationKm;

            if(mcEvent->nearestKnownBaseSurfaceSeparationKm < surfaceDistThresholdKm){
              mcEvent->nearestKnownBaseCluster = clusterInd;
            }
          }

          if(ll < mcEvent->nearestKnownBaseLogLikelihood && mcEvent->nearestKnownBaseSurfaceSeparationKm >= surfaceDistThresholdKm){
            mcEvent->nearestKnownBaseCluster = clusterInd;
          }

          for(int z=0; z < mcEvent->nThresholds; z++){
            if(mcEvent->cluster[z] < 0 && (ll < llEventCuts.at(z) || mcEvent->nearestKnownBaseSurfaceSeparationKm < surfaceDistThresholdKm)){
              clusters[z][clusterInd].sumMcWeights += mcEvent->weight;
              mcEvent->cluster[z] = clusterInd;
            }
          }
        }
      }
    }
    p.inc(eventInd);
  }
}



void Acclaim::Clustering::LogLikelihoodMethod::doClustering(const char* dataGlob, const char* mcGlob, const char* outFileName, bool useAcclaimFiles){

  useAcclaimFiles ? readInSummaries(dataGlob) : readInSampleSummaries(dataGlob, false);
  useAcclaimFiles ? readInSummaries(mcGlob) : readInSampleSummaries(mcGlob, true);
//  useAcclaimFiles ? readInSummaries(dataGlob) : readInTMVATreeSummaries(dataGlob, 0);
//  useAcclaimFiles ? readInSummaries(mcGlob) : readInTMVATreeSummaries(mcGlob, 1);
  std::cout << "Sorting events...";
  std::sort(events.begin(), events.end());
  std::cout << "done" << std::endl;

  const char* fakeArgv[1] = {outFileName};
  TFile* fOut = 0;
  OutputConvention oc(1, const_cast<char**>(fakeArgv));
  fOut = useAcclaimFiles ? oc.makeFile() : new TFile(outFileName, "RECREATE");

  initKDTree();
  // fEventsAlreadyClustered = false;

  if(fUseBaseList){
  
    if (!fEventsAlreadyClustered) {
    
      readInBaseList();
      doBaseEventClustering();
    }
    
    doMcBaseClustering();
  }

  if(!fEventsAlreadyClustered) doEventEventClustering();
  
  doMcEventClustering();

  makeSummaryTrees();

  // testSmallClustersFromPointSource();

  fOut->Write();
  fOut->Close();
  return;


}

void Acclaim::Clustering::LogLikelihoodMethod::testSmallClusters(const char* dataGlob, const char* outFileName, int clusterSizeMin, int clusterSizeMax, int nAttempts){ 

  TFile* fOut = 0;
  fOut = new TFile(outFileName, "RECREATE");
  

  TH2D* h = new TH2D("LLDist", "LLDist", 200, 0, 100, 200, 0, 100);
  h->GetXaxis()->SetTitle("nearest event distance");
  h->GetYaxis()->SetTitle("nearest event LL");
  TH2D* h2 = new TH2D("clusterStuff", "clusterStuff", 100, 0, 100, 100, -.005, .995);
  h2->GetXaxis()->SetTitle("LL threshold");
  h2->GetYaxis()->SetTitle("percent singlets");

  readInSampleSummariesForTesting(dataGlob);
//  readInSummariesForTesting(dataGlob);

//  std::cout << "Is this function called by default? No, it's not." << std::endl;

  for(int i = 0; i < nAttempts; i++)
  {
    Int_t n_in_cluster = tr3->Uniform(clusterSizeMin, clusterSizeMax);
    pickSampleEventsFromList(n_in_cluster);
//    pickEventsFromList(n_in_cluster);
    std::sort(events.begin(), events.end());
    initKDTree();
    doEventEventClustering();
    addToHistograms(h, h2);
    
    for(int j = 0; j < clusters.size(); j++) {
      clusters.at(j).clear();
    }

    events.clear();
  }

  fOut->cd();
  h->Write();
  h2->Write();
  fOut->Write();
  fOut->Close();
  return;
}
