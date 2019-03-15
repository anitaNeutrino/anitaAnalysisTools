/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             A class to cluster. What were you expecting?
***********************************************************************************************************/

#ifndef ACCLAIM_CLUSTERING_H
#define ACCLAIM_CLUSTERING_H

#include "UsefulAdu5Pat.h"
#include "AnitaConventions.h"
#include "BaseList.h"
#include "AnitaEventSummary.h"
#include "TKDTree.h"
#include "TCut.h"
#include "TEntryList.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "AcclaimOpenMP.h"
#include "TRandom3.h"
#include "TChain.h"

class TTree;
class TGraphAntarctica;
class TFile;
class TH2DAntarctica;
class TCanvas;
class TArrowAntarctica;


namespace Acclaim{

  namespace Clustering {

    const Double_t default_sigma_theta = 0.25;
    const Double_t default_sigma_phi = 0.5;
    const Double_t default_var_theta = default_sigma_theta * default_sigma_theta;
    const Double_t default_var_phi = default_sigma_phi * default_sigma_phi;
    const Double_t default_range_easting_northing = 750e3;
    const Double_t default_horizon_distance = 750e3;

    void getAngularResolution(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd, double& sigma_theta, double& sigma_phi);
    void getAngularResolution(double snr, double& sigma_theta, double& sigma_phi);
    TCanvas* drawAngularResolutionModel(double maxSnr);

    void getAngularVariance(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd, double& var_theta, double& var_phi);
    void getAngularVariance(double snr, double& var_theta, double& var_phi);
    TCanvas* drawAngularVarianceModel(double maxSnr);

    template <class T>
    Bool_t isVaguelyNearMcMurdo(const T& t){
      return t.longitude >= 150 && t.longitude < 180 && t.latitude >= -85 && t.latitude < -75;
    }
    template <class T>
    Bool_t inSandbox(const T& t){
      return t.longitude >= -90 && t.longitude < 90;
    }

    /**
     * @class Event
     * @brief Minimum required information about an ANITA event to be clustered
     */
    class Event{
    public:

      //--------------------------------------------------------------------------------
      // determined by reconstruction
      //--------------------------------------------------------------------------------
      UInt_t eventNumber;				/// Event number
      Int_t run;					/// Run
      AnitaPol::AnitaPol_t pol;				/// Polarization
      Int_t peakIndex;					/// Which peak in the map does this represent?

      Double_t centre[3];//!				/// Cartesian coordinates, does not persist in ROOT
      Double_t latitude;				/// latitude
      Double_t longitude;				/// longitude
      Double_t altitude;				/// longitude
      Double_t easting;					/// easting
      Double_t northing;				/// northing
      AnitaEventSummary::PayloadLocation anita;		/// Anita's position

      Double_t theta;					/// reconstructed theta
      Double_t phi;					/// reconstructed phi
      Double_t thetaAdjustmentRequired;			/// the adjustment from traceBackToContinent

      Double_t sigmaTheta;				/// resolution associated with this snr?
      Double_t sigmaPhi;				/// resolution associated with this snr?

      Double_t varTheta;				/// variance associated with this snr?
      Double_t varPhi;				        /// (spherical) variance associated with this snr?

      //--------------------------------------------------------------------------------
      // determined by clustering
      //--------------------------------------------------------------------------------
      Int_t nThresholds;                        
      Int_t* cluster;//[nThresholds]			/// which cluster am I associated with?
      Double_t* dThetaCluster;//[nThresholds]		/// theta distance to cluster
      Double_t* dPhiCluster;//[nThresholds]		/// phi distance to cluster

      Bool_t eventEventClustering;			/// Remove huge clusters near MCM before doing event-to-event clustering
      Double_t nearestKnownBaseLogLikelihood;		/// How far to the nearest known base?
      Double_t nearestKnownBaseSurfaceSeparationKm;	/// How far to the nearest known base?
      Int_t nearestKnownBaseCluster;			/// How far to the nearest known base?
      Int_t nearestKnownBaseClusterSurface;		/// How far to the nearest known base?
      Double_t selfLogLikelihood;			/// If the event is above the continent surface, this may be non-zero
      Double_t nearestEventSurfaceDistanceKm;		/// How far away to the nearest event, in kilometers?
      UInt_t nearestEventSurfaceEventNumber;		/// What's the event number of the nearest surface neighbour?
      Double_t nearestEventSurfaceLogLikelihood;	/// What's the fitted log likelihood to the nearest surface neighbour?
      UInt_t nearestEventSurfaceLLEventNumber;		/// What's the event number of the nearest surface neighbour by LL?

      Int_t antarcticaHistBin;				/// Which global bin in the TH2DAntarctica?
      UsefulAdu5Pat usefulPat; //!			/// Only construct this once
      mutable Bool_t fDebug; //!

      Event(Int_t nT=0);
      Event(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd, Int_t nT=0);
      Event(const Event& event);
      Event& operator=(const Event& event);
      Event(int pol, int peakInd, double peak_phi, double peak_theta, int nT, UInt_t eventNumber, Int_t run,
	    double anita_longitude, double anita_latitude, double anita_altitude, double anita_heading, double coherent_filtered_snr);


      TArrowAntarctica* makeArrowFromAnitaToEvent();
      void setupUsefulPat(bool calculateNow = true);
      void resetClusteringNumbers();
      void deleteArrays();
      void setNThresholds(int n);

      double logLikelihoodFromPoint(double sourceLon, double sourceLat, double sourceAlt, bool addOverHorizonPenalty=false) const;
      template<class T>
      double logLikelihoodFromPoint(const T& point, bool addOverHorizonPenalty=false) const {
	return logLikelihoodFromPoint(point.longitude, point.latitude, point.altitude, addOverHorizonPenalty);
      }
      template<class T>
      inline double cartesianSeparation(const T& event2){
	double d0 = centre[0] - event2.centre[0];
	double d1 = centre[1] - event2.centre[1];
	double d2 = centre[2] - event2.centre[2];
	return TMath::Sqrt(d0*d0 + d1*d1 + d2*d2);
      }

      /** 
       * For sorting based on event number
       * @param other is another event
       * 
       * @return true if this event has a lower eventNumber
       */
      inline bool operator < (const Event& other) const {
        return (eventNumber < other.eventNumber);
      }

      virtual ~Event();
      ClassDef(Event, 14)
    };



    /** 
     * @class McEvent
     * @brief Same as event, but with an energy and a weight
     */
    class McEvent : public Event{
    public:
      Double_t weight;						/// MC information
      Double_t energy;						/// MC information

      McEvent(Int_t nT = 0);
      McEvent(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd, Int_t nT=0);
      McEvent(double weight, double energy, int pol, int peakInd, double peak_phi, double peak_theta, int nT, UInt_t eventNumber, Int_t run,
	      double anita_longitude, double anita_latitude, double anita_altitude, double anita_heading, double coherent_filtered_snr);

      virtual ~McEvent(){;}
	
      ClassDef(McEvent, 1)
    };    





    //--------------------------------------------------------------------------------------------------------
    /**
     * @class Cluster
     * @brief A position on the continent, with which a bunch of events are associated
     */
    class Cluster{
    public:
      Cluster(Int_t i=-1);
      Cluster(const Event& seedEvent, Int_t i=-1);
      Cluster(const BaseList::base& base, Int_t i=-1);

      virtual ~Cluster(){ ;}	

      Double_t centre[3]; //! Does not persist			/// Center in cartesian
	
      Double_t latitude;					/// cluster center latitude
      Double_t longitude;					/// cluster center longitude
      Double_t altitude;					/// cluster center altitude
	
      Int_t numDataEvents;					/// How many data events does this cluster contain?
      Double_t sumMcWeights;					/// How many MC events does this cluster contain?
      Int_t knownBase;						/// Known base == 0, Pseudo-base == 1
      Int_t index;                                              /// Where am I in the cluster array?
      Int_t llEventCutInd; /// which entry in the llEventCut array does this correspond to?
      Double_t llEventCut; /// and what's the value of that cut?

      Int_t antarcticaHistBin; //!				/// Which global bin in the TH2DAntarctica?
      Int_t seedEvent; //!			                /// Which event seeded the cluster?

      void resetClusteringNumbers();

      ClassDef(Cluster, 5)
    };

    
    /**
     * @class LogLikelihoodMethod
     * @brief Workhorse class: groups events into clusters around bases and other events using the log likelihood method
     *
     */
    class LogLikelihoodMethod {
    public:

      static const Int_t SmallClusterSizeThreshold = 100;

      LogLikelihoodMethod();
      virtual ~LogLikelihoodMethod();

      void doClustering(const char* dataGlob, const char* mcGlob, const char* outFileName, bool useAcclaimFiles=true);
      void testSmallClusters(const char* dataGlob, const char* outFileName, int clusterSizeMin=5, int clusterSizeMax=16, int nAttempts=100);
      
      bool getUseBaseList(){return fUseBaseList;}
      void setUseBaseList(bool useBaseList){ // *TOGGLE *GETTER=GetUseBaseList
	fUseBaseList = useBaseList;
      }

      void setCut(TCut cut){ fCut = cut; }
      void setCutHical(bool hc){ fCutHical = hc; } // *TOGGLE *GETTER=GetCutHical
      void setSelfLLMax(double llmax){ fSelfLLMax = llmax; } 
      void setSurfaceDistThresholdKm(Double_t dist){ surfaceDistThresholdKm = dist; }
      void setPercentOfMC(Int_t percent){ fPercentOfMC = percent; }

      bool getDebug(){return fDebug;}
      void setDebug(bool db){fDebug = db;} // *TOGGLE *GETTER=GetDebug

      /** 
       * Add a ThermalChain style cut
       * 
       * @return number of cuts to be applied
       */
      size_t addCut(const TCut c){
	fThermalChainCuts.push_back(c);
	return fThermalChainCuts.size();
      }

      bool fStoreUnclusteredHistograms;
      
    private:
      Double_t getDistSqEventCluster(const Event& event, const Cluster& cluster);
      Double_t getAngDistSqEventCluster(const Event& event, const Cluster& cluster);
      void getDeltaThetaDegDeltaPhiDegEventCluster(const Event& event, const Cluster& cluster, Double_t& deltaThetaDeg, Double_t& deltaPhiDeg);


      Long64_t readInSummaries(const char* summaryGlob);
      Long64_t readInTMVATreeSummaries(const char* summaryGlob, bool isMC);
      Long64_t readInSampleTreeSummaries(const char* summaryGlob);
      void readInSummariesForTesting(const char* summaryGlob);
      void readInSampleSummariesForTesting(const char* summaryGlob);
      void pickEventsFromList(int n_in_cluster);
      void pickSampleEventsFromList(int n_in_cluster);
      size_t addEvent(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd);
      size_t addMcEvent(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd);      
      void assignSingleEventToCloserCluster(Int_t eventInd, Int_t isMC, Cluster& cluster, Int_t z, double llEventCut = -1);
      void readInBaseList();

      TRandom3* tr3;

      void assignMcEventsToClusters();
      void setInitialBaseClusters();
      void makeSummaryTrees();
      void addToHistograms(TH2D* h, TH2D* h2);
      void resetClusters();
      Double_t getSumOfMcWeights();
      void initKDTree();
      Acclaim::Clustering::Event* nextEvent();

      Double_t dMin(const Event* event1, const Event* event2);
      Double_t dSum(const Event* event1, const Event* event2);
      Double_t dAsym(const Event* event1, const Event* event2);
      Double_t dFit(const Event* event1, const Event* event2);

      void testTriangleInequality();
      void testSmallClustersFromPointSource();
      Int_t removeLargeBasesNearMcMurdo();

      void doBaseEventClustering();
      void doEventEventClustering();

      void doMcEventClustering();
      void doMcBaseClustering();
      bool considerBin(const Event& event, Int_t bx, Int_t by, double& easting, double& northing);
      void nearbyEvents2(UInt_t eventInd, std::vector<UInt_t>& nearbyEvents);

      void makeAndWriteNSquaredEventEventHistograms();
      Double_t evalPairLogLikelihoodAtLonLat(const Double_t* params);




      std::vector<const Acclaim::Clustering::Event*> fFitEvent1s; /// First event in the pairwise fit
      std::vector<const Acclaim::Clustering::Event*> fFitEvent2s; /// Second event in the pairwise fit
      Int_t fMaxFitterAttempts; /// How many times should I try if I don't reach a good minimum?
      TCut fCut; /// What to cut on when using tmva trees
      bool fCutHical; /// cut hical? (when using tmva trees)
      bool fSelfLLMax; /// basically a way of excluding above horizon from sample
      TEntryList* fEntryList; /// Entry list that gets filled by cut when using tmva trees

      UInt_t fTestEvent1; /// For debugging
      UInt_t fTestEvent2; /// For debugging
      std::vector<Double_t> fFitEastings;  /// Where the fitter found the potential source could come from
      std::vector<Double_t> fFitNorthings; /// Where the fitter found the potential source could come from
      const Int_t numMcDivisions; /// How much of the mC should I read in?
      Int_t mcDivision; // Which of the MC divisions should I read in? (runs from 0 to numMcDivisions-1)

      std::vector<Double_t> llEventCuts;                        /// Try doing a range of llEventCuts at once...
      Double_t surfaceDistThresholdKm;
      Bool_t fEventsAlreadyClustered;

      std::vector<std::vector<Acclaim::Clustering::Cluster> >clusters;	/// Vector of clusters,
      // std::vector<Acclaim::Clustering::Cluster> clusters;	/// Vector of clusters,      
      std::vector<Acclaim::Clustering::Event> events;		/// Vector of data events
      std::vector<Acclaim::Clustering::McEvent> mcEvents;	/// Vector of Monte Carlo events
      std::vector<std::vector<std::vector<UInt_t> > > fLookupEN;/// Event index lookup
      AntarcticaBackground fMyBackground;


      TKDTreeID* fKDTree;                     /// ROOT's implementation of a KDTree, typedef'd for int/double
      std::vector<Double_t> fEventEastings;
      std::vector<Double_t> fEventNorthings;

      std::vector<TH2DAntarctica*> hUnclusteredEvents;	/// Histograms of events clustered to non-base clusters

      bool fDebug;
      bool fUseBaseList;
      TChain* fChain;
      Int_t fPercentOfMC;
      TGraph* grTestMinimizerWalk;
      TGraph* grTestMinimizerValue;
      TGraph* grTestMinimumPosition;
      std::vector<ROOT::Math::Minimizer*> fMinimizers;
      std::vector<ROOT::Math::Functor> fFunctors;
      Int_t fROOTgErrorIgnoreLevel;
      bool fDrawNewNearbyEventsHistograms;
      bool fReadInBaseList;

      std::vector<TCut> fThermalChainCuts;
    };
  }
}

#endif
