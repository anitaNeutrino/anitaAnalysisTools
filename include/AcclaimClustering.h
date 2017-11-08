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
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "AcclaimOpenMP.h"

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
    void getAngularResolution(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd, double& sigma_theta, double& sigma_phi);
    void getAngularResolution(double snr, double& sigma_theta, double& sigma_phi);
    TCanvas* drawAngularResolutionModel(double maxSnr);

    /**
     * @class Event
     * @brief Minimum required information about an ANITA event to be clustered
     */
    class Event{
    public:

      //--------------------------------------------------------------------------------
      // determined by reconstruction
      //--------------------------------------------------------------------------------	
      UInt_t eventNumber;			/// Event number
      Int_t run;				/// Run
      AnitaPol::AnitaPol_t pol;			/// Polarization
      Int_t peakIndex;				/// Which peak in the map does this represent?

      Double_t centre[3];//!			/// Cartesian coordinates, does not persist in ROOT
      Double_t latitude;			/// latitude
      Double_t longitude;			/// longitude
      Double_t altitude;			/// longitude
      Double_t easting; //!                     /// easting
      Double_t northing;//!                     /// northing
      AnitaEventSummary::PayloadLocation anita;	/// Anita's position

      Double_t theta;				/// reconstructed theta
      Double_t phi;				/// reconstructed phi

      Double_t sigmaTheta;			/// resolution associated with this snr?
      Double_t sigmaPhi;			/// resolution associated with this snr?

      //--------------------------------------------------------------------------------
      // determined by clustering
      //--------------------------------------------------------------------------------
      Double_t dThetaCluster;			/// theta distance to cluster
      Double_t dPhiCluster;			/// phi distance to cluster

      Double_t logLikelihood;			/// log likelihood = -2 * ln (...)
      Int_t cluster;				/// which cluster am I associated with?

      Double_t logLikelihood2;			/// log likelihood to second closest cluster
      Int_t cluster2;              		/// what cluster am I second closest to?
      
      UInt_t numNeighbours;                     /// How many neighbouring events were found within queryRegion?

      UInt_t nearestNeighbourEventNumber;       /// What is the eventNumber of the event am I closest to?
      Double_t nearestNeighbourLogLikelihood;   /// And what is the log likelihood?

      Int_t antarcticaHistBin; //!		/// Which global bin in the TH2DAntarctica?
      mutable UsefulAdu5Pat usefulPat; //!      /// Only construct this once, mutable since not const correct
      Event();
      Event(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd);
      TArrowAntarctica* makeArrowFromAnitaToEvent();
      void setupUsefulPat();

      virtual ~Event(){ ;}
      ClassDef(Event, 4)
    };



    /** 
     * @class McEvent
     * @brief Same as event, but with an energy and a weight
     */
    class McEvent : public Event{
    public:
      Double_t weight;						/// MC information
      Double_t energy;						/// MC information

      McEvent();
      McEvent(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd)
	: Event(sum, pol,  peakInd)
      {
	weight = sum->mc.weight;
	energy = sum->mc.energy;
      }

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
      Cluster();
      // Cluster(const Event& seedEvent);
      Cluster(const BaseList::base& base);

      virtual ~Cluster(){ ;}	

      Double_t centre[3]; //! Does not persist			/// Center in cartesian
	
      Double_t latitude;					/// cluster center latitude
      Double_t longitude;					/// cluster center longitude
      Double_t altitude;					/// cluster center altitude
	
      Int_t numDataEvents;					/// How many data events does this cluster contain?
      Double_t sumMcWeights;					/// How many MC events does this cluster contain?
      Int_t knownBase;						/// Known base == 0, Pseudo-base == 1

      Int_t antarcticaHistBin; //!				/// Which global bin in the TH2DAntarctica?
      Int_t seedEvent;                   //!			/// Which event seeded the cluster?

      ClassDef(Cluster, 4)
    };


    /**
     * @class DistIndex
     * Slightly over-engineered pair, with some default values
     * which allows sorting of eventIndices based on separation
     */
    class DistIndex {
    public:
      DistIndex() : d(-1), i(-1) {;}    
      DistIndex(Int_t eventInd) : d(-1), i(eventInd) {;}
      DistIndex(Double_t distance, Int_t eventInd) : d(distance), i(eventInd) {;}

      friend bool operator<(const DistIndex &di1, const DistIndex &di2){
	return di1.d < di2.d;
      }
      friend bool operator<=(const DistIndex &di1, const DistIndex &di2){
	return di1.d <= di2.d;
      }
      friend bool operator>(const DistIndex &di1, const DistIndex &di2){
	return di1.d > di2.d;
      }
      friend bool operator>=(const DistIndex &di1, const DistIndex &di2){
	return di1.d >= di2.d;
      }
      friend bool operator==(const DistIndex &di1, const DistIndex &di2){
	return di1.d == di2.d;
      }
      friend bool operator!=(const DistIndex &di1, const DistIndex &di2){
	return di1.d != di2.d;
      }
      Double_t d;
      Int_t i;
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

      void doClustering(const char* dataGlob, const char* mcGlob, const char* outFileName);
      void doClusteringDBSCAN(const char* dataGlob, const char* mcGlob, const char* outFileName);
      
      TGraphAntarctica* makeClusterSummaryTGraph(Int_t clusterInd);

      bool getUseBaseList(){return fUseBaseList;}
      void setUseBaseList(bool useBaseList){ // *TOGGLE *GETTER=GetUseBaseList
	fUseBaseList = useBaseList;
      }

      bool getDebug(){return fDebug;}
      void setDebug(bool db){fDebug = db;} // *TOGGLE *GETTER=GetDebug

    private:
      Double_t getDistSqEventCluster(Int_t eventInd, const Acclaim::Clustering::Cluster& cluster);
      // Double_t getDistSqClusterCluster(Int_t clusterInd1, Int_t clusterInd2);
      Double_t getAngDistSqEventCluster(Int_t eventInd, Int_t clusterInd);
      void getDeltaThetaDegDeltaPhiDegEventCluster(Int_t eventInd, Int_t clusterInd, Double_t& deltaThetaDeg, Double_t& deltaPhiDeg);

      Long64_t readInSummaries(const char* summaryGlob);
      size_t addEvent(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd);
      size_t addMcEvent(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd);
      void assignSingleEventToCloserCluster(Int_t pointInd, Int_t isMC, Int_t clusterInd);
      void initializeBaseList();
      void initializeEmptyBaseList();
      Int_t histogramUnclusteredEvents(Int_t& globalMaxBin);
      void recursivelyAddClustersFromData(Int_t minBinContent);
      void assignMcEventsToClusters();
      void assignEventsToBaseClusters();
      void writeAllGraphsAndHists();
      void findClosestEventToClustersOfSizeOne();
      void makeSummaryTrees();
      void resetClusters();
      Double_t getSumOfMcWeights();
      void redoSmallClusters();
      void initKDTree();
      Double_t dMin(const Event& event1, const Event& event2);
      Double_t dSum(const Event& event1, const Event& event2);
      Double_t dAsym(const Event& event1, const Event& event2);
      void testTriangleInequality();
      void sortEventIndices(Int_t eventInd, std::vector<Int_t>& unsorted, std::vector<DistIndex>& sorted);
      bool isCorePointDBSCAN(Event& event);
      void testAngleFindingSpeed();
      
      void DBSCAN();
      void OPTICS();
      const UInt_t fMinPts = 4;
      void rangeQueryEastingNorthing(Int_t eventInd, Double_t range, std::vector<Int_t>& seed);
      void makeAndWriteNSquaredEventEventHistograms();
      Double_t dPoint(const Event& eventInd1, Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt, bool addOverHorizonPenalty=false);
      Double_t evalPairLogLikelihoodAtLonLat(const Double_t* params);
      const Event* fFitEvent1; /// First event in the pairwise fit
      const Event* fFitEvent2; /// Second event in the pairwise fit
      Int_t fMaxFitterAttempts; /// How many times should I try if I don't reach a good minimum?
      Double_t fFitHorizonDistM; ///700e3 metres, distance at which a penalty is added to source location fitting
      Double_t dFit(const Event& event1, const Event& event2);
      UInt_t fTestEvent1; /// For debugging
      UInt_t fTestEvent2; /// For debugging
      Double_t fFitEasting;  /// Where the fitter found the potential source could come from
      Double_t fFitNorthing; /// Where the fitter found the potential source could come from
      TH2DAntarctica* makeUnevenlyBinnedEventHistogram();

      
      Double_t llCut;						/// The cut-off for log-likelihood, which defines the boundary of a cluster
      Double_t maxDistCluster;					/// Only consider adding to a cluster when closer than this value
      Bool_t doneBaseClusterAssignment;				/// Set to true once all read in data events were clustered to bases

      // Int_t numClusters;
      Int_t numCallsToRecursive;

      std::vector<Acclaim::Clustering::Cluster> clusters;	/// Vector of clusters,
      std::vector<Acclaim::Clustering::Event> events;		/// Vector of data events
      std::vector<Acclaim::Clustering::McEvent> mcEvents;	/// Vector of Monte Carlo events
      TKDTreeID* fKDTree;                                       /// ROOT's implementation of a KDTree, typedef'd for int/double

      std::vector<TGraphAntarctica*> grBaseClusterCenters;	/// The locations of the bases
      std::vector<TH2DAntarctica*> hBaseClusteredEvents;	/// Histograms of events clustered to bases
      std::vector<TGraphAntarctica*> grNonBaseClusterCenters;	/// The locations of the non-base clusters
      std::vector<TH2DAntarctica*> hNonBaseClusteredEvents;	/// Histograms of events clustered to non-base clusters
      TH2DAntarctica* hClusters;                                /// Filled with clusters (allows access to the bin of the cluster)
      TH2DAntarctica* hEvents;
      bool fDebug;
      bool fUseBaseList;
      TGraph* grTestMinimizerWalk;

      ROOT::Math::Minimizer* fMinimizer;
      ROOT::Math::Functor fFunctor;
    };


  }
}

#endif
