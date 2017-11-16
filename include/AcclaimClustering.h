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
    const Double_t default_range_easting_northing = 700e3;

    void getAngularResolution(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd, double& sigma_theta, double& sigma_phi);
    void getAngularResolution(double snr, double& sigma_theta, double& sigma_phi);
    TCanvas* drawAngularResolutionModel(double maxSnr);

    template <class T>
    Bool_t isVaguelyNearMcMurdo(const T& t){
      return t.longitude >= 150 && t.longitude < 180 && t.latitude >= -85 && t.latitude < -75;
    }
    template <class T>
    Bool_t inSandbox(const T& t){
      return t.longitude >= 60 && t.longitude < 90 && t.latitude >= -75 && t.latitude < -65;
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
      Int_t nThresholds;
      Int_t* cluster;//[nThresholds]		/// which cluster am I associated with?
      Double_t* dThetaCluster;//[nThresholds]	/// theta distance to cluster
      Double_t* dPhiCluster;//[nThresholds]     /// phi distance to cluster

      UInt_t nearestNeighbourEventNumber;       /// What is the eventNumber of the event am I closest to?
      Double_t nearestNeighbourLogLikelihood;   /// And what is the log likelihood?
      Bool_t eventEventClustering;              /// Remove huge clusters near MCM before doing event-to-event clustering
      Double_t nearestKnownBaseLogLikelihood;   /// How far to the nearest known base?

      Int_t antarcticaHistBin; //!		/// Which global bin in the TH2DAntarctica?
      mutable UsefulAdu5Pat usefulPat; //!      /// Only construct this once, mutable since not const correct

      Event(Int_t nT=0);
      Event(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd, Int_t nT=0);
      Event(const Event& event);

      TArrowAntarctica* makeArrowFromAnitaToEvent();
      void setupUsefulPat();
      void resetClusteringNumbers();
      void deleteArrays();
      void setNThresholds(int n);

      virtual ~Event();
      ClassDef(Event, 8)
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

      void doClustering(const char* dataGlob, const char* mcGlob, const char* outFileName);
      
      TGraphAntarctica* makeClusterSummaryTGraph(Int_t clusterInd);

      bool getUseBaseList(){return fUseBaseList;}
      void setUseBaseList(bool useBaseList){ // *TOGGLE *GETTER=GetUseBaseList
	fUseBaseList = useBaseList;
      }

      bool getDebug(){return fDebug;}
      void setDebug(bool db){fDebug = db;} // *TOGGLE *GETTER=GetDebug

    private:
      Double_t getDistSqEventCluster(const Event& event, const Cluster& cluster);
      Double_t getAngDistSqEventCluster(const Event& event, const Cluster& cluster);
      void getDeltaThetaDegDeltaPhiDegEventCluster(const Event& event, const Cluster& cluster, Double_t& deltaThetaDeg, Double_t& deltaPhiDeg);



      Long64_t readInSummaries(const char* summaryGlob);
      size_t addEvent(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd);
      size_t addMcEvent(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd);      
      void assignSingleEventToCloserCluster(Int_t eventInd, Int_t isMC, Cluster& cluster);

      void initializeBaseList();
      Int_t histogramUnclusteredEvents(Int_t& globalMaxBin);
      void recursivelyAddClustersFromData(Int_t minBinContent);
      void assignMcEventsToClusters();
      void assignEventsToBaseClusters();
      void writeAllGraphsAndHists();
      void makeSummaryTrees();
      void resetClusters();
      Double_t getSumOfMcWeights();
      void redoSmallClusters();
      void initKDTree();
      Double_t dMin(const Event& event1, const Event& event2);
      Double_t dSum(const Event& event1, const Event& event2);
      Double_t dAsym(const Event& event1, const Event& event2);
      void testTriangleInequality();
      void testAngleFindingSpeed();
      Int_t removeLargeBasesNearMcMurdo();

      void doEventEventClustering();      
      void nearbyEvents(Int_t eventInd, std::vector<Int_t>& nearbyEvents, std::vector<double>& nearbyEventLLs, double llRange, double llFitThreshold=-1, double rangeEastingNorthing=default_range_easting_northing);
      void makeAndWriteNSquaredEventEventHistograms();
      Double_t dPoint(const Event& eventInd1, Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt, bool addOverHorizonPenalty=false);
      Double_t evalPairLogLikelihoodAtLonLat(const Double_t* params);
      std::vector<const Event*> fFitEvent1s; /// First event in the pairwise fit
      std::vector<const Event*> fFitEvent2s; /// Second event in the pairwise fit
      Int_t fMaxFitterAttempts; /// How many times should I try if I don't reach a good minimum?
      Double_t fFitHorizonDistM; ///700e3 metres, distance at which a penalty is added to source location fitting
      Double_t dFit(const Event& event1, const Event& event2);
      UInt_t fTestEvent1; /// For debugging
      UInt_t fTestEvent2; /// For debugging
      std::vector<Double_t> fFitEastings;  /// Where the fitter found the potential source could come from
      std::vector<Double_t> fFitNorthings; /// Where the fitter found the potential source could come from
      TH2DAntarctica* makeUnevenlyBinnedEventHistogram();
      std::vector<Double_t> llEventCuts;                        /// Try doing a range of llEventCuts at once...
      // Double_t llEventCut;					/// The cut-off for log-likelihood, which defines the boundary of a cluster
      Double_t llClusterCut;				       	/// The cut-off for log-likelihood, which defines the boundary of a cluster
      Bool_t doneBaseClusterAssignment;				/// Set to true once all read in data events were clustered to bases

      // Int_t numClusters;
      Int_t numCallsToRecursive;

      std::vector<Acclaim::Clustering::Cluster> clusters;	/// Vector of clusters,
      std::vector<Acclaim::Clustering::Event> events;		/// Vector of data events
      std::vector<Acclaim::Clustering::McEvent> mcEvents;	/// Vector of Monte Carlo events
      TKDTreeID* fKDTree;                                       /// ROOT's implementation of a KDTree, typedef'd for int/double
      std::vector<Double_t> fEventEastings;
      std::vector<Double_t> fEventNorthings;

      std::vector<TGraphAntarctica*> grBaseClusterCenters;	/// The locations of the bases
      std::vector<TH2DAntarctica*> hBaseClusteredEvents;	/// Histograms of events clustered to bases
      std::vector<TGraphAntarctica*> grNonBaseClusterCenters;	/// The locations of the non-base clusters
      std::vector<TH2DAntarctica*> hNonBaseClusteredEvents;	/// Histograms of events clustered to non-base clusters
      TH2DAntarctica* hClusters;                                /// Filled with clusters (allows access to the bin of the cluster)
      TH2DAntarctica* hEvents;
      bool fDebug;
      bool fUseBaseList;
      TGraph* grTestMinimizerWalk;

      std::vector<ROOT::Math::Minimizer*> fMinimizers;
      std::vector<ROOT::Math::Functor*> fFunctors;
      Int_t fROOTgErrorIgnoreLevel;
    };


  }
}

#endif
