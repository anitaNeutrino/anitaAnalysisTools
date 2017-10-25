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

class TTree;
class TGraphAntarctica;
class TFile;
class TH2DAntarctica;
class TCanvas;

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

      Int_t antarcticaHistBin; //!		/// Which global bin in the TH2DAntarctica?

      Event();
      Event(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd);
	
      virtual ~Event(){ ;}
      ClassDef(Event, 3)
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
     * @class LogLikelihoodMethod
     * @brief Workhorse class: groups events into clusters around bases and other events using the log likelihood method
     *
     */
    class LogLikelihoodMethod {
    public:

      LogLikelihoodMethod();
      virtual ~LogLikelihoodMethod();

      void doClustering(const char* dataGlob, const char* mcGlob, const char* outFileName);
      
      TGraphAntarctica* makeClusterSummaryTGraph(Int_t clusterInd);

      Int_t getNumClusters(){
	return numClusters;
      }

    private:
      Double_t getDistSqEventCluster(Int_t eventInd, const Acclaim::Clustering::Cluster& cluster);
      // Double_t getDistSqClusterCluster(Int_t clusterInd1, Int_t clusterInd2);
      Double_t getAngDistSqEventCluster(Int_t eventInd, Int_t clusterInd, UsefulAdu5Pat& usefulPat);
      void getDeltaThetaDegDeltaPhiDegEventCluster(Int_t eventInd, Int_t clusterInd, UsefulAdu5Pat& usefulPat, Double_t& deltaThetaDeg, Double_t& deltaPhiDeg);

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

      Double_t llCut;						/// The cut-off for log-likelihood, which defines the boundary of a cluster
      Double_t maxDistCluster;					/// Only consider adding to a cluster when closer than this value
      Bool_t doneBaseClusterAssignment;				/// Set to true once all read in data events were clustered to bases
      Int_t fSmallClusterSizeThreshold;

      Int_t numClusters;
      Int_t numCallsToRecursive;

      std::vector<Acclaim::Clustering::Cluster> clusters;	/// Vector of clusters,
      std::vector<Acclaim::Clustering::Event> events;		/// Vector of data events
      std::vector<Acclaim::Clustering::McEvent> mcEvents;	/// Vector of Monte Carlo events

      std::vector<TGraphAntarctica*> grBaseClusterCenters;	/// The locations of the bases
      std::vector<TH2DAntarctica*> hBaseClusteredEvents;	/// Histograms of events clustered to bases
      std::vector<TGraphAntarctica*> grNonBaseClusterCenters;	/// The locations of the non-base clusters
      std::vector<TH2DAntarctica*> hNonBaseClusteredEvents;	/// Histograms of events clustered to non-base clusters
      TH2DAntarctica* hClusters;                                /// Filled with clusters (allows access to the bin of the cluster)
      bool fDebug;
    };


  }
}

#endif
