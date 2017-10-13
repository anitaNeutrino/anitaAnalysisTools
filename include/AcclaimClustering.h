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

namespace Acclaim{

  namespace Clustering {

    const double default_sigma_theta = 0.25;
    const double default_sigma_phi = 0.5;

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
      Int_t peakIndex;                          /// Which peak in the map does this represent?

      Double_t centre[3];//!			/// Cartesian coordinates, does not persist in ROOT
      Double_t latitude;			/// latitude
      Double_t longitude;			/// longitude
      Double_t altitude;			/// longitude
      AnitaEventSummary::PayloadLocation anita; /// Anita's position

      Double_t theta;				/// reconstructed theta
      Double_t phi;				/// reconstructed phi

      Double_t sigmaTheta;			/// resolution associated with this snr?
      Double_t sigmaPhi;			/// resolution associated with this snr?

      //--------------------------------------------------------------------------------
      // determined by clustering
      //--------------------------------------------------------------------------------
      Double_t dThetaCluster;			/// theta distance to cluster
      Double_t dPhiCluster;			/// phi distance to cluster

      Double_t ll;				/// log likelihood = -2 * ln (...)
      Int_t inCluster;				/// which cluster am I associated with?

      Double_t llSecondBest;			/// log likelihood to second closest cluster
      Int_t secondClosestCluster;		/// what cluster am I second closest to?

      Int_t antarcticaHistBin; //!              /// Which global bin in the TH2DAntarctica?

      Event();
      Event(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd);
	
      virtual ~Event(){ ;}
      ClassDef(Event, 3)
    };



    /** 
     * @class McEvent
     * @brief Holds the position of monte carlo generated neutrino event
     * 
     * Requires a little more information than a data point to keep track of the neutrino energy and weight
     */
    class McEvent : public Event{
    public:
      Double_t weight;			/// MC information
      Double_t energy;			/// MC information

      McEvent();
      McEvent(const AnitaEventSummary* sum, AnitaPol::AnitaPol_t pol, Int_t peakInd)
	: Event(sum, pol,  peakInd)
      {
	weight = sum->mc.weight;
	energy = sum->mc.energy;
      }

      virtual  ~McEvent(){;}
	
      ClassDef(McEvent, 1)
    };    





    //--------------------------------------------------------------------------------------------------------
    /**
     * @class Cluster
     * @brief A position on the continent, with which a bunch of points are associated
     */
    class Cluster{
    public:
      Cluster();
      // Cluster(const Event& seedEvent);
      Cluster(const BaseList::base& base);

      virtual ~Cluster(){ ;}	

      Double_t centre[3]; //! Does not persist  /// Center in cartesian
	
      Double_t latitude; /// cluster center latitude
      Double_t longitude;/// cluster center longitude
      Double_t altitude;/// cluster center altitude
	
      Int_t numDataEvents; /// How many data events does this contain?

      ClassDef(Cluster, 2)
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

      double llCut;
      double maxDistCluster;
      std::vector<Int_t> numIsolatedSmallClusters;
      std::vector<Int_t> numIsolatedSmallBaseClusters;
      Double_t numMcIsolatedSinglets;
      int maxRetestClusterSize;

      std::vector<TGraphAntarctica*> grBaseClusterCenters;
      std::vector<TH2DAntarctica*> hBaseClusteredEvents;

      std::vector<TGraphAntarctica*> grNonBaseClusterCenters;
      std::vector<TH2DAntarctica*> hNonBaseClusteredEvents;
      
      Int_t numIter;
      Int_t numClusters;
      Int_t numCallsToRecursive;

      // The clusters
      std::vector<Cluster> clusters; /// Vector of clusters, 
      std::vector<Event> events; /// Vector of data events
      std::vector<McEvent> mcEvents; /// Vector of Monte Carlo events

      bool doneBaseClusterAssignment; /// Set to true once all read in data events were clustered to bases

    };
  }
}

#endif
