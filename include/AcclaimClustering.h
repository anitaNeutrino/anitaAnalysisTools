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
     * @brief Internal representation of a thing being clustered
     */
    class Event{
    public:

      //--------------------------------------------------------------------------------
      // determined by reconstruction
      //--------------------------------------------------------------------------------	
      UInt_t eventNumber;		/// Event number
      Int_t run;			/// Run
      AnitaPol::AnitaPol_t pol;	/// Polarization
	
      Double_t centre[3];//!		/// Cartesian coordinates, does not persist in ROOT
      Double_t latitude;		/// latitude
      Double_t longitude;		/// longitude
      Double_t altitude;		/// longitude

      Double_t theta;			/// reconstructed theta
      Double_t phi;			/// reconstructed phi
	
      Double_t sigmaTheta;		/// resolution associated with this snr?
      Double_t sigmaPhi;		/// resolution associated with this snr?


      //--------------------------------------------------------------------------------
      // determined by clustering
      //--------------------------------------------------------------------------------
      Double_t dThetaCluster;		/// theta distance to cluster
      Double_t dPhiCluster;		/// phi distance to cluster

      Double_t ll;			/// log likelihood = -2 * ln (...)
      Int_t inCluster;		/// which cluster am I associated with?

      Double_t llSecondBest;		/// log likelihood to second closest cluster
      Int_t secondClosestCluster;	/// what cluster am I second closest to?	

      Event();
      Event(Adu5Pat* pat, Double_t lat=0, Double_t lon=0, Double_t alt=0,
	    Double_t thetaDeg=0, Double_t phiDeg=0,
	    Double_t sigmaTheta = default_sigma_theta,
	    Double_t sigmaPhi = default_sigma_phi,
	    Int_t polIn=AnitaPol::kVertical);
	
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
      Double_t weight;
      Double_t energy;

      explicit McEvent();
      explicit McEvent(Adu5Pat* pat,
		       Double_t lat=0, Double_t lon=0, Double_t alt=0,
		       Double_t thetaDeg = 0, Double_t phiDeg = 0,
		       Double_t sigmaTheta = default_sigma_theta,
		       Double_t sigmaPhi = default_sigma_phi,
		       Int_t polIn=AnitaPol::kVertical,
		       Double_t theWeight=1, Double_t theEnergy=0);
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
      explicit Cluster(const Event& seedEvent);
      explicit Cluster(const BaseList::base& base);

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
      
      size_t addEvent(Adu5Pat* pat, Double_t latitude, Double_t longitude, Double_t altitude, Int_t run, UInt_t eventNumber, Double_t thetaDeg, Double_t phiDeg, Double_t sigmaThetaDeg, Double_t sigmaPhiDeg, AnitaPol::AnitaPol_t pol);
      size_t addMcEvent(Adu5Pat* pat, Double_t latitude, Double_t longitude, Double_t altitude, Int_t run, UInt_t eventNumber, Double_t thetaDeg, Double_t phiDeg, Double_t sigmaThetaDeg, Double_t sigmaPhiDeg, AnitaPol::AnitaPol_t pol, Double_t weight, Double_t energy);
      
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

      Cluster seedCluster(Event& point);
      
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
      std::vector<Cluster> clusters; // only variables relevant to clustering
      std::vector<Int_t> seedEvents; // which points seeded which clusters
      
      // Data
      std::vector<Event> points; // only variables relevant to clustering
      std::vector<Adu5Pat> pats;
      std::vector<UInt_t> eventNumbers; // keep track of these separately
      std::vector<Int_t> runs; // keep track of these separately

      // MC
      std::vector<McEvent> mcEvents;
      std::vector<Adu5Pat> mcPats;
      std::vector<UInt_t> mcEventNumbers;
      std::vector<Int_t> mcRuns;

      // ok since things got very slow all of a sudden I need to speed them up
      // since my interupts in lldb seem to fall inside the UsefulAdu5Pat functions
      // I will cache the values.
      // For each point I need the value for each cluster... use pair<pointInd, clusterInd>
      // the figure of merit is the loglikelihood, with a distance d.

      std::vector<Int_t> pointsInCluster;
      Int_t theMinCluster;
      Double_t sumOfAngularErrorsFromLatLon(const Double_t* latLon);
      std::vector<Int_t> ampBinNumbers;
      std::vector<Int_t> ampBinNumbers2;

      bool doneBaseClusterAssignment; /// Set to true once all read in data events were clustered to bases

    };
  }
}

#endif
