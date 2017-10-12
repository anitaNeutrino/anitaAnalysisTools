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

  class AntarcticaMapPlotter;  

  namespace Clustering {

    //--------------------------------------------------------------------------------------------------------
    /**
     * @class Event
     * @brief Will populate a tree showing the results of the clustering
     * 
     */
    class Event{
    public:
      UInt_t eventNumber;				/// The event number
      Int_t run;					/// The run 
      AnitaPol::AnitaPol_t pol;				/// Polarisation of waveform
      Double_t sourceLat;				/// latitude of position of event
      Double_t sourceLon;				/// longitude of position of event
      Double_t sourceAlt;				/// altitude of posiiton of event
      
      Double_t thetaDeg;				/// reconstruction angle
      Double_t phiDeg;					/// reconstruction angle

      // AnitaEventSummary::PayloadLocation anitaPosition; /// Where is ANITA?

      Double_t sigmaThetaDeg;
      Double_t sigmaPhiDeg;

      Double_t distanceToClusterCentroid;		/// distance (km) from point to centroid
      Double_t distanceToClusterCentroidSecondBest;	/// distance (km) from point to centroid
      
      Double_t minusTwoLogLikelihood;			/// distance in some normalized error units TODO
      Double_t minusTwoLogLikelihoodSecondBest;		/// distance in some normalized error units TODO

      Double_t deltaThetaDeg;				/// angular distance to cluster centre
      Double_t deltaPhiDeg;				/// angular distance to cluster centre

      Int_t inCluster;					/// ID of cluster
      Int_t secondClosestCluster;			/// ID of cluster
      Int_t numPointsWithinMinLL;			/// if > -1 then is a singlet, if > 0 then is a non-isolated singlet

      Int_t isBase;
      Int_t numEventsInCluster;				/// number of events in the cluster containing this event

							/// Double_t clusterPosition[3]; // centroid of cluster cartesian (m)
      Double_t clusterLat;				/// latitude of centroid of cluster
      Double_t clusterLon;				/// longitude of centroid of cluster
      Double_t clusterAlt;				/// altitude of centroid of cluster

      Double_t anitaLat;				/// latitude of centroid of cluster
      Double_t anitaLon;				/// longitude of centroid of cluster
      Double_t anitaAlt;				/// altitude of centroid of cluster

      Int_t numClusters;				/// Total number of clusters (maybe gratuitous)
      Int_t numIterations;				/// number of loop iterations (maybe gratuitous)

      Event(){
	numPointsWithinMinLL = -1;
      }
    };

    /**
     * @class MCEvent
     * @brief Holds additional monte carlo information
     */
    class MCEvent : public Event {
    public:
      /// @todo add constructors
      Double_t weight;
      Double_t energy;
    };



    /**
     * @class LogLikelihoodMethod
     * @brief Workhorse class: groups events into clusters around bases and other events using the log likelihood method
     *
     */
    class LogLikelihoodMethod {

    public:



      



      /**
       * @class Point
       * @brief Internal representation of a thing being clustered
       */
      class Point{
      public:
	Double_t centre[3];		/// Cartesian coordinates

	Double_t latitude;		/// latitude
	Double_t longitude;		/// longitude
	Double_t altitude;		/// longitude

	Double_t thetaDeg;		/// reconstructed theta
	Double_t phiDeg;		/// reconstructed phi
	
	Double_t sigmaThetaDeg;		/// resolution associated with this snr?
	Double_t sigmaPhiDeg;		/// resolution associated with this snr?


	

	Double_t dTheta;		/// theta distance to cluster
	Double_t dPhi;			/// phi distance to cluster
	
	Double_t ll;			/// log likelihood = -2 * ln (...)
	Int_t inCluster;		/// which cluster am I associated with?
	
	Double_t llSecondBest;		/// log likelihood to second closest cluster
	Int_t secondClosestCluster;	/// what cluster am I second closest to?
	
	AnitaPol::AnitaPol_t pol;	/// polarization

	Point();
	Point(Adu5Pat* pat, Double_t lat=0, Double_t lon=0, Double_t alt=0,
	      Double_t thetaDeg=0, Double_t phiDeg=0,
	      Double_t sigmaTheta = 0.25, Double_t sigmaPhi = 0.5,
	      Int_t polIn=AnitaPol::kVertical);
	
	virtual ~Point(){ ;}
	ClassDef(Point, 2)
      };


      

      /** 
       * @class McPoint
       * @brief Holds the position of monte carlo generated neutrino event
       * 
       * Requires a little more information than a data point to keep track of the neutrino energy and weight
       */
      class McPoint : public Point{
      public:
	Double_t weight;
	Double_t energy;

	explicit McPoint();
	explicit McPoint(Adu5Pat* pat,
			 Double_t lat=0, Double_t lon=0, Double_t alt=0,
			 Double_t thetaDeg = 0, Double_t phiDeg = 0,
			 Double_t sigmaTheta = 0.25, Double_t sigmaPhi = 0.5,
			 Int_t polIn=AnitaPol::kVertical,
			 Double_t theWeight=1, Double_t theEnergy=0);
	virtual  ~McPoint(){;}
	
	ClassDef(McPoint, 1)
      };



      //--------------------------------------------------------------------------------------------------------
      /**
       * @class Cluster
       * @brief A position on the continent, with which a bunch of points are associated
       */
      class Cluster{
      public:
	Cluster();
	explicit Cluster(const Point& seedPoint);
	explicit Cluster(const BaseList::base& base);

	virtual ~Cluster(){ ;}	

	Double_t centre[3];
	
	Double_t latitude;
	Double_t longitude;
	Double_t altitude;
	
	Int_t numEvents;
	Double_t totalError;
	Double_t maxDist;
	Int_t numPointsWithinMinLL;

	ClassDef(Cluster, 1)
      };



      

      LogLikelihoodMethod();
      virtual ~LogLikelihoodMethod();

      void doClustering(const char* dataGlob, const char* mcGlob, const char* outFileName);
      
      TGraphAntarctica* makeClusterSummaryTGraph(Int_t clusterInd);
      TTree* makeClusterSummaryTree(TFile* fOut, TFile* fSignalBox);

      Int_t getNumClusters(){
	return numClusters;
      }

    private:

      Long64_t readInSummaries(const char* summaryGlob);
      
      size_t addPoint(Adu5Pat* pat, Double_t latitude, Double_t longitude, Double_t altitude, Int_t run, UInt_t eventNumber, Double_t thetaDeg, Double_t phiDeg, Double_t sigmaThetaDeg, Double_t sigmaPhiDeg, AnitaPol::AnitaPol_t pol);
      size_t addMcPoint(Adu5Pat* pat, Double_t latitude, Double_t longitude, Double_t altitude, Int_t run, UInt_t eventNumber, Double_t thetaDeg, Double_t phiDeg, Double_t sigmaThetaDeg, Double_t sigmaPhiDeg, AnitaPol::AnitaPol_t pol, Double_t weight, Double_t energy);
      
      void assignSinglePointToCloserCluster(Int_t pointInd, Int_t isMC, Int_t clusterInd);
      void initializeBaseList();
      void initializeEmptyBaseList();

      Int_t histogramUnclusteredEvents(Int_t& globalMaxBin);
      void recursivelyAddClusters(Int_t minBinContent);
      void assignMcPointsToClusters();
      void assignEventsToDefaultClusters();
      void findClosestPointToClustersOfSizeOne();

      void resetClusters();
      Double_t getSumOfMcWeights();

      Cluster seedCluster(Point& point);
      
      double llCut;
      double maxDistCluster;
      std::vector<Int_t> numIsolatedSmallClusters;
      std::vector<Int_t> numIsolatedSmallBaseClusters;
      Double_t numMcIsolatedSinglets;
      int maxRetestClusterSize;
      
      std::vector<TH2DAntarctica*> hUnclustereds;
      Int_t numIter;
      Int_t numClusters;
      Int_t numCallsToRecursive;


      // The clusters
      std::vector<Cluster> clusters; // only variables relevant to clustering
      std::vector<Int_t> seedPoints; // which points seeded which clusters
      
      // Data
      std::vector<Point> points; // only variables relevant to clustering
      std::vector<Adu5Pat> pats;
      std::vector<UInt_t> eventNumbers; // keep track of these separately
      std::vector<Int_t> runs; // keep track of these separately

      // MC
      std::vector<McPoint> mcPoints;
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

      bool doneDefaultAssignment;

    };
  }
}

#endif
