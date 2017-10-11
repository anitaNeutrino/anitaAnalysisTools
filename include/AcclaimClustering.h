/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             A class to cluster. What were you expecting?
***********************************************************************************************************/

#ifndef ACCLAIM_CLUSTERING_H
#define ACCLAIM_CLUSTERING_H

#include "UsefulAdu5Pat.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TTree.h"
#include "TFile.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "RootTools.h"
#include "BaseList.h"
#include "ProgressBar.h"
#include "AntarcticaMapPlotter.h"

#include "assert.h"
#include "AnitaConventions.h"

#define nDim 3


namespace Acclaim{

  namespace Clustering {

    //--------------------------------------------------------------------------------------------------------
    /**
     * @class Event
     * @brief Stores the useful subset of AnitaEventSummary information for clustering
     * 
     */
    class Event{
    public:
      UInt_t eventNumber; ///The event number
      Int_t run; /// The run 
      AnitaPol::AnitaPol_t pol; /// Polarisation of waveform
      Double_t eventLat; // latitude of position of event
      Double_t eventLon; // longitude of position of event
      Double_t eventAlt; // altitude of posiiton of event

      Double_t distanceToClusterCentroid; // distance (km) from point to centroid
      Double_t distanceToClusterCentroidSecondBest; // distance (km) from point to centroid
      Double_t minusTwoLogLikelihood; // distance in some normalized error units TODO
      Double_t minusTwoLogLikelihoodSecondBest; // distance in some normalized error units TODO

      Double_t deltaThetaDeg; // angular distance to cluster centre
      Double_t deltaPhiDeg; // angular distance to cluster centre
      Double_t thetaDeg; // reconstruction angle
      Double_t phiDeg; // reconstruction angle

      Int_t inCluster; // ID of cluster
      Int_t secondClosestCluster; // ID of cluster
      //Int_t numPointsWithinMinLL = -1; // if > -1 then is a singlet, if > 0 then is a non-isolated singlet
      Int_t numPointsWithinMinLL; // if > -1 then is a singlet, if > 0 then is a non-isolated singlet

      Int_t isBase;
      Int_t numEventsInCluster; // number of events in the cluster containing this event
      // Double_t clusterPosition[nDim]; // centroid of cluster cartesian (m)
      Double_t clusterLat; // latitude of centroid of cluster
      Double_t clusterLon; // longitude of centroid of cluster
      Double_t clusterAlt; // altitude of centroid of cluster

      Double_t anitaLat; // latitude of centroid of cluster
      Double_t anitaLon; // longitude of centroid of cluster
      Double_t anitaAlt; // altitude of centroid of cluster

      Int_t numClusters; // Total number of clusters (maybe gratuitous)
      Int_t numIterations; // number of loop iterations (maybe gratuitous)

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
      Int_t isMC;
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
       * @brief Minimalistic representation of a thing inside the LogLikelihoodMethod
       */
      class Point{
      public:
	Double_t centre[nDim];
	Double_t latitude;
	Double_t longitude;
	Double_t altitude;
	Double_t thetaDeg;
	Double_t phiDeg;
	Double_t dTheta; // theta distance to cluster
	Double_t dPhi; // phi distance to cluster
	Double_t sigmaThetaDeg; // resolution associated with this snr?
	Double_t sigmaPhiDeg; // resolution associated with this snr?
	Double_t ll; //
	Int_t inCluster; // which cluster am I associated with?
	Double_t llSecondBest;
	Int_t secondClosestCluster;
	AnitaPol::AnitaPol_t pol; // polarization

	Point();	
	Point(Adu5Pat* pat, Double_t lat=0, Double_t lon=0, Double_t alt=0,
	      Double_t sigmaTheta = 0.5, Double_t sigmaPhi = 1, Int_t polIn=AnitaPol::kVertical);
	
	virtual ~Point(){ ;}
	ClassDef(Point, 2)
      };


      

      /** 
       * @class MCPoint
       * @brief Holds the position of monte carlo generated neutrino event
       * 
       * Requires a little more information than a data point to keep track of the neutrino energy and weight
       */
      class MCPoint : public Point{
      public:
	Double_t weight;
	Double_t energy;
	explicit MCPoint();
	explicit MCPoint(Adu5Pat* pat,
			 Double_t lat=0, Double_t lon=0, Double_t alt=0,
			 Double_t sigmaTheta = 0.5, Double_t sigmaPhi = 1,
			 Int_t polIn=AnitaPol::kVertical,
			 Double_t theWeight=1, Double_t theEnergy=0);
	ClassDef(MCPoint, 1)
      };


      //--------------------------------------------------------------------------------------------------------
      /**
       * @class Cluster
       * @brief Where the events are clustered
       */
      class Cluster{
      public:
	Cluster();
	explicit Cluster(const Point& seedPoint);
	explicit Cluster(const BaseList::base& base);

	Double_t centre[nDim];
	Double_t latitude;
	Double_t longitude;
	Double_t altitude;
	Int_t numEvents;
	Double_t totalError;
	Double_t maxDist;
	Int_t numPointsWithinMinLL;

	virtual ~Cluster(){ ;}
	ClassDef(Cluster, 1)
      };



      

      LogLikelihoodMethod(Int_t nClusters, Int_t numIterations, Int_t approxNumPoints=0);
      void doClustering(const char* dataGlob, const char* mcGlob);
      

      void initializeBaseList();
      void initializeEmptyBaseList();

      TGraph* makeClusterSummaryTGraph(Int_t clusterInd);
      TTree* makeClusterSummaryTree(TFile* fOut, TFile* fSignalBox);

      Int_t getNumClusters(){
	return numClusters;
      }

      Int_t histogramUnclusteredEvents(Int_t& globalMaxBin);
      void recursivelyAddClusters(Int_t minBinContent);
      void assignMCPointsToClusters();
      void assignEventsToDefaultClusters();
      void findClosestPointToClustersOfSizeOne();
      double llCut;
      double maxDistCluster;
      std::vector<Int_t> numIsolatedSmallClusters;
      std::vector<Int_t> numIsolatedSmallBaseClusters;
      Double_t numMcIsolatedSinglets;

      void resetClusters();
      Double_t getSumOfMcWeights();
      int maxRetestClusterSize;

    private:

      Long64_t readInSummaries(const char* summaryGlob);
      size_t addPoint(Adu5Pat* pat, Double_t latitude, Double_t longitude, Double_t altitude, Int_t run, UInt_t eventNumber, Double_t sigmaThetaDeg, Double_t sigmaPhiDeg, AnitaPol::AnitaPol_t pol);
      size_t addMCPoint(Adu5Pat* pat, Double_t latitude, Double_t longitude, Double_t altitude, Int_t run, UInt_t eventNumber, Double_t sigmaThetaDeg, Double_t sigmaPhiDeg, AnitaPol::AnitaPol_t pol, Double_t weight, Double_t energy);
      

      void assignSinglePointToCloserCluster(Int_t pointInd, Int_t isMC, Int_t clusterInd);
      AntarcticaMapPlotter* amp;
      Cluster seedCluster(Point& point);
      std::vector<TH2D*> hUnclustereds;
      Int_t numIter;
      Int_t numClusters;
      Int_t numCallsToRecursive;

      std::vector<Point> points; // only variables relevant to clustering
      std::vector<MCPoint> mcPoints; // only variables relevant to clustering
      std::vector<Adu5Pat*> pats;
      std::vector<Adu5Pat*> mcpats;
      std::vector<Cluster> clusters; // only variables relevant to clustering
      std::vector<Int_t> seedPoints; // which points seeded which clusters
      std::vector<UInt_t> eventNumbers; // keep track of these separately
      std::vector<UInt_t> mceventNumbers; // keep track of these separately
      std::vector<Int_t> runs; // keep track of these separately
      std::vector<Int_t> mcruns; // keep track of these separately

      Bool_t initialized;
      Double_t minimalImprovement;


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
