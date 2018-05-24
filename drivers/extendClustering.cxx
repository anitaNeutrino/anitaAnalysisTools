#include "AcclaimClustering.h"
#include "OutputConvention.h"
#include "DrawStrings.h"
#include <iostream>

using namespace Acclaim;

int main(int argc, char* argv[]){

  if(argc!=3){
    std::cerr << argv[0] << " 'alreadyClusteredOutput' 'extraDataGlob'" << std::endl;
    return 1;
  }
  const char* outFileName = argv[0];
  const char* dataGlob = argv[1];
  const char* extendedGlob = argv[2];

  Clustering::LogLikelihoodMethod clusterer;
  // clusterer.setDebug(true);
  clusterer.fStoreUnclusteredHistograms = false;
  TCut subThreshold = "eventNumber==60832108";
  // clusterer.addCut(!ThermalTree::isAboveHorizontal + ThermalTree::passAllQualityCuts + ThermalTree::isNotTaggedAsPulser + ThermalTree::fisherCut + !ThermalTree::closeToHiCal + ThermalTree::closeToMC);
  clusterer.addCut(subThreshold);
  clusterer.doClustering(dataGlob, extendedGlob, outFileName);

  return 0;
}
