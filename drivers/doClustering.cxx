#include "AcclaimClustering.h"
#include "OutputConvention.h"
#include "DrawStrings.h"
#include <iostream>

using namespace Acclaim;

int main(int argc, char* argv[]){

  if(!(argc == 2 || argc==3)){
    std::cerr << argv[0] << " 'signalAndBackgroundGlob' " << std::endl;
    std::cerr << argv[0] << " 'signalGlob' 'backGroundGlob'" << std::endl;
    return 1;
  }
  const char* outFileName = argv[0];
  const char* dataGlob = argv[1];
  const char* mcGlob = argc >= 2 ? argv[2] : NULL;

  Clustering::LogLikelihoodMethod clusterer;
  // clusterer.setDebug(true);
  clusterer.fStoreUnclusteredHistograms = false;
  clusterer.addCut(!ThermalTree::isAboveHorizontal + ThermalTree::passAllQualityCuts + ThermalTree::isNotTaggedAsPulser + ThermalTree::fisherCut + !ThermalTree::closeToHiCal + ThermalTree::closeToMC);

  clusterer.doClustering(dataGlob, mcGlob, outFileName);

  return 0;
}
