#include "CutOptimizer.h"
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
  const char* signalGlob = argv[1];
  const char* backgroundGlob = argc >= 2 ? argv[2] : NULL;

  CutOptimizer co(signalGlob, backgroundGlob);


  std::vector<const TCut *> signalSelection;
  

  signalSelection.push_back(&ThermalTree::closeToMC);
  
  std::vector<const TCut *> backgroundSelection;

  backgroundSelection.push_back(&ThermalTree::isAboveHorizontal);
  backgroundSelection.push_back(&ThermalTree::anita3QuietTime);
  backgroundSelection.push_back(&ThermalTree::isNotTaggedAsPulser);
  backgroundSelection.push_back(&ThermalTree::passAllQualityCuts);


  std::vector<TString> variables;

  variables.push_back("coherent_filtered_fracPowerWindowGradient/deconvolved_filtered_fracPowerWindowGradient");
  variables.push_back("deconvolved_filtered_fracPowerWindowGradient");
  // variables.push_back("coherent_filtered_fracPowerWindowGradient");
  // variables.push_back("deconvolved_filtered_fracPowerWindowGradient");

  variables.push_back("deconvolved_filtered_impulsivityMeasure/coherent_filtered_impulsivityMeasure");
  variables.push_back("deconvolved_filtered_impulsivityMeasure");
  // variables.push_back("coherent_filtered_impulsivityMeasure");
  // variables.push_back("deconvolved_filtered_impulsivityMeasure");
  variables.push_back("coherent_filtered_peakHilbert");
  variables.push_back("deconvolved_filtered_peakHilbert");
  variables.push_back("peak_value");
  
  // co.setDebug(true);
  co.optimize(signalSelection, backgroundSelection, variables, outFileName);

  // gSystem->Exit(0);
  return 0;
}
