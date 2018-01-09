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

  const TCut realSnr = "!TMath::IsNaN(coherent_filtered_snr) && !TMath::IsNaN(deconvolved_filtered_snr) && TMath::Finite(coherent_filtered_snr) && TMath::Finite(deconvolved_filtered_snr)";  

  std::vector<const TCut *> signalSelection;
  
  const TCut closeToMC = "fabs(FFTtools::wrap(peak_phi - mc_phi)) < 5 && fabs(peak_theta + mc_theta) < 3.5";
  signalSelection.push_back(&closeToMC);
  
  std::vector<const TCut *> backgroundSelection;
  const TCut isAboveHorizontal = "peak_theta > 0";
  const TCut notPulser = "flags_pulser == 0";
  const TCut isRF = "flags_isRF == 0";
  const TCut notShortWaveform = "flags_isVarner2 == 0";
  
  backgroundSelection.push_back(&isAboveHorizontal);
  backgroundSelection.push_back(&SumTree::anita3QuietTime);
  backgroundSelection.push_back(&notPulser);
  backgroundSelection.push_back(&realSnr);
  backgroundSelection.push_back(&notShortWaveform);
  
  std::vector<TString> variables;

  variables.push_back("coherent_filtered_fracPowerWindowGradient");
  variables.push_back("deconvolved_filtered_fracPowerWindowGradient");
  variables.push_back("coherent_filtered_impulsivityMeasure");
  variables.push_back("deconvolved_filtered_impulsivityMeasure");    
  variables.push_back("coherent_filtered_peakHilbert");
  variables.push_back("deconvolved_filtered_peakHilbert");
  
  // co.setDebug(true);
  co.optimize(signalSelection, backgroundSelection, variables, outFileName);

  // gSystem->Exit(0);
  return 0;
}
