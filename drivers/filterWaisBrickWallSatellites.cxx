#include "AnalysisFlow.h"
#include "BasicFilters.h"
#include "AcclaimFilters.h"

using namespace Acclaim;

int main(int argc, char* argv[]){

  AnitaVersion::set(3);
  int run = argc > 1 ? atoi(argv[1]) : 352;

  std::map<TString, FilterStrategy*> filterStrats;
  bool saveOutput = true;  
  Filters::appendFilterStrategies(filterStrats, saveOutput);

  FilterStrategy* stupidNotchStrat = Filters::findStrategy(filterStrats, "BrickWallSatellites");
  if(!stupidNotchStrat){ 
    std::cerr << "Well, this script is pointless... I give up." << std::endl;
    return 1;
  }
  AnalysisFlow analysisSimpleNotchFilters(argv[0], run, AnalysisFlow::kWaisPulser, stupidNotchStrat);
  analysisSimpleNotchFilters.doAnalysis();
    
  return 0;
}
