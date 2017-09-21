#include "AnalysisFlow.h"
#include "BasicFilters.h"
#include "AcclaimFilters.h"
#include "AcclaimCmdLineArgs.h"

using namespace Acclaim;

int main(int argc, char* argv[]){

  Acclaim::CmdLineArgs args(argc, argv);  

  std::map<TString, FilterStrategy*> filterStrats;
  bool saveOutput = false;
  Filters::appendFilterStrategies(filterStrats, saveOutput);
  
  FilterStrategy* strat = Filters::findStrategy(filterStrats, "RayleighFilter");
  if(!strat){ 
    std::cerr << "Well, this script is pointless... I give up." << std::endl;
    return 1;
  }
  AnalysisFlow analysis(&args, strat);
  analysis.doAnalysis();
    
  return 0;
}
