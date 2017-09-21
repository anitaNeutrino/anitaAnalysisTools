#include "AnalysisFlow.h"
#include "BasicFilters.h"
#include "AcclaimFilters.h"
#include "UCFilters.h"
#include "AcclaimCmdLineArgs.h"

using namespace Acclaim;

int main(int argc, char* argv[]){

  Acclaim::CmdLineArgs args(argc, argv);

  std::map<TString, FilterStrategy*> filterStrats;
  bool saveOutput = false;
  Filters::appendFilterStrategies(filterStrats, saveOutput);
  
  // FilterStrategy* strat = Filters::findStrategy(filterStrats, "RayleighFilter");
  FilterStrategy* strat = new FilterStrategy();
  UCorrelator::fillStrategyWithKey(strat,Acclaim::Filters::getCosminsFavouriteSineSubName());
  
  AnalysisFlow analysis(&args, strat);
  analysis.doAnalysis();
    
  return 0;
}
