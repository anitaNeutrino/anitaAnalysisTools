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
  
  FilterStrategy* strat = new FilterStrategy();
  UCorrelator::fillStrategyWithKey(strat,Acclaim::Filters::getCosminsFavouriteSineSubName());

  FilterStrategy* stupidNotchStrat = Filters::findStrategy(filterStrats, "BrickWallSatellites");

  for(unsigned i=0; i < stupidNotchStrat->nOperations(); i++){
    strat->addOperation((FilterOperation*)stupidNotchStrat->getOperation(i));
  }
  
  AnalysisFlow analysis(&args, strat);
  analysis.doAnalysis();
    
  return 0;
}
