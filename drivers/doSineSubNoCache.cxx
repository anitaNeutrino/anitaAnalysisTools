#include "AnalysisFlow.h"
#include "BasicFilters.h"
#include "AcclaimFilters.h"
#include "UCFilters.h"
#include "AcclaimCmdLineArgs.h"

using namespace Acclaim;

int main(int argc, char* argv[]){

  Acclaim::CmdLineArgs args(argc, argv);

  bool useCache = false;
  UCorrelator::SineSubtractFilter::setUseCache(useCache);
  
  FilterStrategy* strat = new FilterStrategy();
  ALFAFilter* alfaFilter = new ALFAFilter(Filters::Bands::alfaLowPassGHz);
  strat->addOperation(alfaFilter);
  UCorrelator::fillStrategyWithKey(strat, Acclaim::Filters::getCosminsFavouriteSineSubName());
  
  AnalysisFlow analysis(&args, strat);
  analysis.doAnalysis();
    
  return 0;
}
