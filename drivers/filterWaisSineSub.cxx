#include "AnalysisFlow.h"
#include "BasicFilters.h"
#include "AcclaimFilters.h"
#include "UCFilters.h"

using namespace Acclaim;

int main(int argc, char* argv[]){

  AnitaVersion::set(3);
  int run = argc > 1 ? atoi(argv[1]) : 352;

  std::map<TString, FilterStrategy*> filterStrats;
  bool saveOutput = false;
  Filters::appendFilterStrategies(filterStrats, saveOutput);
  
  // FilterStrategy* strat = Filters::findStrategy(filterStrats, "RayleighFilter");
  FilterStrategy* strat = new FilterStrategy();
  UCorrelator::fillStrategyWithKey(strat,"sinsub_10_3_ad_2");  

  for(unsigned i=0; i < strat->nOperations(); i++){
    std::cout << i << "\t" << strat->getOperation(i)->description() << std::endl;
  }
  // return 0;
  
  AnalysisFlow analysis(argv[0], run, AnalysisFlow::kWaisPulser, strat, AnitaDataset::kDefault);
  analysis.doAnalysis();
    
  return 0;
}
