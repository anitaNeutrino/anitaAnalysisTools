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

  FilterStrategy* strat = new FilterStrategy();

  AnalysisFlow analysis(argv[0], run, AnalysisFlow::kAll, strat, AnitaDataset::kNoBlinding);
  analysis.doAnalysis();

  return 0;
}
