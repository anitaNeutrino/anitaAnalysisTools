#include "AnalysisFlow.h"
#include "BasicFilters.h"
#include "AcclaimFilters.h"
#include "UCFilters.h"
#include "AcclaimCmdLineArgs.h"
#include <memory>

using namespace Acclaim;

int main(int argc, char* argv[]){

  auto geom = AnitaGeomTool::Instance();
  geom->usePhotogrammetryNumbers(true);
  
  auto calib = AnitaEventCalibrator::Instance();
  for(auto& channelsPerSurf : calib->relativePhaseCenterToAmpaDelays){
    for(auto& chanExtraDt : channelsPerSurf){
      chanExtraDt = 0;
    }
  }

  
  Acclaim::CmdLineArgs args(argc, argv);

  auto strat = std::make_shared<FilterStrategy>();
  // ALFAFilter* alfaFilter = new ALFAFilter(Filters::Bands::alfaLowPassGHz);
  // strat->addOperation(alfaFilter);
  // UCorrelator::fillStrategyWithKey(strat, Acclaim::Filters::getCosminsFavouriteSineSubName());
  
  AnalysisFlow analysis(&args, strat.get());
  analysis.doAnalysis();

  return 0;
}
