#include "NoiseMonitor.h"
#include "BasicFilters.h"
#include "AcclaimFilters.h"
#include "UCFilters.h"
#include "AnitaDataset.h"
#include "RawAnitaHeader.h"
#include "AcclaimCmdLineArgs.h"

using namespace Acclaim;

int main(int argc, char* argv[]){

  Acclaim::CmdLineArgs args(argc, argv);
  UCorrelator::SineSubtractFilter::setUseCache(true);

  FilterStrategy* strat = new FilterStrategy();
  ALFAFilter* alfaFilter = new ALFAFilter(Filters::Bands::alfaLowPassGHz);
  strat->addOperation(alfaFilter);
  UCorrelator::fillStrategyWithKey(strat, Acclaim::Filters::getCosminsFavouriteSineSubName());  

  // return 0;
  NoiseMonitor nm(strat);
  AnitaDataset d(args.run);

  // this should force the generation
  std::cout << nm.getRMS(AnitaPol::kHorizontal, 0, d.header()->realTime) << std::endl;
  
  return 0;
}
