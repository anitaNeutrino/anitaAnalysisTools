#include "AnalysisFlow.h"
#include "BasicFilters.h"
#include "AcclaimFilters.h"
#include "UCFilters.h"
#include "AcclaimCmdLineArgs.h"

using namespace Acclaim;

int main(int argc, char* argv[]){

  Acclaim::CmdLineArgs args(argc, argv);
  AnalysisFlow analysis(&args);
  analysis.doAnalysis();
    
  return 0;
}
