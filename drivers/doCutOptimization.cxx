#include "CutOptimizer.h"

using namespace Acclaim;

int main(int argc, char* argv[]){
  (void) argc;

  const TString outFileName = argv[0];
  const TString signalFileName = "filterWaisAcclaim_352_2017-06-20_18-48-08.root";
  const TString backgroundFileName = "filterDecimatedAcclaim_390_08_2017-05-18_22-26-28.root";
  
  CutOptimizer co(outFileName, signalFileName, backgroundFileName);
  co.optimize();

  return 0;
}
