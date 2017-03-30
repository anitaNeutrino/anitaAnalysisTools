#include "AnalysisFlow.h"


using namespace Acclaim;

int main(int argc, char* argv[]){

  AnitaVersion::set(3);
  int run = argc > 1 ? atoi(argv[1]) : 352;  
  AnalysisFlow analysis("analyzeWais", run, AnalysisFlow::kWaisPulser);
  analysis.doAnalysis();
  

}
