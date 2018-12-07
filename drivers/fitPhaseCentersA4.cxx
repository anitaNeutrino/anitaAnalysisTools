#include "AcclaimPhaseCenterFitter.h"

int main(int argc, char* argv[]){

  Acclaim::PhaseCenterFitter fitter("AcclaimCorrelationSummary_*.root");
  fitter.fit();
  return 0;
}
