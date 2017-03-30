/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             Just for fun.
*************************************************************************************************************** */


#include "OutputConvention.h"
#include <unistd.h>


using namespace Acclaim;

int main(int argc, char* argv[]){
  

  OutputConvention oc(argc, argv);

  TString outFileName = oc.getOutputFileName();

  TFile* f = new TFile(outFileName, "recreate");
  f->Close();
  
  
  return 0;
}
