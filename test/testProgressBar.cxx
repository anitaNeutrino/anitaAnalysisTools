/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             Just for fun.
*************************************************************************************************************** */

#include "ProgressBar.h"
#include <unistd.h>


using namespace Acclaim;

int main(){
  
  ProgressBar progress(100);
  for(int i=0; i<100; i++){
    progress++;
  }

  Int_t n2 = 11024;
  ProgressBar progress2(n2);
  for(int i=0; i<n2; i++){
    progress2++;
    usleep(20000);
    // std::cerr << "some text output" << std::endl;
  }


  return 0;
}
