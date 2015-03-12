/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             Just for fun.
*************************************************************************************************************** */

#include "ProgressBar.h"
#include <unistd.h>

int main(){
  
  ProgressBar progress(100);
  for(int i=0; i<100; i++){
    progress++;
  }

  ProgressBar progress2(10);
  for(int i=0; i<10; i++){
    progress2++;
    usleep(200000);
    std::cerr << "some text output" << std::endl;
  }


  return 0;
}
