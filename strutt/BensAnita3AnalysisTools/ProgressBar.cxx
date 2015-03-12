/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             Just a bit of fun.
*************************************************************************************************************** */


#include "ProgressBar.h"

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

ClassImp(ProgressBar);

ProgressBar::ProgressBar(){
  std::cerr << "Assuming 100 events in ProgressBar" << std::endl;
  maxEntry = 100;
  counter = 0;
  percentage = 0;
}

ProgressBar::ProgressBar(Long64_t maxEntryInit){
  maxEntry = maxEntryInit;
  counter = 0;
  percentage = 0;
}

void ProgressBar::operator++(int){

  if(percentage>=100) return;
  
  counter++;
  double ratio = double(counter)/maxEntry;

  if(ratio*100 > percentage){

    while(ratio*100 > percentage){
      percentage++;
    }

    // printf("\n\033[F\033[J");
    // std::cout << std::endl;
    // printf("\033[F\033[J");
    printf("\r");

    // Show the percentage complete.
    printf(ANSI_COLOR_RED);
    printf("%3u%%", (UInt_t)(percentage) );
    printf(ANSI_COLOR_RESET); 
    printf(" [");

    // Show the load bar.
    printf(ANSI_COLOR_BLUE);
    for (UInt_t i=0; i<percentage; i++){
      printf("=");
    }
    printf(ANSI_COLOR_RESET);
 
    for (UInt_t i=percentage; i<100; i++){
      printf(" ");
    }
 
    printf("]");
  }

  if(percentage>=100) printf("\n");
  
  return;

}
