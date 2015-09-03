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


/*!
  \brief Default constructor - don't use this
*/
ProgressBar::ProgressBar(){
  std::cerr << "Assuming 100 events in ProgressBar" << std::endl;
  maxEntry = 100;
  counter = 0;
  percentage = 0;
  watch.Start(kTRUE);
}


/*!
  \brief Useful constructor - do use this one.
  \param maxEntryInit is the number of events you want to loop over
*/
ProgressBar::ProgressBar(Long64_t maxEntryInit){
  maxEntry = maxEntryInit;
  counter = 0;
  percentage = 0;
  watch.Start(kTRUE);
}


/*!
  \brief Increment operator, use when you have completed one iteration of the main loop
*/
void ProgressBar::operator++(int){

  if(percentage>=100) return;
  
  // Stops the watch
  Int_t seconds = Int_t(watch.RealTime());
  Int_t hours = seconds / 3600;
  hours = hours < 0 ? 0 : hours;
  seconds = seconds - hours * 3600;
  Int_t mins = seconds / 60;
  mins = mins < 0 ? 0 : mins;
  
  seconds = seconds - mins * 60;
  
  counter++;
  double ratio = double(counter)/maxEntry;
  // std::cout << ratio << "\t" << counter << "\t" << maxEntry << std::endl;

  if(ratio*100 > percentage){

    while(ratio*100 > percentage){
      percentage++;
    }

    // fprintf(stderr, "\n\033[F\033[J");
    // std::cout << std::endl;
    // fprintf(stderr, "\033[F\033[J");
    fprintf(stderr, "\r");

    // Show the percentage complete.
    fprintf(stderr, ANSI_COLOR_RED);
    fprintf(stderr, "%3u%%", (UInt_t)(percentage) );
    fprintf(stderr, ANSI_COLOR_RESET); 
    fprintf(stderr, " [");

    // Show the load bar.
    fprintf(stderr, ANSI_COLOR_BLUE);
    fprintf(stderr, "%02d:%02d:%02d", hours, mins, seconds);
    for (UInt_t i=8; i<percentage; i++){
      fprintf(stderr, "=");
    }
    fprintf(stderr, ANSI_COLOR_RESET);
 
    Int_t startSpace = percentage > 8 ? percentage : 8;
    for (Int_t i=startSpace; i<100; i++){
      fprintf(stderr, " ");
    }
 
    fprintf(stderr, "]");
  }

  if(percentage>=100) fprintf(stderr, "\n");
  watch.Start(kFALSE);
  return;

}

/*!
  \brief For debugging, prints state of internal variables
*/
void ProgressBar::status(){
  std::cout << percentage << "\t" << counter << "\t" << maxEntry << std::endl;
}
