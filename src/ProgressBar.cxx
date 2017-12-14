#include "ProgressBar.h"

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"



Int_t Acclaim::ProgressBar::progState = 0;

//---------------------------------------------------------------------------------------------------------
/**
 * @brief Default constructor - don't use this
*/
Acclaim::ProgressBar::ProgressBar(){
  std::cerr << "Assuming 100 events in ProgressBar" << std::endl;
  fMaxEntry = 100;
  fCounter = 0;
  fPercentage = 0;
  fWatch.Start(kTRUE);
  fSetHandler = 0;
  fNumBreakTries = 0;
  fLastProgState = 0;
}






//---------------------------------------------------------------------------------------------------------
/**
 * @brief Useful constructor - do use this one.
 * @param maxEntryInit is the number of events you want to loop over
*/
Acclaim::ProgressBar::ProgressBar(Long64_t maxEntryInit){
  fMaxEntry = maxEntryInit > 1 ? maxEntryInit - 1 : 1;
  fCounter = 0;
  fPercentage = 0;
  fWatch.Start(kTRUE);
  fSetHandler = 0;
  fNumBreakTries = 0;
  fLastProgState = 0;  
}






//---------------------------------------------------------------------------------------------------------
/**
 * @brief Increment operator, use when you have completed one iteration of the main loop
*/
void Acclaim::ProgressBar::operator++(int){

  if(fPercentage>=100) return;
  
  // Stops the fWatch
  Int_t seconds = Int_t(fWatch.RealTime());
  Int_t hours = seconds / 3600;
  hours = hours < 0 ? 0 : hours;
  seconds = seconds - hours * 3600;
  Int_t mins = seconds / 60;
  mins = mins < 0 ? 0 : mins;
  
  seconds = seconds - mins * 60;
  
  fCounter++;
  double ratio = double(fCounter)/fMaxEntry;
  // std::cout << ratio << "\t" << fCounter << "\t" << maxEntry << std::endl;

  if(ratio*100 > fPercentage){

    while(ratio*100 > fPercentage){
      fPercentage++;
    }

    // fprintf(stderr, "\n\033[F\033[J");
    // std::cout << std::endl;
    // fprintf(stderr, "\033[F\033[J");
    fprintf(stderr, "\r");

    // Show the fPercentage complete.
    fprintf(stderr, ANSI_COLOR_RED);
    fprintf(stderr, "%3u%%", (UInt_t)(fPercentage) );
    fprintf(stderr, ANSI_COLOR_RESET); 
    fprintf(stderr, " [");

    // Show the load bar.
    fprintf(stderr, ANSI_COLOR_BLUE);
    fprintf(stderr, "%02d:%02d:%02d", hours, mins, seconds);
    for (UInt_t i=8; i<fPercentage; i++){
      fprintf(stderr, "=");
    }
    fprintf(stderr, ANSI_COLOR_RESET);
 
    Int_t startSpace = fPercentage > 8 ? fPercentage : 8;
    for (Int_t i=startSpace; i<100; i++){
      fprintf(stderr, " ");
    }
 
    fprintf(stderr, "]");
  }

  if(fPercentage>=100) fprintf(stderr, "\n");
  fWatch.Start(kFALSE);
  return;

}



//---------------------------------------------------------------------------------------------------------
/** 
 * @brief New primary function to move through main for loop in analysis program
 *
 * This function wraps the Long64_t implementation of inc
 * 
 * @param entry is a reference to the loop varible, assumed to start at 0.
 * @param numEntries is the maximum entry, assumed that the loop condition is entry < numEntries
 */
void Acclaim::ProgressBar::inc(UInt_t& entry, Long64_t numEntries){
  Long64_t entryL64 = entry;
  inc(entryL64, numEntries);
  entry = entryL64;
}


//---------------------------------------------------------------------------------------------------------
/** 
 * @brief New primary function to move through main for loop in analysis program
 *
 * This function wraps the Long64_t implementation of inc
 * 
 * @param entry is a reference to the loop varible, assumed to start at 0.
 * @param numEntries is the maximum entry, assumed that the loop condition is entry < numEntries
 */
void Acclaim::ProgressBar::inc(Int_t& entry, Long64_t numEntries){
  Long64_t entryL64 = entry;
  inc(entryL64, numEntries);
  entry = entryL64;
}



//---------------------------------------------------------------------------------------------------------
/** 
 * @brief New primary function to move through main for loop in analysis program
 *
 * This now looks for SIGINT signals and alters the main loop variable, entry = numEntries.
 * 
 * @param entry is a reference to the loop varible, assumed to start at 0.
 * @param numEntries is the maximum entry, assumed that the loop condition is entry < numEntries
 */
void Acclaim::ProgressBar::inc(Long64_t& entry, Long64_t numEntries){
  numEntries = numEntries < 0 ? fMaxEntry : numEntries;

  if(fSetHandler==0){
    signal (SIGINT, Acclaim::ProgressBar::mainLoopSigintHandle);
    fSetHandler = 1;
  }

  if(fSetHandler==1 && fLastProgState!=0){
    if(fNumBreakTries==0)
    {
      std::cerr << "Program with ProgressBar received SIGINT, will try and exit main loop gracefully. " << std::endl;
      entry=numEntries;
    }
    else
    {
      std::cerr << "Unable to exit main loop gracefully, raising SIGINT properly." << std::endl;
      signal(SIGINT, SIG_DFL); // now sigint points to the default handler again
      raise(SIGINT); // raise sigint
    }
    fNumBreakTries++;
  }

  int diff = entry - fCounter;
  for(int i=0; i < diff; i++){
    // std::cout << fCounter << "\t" << entry << std::endl;
    (*this)++;
  }
  fLastProgState = progState;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief For debugging, prints state of internal variables
*/
void Acclaim::ProgressBar::status(){
  std::cout << fPercentage << "\t" << fCounter << "\t" << fMaxEntry << std::endl;
}






//---------------------------------------------------------------------------------------------------------
/** 
 * @brief Custom handler, sets variable when signal is received.
 * 
 * @param param is the signal, I presume
 */
void Acclaim::ProgressBar::mainLoopSigintHandle(int param){
  progState = param;
}
