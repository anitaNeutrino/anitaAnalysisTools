/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             Just a bit of fun while waiting for jobs to run.
*************************************************************************************************************** */



#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#include "TObject.h"
#include "TStopwatch.h"
#include <iostream>
#include <signal.h>

namespace Acclaim
{

  /** @class ProgressBar
   * @brief Prints a progress bar and timer to stderr
   */
  class ProgressBar{

  public:
    ProgressBar();
    ProgressBar(Long64_t maxEntry);

    void operator++(int);
    void status();

    void inc(Long64_t& entry, Long64_t numEntries=-1);
    void inc(UInt_t& entry,   Long64_t numEntries=-1);
    void inc(Int_t& entry,    Long64_t numEntries=-1);

    static void mainLoopSigintHandle(int param);
    static int progState;
  
  private:
    Long64_t fMaxEntry;		/// Number of events you will loop over
    Long64_t fCounter;		/// Number of loops completed
    UInt_t fPercentage;		/// Percentage to print to the screen
    TStopwatch fWatch;		/// ROOT's stopwatch class, used to time the progress since object construction
    Int_t fSetHandler;		/// Have we set the interupt signal handler?
    Int_t fNumBreakTries;	/// Did I manage to break of the the loop, not if this number goes above 1
    Int_t fLastProgState;	/// Was the program state OK the last time I incremented
  };
}

#endif //PROGRESSBAR_H
