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

    void inc(Long64_t& entry, Long64_t numEntries);
    void inc(UInt_t& entry, Long64_t numEntries);
    void inc(Int_t& entry, Long64_t numEntries);

    static void mainLoopSigintHandle(int param);
    static int progState;
  
  private:
    Long64_t maxEntry; //!< Number of events you will loop over
    Long64_t counter; //!< Number of loops completed
    UInt_t percentage; //!< Percentage to print to the screen
    TStopwatch watch; //!< ROOT's stopwatch class, used to time the progress since object construction
    Int_t setHandler; //!< Have we set the interupt signal handler?
    Int_t numBreakTries; ///!< Did I manage to break of the the loop, not if this number goes above 1?
  };
}

#endif //PROGRESSBAR_H
