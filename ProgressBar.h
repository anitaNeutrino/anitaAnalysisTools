/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             Just a bit of fun while waiting for jobs to run.
*************************************************************************************************************** */



#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#include <TObject.h>
#include <TStopwatch.h>
#include <iostream>

class ProgressBar : public TObject{

public:
  ProgressBar();
  ProgressBar(Long64_t maxEntry);

  void operator++(int);
  void status();

private:
  Long64_t maxEntry; ///< Number of events you will loop over
  Long64_t counter; ///< Number of loops completed
  UInt_t percentage; ///< Percentage to print to the screen
  TStopwatch watch; ///< ROOT's stopwatch class, used to time the progress since object construction

  ClassDef(ProgressBar, 0);

};


#endif //PROGRESSBAR_H
