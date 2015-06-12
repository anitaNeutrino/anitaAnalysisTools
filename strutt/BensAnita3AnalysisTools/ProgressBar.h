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

  /* handles everything */
  void operator++(int);

  /* For debugging */
  void status();

private:
  Long64_t maxEntry;
  Long64_t counter;    
  UInt_t percentage;
  TStopwatch watch;

  ClassDef(ProgressBar, 0);

};


#endif //PROGRESSBAR_H
