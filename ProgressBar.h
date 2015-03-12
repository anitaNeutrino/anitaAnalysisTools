/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             Just a bit of fun while waiting for jobs to run.
*************************************************************************************************************** */



#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#include <TObject.h>
#include <iostream>

class ProgressBar : public TObject{

public:
  ProgressBar();
  ProgressBar(Long64_t maxEntry);

  /* handles everything */
  void operator++(int);

private:
  Long64_t maxEntry;
  Long64_t counter;    
  UInt_t percentage;

  ClassDef(ProgressBar, 0);

};


#endif //PROGRESSBAR_H
