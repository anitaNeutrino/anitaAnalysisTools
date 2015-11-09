/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             A class to manage output file names from environment variables, exec arguments and time.
*************************************************************************************************************** */

#ifndef OUTPUTCONVENTION_H
#define OUTPUTCONVENTION_H


#include "TDatime.h"
#include "TString.h"
#include "TFile.h"
#include <stdlib.h>

class OutputConvention{

public:

  OutputConvention(int argcIn, char* argvIn[]);
  ~OutputConvention();

  TString getOutputFileName();

  
private:
  int argc;
  char** argv;
  TDatime dateTime;
  TString outputDir;
  TString dateTimeSuffix;
  TString outFileName;
  
  TString getDateTimeSuffix();
  TString getOutputDir();

  
  
};




#endif
