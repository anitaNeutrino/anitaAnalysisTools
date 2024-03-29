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
#include "TChain.h"

#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>
#include <iostream>
#include <algorithm>



namespace Acclaim
{

  /**
   * @class OutputConvention
   * @brief A class to systematically name files produced by my analysis programs.
   *
   * Uses the program name, arguments, date, and time.
   */
  class OutputConvention{

  public:

    OutputConvention(int argcIn, char* argvIn[]);

    TString getOutputFileName(TString ext="");
    TString getOutputDir();
    TFile* makeFile(); //!< Create the new output file with proper name

    static TFile* getFile(TString fileNameWithWildcards); //!< Opens matching file with the most recent suffix

  private:
    int argc; /// The argc from the main program
    char** argv; /// The argv from the main program
    TDatime dateTime; /// The dateTime type containing the date/time.
    TString outputDir; /// TString contining the output directory.
    TString dateTimeSuffix; /// TString contining the fileName suffix extracted from the date/time.
    TString outFileName; /// The output total file output name.
    TString getDateTimeSuffix();


  };


}

#endif
