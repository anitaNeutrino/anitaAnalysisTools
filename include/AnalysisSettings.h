/* -*- C++ -*-.***************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             Never again will a forgetful analyst (me) not be sure what the state of the analysis software was when he ran a script. This class should track all the configurable numbers in the analysis and automatically write them in an easily accessible way into an output file (though that may be handled elsewhere).
             The mantra here is no arbitrary choices in code, that don't get set from this class.
*****************************************************************************************************/

#ifndef ACCLAIM_ANALYSIS_SETTINGS_H
#define ACCLAIM_ANALYSIS_SETTINGS_H

#include <iostream>
#include "TFile.h"
#include "TString.h"
#include "TXMLEngine.h"


namespace Acclaim {

class AnalysisSettings{

 public:
  AnalysisSettings(const char* fileName = "example.conf");

 protected:
  Bool_t stringIsKeyValuePair(const TString& commentStrippedLine, TString& key, TString& value);
  Bool_t stringIsSection(const TString& commentStrippedLine, TString& secName);
  Bool_t stringIsAlphaNumeric(const TString& commentStrippedLine);
  
  void parseSettingsFile(const char* fileName);
  
  void handleKeyValue(const TString& key, const TString& value);
  void handleSection(const TString& section);
  
  TString settingsFileAsString;
};


}


#endif
