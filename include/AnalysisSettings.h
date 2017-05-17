/* -*- C++ -*-.***************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             Never again will a forgetful analyst (me) not be sure what the state of the analysis software was when he ran a script.
             This class should read, set and store all the configurable numbers in the analysis.
             The mantra here is no arbitrary or semi-arbitrary choices should be hard coded. They should be here.
*****************************************************************************************************/

#ifndef ACCLAIM_ANALYSIS_SETTINGS_H
#define ACCLAIM_ANALYSIS_SETTINGS_H

#include <iostream>
#include "TString.h"
#include "TFile.h"
#include <map>

namespace Acclaim {

class AnalysisSettings{

  // nested maps
  typedef std::map<TString, TString> VariableMap_t;
  typedef std::map<TString, VariableMap_t*> SectionMap_t;

 public:
  AnalysisSettings(const char* fName = NULL);
  // AnalysisSettings(const char* fName = "AcclaimSettings.conf");  
  void apply(TObject* obj);
  void write(TFile* f);
  void print();
  
 protected:
  Bool_t stringIsKeyValuePair(const TString& commentStrippedLine, TString& key, TString& value);
  Bool_t stringIsSection(const TString& commentStrippedLine, TString& secName);
  Bool_t stringIsAlphaNumeric(const TString& commentStrippedLine);
  
  void handleKeyValue(VariableMap_t* variableMap, const TString& key, const TString& value);
  VariableMap_t* handleSection(const TString& section);

  void parseSettingsFile();
  Bool_t tryFile(const char* fName);
  void findFile();

  TString fileName;
  TString parsedFileCopy;
  SectionMap_t sectionMap;
};


}


#endif
