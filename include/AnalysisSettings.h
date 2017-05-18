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

// Macro used to define an analysis SettingVariable and create the associate getter/setter functions
//
// This is evil, but a lesser of two evils I think.
// It ensures that the functions always match the naming conventions expected by ROOT.
// Using ANALYSIS_SETTING(Int_t, Example) in a class definition defines the following
// Int_t fExample;
// Int_t GetExample();
// void SetExample(Int_t)
//
// NOTE: Does NOT work with string-like types. TString, const char*, Option_t* etc all cause problems!
//       Known working types are Int_t, Bool_t, Double_t, stick to those.
//       Current work around is define your string option in an array and pick with an index.

#define ANALYSIS_SETTING(var_type, SettingVariable) \
    protected: \
       var_type f##SettingVariable; \
    public: \
       var_type Get##SettingVariable() const \
       {\
          return f##SettingVariable; \
       }\
       void Set##SettingVariable(var_type val) \
       {\
          f##SettingVariable = val; \
       }



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
