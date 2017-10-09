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

// This one is useful for converting an enum in a file to Int_t
#define ENUM_ANALYSIS_SETTING(enum_type, SettingVariable) \
    protected:                                    \
       enum_type f##SettingVariable;              \
    public:                                       \
       Int_t Get##SettingVariable() const         \
       {                                          \
	 return (Int_t)f##SettingVariable;        \
       }                                          \
       enum_type GetEnum##SettingVariable() const \
       {                                          \
	 return f##SettingVariable;		  \
       }                                          \
       void Set##SettingVariable(Int_t val)       \
       {                                          \
	 f##SettingVariable = (enum_type)val;	  \
       }                                          \
       void Set##SettingVariable(enum_type val)   \
       {                                          \
	 f##SettingVariable = val;                \
       }



namespace Acclaim {

class AnalysisSettings{

  // nested maps
  typedef std::map<TString, TString> VariableMap_t;
  typedef std::map<TString, VariableMap_t*> SectionMap_t;

 public:
  AnalysisSettings(const char* fName = NULL);
  void apply(TObject* obj) const;
  void write(TFile* f) const;
  void print() const;

  Bool_t getSetting(const char* settingName, Bool_t&   settingVal) const;
  Bool_t getSetting(const char* settingName, Int_t&    settingVal) const;
  Bool_t getSetting(const char* settingName, Double_t& settingVal) const;
  
 protected:
  Bool_t stringIsKeyValuePair(const TString& commentStrippedLine, TString& key, TString& value) const;
  Bool_t stringIsSection(const TString& commentStrippedLine, TString& secName) const;
  Bool_t stringIsAlphaNumeric(const TString& commentStrippedLine) const;
  
  void handleKeyValue(VariableMap_t* variableMap, const TString& key, const TString& value);
  VariableMap_t* handleSection(const TString& section);

  Bool_t getSetting(const char* settingName, TString& settingVal) const;
  

  void parseSettingsFile();
  Bool_t tryFile(const char* fName);
  void findFile();

  TString fileName;
  TString parsedFileCopy;
  SectionMap_t sectionMap;
};


}


#endif
