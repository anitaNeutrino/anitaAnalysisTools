#include "AnalysisSettings.h"
#include <fstream>
#include "TObjString.h"
#include "TString.h"
#include "TObjArray.h"
#include "TRegexp.h"


Bool_t Acclaim::AnalysisSettings::stringIsKeyValuePair(const TString& commentStrippedLine, TString& key, TString& value){
  key = "";
  value = "";
  
  Bool_t isKeyValuePair = false;

  TObjArray* tokens = commentStrippedLine.Tokenize("=");
  int nTokens = tokens->GetEntries();

  // If there are two tokens, then there was one delimeter so we're good
  if(nTokens == 2){
    key = ((TObjString*) tokens->At(0))->GetString();
    value = ((TObjString*) tokens->At(1))->GetString();
    isKeyValuePair = true;
  }
  delete tokens;
  return isKeyValuePair;
}


Bool_t Acclaim::AnalysisSettings::stringIsAlphaNumeric(const TString& str){
  // really this should be "stringIsNotWhiteSpace", but that's hard to figure out how to do
  TRegexp reggie("[a-zA-Z0-9]");
  Ssiz_t len = str.Length();

  Bool_t isAlphaNumeric = reggie.Index(str, &len) != -1;
  return isAlphaNumeric;
}



Bool_t Acclaim::AnalysisSettings::stringIsSection(const TString& commentStrippedLine, TString& secName){

  secName = "";
  Bool_t isSection = false;
  Ssiz_t len = commentStrippedLine.Length();

  if(len > 0){
    TString firstChar = TString(commentStrippedLine(0,1));
    if(firstChar.Contains("[") && commentStrippedLine.Contains("]")){
      
      TObjArray* tokens = TString(commentStrippedLine(1, len)).Tokenize("]");
      int nTokens = tokens->GetEntries();
      TString secNameMaybe = ((TObjString*) tokens->At(0))->GetString();

      if(nTokens<=2){
        Bool_t trailingGood = true;
        if(nTokens==2){
          TString trailing = ((TObjString*) tokens->At(1))->GetString();
          if(stringIsAlphaNumeric(trailing)){
            trailingGood = false;
          }
        }
        if(trailingGood){
          secName = secNameMaybe;
          isSection = true;
        }
      }
      delete tokens;
    }
  }  
  return isSection;
}


void Acclaim::AnalysisSettings::handleSection(const TString& section){
  std::cout << "I found a section: " << section << std::endl;
}

void Acclaim::AnalysisSettings::handleKeyValue(const TString& key, const TString& value){
  std::cout << "I found a key/value pair: " << key << "\t" << value << std::endl;
}

void Acclaim::AnalysisSettings::parseSettingsFile(const char* fileName){

  std::ifstream settingsFile(fileName);

  if(!settingsFile.is_open()){
    std::cerr << "Oops" << std::endl;
  }
  else{
    int lineNum = 1;
    while(!settingsFile.eof()){

      std::string thisLine;
      std::getline(settingsFile, thisLine);

      settingsFileAsString.Append(thisLine);
      settingsFileAsString.Append("\n");

      // First cut out the comment, which is all characters after the first #, including the #
      std::size_t found = thisLine.find("#");

      // Here we switch to TString because of its lovely tokenization methods
      TString thisLineCommentsRemoved(thisLine.substr(0, found));

      TString section;
      Bool_t isSection = stringIsSection(thisLineCommentsRemoved, section);


      if(isSection){
        handleSection(section);
        // take action
      }
      else{
        TString key, value;
        Bool_t isKeyValuePair = stringIsKeyValuePair(thisLineCommentsRemoved, key, value);        
        if(isKeyValuePair){
          handleKeyValue(key, value);
        }
        else{
          Bool_t isAlphaNumeric = stringIsAlphaNumeric(thisLineCommentsRemoved);
          if(isAlphaNumeric){
            std::cerr << "Warning in " << __PRETTY_FUNCTION__ 
                      << ". I couldn't parse line " << lineNum 
                      << " in " << fileName << ". It said: " << std::endl;
            std::cerr << thisLine << std::endl;
          }
          else{
            // then I think it was whitespace...
          }
        }
      }
      
      // std::cout << settingsFileAsString << std::endl;
      lineNum++;
    }
    
  }
}

Acclaim::AnalysisSettings::AnalysisSettings(const char* fileName){

  parseSettingsFile(fileName);
}


