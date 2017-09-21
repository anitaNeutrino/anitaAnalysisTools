#include "AnalysisSettings.h"
#include <fstream>
#include "TObjString.h"
#include "TString.h"
#include "TObjArray.h"
#include "TRegexp.h"
#include "TClass.h"
#include "TDataMember.h"
#include "TMethodCall.h"
#include <cstdlib>
#include <stdexcept>


const char* default_filename = "AcclaimSettings.conf";

Acclaim::AnalysisSettings::AnalysisSettings(const char* fName) : fileName(fName) {

  parseSettingsFile();
}


void Acclaim::AnalysisSettings::write(TFile* f){

  if(!f){
    std::cerr << "Error in " << __PRETTY_FUNCTION__
              << ", got NULL pointer to TFile. Cannot save AnalysisSettings." << std::endl;
    return;
  }

  if(!f->IsWritable()){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", TFile " << f->GetName()
              << " is not writable. Cannot save AnalysisSettings." << std::endl;
    return;
  }

  if(!f->IsOpen()){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", TFile " << f->GetName()
              << " is not open. Cannot save AnalysisSettings." << std::endl;
    return;
  }
  
  f->cd();

  TObjArray* tokens = fileName.Tokenize("/");
  TString shortFileName = ((TObjString*) tokens->Last())->GetString();  

  TNamed* tn = new TNamed(shortFileName, parsedFileCopy);
  tn->Write();
  delete tn;
  
}


void Acclaim::AnalysisSettings::print(){
  // For debugging
  std::cout << __PRETTY_FUNCTION__ << std::endl;  

  SectionMap_t::iterator sit = sectionMap.begin();
  for(; sit!= sectionMap.end(); ++sit){
    VariableMap_t* vm = sit->second;
    std::cout << "[" << sit->first << "]" << std::endl;

    VariableMap_t::iterator vit = vm->begin();
    for(; vit!=vm->end(); ++vit){
      std::cout << vit->first << "=" << vit->second << std::endl;
    }

    std::cout << std::endl;
  }
}




void Acclaim::AnalysisSettings::apply(TObject* obj){

  if(obj==NULL){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", got NULL TObject pointer. Doing nothing." << std::endl;
    return;
  }

  TString className = obj->ClassName();

  SectionMap_t::iterator sit = sectionMap.find(className);
  if(sit==sectionMap.end()){
    std::cerr << "Fatal error in " << __PRETTY_FUNCTION__ << ", found no settings for " << className
              << "." << std::endl;
    throw std::runtime_error("Missing settings");
  }
  
  VariableMap_t* varMap = sit->second;
  TClass* cl = obj->IsA();

  VariableMap_t::iterator vit = varMap->begin();
  for(; vit!= varMap->end(); ++vit){
    TDataMember *dm = cl->GetDataMember(vit->first);

    if(!dm){
      std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", unable to find data member for "
                << className << "." << vit->first << "." << std::endl;
      throw std::runtime_error("Missing data member");
    }
    else{
      TMethodCall *sm = dm->SetterMethod(cl);
      if(!sm){
        std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", unable to find setter method for "
                  << className << "." << vit->first << ". Skipping this setting." << std::endl;
        throw std::runtime_error("Missing Setter function");        
      }
      else{
        sm->Execute(obj, vit->second);
      }
    }
  }
}



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


Acclaim::AnalysisSettings::VariableMap_t* Acclaim::AnalysisSettings::handleSection(const TString& section){  
  
  // std::cout << "I found a section: " << section << std::endl;

  SectionMap_t::iterator it = sectionMap.find(section);
  VariableMap_t* variableMap = NULL;
  if(it!=sectionMap.end()){
    variableMap = it->second;
  }
  else{
    variableMap = new VariableMap_t();
    sectionMap[section] = variableMap;
  }
  return variableMap;
}


void Acclaim::AnalysisSettings::handleKeyValue(VariableMap_t* variableMap, const TString& key, const TString& value){

  if(variableMap==NULL){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << " found variable " << key << " with value " << value << " but has no section! Variables must have a section. Skipping this entry." << std::endl;
    return;
  }

  VariableMap_t::iterator it = variableMap->find(key);
  if(it!=variableMap->end()){
    // something funky going on here...
    TString sectionName;
    SectionMap_t::iterator it2 = sectionMap.begin();    
    for(; it2 != sectionMap.end(); ++it2){
      if(variableMap==it2->second){
        sectionName = it2->first;
        break;
      }
    }
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", while parsing section " << sectionName \
              << "I found a duplicate entry for " << key << std::endl;
    std::cerr << "The original value is " << it->second << ", the duplicated value is " << value << std::endl;
    std::cerr << "Ignoring the second value." << std::endl;
  }
  else{
    (*variableMap)[key] = value;
  }
}


Bool_t Acclaim::AnalysisSettings::tryFile(const char* fName){
  Bool_t openedFile = false;

  TString n(fName);
  if(n.Length() > 0 && n[n.Length()-1] == '/'){
    std::cerr << "Info in " << __PRETTY_FUNCTION__ << ", interpreting forward slash terminated input as directory, appending " << default_filename << std::endl;
    n += TString::Format("%s", default_filename);
  }

  std::cout << "Info in " << __PRETTY_FUNCTION__ << ", trying "<< n.Data() << std::endl;  
  
  std::ifstream a(n.Data());
  if(a.is_open()){
    openedFile = true;
    fileName = n.Data();
  }
  return openedFile;
}

void Acclaim::AnalysisSettings::findFile(){

  if(fileName.Length()!=0){
    // try to open current file name
    if(tryFile(fileName)){      
      return;
    }
  }

  // that didn't work so try a couple more things ...

  // default file name in the present working directory
  const char* fName = default_filename;
  if(tryFile(fName)){
    fileName = TString(fName);
    return;
  }

  const char* anitaUtilInstallDir = getenv("ANITA_UTIL_INSTALL_DIR");
  if(!anitaUtilInstallDir){
    // you've got some major issues...
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", ANITA_UTIL_INSTALL_DIR not set." << std::endl;
  }

  TString fName2 = TString::Format("%s/share/Acclaim/%s", anitaUtilInstallDir, default_filename);
  if(tryFile(fName2.Data())){
    fileName = TString(fName2);
    return;
  }

  // if you get here, then there's some problems...
  std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", couldn't find any settings files! Giving up." << std::endl;
  exit(1);
}



void Acclaim::AnalysisSettings::parseSettingsFile(){

  findFile();
  std::cout << "After findFile(), fileName = " << fileName << std::endl;
  
  std::ifstream settingsFile(fileName);

  if(!settingsFile.is_open()){
    // if you get here, then there's some problems...
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", couldn't find settings file " << fileName
              << ". Giving up." << std::endl;
    throw "Acclaim::AnalysisSettings Error";
  }
  else{
    int lineNum = 1;

    VariableMap_t* varMap = NULL;
    
    while(!settingsFile.eof()){

      std::string thisLine;
      std::getline(settingsFile, thisLine);

      parsedFileCopy.Append(thisLine);
      parsedFileCopy.Append("\n");

      // First cut out the comment, which is all characters after the first #, including the #
      std::size_t found = thisLine.find("#");

      // Here we switch to TString because of its lovely tokenization methods
      TString thisLineCommentsRemoved(thisLine.substr(0, found));

      TString section;
      Bool_t isSection = stringIsSection(thisLineCommentsRemoved, section);

      if(isSection){
        varMap = handleSection(section);
      }
      else{
        TString key, value;
        Bool_t isKeyValuePair = stringIsKeyValuePair(thisLineCommentsRemoved, key, value);        
        if(isKeyValuePair){
          handleKeyValue(varMap, key, value);
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
      
      lineNum++;
    }
    
  }
}



