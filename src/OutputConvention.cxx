#include "OutputConvention.h"





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Constructor
 *
 * @param argcIn should be the main executable's argc value.
 * @param argvIn should be the main executable's argv value.
 */
OutputConvention::OutputConvention(int argcIn, char* argvIn[]){
  dateTime = TDatime();
  argc = argcIn;
  argv = argvIn;
  outFileName = "";
  dateTimeSuffix = "";
  outputDir = "";
  // subDir = "";
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Get the name of the output file from the program name (and system time).
 *
 * @param ext is an optional parameter file extension you want (e.g. .txt, .root, .csv), if none is given then .root is used.
 * @return the output file name
 */
TString OutputConvention::getOutputFileName(TString ext){

  if(outFileName==""){

    // Directory
    outFileName = getOutputDir();

    // Excutable name
    TString theArgv0 = TString::Format("%s", argv[0]);
    TObjArray* tkns = theArgv0.Tokenize("/");
    if(tkns->GetSize() > 0){
      TObjString* lastPartOfName = (TObjString*) tkns->At(tkns->GetLast());
      theArgv0 = lastPartOfName->GetString();
    }
    
    outFileName += theArgv0;

    // Excutable args
    if(argc > 1){
      for(int argInd=1; argInd < argc; argInd++){
	outFileName += TString::Format("_%s", argv[argInd]);
      }
    }
    
    // Date and time of running executable 
    outFileName += getDateTimeSuffix();

    if(ext==""){
      // ROOT suffix      
      outFileName += ".root";
    }
    else{
      if(strncmp(ext.Data(), ".", 1)!=0){
	ext = TString::Format(".%s", ext.Data());
      }
      outFileName += ext;
    }
  }

  
  return outFileName;

  
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Get the file suffix from the system time.
 *
 * @return the suffix for the file name, based on the date and time.
 */
TString OutputConvention::getDateTimeSuffix(){
  if(dateTimeSuffix==""){
    // dateTimeSuffix = TString::Format("_%d-%d-%d_%d-%d-%d",
    // 				     dateTime.GetYear(),
    // 				     dateTime.GetMonth(),
    // 				     dateTime.GetDay(),
    // 				     dateTime.GetHour(),
    // 				     dateTime.GetMinute(),
    // 				     dateTime.GetSecond()
    // 				     );
    dateTimeSuffix = TString::Format("_%d", dateTime.GetYear());

    Int_t m = dateTime.GetMonth();
    if(m < 10){
      dateTimeSuffix += TString::Format("-0%d", m);
    }
    else{
      dateTimeSuffix += TString::Format("-%d", m);
    }

    Int_t d = dateTime.GetDay();
    if(d < 10){
      dateTimeSuffix += TString::Format("-0%d", d);
    }
    else{
      dateTimeSuffix += TString::Format("-%d", d);
    }

    Int_t h = dateTime.GetHour();
    if(h < 10){
      dateTimeSuffix += TString::Format("_0%d", h);
    }
    else{
      dateTimeSuffix += TString::Format("_%d", h);
    }

    Int_t m2 = dateTime.GetMinute();
    if(m2 < 10){
      dateTimeSuffix += TString::Format("-0%d", m2);
    }
    else{
      dateTimeSuffix += TString::Format("-%d", m2);
    }

    Int_t s = dateTime.GetSecond();
    if(s < 10){
      dateTimeSuffix += TString::Format("-0%d", s);
    }
    else{
      dateTimeSuffix += TString::Format("-%d", s);
    }
  }
  return dateTimeSuffix;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Looks for an environment variable called OUTPUT_DIR and if it exists, sets it as the output dir.
 *
 * @return The output dir.
 */
TString OutputConvention::getOutputDir(){
  
  outputDir = "";
  const char* outputDirPoss = getenv("OUTPUT_DIR");
  if(outputDirPoss!=NULL){
    outputDir += TString::Format("%s/", outputDirPoss); // Add trailing forward slash...
  }

  return outputDir;
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Opens a TFile matching a fileName with wildcards. If multiple matches gets the "greatest" TString, which hopefully corresponds to the file with the latest date suffix.
 *
 * @param fileNameWithWildcards is the name of the file (with wildcards) that you wish to open.
 * @return NULL if no matches, the opened file if there is a match.
 * 
 * Sorts all matching files into increasing order of fileName.
 * If the files have my standard date suffix attached to them, then this should correspond to the most recent file.
 */
TFile* OutputConvention::getFile(TString fileNameWithWildcards){
  // This might be my favourite little bit of stand alone code.
  
  TFile* theFile = NULL;

  TChain* tempChain = new TChain("tempChain");
  tempChain->Add(fileNameWithWildcards);
  
  TObjArray* fileList = tempChain->GetListOfFiles();

  const int numFiles = fileList->GetEntries();
  if(numFiles > 0){
    // std::cerr << "numFiles = " << numFiles << std::endl;
    std::vector<TString> fileNames(numFiles, "");;  

    for(int fileInd=0; fileInd < numFiles; fileInd++){
      fileNames.at(fileInd) = TString::Format("%s", fileList->At(fileInd)->GetTitle());
    }
    std::sort(fileNames.begin(), fileNames.end(), std::greater<TString>());
    for(int fileInd=0; fileInd < numFiles; fileInd++){
      // std::cerr << fileInd << "\t" << fileNames.at(fileInd) << std::endl;
    }
    theFile = TFile::Open(fileNames.at(0));    
  }
  delete tempChain;

  return theFile;

}

