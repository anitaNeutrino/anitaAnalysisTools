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
 * @return the output file name
 */
TString OutputConvention::getOutputFileName(){

  if(outFileName==""){

    // Directory
    outFileName = getOutputDir();

    // Excutable name
    outFileName += TString::Format("%sPlots", argv[0]);

    // Excutable args
    if(argc > 1){
      for(int argInd=1; argInd < argc; argInd++){
	outFileName += TString::Format("_%s", argv[argInd]);
      }
    }
    
    // Date and time of running executable 
    outFileName += getDateTimeSuffix();

    // ROOT suffix
    outFileName += ".root";
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

  // TString outputDirPlusSubDir = outputDir;
  // if(subDir!=""){
  //   outputDirPlusSubDir += subDir + "/";
  // }

  // // Try to make output dir plus subDir
  // int status = mkdir(outputDirPlusSubDir.Data() , S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  // if(status!=0 && errno!=EEXIST){
  //   // If that fails try to make output dir only
  //   status = mkdir(outputDir.Data() , S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  //   if(status!=0 && errno!=EEXIST){
  //     // If that fails then give up and set output dir to nothing (i.e. the pwd)
  //     outputDir = "";
  //   }
  // }
  // else{ // If subdir was fine then we can return it.
  //   outputDir = outputDirPlusSubDir;
  // }

  return outputDir;
}
