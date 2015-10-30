#include "OutputConvention.h"

OutputConvention::OutputConvention(int argcIn, char* argvIn[]){
  dateTime = TDatime();
  argc = argcIn;
  argv = argvIn;
  outFileName = "";
  dateTimeSuffix = "";
  outputDir = "";  
}

OutputConvention::~OutputConvention(){
  
}


TString OutputConvention::getOutputFileName(){

  if(outFileName==""){

    // Directory
    outFileName = getOutputDir();

    // Excutable name
    outFileName += TString::Format("%sPlots", argv[0]);

    // Excutable args
    for(int argInd=1; argInd < argc; argInd++){
      outFileName += TString::Format("_%s", argv[argInd]);
    }

    // Date and time of running executable 
    outFileName += getDateTimeSuffix();

    // ROOT suffix
    outFileName += ".root";
  }

  
  return outFileName;

  
}

TString OutputConvention::getDateTimeSuffix(){
  if(dateTimeSuffix==""){
    dateTimeSuffix = TString::Format("_%d-%d-%d_%d-%d-%d",
				     dateTime.GetYear(),
				     dateTime.GetMonth(),
				     dateTime.GetDay(),
				     dateTime.GetHour(),
				     dateTime.GetMinute(),
				     dateTime.GetSecond()
				     );
  }
  return dateTimeSuffix;
}


TString OutputConvention::getOutputDir(){
  
  outputDir = "";
  const char* outputDirPoss = getenv("OUTPUT_DIR");
  if(outputDirPoss!=NULL){
    outputDir += TString::Format("%s/", outputDirPoss); // Add trailing forward slash...
  }  
  return outputDir;
}
