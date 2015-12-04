/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             What I really want is a FancyFFTs constructor or destructor.
	     But since I made it a static class I need another class which has a constructor.
	     Which I can make a static instance of inside the FancyFFT class.
*************************************************************************************************************** */

#include "FancyFFTsWisdomManager.h"

FancyFFTsWisdomManager::FancyFFTsWisdomManager(){

  // std::cout << "FancyFFTsWisdomManager::FancyFFTsWisdomManager()" << std::endl;  
  const char* anitaUtilEnv = "ANITA_UTIL_INSTALL_DIR";
  const char* wisdomDirTemp = getenv(anitaUtilEnv);
  char* wisdomDirTemp2 = new char[FILENAME_MAX];
  sprintf(wisdomDirTemp2, "%s/etc", wisdomDirTemp);
  wisdomDir = wisdomDirTemp2;
  FILE* wisdomFile = NULL;
  
  if(wisdomDir==NULL){// no environment variable
    std::cerr << "Warning in FancyFFTsWisdomManager! Please set environment variable " 
	      <<  anitaUtilEnv << " so I know where to save fftw3 wisdom." << std::endl;
  }
  else{
    char wisdomFileName[FILENAME_MAX];
    sprintf(wisdomFileName, "%s/fftw.wisdom", wisdomDir);

    /* ... sigh... is more backwards compatible but doing even more c style*/
    wisdomFile = fopen(wisdomFileName, "r");
    if(wisdomFile!=NULL){
      // wisdomImportSuccess = fftw_import_wisdom_from_filename(wisdomFileName);
      fftw_import_wisdom_from_file(wisdomFile);
      fclose(wisdomFile);
    }
  }

  if(wisdomFile==NULL){//then no file
    std::cerr << "Warning in FancyFFTsWisdomManager! Failed to import fftw wisdom, "
	      << "creating plans may take some time." << std::endl;
  }
  
}


FancyFFTsWisdomManager::~FancyFFTsWisdomManager(){

  // std::cout << "FancyFFTsWisdomManager::~FancyFFTsWisdomManager()" << std::endl;  
  if(wisdomDir != NULL){
    char wisdomFileName[FILENAME_MAX];
    sprintf(wisdomFileName, "%s/fftw.wisdom", wisdomDir);

    /* ... sigh... is more backwards compatible but doing even more c style*/
    FILE* wisdomFile = fopen(wisdomFileName, "w");
    if(wisdomFile!=NULL){
      fftw_export_wisdom_to_file(wisdomFile);
      fclose(wisdomFile);
    }
    else{
      std::cerr << "Warning in FancyFFTsWisdomManager! Failed to export fftw wisdom." << std::endl;
    }
  }
}