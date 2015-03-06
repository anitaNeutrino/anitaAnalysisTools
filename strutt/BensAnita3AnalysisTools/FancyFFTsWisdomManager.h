/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             What I really want is a FancyFFTs constructor or destructor.
	     But since I made it a static class I need another class which has a constructor.
	     Which I can make a static instance of inside the FancyFFT class.
*************************************************************************************************************** */


#ifndef FANCYFFTS_WISDOM_MANAGER_H
#define FANCYFFTS_WISDOM_MANAGER_H

#include <iostream>
#include <fftw3.h>
#include <TObject.h>
#include <TSystem.h>

class FancyFFTsWisdomManager : public TObject{
public:
  FancyFFTsWisdomManager();
  ~FancyFFTsWisdomManager();
  
private:
  int wisdomImportSuccess;
  const char* wisdomDir;

  ClassDef(FancyFFTsWisdomManager, 0);
};


#endif //FANCYFFTS_WISDOM_MANAGER_H
