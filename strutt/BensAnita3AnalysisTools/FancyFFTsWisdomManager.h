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

#ifdef __CINT__
#ifdef FFTW_64_BIT // Hack for Hawaii install of FFTW
typedef struct {char a[16];} __float128; /* 16 chars have the same size as one __float128 */
#endif
#endif

#include <iostream>
#include <fftw3.h>
#include "TObject.h"
#include "TSystem.h"
#include <stdlib.h>

class FancyFFTsWisdomManager{
public:
  FancyFFTsWisdomManager();
  ~FancyFFTsWisdomManager();
  
private:
  const char* wisdomDir;

};


#endif //FANCYFFTS_WISDOM_MANAGER_H
