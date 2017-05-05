/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             My ANITA-3 data quality cuts, to remove bad events that mess up things like averages
***********************************************************************************************************/

#ifndef QUALITYCUTS_H
#define QUALITYCUTS_H

#include "TObject.h"
#include "TString.h"

#include "FilteredAnitaEvent.h"
#include "AnalysisWaveform.h"
#include "RootTools.h"

#include "AnitaEventSummary.h"

namespace Acclaim{

  /**
   * Event selection cut
   */
  class QualityCut : public TObject {

  public:
    static Bool_t applyAll(const UsefulAnitaEvent* usefulEvent, AnitaEventSummary* sum=NULL); // static utility function, applies all cuts defined here

    
    QualityCut(){;}
    virtual ~QualityCut(){;}
    // virtual void apply(FilteredAnitaEvent* fEv) = 0;
    virtual void apply(const UsefulAnitaEvent* useful, AnitaEventSummary* sum = NULL) = 0;    
    TString description;
    Bool_t eventPassesCut;
    ClassDef(QualityCut, 1)
  };



  class SurfSaturationCut : public QualityCut {
    ClassDef(Acclaim::SurfSaturationCut, 1);

  protected:
    double maxLimit;
    double minLimit;
    double asymLimit;

    double maxVolts;
    int maxVoltsAnt;
    AnitaPol::AnitaPol_t maxVoltsPol;

    double minVolts;
    int minVoltsAnt;
    AnitaPol::AnitaPol_t minVoltsPol;
    
    double asymVolts;
    int asymVoltsAnt;
    AnitaPol::AnitaPol_t asymVoltsPol;

  public:
    SurfSaturationCut();
    virtual ~SurfSaturationCut(){;}    
    // virtual void apply(FilteredAnitaEvent* fEv);
    virtual void apply(const UsefulAnitaEvent* useful, AnitaEventSummary* sum = NULL);    
  };




  

  class SelfTriggeredBlastCut : public QualityCut {
    ClassDef(Acclaim::SelfTriggeredBlastCut, 1);    
    
  protected:
    double ratioCutHigh;
    double ratioCutLow;


    double maxRatio;
    int maxRatioPhi;
    AnitaPol::AnitaPol_t maxRatioPol;

  public:
    SelfTriggeredBlastCut();   
    virtual ~SelfTriggeredBlastCut(){;}
    // virtual void apply(FilteredAnitaEvent* fEv);
    virtual void apply(const UsefulAnitaEvent* useful, AnitaEventSummary* sum = NULL);
  };
  
  
  class NumPointsCut : public QualityCut {
    ClassDef(Acclaim::NumPointsCut, 1);    
    
  protected:
    int numPointsCutLow; // presumably don't need a cut on the other side of this...
    int minNumPoints;

  public:
    NumPointsCut();
    virtual ~NumPointsCut(){;}
    // virtual void apply(FilteredAnitaEvent* fEv);
    virtual void apply(const UsefulAnitaEvent* useful, AnitaEventSummary* sum = NULL);
  };
  
}

#endif //ANALYSISCUTS
