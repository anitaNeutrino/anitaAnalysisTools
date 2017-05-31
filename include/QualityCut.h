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
   * @class QualityCut
   * @brief Base class from which all my quality cuts inherit
   *
   * Contains a description, boolian to indicate pass/fail and apply function which must be implemented in descendents
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



  /** 
   * @class SurfSaturationCut
   * @brief Tries to find unphysical spikes in a waveform that are characteristic of a kind of digitiser corruption
   *
   * Currently comparses the maximum, minimum voltages and the difference in their magnitude to try to find spikes
   */
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




  
  /** 
   * @class PayloadBlastCut
   * @brief Tries to determine whether an event is a payload blast
   *
   * Payload blasts tend to have lots of power in the bottom ring, but not the top.
   * The cut comparse the peak-to-peak voltages in the top and bottom rings in matching phi-sectors
   * 
   */
  class PayloadBlastCut : public QualityCut {
    ClassDef(Acclaim::PayloadBlastCut, 1);    
    
  protected:
    double ratioCutHigh;
    double ratioCutLow;


    double maxRatio;
    int maxRatioPhi;
    AnitaPol::AnitaPol_t maxRatioPol;

  public:
    PayloadBlastCut();   
    virtual ~PayloadBlastCut(){;}
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
