/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             My ANITA-3 data quality cuts, to remove bad events
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
   * Some of my filters use rolling averages so "bad events" could ruin them.
   * Therefore apply() acts on a UsefulAnitaEvent rather than a FilteredAnitaEvent so that we can decide whether or not we
   * want to exponse the a FilterStrategy with an internal state to a bad event.
   * 
   */
  class QualityCut : public TObject {

  public:
    static Bool_t applyAll(const UsefulAnitaEvent* usefulEvent, AnitaEventSummary* sum=NULL); // static utility function, applies all cuts defined here
    static Bool_t passedAll(const AnitaEventSummary* sum, bool describe = false); // static utility function, reads flags in AnitaEventSummary, returns true if event passed
    

    
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
  


  /** 
   * @class NumPointsCut
   * @brief Removes events which have short waveforms
   *
   * Very simple, just counts the number of points in each chanel of a UsefulAnitaEvent
   * 
   */
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
