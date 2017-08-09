/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             Strictly templated cut functions for the AnalysisPlot class
***********************************************************************************************************/

#ifndef ANALYSIS_CUTS_H
#define ANALYSIS_CUTS_H

#include "TString.h"
#include "AnitaConventions.h"

class AnitaEventSummary;

namespace Acclaim
{


/**
 * @class AnalysisCut is a class to perform a simple cut using the results stored in an AnitaEventSummary
 *
 * The base AnalysisCut class is designed to be inherited from.
 * It contains a purely virtual apply(const AnitaEventSummary* sum) const function.
 * The apply function is then overloaded in the derived classes where the actual cuts are implemented.
 * To work properly with the AnalysisPlot class, apply must return an integer from 0 to fMaxRetVal-1 (inclusive).
 * In addition to the expected return values, the name/title members are also used in the AnalysisPlot for histogram names and titles.
 */

  class AnalysisCut {
   public:
    AnalysisCut(const char* name, const char* title, int mrv);
    virtual int apply(const AnitaEventSummary* sum) const = 0; // to be overloaded with actual cut
    inline int getMaximumReturnValue() const {return fMaxRetVal;}
    inline const char* getName() const {return fName.Data();}
    inline const char* getTitle() const {return fTitle.Data();}
   protected:
    TString fName;
    TString fTitle;
    int fMaxRetVal;
  };

  class IsAboveHorizontal : public AnalysisCut
  {
   public:
    IsAboveHorizontal() : AnalysisCut("isAboveHorizontal", "Above Horizontal", 2) {;}
    virtual int apply(const AnitaEventSummary* sum) const; /// Returns false(0) or true(1)
  };

  class IsTaggedAsWaisPulser : public AnalysisCut
  {
   public:
    IsTaggedAsWaisPulser() : AnalysisCut("isTaggedAsWaisPulser", "Tagged As WAIS Pulser", 2) {;}
    virtual int apply(const AnitaEventSummary* sum) const; /// Returns false(0) or true(1)
  };

  class IsTaggedAsLDBPulser : public AnalysisCut
  {
   public:
    IsTaggedAsLDBPulser() : AnalysisCut("isTaggedAsLDBPulser", "Tagged As LDB Pulser", 2) {;}
    virtual int apply(const AnitaEventSummary* sum) const; /// Returns false(0) or true(1)
  };

  class HigherPol : public AnalysisCut
  {
   public:
    HigherPol() : AnalysisCut("higherPol", "Polarization", AnitaPol::kNotAPol) {;}
    virtual int apply(const AnitaEventSummary* sum) const; /// Returns 0 for HPol, 1 for VPol
  };

  class HasSourceLocation : public AnalysisCut
  {
   public:
    HasSourceLocation() : AnalysisCut("hasSourceLocation", "Has Source Location", 2) {;}
    virtual int apply(const AnitaEventSummary* sum) const; /// Returns false(0) or true(1)
  };

  class IsOnContinent : public AnalysisCut
  {
   public:
    IsOnContinent() : AnalysisCut("IsOnContinent", "Reconstructs to Land", 2) {;}
    virtual int apply(const AnitaEventSummary* sum) const; /// Returns false(0) or true(1)
  };

  class IsTaggedAsPayloadBlast : public AnalysisCut
  {
   public:
    IsTaggedAsPayloadBlast() : AnalysisCut("isTaggedAsPayloadBlast", "Tagged as Payload Blast", 2) {;}
    virtual int apply(const AnitaEventSummary* sum) const; /// Returns false(0) or true(1)
  };

  class IsWithin20DegreesOfSunInPhi : public AnalysisCut{
   public:
    IsWithin20DegreesOfSunInPhi() : AnalysisCut("isWithin20DegreesOfSunInPhi", "|#delta#phi_{sun}| < 20", 2) {;}
    virtual int apply(const AnitaEventSummary* sum) const; /// Returns false(0) or true(1)
  };

  class IsGood : public AnalysisCut{
   public:
    IsGood() : AnalysisCut("isGood", "Is Good", 2) {;}
    virtual int apply(const AnitaEventSummary* sum) const; /// Returns false(0) or true(1)
  };


  class GoodGPS : public AnalysisCut{
   public:
    GoodGPS() : AnalysisCut("goodGPS", "Good GPS", 2) {;}
    virtual int apply(const AnitaEventSummary* sum) const; /// Returns false(0) or true(1)
  };



  // const globals so you don't need to instantiate these yourself
  namespace AnalysisCuts{
    const IsAboveHorizontal isAboveHorizontal;
    const IsTaggedAsWaisPulser isTaggedAsWaisPulser;
    const HigherPol higherPol;
    const HasSourceLocation hasSourceLocation;
    const IsOnContinent isOnContinent;
    const IsTaggedAsPayloadBlast isTaggedAsPayloadBlast;
    const IsWithin20DegreesOfSunInPhi isWithin20DegreesOfSunInPhi;
    const IsGood isGood;
    const GoodGPS goodGPS;
  }
}

#endif
