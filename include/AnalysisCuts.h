/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             Strictly templated cut functions for the AnalysisPlot class
***********************************************************************************************************/

#ifndef ANALYSIS_CUTS_H
#define ANALYSIS_CUTS_H

class AnitaEventSummary;

namespace Acclaim
{


/** 
 * @namespace AnalysisCuts is a set of cut functions designed for use with the AnalysisPlot class
 * 
 * In order to be of use for the AnalysisPlot class each function must:
 *    Return an int, that runs from 0 to some maximum value (that is hopefully small).
 *    Take a single argument, a pointer to a const AnitaEventSummary.
 *    In the case that a NULL pointer is passed, return the maximum possible return value.
 */

  namespace AnalysisCuts {


  int isAboveHorizontal(const AnitaEventSummary* sum); /// /// Returns false(0) or true(1)
  int isTaggedAsWaisPulser(const AnitaEventSummary* sum); /// /// Returns false(0) or true(1)
  int higherPol(const AnitaEventSummary* sum); /// Returns 0 for HPol, 1 for VPol 
  int hasSourceLocation(const AnitaEventSummary* sum); /// Returns false(0) or true(1)
  int isOnContinent(const AnitaEventSummary* sum); /// Returns false(0) or true(1)



}
}

#endif
