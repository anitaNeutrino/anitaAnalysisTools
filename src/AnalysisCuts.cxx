#include "AnalysisCuts.h"
#include "AnitaEventSummary.h"
#include "RampdemReader.h"

/**
 * Generic constructor, assigns name title and maximum return value
 *
 * @param name is the cut name
 * @param title is the cut title
 * @param mrv is the maximum value it is possible for apply() to return
 */
Acclaim::AnalysisCut::AnalysisCut(const char* name, const char* title, int mrv)
    : fName(name), fTitle(title), fMaxRetVal(mrv)
{
  // just assign name, title, and maximum return value
}



/**
 * Checks if the event was above horizontal
 *
 * @param sum is the AnitaEventSummary
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::IsAboveHorizontal::apply(const AnitaEventSummary* sum) const
{
  return sum->higherPeak().theta > 0 ? 1 : 0;
}




/**
 * Was the event tagged as a WAIS pulser (from timing info)
 *
 * @param sum is the AnitaEventSummary
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::IsTaggedAsWaisPulser::apply(const AnitaEventSummary* sum) const
{
  return sum->flags.pulser == AnitaEventSummary::EventFlags::WAIS;
}


/**
 * Get the higher polarisation, wraps the higherPeakPol() function
 *
 * @param sum is the AnitaEventSummary
 *
 * @return 0 for HPol, 1 for VPol
 */
int Acclaim::HigherPol::apply(const AnitaEventSummary* sum) const
{
  return sum->higherPeakPol();
}

/**
 * Does the event reconstruct to somewhere on the Earth (land or sea)?
 *
 * @param sum is the AnitaEventSummary
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::HasSourceLocation::apply(const AnitaEventSummary* sum) const
{
  bool didReconstruct = (sum->higherPeak().latitude < -900 || sum->higherPeak().theta_adjustment_needed > 0) ? false : true;
  return didReconstruct;
}

/**
 * Does the event reconstruct to somewhere on the Antarctic land mass?
 * Wraps RampdemReader::isOnContinent(lon, lat)
 *
 * @param sum is the AnitaEventSummary
 *
 * @return 1 if true, 0 if false
 */
int Acclaim::IsOnContinent::apply(const AnitaEventSummary* sum) const
{
  return RampdemReader::isOnContinent(sum->higherPeak().longitude, sum->higherPeak().latitude);
}
