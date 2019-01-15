#include "AcclaimCorrelationSummary.h"

ClassImp(Acclaim::CorrelationPair)
ClassImp(Acclaim::CorrelationSummary)


Acclaim::CorrelationPair::CorrelationPair()
: pol(-1),
  ant1(-1),
  ant2(-1),
  dt(uninitialized),
  dt_expected(uninitialized),
  correlation(uninitialized)
{;}

Acclaim::CorrelationPair::~CorrelationPair(){;}

Acclaim::CorrelationPair::CorrelationPair(int a1, int a2, double t, double rho,
					  int pol, double phi, double theta, UInt_t eventNum)
	: pol(pol),
	  ant1(a1),
	  ant2(a2),
	  dt(t),
	  dt_expected(uninitialized),
	  correlation(rho)
	  // phiDeg(phi),
	  // thetaDeg(theta),
	  // eventNumber(eventNum)
{

}




Acclaim::CorrelationSummary::CorrelationSummary(){;}

Acclaim::CorrelationSummary::~CorrelationSummary(){;}

Acclaim::CorrelationSummary::CorrelationSummary(AnitaPol::AnitaPol_t pol, UInt_t eventNumber,  double phiDeg, double thetaDeg, int phiSector, Adu5Pat* pat)
      : fPol(static_cast<int>(pol)),
	fEventNumber(eventNumber),
	fPhiSector(phiSector),
	fPhiDeg(phiDeg),
	fThetaDeg(thetaDeg),
	fPat(*pat)
    {}
