#include "AcclaimCorrelationSummary.h"

ClassImp(Acclaim::CorrelationPair)

Acclaim::CorrelationPair::CorrelationPair(){;}
Acclaim::CorrelationPair::~CorrelationPair(){;}

Acclaim::CorrelationPair::CorrelationPair(int a1, int a2, double t, double rho,
					  int pol, double phi, double theta, UInt_t eventNum)
	: pol(pol),
	  ant1(a1),
	  ant2(a2),
	  dt(t),
	  correlation(rho),
	  phiDeg(phi),
	  thetaDeg(theta),
	  eventNumber(eventNum)
{}
