/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             Finds the maximum value of cross-correlations.
***********************************************************************************************************/

#ifndef ACCLAIM_CORRELATION_SUMMARY_H
#define ACCLAIM_CORRELATION_SUMMARY_H

#include "TObject.h"
#include "AnitaConventions.h"
#include "Adu5Pat.h"

namespace Acclaim
{

  /**
   * @class CorrelationPair
   * @brief Just a struct to hold the correlations for an individual antenna pair
   * 
   * Should have enough info for a TTree::Draw to be useful
   */

  class CorrelationPair {
  public:
    CorrelationPair();
    CorrelationPair(int a1, int a2, double t, double rho,
		    int pol, double phi, double theta, UInt_t eventNum);
    virtual ~CorrelationPair();
    int pol;
    int ant1;
    int ant2;
    double dt;
    double correlation;
    // double phiDeg;
    // double thetaDeg;
    // UInt_t eventNumber;

    ClassDef(CorrelationPair,  2);
  };

  

  /**
   * @class CorrelationSummary
   * @brief For calibration, holds a set of dts.
   *
   */
  class CorrelationSummary {

  public:
    CorrelationSummary();

    CorrelationSummary(AnitaPol::AnitaPol_t pol, UInt_t eventNumber,  double phiDeg, double thetaDeg, int phiSector, Adu5Pat* pat);
    
    virtual ~CorrelationSummary();
    
    inline std::size_t add(int ant1, int ant2, double dt, double rho){
      fPairs.emplace_back(CorrelationPair(ant1, ant2, dt, rho, fPol, fPhiDeg, fThetaDeg, fEventNumber));
      return N();
    }

    inline std::size_t N() const {
      return fPairs.size();
    }
    

    inline const CorrelationPair& get(int i) const {
      return fPairs.at(i);
    }

    
    Int_t fPol = -1;
    UInt_t fEventNumber = 0;
    int fPhiSector = -1;
    double fPhiDeg = -999;
    double fThetaDeg = -999;
    Adu5Pat fPat;
    std::vector<CorrelationPair> fPairs;

    ClassDef(CorrelationSummary, 2);

  };
}
#endif
