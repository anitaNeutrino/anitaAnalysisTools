


#ifndef ACCLAIM_THERMAL_CHAIN_H
#define ACCLAIM_THERMAL_CHAIN_H

#include "TCut.h"
#include "AnitaConventions.h"
#include "Adu5Pat.h"

class TChain;
class TEntryList;

namespace Acclaim {


  /**
   * @class ThermalChain
   * @brief A class to handle loading Thermal Trees into chains and applying cuts
   * 
   * Because the ROOT TEntryList interface to TChains is quite annoying
   */

  class ThermalChain {

  public:
    ThermalChain(const char* glob, const char* treeName="thermalTree");
    virtual ~ThermalChain();
    void addCut(const TCut& cut);
    void addCut(const char* cut);
    void setCut(const TCut& cut);
    void setCut(const char* cut);

    TChain* getChain() const {return fChain;}
    Long64_t N() const;

    Long64_t getEntry(Long64_t entry);
    bool GetUseProof(){return fUseProof;}
    void SetUseProof(bool useProof=true);


    Int_t run;
    UInt_t eventNumber;
    UInt_t realTime;
    AnitaPol::AnitaPol_t pol;
    Int_t peakInd;
    Float_t peak_phi;
    Float_t peak_theta;
    Float_t anita_longitude;
    Float_t anita_latitude;
    Float_t anita_altitude;
    Float_t anita_heading;
    Float_t coherent_filtered_snr;
    Float_t weight;

    Adu5Pat pat();
    
  private:
    TChain* fChain;
    TCut fCut;
    mutable bool fEntryListDirty;
    mutable TEntryList* fEntryList;
    bool fUseProof;
    void makeSelection() const;

    void setBranches();
    Int_t eventNumberInt;
    Int_t realTimeInt;
    Float_t polFloat;
    Float_t peakIndFloat;
  };
  


}



#endif
