


#ifndef ACCLAIM_THERMAL_CHAIN_H
#define ACCLAIM_THERMAL_CHAIN_H

#include "TCut.h"

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
    
  private:
    TChain* fChain;
    TCut fCut;
    mutable bool fEntryListDirty;
    mutable TEntryList* fEntryList;

    void makeSelection() const;
  };
  


}



#endif
