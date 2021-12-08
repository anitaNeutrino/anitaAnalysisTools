


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
   * Because the ROOT TEntryList interface to TChains is quite annoying, and I only want to type it once.
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
    Long64_t getEvent(UInt_t eventNumber);
    bool GetUseProof(){return fUseProof;}
    void SetUseProof(bool useProof=true);



    /**
     * From the thermal chain
     */
    
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
    Float_t deconvolved_filtered_snr;
    Float_t weight;
    Float_t mc_energy;
    Float_t peak_value;
    Float_t coherent_filtered_peakHilbert;
    Float_t deconvolved_filtered_peakHilbert;
    Float_t coherent_filtered_impulsivityMeasure;
    Float_t deconvolved_filtered_impulsivityMeasure;
    Float_t coherent_filtered_fracPowerWindowGradient;
    Float_t deconvolved_filtered_fracPowerWindowGradient;    


    /**
     * From the hical chain
     */
    Int_t duringHiCal;
    Double_t hiCalPhi;
    Double_t hiCalTheta;
    UInt_t eventNumber2;
    Int_t run2;



    /**
     * From the surface chain
     */
    Double_t longitude;
    Double_t latitude;
    Double_t altitude;
    Double_t thetaAdjustmentRequired;
    Int_t onContinent;
    Int_t onIceShelf;
    Double_t iceThickness;
    UInt_t eventNumber3;
    Int_t run3;


    Adu5Pat pat();
    Double_t fisherDiscriminant();
    
  private:
    TChain* fChain;
    TChain* fFriendChain1;
    TChain* fFriendChain2;

    TCut fCut;
    mutable bool fEntryListDirty;
    mutable TEntryList* fEntryList;
    bool fUseProof;
    void makeSelection() const;
    void doTypeConversions();
    
    Int_t fAnitaVersion;

    void setBranches();
    Int_t eventNumberInt;
    Int_t realTimeInt;
    Float_t polFloat;
    Float_t peakIndFloat;
  };
  


}



#endif
