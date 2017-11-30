#include "TObject.h"
#include "TString.h"
#include "TChain.h"
#include "AnitaConventions.h"
#include "TGraphAntarctica.h"

// This is a guess at the version number, if this doesn't work for you
// feel free to try harder to track down the change
#if ROOT_VERSION_CODE < ROOT_VERSION(6,0,8)
#define TCHAIN_NENTRIES_DEFAULT TChain::kBigNumber
#else
#define TCHAIN_NENTRIES_DEFAULT TChain::kMaxEntries
#endif   

class AnitaEventSummary;
class TH2D;
class TH2DAntarctica;
class TProfile2DAntarctica;
class TProof;

namespace Acclaim
{

  class AnalysisPlot;
  class AnalysisProf;

  class SummarySet {
  public:
    SummarySet(const char* pathToSummaryFiles, const char* treeName = "sumTree", const char* summaryBranchName = "sum", bool useProof = false);
    virtual ~SummarySet();

    Long64_t N() const {return fN;}
    AnitaEventSummary * summary() const { return fSum;}

    Long64_t getEntry(Long64_t entry);
    Long64_t getEvent(UInt_t eventNumber);
    Long64_t first(){return getEntry(0);}
    Long64_t last(){return fN > 0 ? getEntry(fN-1) : -1;}

    UInt_t getFirstTime(){return fFirstTime;}
    UInt_t getLastTime(){return fLastTime;}
    UInt_t getFirstEventNumber(){return fFirstEventNumber;}
    UInt_t getLastEventNumber(){return fLastEventNumber;}

    AnalysisProf* bookTimeAnalysisProf(const char* name, const char* title, int nx, int ny, double yMin, double yMax);
    AnalysisProf* bookEventNumberAnalysisProf(const char* name, const char* title, int nx, int ny, double yMin, double yMax);

    AnalysisPlot* bookTimeAnalysisPlot(const char* name, const char* title, int nx, int ny, double yMin, double yMax);
    AnalysisPlot* bookEventNumberAnalysisPlot(const char* name, const char* title, int nx, int ny, double yMin, double yMax);

    TH2D* bookTimeHistogram(const char* name, const char* title, int nx, int ny, double yMin, double yMax);
    TH2D* bookEventNumberHistogram(const char* name, const char* title, int nx, int ny, double yMin, double yMax);


    Long64_t Process(TSelector* selector, Option_t* option = "", Long64_t nentries = TCHAIN_NENTRIES_DEFAULT, Long64_t firstentry = 0);

    TProfile2DAntarctica* makeAntarcticaProf(AnitaPol::AnitaPol_t pol=AnitaPol::kNotAPol, const char* name="", const char* title="", Int_t nx=-1, Int_t ny=-1);
    TH2DAntarctica*       makeAntarcticaHist(AnitaPol::AnitaPol_t pol=AnitaPol::kNotAPol, const char* name="", const char* title="", Int_t nx=-1, Int_t ny=-1);  


    Double_t getTotalSize() const;
    TGraphAntarctica* makePayloadLocationGraph(int stride = TGraphAntarctica::defaultGpsTreeStride);


    // If I keep adding these wrapper functions at some point it might make sense to inherit from TChain...
    TChain* getChain(){return fChain;}
    Long64_t Draw(const char* varexp, const TCut &selection, Option_t *option = "", Long64_t nentries = TCHAIN_NENTRIES_DEFAULT, Long64_t firstentry = 0);
    Long64_t Draw(const char* varexp,const char* selection="", Option_t* option = "", Long64_t nentries = TCHAIN_NENTRIES_DEFAULT, Long64_t firstentry = 0);

    /** 
     * Returns the approximate mean size of each SummaryTree entry
     * 
     * @return effective AnitaEventSummary branch size
     */
    Double_t getBytesPerEvent() const { 
      return N() > 0 ? getTotalSize()/N() : 0;
    }
    void SetUseProof(bool useProof=true) {fUseProof = useProof;}
    Bool_t GetUseProof() {return fUseProof;}



  protected:

    void init();
    void initProof();
    void renameProofCanvas(const char* varexp);

    TString fPathToSummaryFiles;	/// The glob passed to the TChain
    TString fTreeName;			/// The name of the TTrees, default is "sumTree"
    TString fSummaryBranchName;		/// Branch name of the AnitaEventSummary, default is "sum"

    TChain* fChain;			/// The chain of TTrees containing the AnitaEventSummary
    Long64_t fN;			/// The number of entries in the AnitaEventSummary chain, access with SummarySet::N()
    AnitaEventSummary* fSum;		/// Pointer to the current entry in the chain, access with SummarySet::summary()

    UInt_t fFirstTime;			/// The realTime of the first entry in the summary chain, useful for booking histograms
    UInt_t fFirstEventNumber;		/// The eventNumber of the first entry in the summary chain, useful for booking histograms
    UInt_t fLastTime;			/// The realTime of the last entry in the summary chain, useful for booking histograms
    UInt_t fLastEventNumber;		/// The eventNumber of the last entry in the summary chain, useful for booking histograms
    Bool_t fUseProof;			/// Switch on the Parallel ROOT Facility, for speedy histogram plotting
    TProof* fProof;			/// Pointer to the PROOF session
    Bool_t fBuiltIndex;                 /// Built chain index?
  };
}


