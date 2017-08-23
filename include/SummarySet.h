#include "TObject.h"
#include "TString.h"
#include "TChain.h"

// This is a guess at the version number, if this doesn't work for you
// feel free to try harder to track down the change
#if ROOT_VERSION_CODE < ROOT_VERSION(6,0,8)
#define TCHAIN_NENTRIES_DEFAULT TChain::kBigNumber
#else
#define TCHAIN_NENTRIES_DEFAULT TChain::kMaxEntries
#endif   

class AnitaEventSummary;
class TH2D;

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


   // If I keep adding these wrapper functions at some point it might make sense to inherit from TChain...
   TChain* getChain(){return fChain;}
   Long64_t Draw(const char* varexp, const TCut &selection, Option_t *option = "", Long64_t nentries = TCHAIN_NENTRIES_DEFAULT, Long64_t firstentry = 0);
   Long64_t Draw(const char* varexp,const char* selection="", Option_t* option = "", Long64_t nentries = TCHAIN_NENTRIES_DEFAULT, Long64_t firstentry = 0);

   void SetUseProof(bool useProof=true) {fUseProof = useProof;}
   Bool_t GetUseProof() {return fUseProof;}

  protected:

   void init();
   void initProof();
   TString fPathToSummaryFiles;
   TString fTreeName;
   TString fSummaryBranchName;

   TChain* fChain;
   Long64_t fN;
   AnitaEventSummary* fSum;

   // useful data about the range of the summary data set
   UInt_t fFirstTime;
   UInt_t fFirstEventNumber;
   UInt_t fLastTime;
   UInt_t fLastEventNumber;
   Bool_t fUseProof;
 };
}


