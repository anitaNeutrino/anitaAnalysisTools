#include "TObject.h"

class AnitaEventSummary;
class TChain;

namespace Acclaim
{

 class SummarySet {
  public:
   SummarySet(const char* pathToSummaryFiles, const char* treeName = "sumTree", const char* summaryBranchName = "sum");
   virtual ~SummarySet();

   Long64_t N() const {return fN;}
   AnitaEventSummary * summary() const { return fSum;}

   Long64_t getEntry(Long64_t entry);
   Long64_t first(){return getEntry(0);}
   Long64_t last(){return fN > 0 ? getEntry(fN-1) : -1;}

  protected:
   TChain* fChain;
   Long64_t fN;
   AnitaEventSummary* fSum;

 };
}


