#include "TObject.h"

class AnitaEventSummary;
class TChain;
class TH2D;

namespace Acclaim
{

 class AnalysisPlot;

 class SummarySet {
  public:
   SummarySet(const char* pathToSummaryFiles, const char* treeName = "sumTree", const char* summaryBranchName = "sum");
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

   AnalysisPlot* bookTimeAnalysisPlot(const char* name, const char* title, int nx, int ny, double yMin, double yMax);
   AnalysisPlot* bookEventNumberAnalysisPlot(const char* name, const char* title, int nx, int ny, double yMin, double yMax);

   TH2D* bookTimeHistogram(const char* name, const char* title, int nx, int ny, double yMin, double yMax);
   TH2D* bookEventNumberHistogram(const char* name, const char* title, int nx, int ny, double yMin, double yMax);

  protected:

   void init(const char* pathToSummaryFiles, const char* treeName, const char* summaryBranchName);

   TChain* fChain;
   Long64_t fN;
   AnitaEventSummary* fSum;

   // useful data about the range of the summary data set
   UInt_t fFirstTime;
   UInt_t fFirstEventNumber;
   UInt_t fLastTime;
   UInt_t fLastEventNumber;

 };
}


