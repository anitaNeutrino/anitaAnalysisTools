#include "SummarySet.h"
#include "ProgressBar.h"
#include "OutputConvention.h"

#include "AnitaEventSummary.h"

#include "TH2DAntarctica.h"
#include "TGraphAntarctica.h"

#include <iostream>


int main(int argc, char* argv[]){


  if(!(argc == 2)){
    std::cerr << argv[0] << " 'AnitaEventSummaryGlob' " << std::endl;
    return 1;
  }
  const char* summaryGlob = argv[1];


  Acclaim::SummarySet ss(summaryGlob);
  Long64_t n = ss.N();

  if(n > 0){

    Acclaim::OutputConvention oc(1, argv); // ignore command line args...
    TFile* fOut = oc.makeFile();
    

    TH2DAntarctica* hEvents[AnitaPol::kNotAPol];
    TProfile2DAntarctica* profEvents[AnitaPol::kNotAPol];
    TGraphAntarctica* grFlightPath = new TGraphAntarctica();
    grFlightPath->SetName("grFlightPath");
    
    for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
      TString name = "Event";
      name.Append(AnitaPol::polAsChar(pol));

      TString title = "";
      title.Append(AnitaPol::polAsChar(pol));
      title += "Pol";
      
      hEvents[pol] = new TH2DAntarctica(TString("h") + name, title + " event count");
      profEvents[pol] = new TProfile2DAntarctica(TString("prof") + name, title + " map peak profile");
    }
    
    Acclaim::ProgressBar p(n);
    for(Long64_t entry=0; entry < n; entry++){
      ss.getEntry(entry);
      AnitaEventSummary* sum = ss.summary();

      for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
        AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
        for(int i=0; i < sum->nPeaks[pol]; i++){
          hEvents[pol]->Fill(sum->peak[pol][i].longitude, sum->peak[pol][i].latitude);
          profEvents[pol]->Fill(sum->peak[pol][i].longitude, sum->peak[pol][i].latitude, sum->peak[pol][i].value);
        }

      }
      if((entry % TGraphAntarctica::defaultGpsTreeStride) == 0){
        grFlightPath->SetPoint(grFlightPath->GetN(), sum->anitaLocation.longitude, sum->anitaLocation.latitude);
      }      
      p.inc(entry, n);
    }



    
    grFlightPath->Write();
    delete grFlightPath;

    for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
      hEvents[pol]->Write();
      delete hEvents[pol];
      
      profEvents[pol]->Write();
      delete profEvents[pol];
    }    

    
    fOut->Write();
    fOut->Close();    
  }

  return 0;
}

