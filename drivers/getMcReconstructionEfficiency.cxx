#include "SummarySet.h"
#include "OutputConvention.h"
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"

int main(int argc, char* argv[]){

  if(!(argc == 2)){
    std::cerr << argv[0] << " 'mcGlob' " << std::endl;
    return 1;
  }

  
  Acclaim::SummarySet ssMc(argv[1]);
  ssMc.SetUseProof(1);

  Acclaim::OutputConvention oc(1, &argv[0]);
  TFile* fOut = oc.makeFile();

  const char* weight = "sum.weight()";

  
  ssMc.Draw("sum.trainingPeak().dPhiMC()>>hTrainingPeakDPhi(1024, -6, 6)", weight);
  TH1F* hTrainingPeakDPhi = (TH1F*) gROOT->FindObject("hTrainingPeakDPhi");
  if(hTrainingPeakDPhi) {
    hTrainingPeakDPhi->Write();
    // delete hTrainingPeakDPhi;
  }
  
  ssMc.Draw("sum.trainingPeak().dThetaMC()>>hTrainingPeakDTheta(1024, -6, 6)", weight);
  TH1F* hTrainingPeakDTheta = (TH1F*) gROOT->FindObject("hTrainingPeakDTheta");
  if(hTrainingPeakDTheta) {
    hTrainingPeakDTheta->Write();
    // delete hTrainingPeakDTheta;
  }

  ssMc.Draw("sum.highestPeak().dPhiMC()>>hHighestPeakDPhi(1024, -6, 6)", weight);
  TH1F* hHighestPeakDPhi = (TH1F*) gROOT->FindObject("hHighestPeakDPhi");
  if(hHighestPeakDPhi) {
    hHighestPeakDPhi->Write();
    // delete hHighestPeakDPhi;
  }

  ssMc.Draw("sum.highestPeak().dThetaMC()>>hHighestPeakDTheta(1024, -6, 6)", weight);
  TH1F* hHighestPeakDTheta = (TH1F*) gROOT->FindObject("hHighestPeakDTheta");
  if(hHighestPeakDTheta) {
    hHighestPeakDTheta->Write();
    // delete hHighestPeakDTheta;
  }
  
  ssMc.Draw("sum.mostImpulsivePeak(0).dPhiMC()>>hMostImpulsive0DPhi(1024, -6, 6)", weight);
  TH1F* hMostImpulsive0DPhi = (TH1F*) gROOT->FindObject("hMostImpulsive0DPhi");
  if(hMostImpulsive0DPhi) {
    hMostImpulsive0DPhi->Write();
    // delete hMostImpulsive0DPhi;
  }

  ssMc.Draw("sum.mostImpulsivePeak(0).dThetaMC()>>hMostImpulsive0DTheta(1024, -6, 6)", weight);
  TH1F* hMostImpulsive0DTheta = (TH1F*) gROOT->FindObject("hMostImpulsive0DTheta");
  if(hMostImpulsive0DTheta) {
    hMostImpulsive0DTheta->Write();
    // delete hMostImpulsive0DTheta;
  }
  
  ssMc.Draw("sum.mostImpulsivePeak(1).dPhiMC()>>hMostImpulsive1DPhi(1024, -6, 6)", weight);
  TH1F* hMostImpulsive1DPhi = (TH1F*) gROOT->FindObject("hMostImpulsive1DPhi");
  if(hMostImpulsive1DPhi) {
    hMostImpulsive1DPhi->Write();
    // delete hMostImpulsive1DPhi;
  }

  ssMc.Draw("sum.mostImpulsivePeak(1).dThetaMC()>>hMostImpulsive1DTheta(1024, -6, 6)", weight);
  TH1F* hMostImpulsive1DTheta = (TH1F*) gROOT->FindObject("hMostImpulsive1DTheta");
  if(hMostImpulsive1DTheta) {
    hMostImpulsive1DTheta->Write();
    // delete hMostImpulsive1DTheta;
  }  
  

  fOut->Write();
  fOut->Close();

  gSystem->Exit(0);
  return 0;
}
