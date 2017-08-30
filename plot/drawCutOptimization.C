#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "CutOptimizer.h"
#include <iostream>

void drawFisherPlot(TFile* f){
  TH1D* hSigInt = (TH1D*) f->Get("hSigInt");
  TH1D* hBackInt = (TH1D*) f->Get("hBackInt");

  hSigInt->SetLineColor(kRed);
  hBackInt->SetLineColor(kBlue);

  auto c1 = new TCanvas();
  hSigInt->Draw();
  hBackInt->Draw("same");

  hSigInt->SetTitle("Fisher Discriminant;Fisher Score;Fraction of events above value");
  hSigInt->GetXaxis()->SetRangeUser(-15, 25);
  hSigInt->GetYaxis()->SetRangeUser(1e-09, 2);
  auto l1 = new TLegend(0.65, 0.75, 0.99, 0.99);
  l1->AddEntry(hSigInt, "MC", "l");
  l1->AddEntry(hBackInt, "Upward pointing thermal", "l");
  l1->Draw();
  
  c1->SetLogy(1);
}


void drawEfficiencies(TFile* f, bool makeFisher){

  TList* l = f->GetListOfKeys();
  std::vector<TEfficiency*> effSnr;
  std::vector<TEfficiency*> effEnergy;

  for(int i=0; i < l->GetEntries(); i++){
    TString name = l->At(i)->GetName();

    bool effPrefix = name.Contains("eff_");
    if(effPrefix && name.Contains("_vs_SNR")){
      effSnr.push_back((TEfficiency*)f->Get(name));
    }
    else if(effPrefix && name.Contains("_vs_Energy")){
      effEnergy.push_back((TEfficiency*)f->Get(name));
    }
  }

  const int nC = 7;
  EColor colors[nC] = {kRed, kBlue, kGreen, kYellow, kMagenta, kOrange, kSpring};
  
  auto c1 = new TCanvas();
  auto l1 = new TLegend(0.8, 0.8, 1, 1);
  for(unsigned i=0; i < effSnr.size(); i++){
    const char* opt = i == 0 ? "" : "same";
    effSnr.at(i)->Draw(opt);
    effSnr.at(i)->SetLineColor(colors[i%nC]);
    l1->AddEntry(effSnr.at(i), effSnr.at(i)->GetName(), "l");
  }
  l1->Draw();

  auto c2 = new TCanvas();
  auto l2 = new TLegend(0.8, 0.8, 1, 1);
  
  for(unsigned i=0; i < effEnergy.size(); i++){
    const char* opt = i == 0 ? "" : "same";
    effEnergy.at(i)->Draw(opt);
    effEnergy.at(i)->SetLineColor(colors[i%nC]);
    l2->AddEntry(effSnr.at(i), effEnergy.at(i)->GetName(), "l");

    if(i==0){
      effEnergy.at(i)->SetTitle("Efficiency of pre-thermal cuts; log10(Energy) eV; Efficiency");
    }
  }  
  l2->Draw();
  
  if(makeFisher){

    Acclaim::CutOptimizer::FisherResult* fr = (Acclaim::CutOptimizer::FisherResult*) f->Get("FisherResult");
    TString cutForm = fr->getFisherFormula();
    const TH1* h_energy = effEnergy.at(0)->GetPassedHistogram();
    double xLow = h_energy->GetBinLowEdge(1);
    double xHigh = h_energy->GetBinLowEdge(h_energy->GetNbinsX());
    const int nx = h_energy->GetNbinsX();
    TString name = "eff_fisher_vs_Energy";
    TTree* signalTree = (TTree*) f->Get("signalTree");


    TString namePassed = name + "_passed";
    TH1D* hPassed = new TH1D(namePassed, namePassed, nx, xLow, xHigh);

    TString passedCommand = "TMath::Log10(mc_energy)>>" + namePassed;
    TString cutCommand = "weight*(" + cutForm + "> 7.5)";
    signalTree->Draw(passedCommand, cutCommand, "goff");
    
    TString nameTotal = name + "_total";
    TH1D* hTotal = new TH1D(nameTotal, nameTotal, nx, xLow, xHigh);

    std::cout << hPassed->Integral() << "\t" << hTotal->Integral() << std::endl;

    TString totalCommand = "TMath::Log10(mc_energy)>>" + nameTotal;
    
                                                                      signalTree->Draw(totalCommand, "weight", "goff");
    
    TEfficiency* effE = new TEfficiency(*hPassed, *hTotal);
effE->SetName(name);
                                                effE->SetTitle("Fisher effiency?");
                                        effE->SetTitle("Thermal cut efficiency at Fisher Score of 7.5; log10(Energy) eV; Efficiency");
                                        
    auto c1 = new TCanvas();
    effE->Draw();
  }
  
  
}

void drawCutOptimization(const char* fileName){

  TFile* f = TFile::Open(fileName);
  drawFisherPlot(f);
  drawEfficiencies(f, true);  

  
}
