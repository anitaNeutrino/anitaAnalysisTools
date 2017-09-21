#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "CutOptimizer.h"
#include <iostream>
#include "TF1.h"
#include "TGraph.h"
#include "TKey.h"
#include "RootTools.h"
#include "TH2D.h"

const int numRfTriggersA3 = 77e6;
const double numDesiredBackground = 0.5; //ish
const double backgroundAcceptance = numDesiredBackground/numRfTriggersA3;
double fisherCutVal = 999;

Acclaim::CutOptimizer::FisherResult* fr = NULL;


double getLowestNonZeroBin(TH1D* h){
  double min = DBL_MAX;
  for(int i=1; i <= h->GetNbinsX(); i++){
    double y = h->GetBinContent(i);
    if(y > 0 && y < min){
      min = y;
    }
  }
  return min;
}


void overlayOneDimDists(TFile* f){
  TList* l = f->GetListOfKeys();

  std::set<TString> suffixes;
  const TString sigPref = "hsignalTree_";
  const TString backPref = "hbackgroundTree_";  
  
  for(int i=0; i < l->GetEntries(); i++){

    TKey* k = dynamic_cast<TKey*>(l->At(i));
    if(k && strlen(k->GetName()) > 0){
      // std::cerr << k->GetName() << std::endl;

      TString name = k->GetName();
      if(name.Contains(sigPref) || name.Contains(backPref)){
        name.ReplaceAll(sigPref, "");
        name.ReplaceAll(backPref, "");

        suffixes.insert(name);
      }      
    }    
  }

  std::set<TString>::iterator it = suffixes.begin();
  for(; it != suffixes.end(); ++it){
    TString suf = *it;
    TH2D* h2s = (TH2D*) f->Get(sigPref + *it);
    TH1D* hs = h2s->ProjectionY();
    TH2D* h2b = (TH2D*) f->Get(backPref + *it);
    TH1D* hb = h2b->ProjectionY();

    hs->Scale(1./hs->Integral());
    hb->Scale(1./hb->Integral());

    auto c1 = new TCanvas();
    c1->SetLogy(1);
    hs->SetLineColor(kRed);
    hb->SetLineColor(kBlue);

    hb->Draw("hist");
    hs->Draw("histsame");

    auto l1 = new TLegend(0.8, 0.8, 1, 1);
    l1->AddEntry(hs, "MC neutrinos", "l");
    l1->AddEntry(hb, "Upward pointing RF", "l");
    l1->Draw();

    double max = TMath::Max(hs->GetBinContent(hs->GetMaximumBin()),
                            hb->GetBinContent(hb->GetMaximumBin()));

    double min = TMath::Min(getLowestNonZeroBin(hs),
                            getLowestNonZeroBin(hb));

    // std::cerr << *it << "\t" << max << std::endl;
    double minFact = 0.1; //0.9* - 0.8*double(c1->GetLogy());
    hs->SetMaximum(max*1.1);
    hs->SetMinimum(min*minFact);
    hb->SetMaximum(max*1.1);
    hb->SetMinimum(min*minFact);
  }
}




void drawFisherPlot(TFile* f){
  TH1D* hSigInt = (TH1D*) f->Get("hSigInt");
  TH1D* hBackInt = (TH1D*) f->Get("hBackInt");

  int fitStartBin = hBackInt->FindLastBinAbove(1e-2);
  int fitEndBin = hBackInt->FindLastBinAbove(1e-7);
  double fitStart = hBackInt->GetXaxis()->GetBinLowEdge(fitStartBin);
  double fitEnd = hBackInt->GetXaxis()->GetBinUpEdge(fitEndBin);
  double xLow = hBackInt->GetXaxis()->GetBinLowEdge(0);
  double xHigh = hBackInt->GetXaxis()->GetBinUpEdge(hBackInt->GetNbinsX());

  std::cerr << fitStart << "\t" << fitEnd << std::endl;
  // TF1* fBackExp = new TF1("fBackExp", "[0]*exp(-[1]*x)", xLow, xHigh);
  TF1* fBackExp = new TF1("fBackExp", "[0]*exp(-[1]*x - [2]*x*x - [3]*x*x*x + [4]*x*x*x*x)", xLow, xHigh);

  auto c1 = new TCanvas();
  hSigInt->Draw();
  hBackInt->Draw("same");
  hBackInt->Fit(fBackExp, "", "", fitStart, fitEnd);

  hSigInt->SetLineColor(kRed);
  hBackInt->SetLineColor(kBlue);

  // TF1* fBackExpSolver = new TF1("fBackExpSolver", TString::Format("pow([0]*exp(-[1]*x) - %lf, 2)", backgroundAcceptance),
  //                               hBackInt->GetXaxis()->GetBinLowEdge(0),
  //                               hBackInt->GetXaxis()->GetBinUpEdge(hBackInt->GetNbinsX()));

  std::cerr << "With " << numRfTriggersA3 << " RF triggers, and only wanting " << numDesiredBackground << " to pass cuts..." << std::endl;
  std::cerr << "I have a background acceptance of " << backgroundAcceptance << std::endl;
  double requiredFisherScore = fBackExp->GetX(backgroundAcceptance, 0.001, 100);
  std::cerr << "Extrapolating my fit to that acceptance, I get a fisher score at " << requiredFisherScore << std::endl;
  fisherCutVal = requiredFisherScore; // set for other plots
  double thermalCutSignalEfficiency = hSigInt->Interpolate(requiredFisherScore);
  std::cerr << "My signal efficiency at this value is " << thermalCutSignalEfficiency << std::endl;

  fBackExp->SetRange(fitStart, requiredFisherScore);
  
  TGraph* grLine = new TGraph();
  const int nInterpPoints = 50;
  for(int i=0; i < nInterpPoints; i++)
  {
    double x = fitEnd + (requiredFisherScore - fitEnd)*(double(i)/(nInterpPoints-1));
    double y = fBackExp->Eval(x);
    grLine->SetPoint(i, x, y);
  }
  grLine->SetPoint(grLine->GetN(), requiredFisherScore, backgroundAcceptance);
  grLine->SetPoint(grLine->GetN(), requiredFisherScore, thermalCutSignalEfficiency);
  grLine->SetLineColor(kMagenta);
  grLine->SetLineStyle(3);
  grLine->Draw("lsame");

  TH1D* hWaisInt = NULL;
  bool doWais = true;
  if(doWais){
    TTree* waisTree = (TTree*) f->Get("waisTree");

    if(waisTree)    
    {
      int nx = hSigInt->GetNbinsX();
      TH1D* hWais = new TH1D("hWais", "hWais", nx, hSigInt->GetXaxis()->GetBinLowEdge(1), hSigInt->GetXaxis()->GetBinUpEdge(nx));
      hWaisInt = new TH1D("hWaisInt", "hWaisInt", nx, hSigInt->GetXaxis()->GetBinLowEdge(1), hSigInt->GetXaxis()->GetBinUpEdge(nx));
      hWaisInt->SetLineColor(kCyan);
      TString drawWais = fr->getFisherFormula() + ">>hWais";
      waisTree->Draw(drawWais, "", "goff");
      hWais->Scale(1./hWais->Integral());
      double cumulativeWais = 1; //hSignal->Integral();
      for(int bx=1; bx <= nx; bx++){
        hWaisInt->SetBinContent(bx,  cumulativeWais);

        cumulativeWais -= hWais->GetBinContent(bx);
      }
      hWaisInt->Draw("same");

      double thermalCutWaisEfficiency = hWaisInt->Interpolate(requiredFisherScore);
      
      std::cerr <<  "By the way, the WAIS pulse efficiency is " << thermalCutWaisEfficiency << std::endl;
    }
  }
  

  hSigInt->SetTitle("Fisher Discriminant;Fisher Score;Fraction of events above value");
  hSigInt->GetXaxis()->SetRangeUser(-15, 25);
  hSigInt->GetYaxis()->SetRangeUser(1e-09, 2);
  auto l1 = new TLegend(0.65, 0.75, 0.99, 0.99);
  l1->AddEntry(hSigInt, "MC", "l");
  if(hWaisInt){
    l1->AddEntry(hWaisInt, "WAIS pulses", "l");
  }
  l1->AddEntry(hBackInt, "Above horizontal RF triggers", "l");
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
  
  auto cSNR = new TCanvas();
  auto lSNR = new TLegend(0.8, 0.8, 1, 1);
  for(unsigned i=0; i < effSnr.size(); i++){
    const char* opt = i == 0 ? "" : "same";
    effSnr.at(i)->Draw(opt);
    effSnr.at(i)->SetLineColor(colors[i%nC]);
    TString n = effSnr.at(i)->GetName();
    n.ReplaceAll("_vs_SNR", "").ReplaceAll("eff_", "");
    lSNR->AddEntry(effSnr.at(i), n, "l");
  }
  lSNR->Draw();
  
  auto cEnergy = new TCanvas();
  auto lEnergy = new TLegend(0.8, 0.8, 1, 1);
  
  for(unsigned i=0; i < effEnergy.size(); i++){
    const char* opt = i == 0 ? "" : "same";
    effEnergy.at(i)->Draw(opt);
    effEnergy.at(i)->SetLineColor(colors[i%nC]);
    TString n = effEnergy.at(i)->GetName();
    n.ReplaceAll("_vs_Energy", "").ReplaceAll("eff_", "");

    lEnergy->AddEntry(effSnr.at(i), n, "l");

    if(i==0){
      effEnergy.at(i)->SetTitle("Current efficiency of pre-clustering cuts (independent); log10(Energy) eV; Efficiency");
    }
  }  

  if(makeFisher && fisherCutVal < 999){

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
    TString cutCommand = "weight*(" + cutForm + TString::Format("> %lf)", fisherCutVal);
    signalTree->Draw(passedCommand, cutCommand, "goff");

    TString nameTotal = name + "_total";
    TH1D* hTotal = new TH1D(nameTotal, nameTotal, nx, xLow, xHigh);

    std::cout << hPassed->Integral() << "\t" << hTotal->Integral() << std::endl;

    TString totalCommand = "TMath::Log10(mc_energy)>>" + nameTotal;

    signalTree->Draw(totalCommand, "weight", "goff");

    TEfficiency* effE = new TEfficiency(*hPassed, *hTotal);
    effE->SetName(name);
    effE->SetTitle("Fisher effiency?");
    effE->SetTitle(TString::Format("Thermal cut efficiency at Fisher Score of %4.2lf; log10(Energy) eV; Efficiency", fisherCutVal));

    // cSNR->cd();
    cEnergy->cd();
    lEnergy->AddEntry(effE, "Thermal cut", "l");
    effE->Draw("same");



    // TString namePassed = name + "_passed";
    // TH1D* hPassed = new TH1D(namePassed, namePassed, nx, xLow, xHigh);

    // TString passedCommand = "TMath::Log10(mc_energy)>>" + namePassed;
    // TString cutCommand = "weight*(" + cutForm + TString::Format("> %lf)", fisherCutVal);
    // signalTree->Draw(passedCommand, cutCommand, "goff");

    // TString nameTotal = name + "_total";
    // TH1D* hTotal = new TH1D(nameTotal, nameTotal, nx, xLow, xHigh);

    // std::cout << hPassed->Integral() << "\t" << hTotal->Integral() << std::endl;

    // TString totalCommand = "TMath::Log10(mc_energy)>>" + nameTotal;

    // signalTree->Draw(totalCommand, "weight", "goff");

    // TEfficiency* effE = new TEfficiency(*hPassed, *hTotal);
    // effE->SetName(name);
    // effE->SetTitle("Fisher effiency?");
    // effE->SetTitle(TString::Format("Thermal cut efficiency at Fisher Score of %4.2lf; SNR; Efficiency", fisherCutVal));
    

    // cSNR->cd();
    // effSNR->Draw("same")
  }
  cEnergy->cd();
  lEnergy->Draw();
  
}


void findReasonsForOutlying(TFile* f, double fisherCut, std::vector<TString>& treeNames){
  Acclaim::CutOptimizer::FisherResult* fr = (Acclaim::CutOptimizer::FisherResult*) f->Get("FisherResult");
  TString fisherForm = fr->getFisherFormula();

  std::vector<TString> fisherComponents; 
  Acclaim::RootTools::tokenize(fisherComponents,  fisherForm.Data(), std::vector<const char*>{"+(",  ")"});


  TString scanCommand = "run:eventNumber";
  for(UInt_t i=0; i < fisherComponents.size(); i++){
    scanCommand += ":" + fisherComponents.at(i);
  }

  TString scanCut  = fisherForm + TString::Format(" > %lf", fisherCut);

  for(auto& treeName : treeNames ){

    TTree* signalTree = (TTree*) f->Get(treeName);
    if(!signalTree){
      std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", couldn't find " << treeName << std::endl;
    }
    else{
      TString histName = "h" + treeName + "Components";    
      TH2D* hSignalComponents = new TH2D(histName, "Components of Fisher Discriminant for " + treeName,
                                         fisherComponents.size(), 0,  fisherComponents.size(),
                                         1024, -20, 20);

    
      // signalTree->Scan(scanCommand, scanCut);
      for(UInt_t i=0; i < fisherComponents.size(); i++){
        TH2D hTemp("hTemp", "hTemp",
                   hSignalComponents->GetNbinsX(), hSignalComponents->GetXaxis()->GetBinLowEdge(1), hSignalComponents->GetXaxis()->GetBinUpEdge(hSignalComponents->GetNbinsX()),
                   hSignalComponents->GetNbinsY(), hSignalComponents->GetYaxis()->GetBinLowEdge(1), hSignalComponents->GetYaxis()->GetBinUpEdge(hSignalComponents->GetNbinsY()));
        signalTree->Draw(fisherComponents[i] + TString::Format(":%d>>%s", i, hTemp.GetName()), "weight", "goff");      
        hSignalComponents->Add(&hTemp);

        std::cerr << "done " << i << " of " << fisherComponents.size() << std::endl;
      }
      for(UInt_t i=0; i < fisherComponents.size(); i++){
        std::vector<TString> binLabel;
        Acclaim::RootTools::tokenize(binLabel, fisherComponents.at(i), "*");
        hSignalComponents->GetXaxis()->SetBinLabel(i+1, binLabel.size() > 1 ? binLabel.at(1) : "Constant");
      }
      auto cs = new TCanvas();
      cs->SetBottomMargin(0.2);
      cs->SetRightMargin(0.2);    
      hSignalComponents->Draw("colz");
      cs->SetLogz(1);
    }    
  }

  
  
}

void drawCutOptimization(const char* fileName){

  TFile* f = TFile::Open(fileName);
  fr = (Acclaim::CutOptimizer::FisherResult*) f->Get("FisherResult");
  std::vector<TString> treeNames = {"signalTree", "backgroundTree"};//, "waisTree", "blastTree"};
  findReasonsForOutlying(f, 4, treeNames);
  drawFisherPlot(f);
  drawEfficiencies(f, true);
  // overlayOneDimDists(f);
  
}
