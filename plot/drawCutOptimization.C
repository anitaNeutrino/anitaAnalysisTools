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

const int numPreFisherEvents = 20760360; //2527162.0; //8298605; //77e6;
// const int numPreFisherEvents = 78e6; //8298605; //77e6;
const double numDesiredBackground = 0.5; //ish
const double backgroundAcceptance = numDesiredBackground/numPreFisherEvents;
double fisherCutVal = 999;

Acclaim::CutOptimizer::FisherResult* fr = NULL;

double backgroundOutlierLimit = 999;
std::vector<UInt_t> outlierEventNumbers;


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





void findBackgroundOutliers(TFile* f, int nOutApprox = 20){
  
  if(backgroundOutlierLimit >= 999){
    TTree* t = (TTree*) f->Get("backgroundTree");
    TH1D* hBack = (TH1D*) f->Get("hbackgroundTreeFisher_py");

    if(t && hBack){
      const int nb = t->GetEntries();
    
      const int nx = hBack->GetNbinsX();
      double reverseIntegral=0;
      for(int bx=nx; bx > 0; bx--){
        reverseIntegral += hBack->GetBinContent(bx);

        if(reverseIntegral*nb < nOutApprox){
          backgroundOutlierLimit = hBack->GetBinLowEdge(bx);
        }
      }

      std::cerr << "Background outlier values are = " << backgroundOutlierLimit << std::endl;    
      TString fc = fr->getFisherFormula();
      std::cout << fc << std::endl;
      t->Draw("run:eventNumber:" + fc, fc + TString::Format(" > %lf", backgroundOutlierLimit), "goff");
      const int n = t->GetSelectedRows();
      double* runs = t->GetV1();
      double* events = t->GetV2();
      double* fs = t->GetV3();
      outlierEventNumbers.reserve(n);
      for(int i=0; i < n; i++){
        outlierEventNumbers.push_back(events[i]);
        std::cout << runs[i] << "\t" << outlierEventNumbers[i] << "\t" << fs[i] << std::endl;
      }
    }
  }
}



TCanvas* drawWithOutliers(TFile* f, TTree* thisTree, const char* twoDimFormula){

  findBackgroundOutliers(f);
  
  TCanvas* c2 = new TCanvas();

  TString fisherForm = fr->getFisherFormula();
  // thisTree->Draw("trainingDeconvolvedFiltered_peakHilbert:trainingPeak_value", "", "colz");
  thisTree->Draw(twoDimFormula, "", "colz");  
      
  // thisTree->Draw("trainingDeconvolvedFiltered_peakHilbert:trainingPeak_value", TString::Format("%s > %lf", fisherForm.Data(), backgroundOutlierLimit), "goff");
  thisTree->Draw(twoDimFormula, TString::Format("%s > %lf", fisherForm.Data(), backgroundOutlierLimit), "goff");  
  TGraph* gr = new TGraph(thisTree->GetSelectedRows(), thisTree->GetV2(), thisTree->GetV1());
  gr->SetMarkerColor(kMagenta);
  gr->SetMarkerStyle(8);
  gr->Draw("psame");
  c2->SetLogz(1);
  return c2;
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

  // int fitStartBin = hBackInt->FindLastBinAbove(1e-2);
  // int fitEndBin = hBackInt->FindLastBinAbove(1e-7);
  // double fitStart = hBackInt->GetXaxis()->GetBinLowEdge(fitStartBin);
  // double fitEnd = hBackInt->GetXaxis()->GetBinUpEdge(fitEndBin);
  // double xLow = hBackInt->GetXaxis()->GetBinLowEdge(0);
  // double xHigh = hBackInt->GetXaxis()->GetBinUpEdge(hBackInt->GetNbinsX());

  // TTree* backgroundTree = (TTree*) f->Get("backgroundTree");
  // const double nPassPreThermalCut = backgroundTree->GetEntries();
  // std::cerr << "There were " << nPassPreThermalCut << " events passing pre-thermal cuts" << std::endl;
  std::cerr << fitStart << "\t" << fitEnd << std::endl;
  // TF1* fBackExp = new TF1("fBackExp", "[0]*exp(-[1]*x)", xLow, xHigh);
  TF1* fBackExp = new TF1("fBackExp", "[0]*exp(-[1]*x - [2]*x*x - [3]*x*x*x + [4]*x*x*x*x)", xLow, xHigh);

  auto c1 = new TCanvas();
  hSigInt->Draw();
  hBackInt->Draw("same");
  hBackInt->Fit(fBackExp, "", "", fitStart, fitEnd);

  hSigInt->SetLineColor(kRed);
  hBackInt->SetLineColor(kBlue);

  std::cerr << "With " << numPreFisherEvents << " RF triggers, and only wanting " << numDesiredBackground << " to pass cuts..." << std::endl;
  std::cerr << "I have a background acceptance of " << backgroundAcceptance << std::endl;
  double requiredFisherScore = 1e9;
  for(int bx=1; bx < hBackInt->GetNbinsX(); bx++){
    double val = hBackInt->GetBinContent(bx);
    // std::cout << val << std::endl;
    if(val < backgroundAcceptance){
      requiredFisherScore = hBackInt->GetXaxis()->GetBinUpEdge(bx);
      std::cout << val << "\t" << hBackInt->GetBinContent(bx-1) << std::endl;
      break;
    }
  }

  // double requiredFisherScore = fBackExp->GetX(backgroundAcceptance, 0.001, 100);  
  std::cerr << "Extrapolating my fit to that acceptance, I get a fisher score at " << requiredFisherScore << std::endl;
  fisherCutVal = requiredFisherScore; // set for other plots
  double thermalCutSignalEfficiency = hSigInt->Interpolate(requiredFisherScore);
  std::cerr << "My signal efficiency at this value is " << thermalCutSignalEfficiency << std::endl;
}


// below horizontal 20760360




void makeEfficiencies(TFile* f){

  TTree* ts = (TTree*) f->Get("signalAnalysisCutTree");
  TTree* tb = (TTree*) f->Get("backgroundAnalysisCutTree");

  const int nMcDef = 2;
  const char* mcDefaultBranches[nMcDef] = {"weight", "closeToMC"};
  
  const int nC = 13;
  const char* cuts[nC] = {"isRfTrigger", "isNotTaggedAsPulser", "npbc0A", "npbc0B", "npbc1",
			  "npbc2", "npbc3", "smallDeltaRough", "goodGPS",
			  "reasonableHilbertPeakTimeShiftAfterDedispersion",
			  "higherHilbertPeakAfterDedispersion",
			  "higherImpulsivityMeasureAfterDedispersion",
			  "lowerFracPowerWindowGradientAfterDedispersion"};
  // ts->Draw("");


  // isRfTrigger && isNotTaggedAsPulser && npbc0A && npbc0B && npbc1 && npbc2 && npbc3 && smallDeltaRough && goodGPS && reasonableHilbertPeakTimeShiftAfterDedispersion && higherHilbertPeakAfterDedispersion && higherImpulsivityMeasureAfterDedispersion && lowerFracPowerWindowGradientAfterDedispersion
  
}



void plotComponents(TFile* f, double fisherCut, std::vector<TString>& treeNames){
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

    TTree* thisTree = (TTree*) f->Get(treeName);
    if(!thisTree){
      std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", couldn't find " << treeName << std::endl;
    }
    else{
      const double fisherPlotRange = 30;
      TString histName = "h" + treeName + "Components";    
      TH2D* hSignalComponents = new TH2D(histName, "Components of Fisher Discriminant for " + treeName,
                                         fisherComponents.size(), 0,  fisherComponents.size(),
                                         1024, -fisherPlotRange, fisherPlotRange);

    
      // thisTree->Scan(scanCommand, scanCut);
      for(UInt_t i=0; i < fisherComponents.size(); i++){
        TH2D hTemp("hTemp", "hTemp",
                   hSignalComponents->GetNbinsX(), hSignalComponents->GetXaxis()->GetBinLowEdge(1), hSignalComponents->GetXaxis()->GetBinUpEdge(hSignalComponents->GetNbinsX()),
                   hSignalComponents->GetNbinsY(), hSignalComponents->GetYaxis()->GetBinLowEdge(1), hSignalComponents->GetYaxis()->GetBinUpEdge(hSignalComponents->GetNbinsY()));
        thisTree->Draw(fisherComponents[i] + TString::Format(":%d>>%s", i, hTemp.GetName()), "weight", "goff");      
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

      
      // drawWithOutliers(f, thisTree, "trainingDeconvolvedFiltered_peakHilbert:trainingPeak_value");
      drawWithOutliers(f, thisTree, "trainingDeconvolved_impulsivityMeasure:trainingDeconvolved_fracPowerWindowGradient");
      // drawWithOutliers(f, thisTree, "trainingDeconvolved_impulsivityMeasure:trainingDeconvolved_fracPowerWindowGradient");
    }
  }
}




void drawCutOptimization(const char* fileName){

  TFile* f = TFile::Open(fileName);
  fr = (Acclaim::CutOptimizer::FisherResult*) f->Get("FisherResult");  
  std::vector<TString> treeNames = {"signalCuts", "backgroundCuts"};//, "waisTree", "blastTree"};
  drawFisherPlot(f);
  // drawEfficiencies(f, true);
  // findBackgroundOutliers(f, 50);
  // plotComponents(f, 4, treeNames);
  
  
  
  overlayOneDimDists(f);
  
}
