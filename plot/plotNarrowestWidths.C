#include "AnitaEventSummary.h"
#include "TGraph.h"
#include "TH2D.h"
#include "AnitaConventions.h"
#include "SummarySet.h"
#include "ProgressBar.h"
#include "TCanvas.h"

void plotNarrowestWidths(){

  // Acclaim::SummarySet ss("~/ANITA/anita3Analysis/sineSubPlusNotches/2017_09_13/filterDecimatedSineSubPlusBrickWall_*.root");
  // Acclaim::SummarySet ssMc("~/ANITA/anita3Analysis/sineSubPlusNotches/2017_09_13/filterAllSineSubPlusBrickWall_*.root");
  // Acclaim::SummarySet ss("~/ANITA/anita3Analysis/sineSub/2017_09_14/filterDecimatedSineSub_*.root");
  // Acclaim::SummarySet ssMc("~/ANITA/anita3Analysis/sineSub/2017_09_14/filterAllSineSub_*.root");
  Acclaim::SummarySet ss("~/ANITA/anita3Analysis/sineSub/2017_09_21/doSineSub_*.root");
  Acclaim::SummarySet ssMc("~/ANITA/anita3Analysis/sineSub/2017_09_23/doSineSub_*.root");
  
  TH2D* hS2 = new TH2D("hS2", "hS2", 6, -5, 55, 1024, 0, 100);
  TH2D* hB2 = new TH2D("hB2", "hB2", 6, -5, 55, 1024, 0, 100);  

  TH2D* hS3 = new TH2D("hS3", "hS3", 1024, 0, 100, 1024, 0, 100);
  TH2D* hB3 = new TH2D("hB3", "hB3", 1024, 0, 100, 1024, 0, 100);


  TH2D* hS4 = new TH2D("hS4", "hS4", 1024, 0, 100, 1024, 0, 1);
  TH2D* hB4 = new TH2D("hB4", "hB4", 1024, 0, 100, 1024, 0, 1);
  
  TH1D* hS = new TH1D("hS", "hS", 1024, 0, 100);
  TH1D* hB = new TH1D("hB", "hB", 1024, 0, 100);
  hS->SetLineColor(kRed);
  hB->SetLineColor(kBlue);
  TH1D* hSc = new TH1D("hSc", "hSc", 1024, 0, 100);
  TH1D* hBc = new TH1D("hBc", "hBc", 1024, 0, 100);
  hSc->SetLineColor(kRed);
  hBc->SetLineColor(kBlue);

  TH1D* hSc2 = new TH1D("hSc2", "hSc2", 1024, 0, 100);
  TH1D* hBc2 = new TH1D("hBc2", "hBc2", 1024, 0, 100);
  hSc2->SetLineColor(kRed);
  hBc2->SetLineColor(kBlue);

  Acclaim::SummarySet* sss[2] = {&ssMc, &ss};
  
  auto c1 = new TCanvas();

  bool firstGraph = true;
  for(int j=0; j < 2; j++){

    const Long64_t n = 100; //ss.N();
    Acclaim::ProgressBar p(n);
  
    for(Long64_t entry=0; entry < n; entry++){
      sss[j]->getEntry(entry);

      AnitaEventSummary* sum = sss[j]->summary();

      if(true){

        bool graphThisEntry = entry < 100;
    
        TGraph* gr = graphThisEntry ? new TGraph() : NULL;
        
        if(j==0 && sum->trainingPeak().closeToMC(3, 3)){
          for(int i=0; i < 5; i++){
            double powFrac = 10 + i*10;
            hS2->Fill(powFrac, sum->trainingDeconvolvedFiltered().fracPowerWindowBegins[i], sum->weight());
            if(graphThisEntry && gr){
              gr->SetPoint(i, powFrac, sum->trainingDeconvolvedFiltered().fracPowerWindowBegins[i]);
              gr->SetLineColor(kRed);
            }                
          }
          hS->Fill(sum->trainingDeconvolvedFiltered().fracPowerWindowGradient(), sum->weight());
          hSc->Fill(sum->trainingDeconvolvedFiltered().fracPowerWindowIntercept(), sum->weight());
          hSc2->Fill(sum->trainingDeconvolvedFiltered().fracPowerWindowChisquare(), sum->weight());
          hS3->Fill(sum->trainingDeconvolvedFiltered().fracPowerWindowGradient(), sum->trainingDeconvolvedFiltered().fracPowerWindowChisquare(), sum->weight());
          hS4->Fill(sum->trainingDeconvolvedFiltered().fracPowerWindowGradient(), sum->trainingDeconvolvedFiltered().impulsivityMeasure, sum->weight());            

        }
        else if(j == 1 && sum->trainingPeak().theta > 0 && sum->flags.isGood == 1 && sum->flags.pulser == 0){
          for(int i=0; i < 5; i++){
            double powFrac = 10 + i*10;
            hB2->Fill(powFrac, sum->trainingDeconvolvedFiltered().fracPowerWindowBegins[i], sum->weight());
            if(graphThisEntry && gr){
              gr->SetPoint(i, powFrac, sum->trainingDeconvolvedFiltered().fracPowerWindowBegins[i]);
              gr->SetLineColor(kBlue);
            }                
          }
          hB->Fill(sum->trainingDeconvolvedFiltered().fracPowerWindowGradient(), sum->weight());
          hBc->Fill(sum->trainingDeconvolvedFiltered().fracPowerWindowIntercept(), sum->weight());
          hBc2->Fill(sum->trainingDeconvolvedFiltered().fracPowerWindowChisquare(), sum->weight());
          hB3->Fill(sum->trainingDeconvolvedFiltered().fracPowerWindowGradient(), sum->trainingDeconvolvedFiltered().fracPowerWindowChisquare(), sum->weight());
          hB4->Fill(sum->trainingDeconvolvedFiltered().fracPowerWindowGradient(), sum->trainingDeconvolvedFiltered().impulsivityMeasure, sum->weight());            

          // if(sum->trainingDeconvolvedFiltered().fracPowerWindowGradient() < 20){
          //   std::cerr << std::endl << sum->run << "\t" << sum->eventNumber << "\t" << polInd << "\t" << peakInd << "\t" << sum->peak[polInd][peakInd].dPhiMC() << "\t" << sum->peak[polInd][peakInd].dThetaMC() << std::endl;
          // }
        }
    
        // std::cerr << graphThisEntry << "\t" << gr << std::endl;
        if(graphThisEntry && gr){
          if(gr->GetN() > 0){
            c1->cd();
            TString opt = firstGraph ? "a " : "same ";
            opt += "l pmc lmc";
            if(firstGraph){
              gr->SetTitle("Narrowest window containing varying fractions of the total power; Percentage of total power; Narrowest window containing power (ns)");
            }
            gr->Draw(opt);
            gr->SetMinimum(0);
            gr->SetMaximum(100);
            firstGraph = false;
          }
        }
      }
      p.inc(entry, n);
    }
  }
  // auto c2 = new TCanvas();  
  // hS2->Draw("colz");
  // c2->SetLogz(1);
  // auto c3 = new TCanvas();  
  // hB2->Draw("colz");
  // c3->SetLogz(1);

  // auto c4 = new TCanvas();
  // hB->Draw();
  // hS->Draw("same");

  // auto c5 = new TCanvas();
  // hBc->Draw();
  // hSc->Draw("same");

  // auto c6 = new TCanvas();
  // hBc2->Draw();
  // hSc2->Draw("same");

  // auto c7 = new TCanvas();
  // hS3->Draw("colz");

  // auto c8 = new TCanvas();
  // hB3->Draw("colz");

  // auto c9 = new TCanvas();
  // hS4->SetLineColor(kRed);
  // hS4->Draw("box");
  // hS4->SetTitle("; gradient; impulsivity");
  // hB4->SetLineColor(kBlue);
  // hB4->Draw("box same");  
  
  
  
}
