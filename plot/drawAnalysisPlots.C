#include "AnitaEventSummary.h"
#include "RootTools.h"

#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TCut.h"
#include "TStyle.h"

#include "TH1D.h"
#include <iostream>
#include "ProgressBar.h"
#include "TFitResult.h"
#include "TF1.h"

using namespace Acclaim;

TChain* waisChain = NULL;
TChain* tenPercentChain = NULL;

const TCut isWaisPulser = "flags.pulser == 1";
const TCut isLDBPulser = "flags.pulser == 2";
const TCut notPulser = "flags.pulser == 0";
const TCut closeToWais = "TMath::Abs(Acclaim::RootTools::getDeltaAngleDeg(sum.higherPeak().phi, wais.phi)) < 5";
const TCut closeToSun = "TMath::Abs(sum.dPhiSun()) < 7 && TMath::Abs(sum.dThetaSun()) < 7";

void drawWaisPlots();
void drawTenPercentPlots();

void drawAnalysisPlots(){

  gStyle->SetOptStat("mre");
  gStyle->SetPalette(kRainBow);
  

  waisChain = new TChain("sumTree");
  waisChain->Add("~/ANITA/anita3Analysis/sineSub/2017_07_28/filterWaisAndNonRfSineSub_*");
  // waisChain->Add("~/ANITA/anita3Analysis/sineSub/filterWaisAndNonRfSineSub_352_2017-07-26_12-10-50.root");  
  // waisChain->Add("~/ANITA/anita3Analysis/sineSub/2017_07_23/filterWaisAndNonRfSineSub_*.root");
  std::cout  << waisChain->GetEntries() << " wais chain entries." << std::endl;


  tenPercentChain = new TChain("sumTree");
  // tenPercentChain->Add("~/ANITA/anita3Analysis/sineSub/2017_07_23/filterDecimatedSineSub_*.root");
  // tenPercentChain->Add("~/ANITA/anita3Analysis/sineSub/2017_07_24/filterDecimatedSineSub_*.root");  
  tenPercentChain->Add("~/ANITA/anita3Analysis/sineSub/2017_07_28/filterDecimatedSineSub_*.root");  
  std::cout  << tenPercentChain->GetEntries() << " decimated chain entries." << std::endl;

   drawWaisPlots();
  // drawTenPercentPlots();
}



void drawTenPercentPlots(){

  auto c3 = new TCanvas();
  auto h3 = new TH2D("hRotCrossCorr", "", 512, 0, 1, 512, 0, 1024);
  tenPercentChain->Draw("sum.higherCoherent().peakHilbert:sum.higherPeak().value>>hRotCrossCorr", "run < 200" + isLDBPulser, "colz");
  c3->SetLogz(1);  

  auto c4 = new TCanvas();
  auto h4 = new TH2D("hRotCrossCorrNot", "", 512, 0, 1, 512, 0, 1024);
  tenPercentChain->Draw("sum.higherCoherent().peakHilbert:sum.higherPeak().value>>hRotCrossCorrNot", "run < 200" + !isLDBPulser, "colz");
  c4->SetLogz(1);  
  return;
  
  
  
  auto c1 = new TCanvas();
  ProgressBar p1(1);
  auto h1 = new TH2D("hHigherPeakVsTime", "", 1024, RootTools::a3StartTime, RootTools::a3EndTime, 128, 0, 1);  
  tenPercentChain->Draw("sum.higherPeak().value:realTime>>hHigherPeakVsTime", "", "colz");
  p1++;
  c1->SetLogz(1);
  c1->Update();

  // removing pulsers
  auto c2 = new TCanvas();
  ProgressBar p2(1);
  auto h2 = new TH2D("hHigherPeakVsTimeNotPulsers", "", 1024, RootTools::a3StartTime, RootTools::a3EndTime, 128, 0, 1);    
  tenPercentChain->Draw("sum.higherPeak().value:realTime>>hHigherPeakVsTimeNotPulser", notPulser, "colz");
  p2++;
  c2->SetLogz(1);
  c2->Update();

                                              
}



void drawWaisPlots(){

    
  // Quickly draw all the plots I need to check the analysis is on track

  // Important WAIS plots

  TH1D* hSunTheta = new TH1D("hSunTheta", "#delta#theta_{sun}", 128, -7, 7);
  TH1D* hSunPhi = new TH1D("hSunPhi", "#delta#phi_{sun}", 128, -7, 7);
  hSunTheta->SetLineColor(kBlue);
  hSunPhi->SetLineColor(kRed);
  
  waisChain->Draw("sum.dThetaSun()>>hSunTheta", notPulser && closeToSun, "goff");
  waisChain->Draw("sum.dPhiSun()>>hSunPhi", notPulser && closeToSun, "goff");
  
  

  const int n0 = 2;
  TH1D* hSuns[n0] = {hSunTheta, hSunPhi};
  auto c0 = Acclaim::RootTools::drawHistsWithStatsBoxes(2, hSuns, "", "mre");
  c0->SetLogy(1);

  hSunTheta->Fit("gaus", "Q");
  auto f0 = (TF1*) hSunTheta->FindObject("gaus");
  if(f0){
    f0->SetLineColor(kBlue);
  }
  
  hSunPhi->Fit("gaus", "Q");
  auto l0 = c0->BuildLegend();
  
  hSunTheta->SetTitle("Sun resolution with sine subtraction; #delta#theta or #delta#phi (Degrees); Events per bin");
  // return;
  
  // waisChain->Draw("Acclaim::RootTools::getDeltaAngleDeg(sum.higherPeak().phi, wais.phi)", isWaisPulser+closeToWais);
  auto h1 = new TH1D("hDecoMinusCohFilt", "deconvolved - coherent_filtered", 128, -50, 50);
  auto h2 = new TH1D("hDecoFiltMinusCohFilt", "deconvolved_filtered - coherent_filtered", 128, -50, 50);
  auto h3 = new TH1D("hCohMinusCohFilt", "coherent - coherent_filtered", 128, -50, 50);
  waisChain->Draw("sum.higherDeconvolved().peakHilbert-sum.higherCoherentFiltered().peakHilbert>>hDecoMinusCohFilt", isWaisPulser+closeToWais, "goff");
  waisChain->Draw("sum.higherDeconvolvedFiltered().peakHilbert-sum.higherCoherentFiltered().peakHilbert>>hDecoFiltMinusCohFilt", isWaisPulser+closeToWais, "goff");  
  waisChain->Draw("sum.higherCoherent().peakHilbert-sum.higherCoherentFiltered().peakHilbert>>hCohMinusCohFilt", isWaisPulser+closeToWais, "goff");

  const int n1 = 3;
  TH1D* h1s[n1] = {h3, h1, h2};
  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlue);  
  h3->SetLineColor(kBlack);

  auto c1 = Acclaim::RootTools::drawHistsWithStatsBoxes(n1, h1s, "", "mre");

  c1->BuildLegend();
  h3->SetTitle("Different in peakHilbert values between the different coherent waveforms; #delta peakHilbert (mV); Events per bin");

  return;

  // waisChain->Draw("sum.higherCoherent().snr-sum.higherDeconvolved().snr", isWaisPulser + closeToWais, "colz");

  // waisChain->Scan("run:eventNumber:sum.higherCoherent().snr-sum.higherDeconvolved().snr", "run==352" + isWaisPulser + closeToWais);

  
  TH1D* hWaisTheta = new TH1D("hWaisTheta", "#delta#theta_{wais}", 128, -7, 7);
  TH1D* hWaisPhi = new TH1D("hWaisPhi", "#delta#phi_{wais}", 128, -7, 7);
  hWaisTheta->SetLineColor(kBlue);
  hWaisPhi->SetLineColor(kRed);
  
  waisChain->Draw("sum.dThetaWais()>>hWaisTheta", isWaisPulser && closeToWais, "goff");
  waisChain->Draw("sum.dPhiWais()>>hWaisPhi", isWaisPulser && closeToWais, "goff");

  const int n2 = 2;
  TH1D* hWaiss[n2] = {hWaisTheta, hWaisPhi};
  auto c2 = Acclaim::RootTools::drawHistsWithStatsBoxes(2, hWaiss, "", "mre");

  hWaisTheta->Fit("gaus", "Q");
  auto f = (TF1*) hWaisTheta->FindObject("gaus");
  f->SetLineColor(kBlue);
  hWaisPhi->Fit("gaus", "Q");
  
  auto l2 = c2->BuildLegend();

  hWaisTheta->SetTitle("WAIS resolution with sine subtraction; #delta#theta or #delta#phi (Degrees); Events per bin");
  
  c2->SetLogy(1);
  
  
  
}
