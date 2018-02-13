#include "DrawStrings.h"
#include "RootTools.h"


using namespace Acclaim;

TH2D* plotParts(TChain* c, bool proof, const TString& drawString, const TCut& cuts, const TString& name, int nx, double xlow, double xup){
  std::vector<TString> parts;
  RootTools::tokenize(parts,  drawString, "+");

  const int nVars= parts.size();
  TH2D* h2= new TH2D(name, name, nx, xlow, xup,  nVars, 0, nVars);
  int by = 1;
  for(const auto& part : parts){
    std::cout << "Histogramming " << part.Data() << std::endl;
    TString tempHistName = name + TString::Format("_tempPart%d", by);
    TString histBit = TString(">>") + tempHistName + TString::Format("(%d,%lf,%lf)", nx, xlow, xup);
    TString drawCommand = part + histBit;
    c->Draw(drawCommand, cuts, "goff");

    TH1F* h = NULL;
    if(proof){
      h = (TH1F*) gProof->GetOutputList()->FindObject(tempHistName);
      gProof->GetOutputList()->Clear();
    }
    else{
      h = (TH1F*) gROOT->FindObject(tempHistName);
    }

    for(int bx = 1; bx <= nx; bx++){
      h2->SetBinContent(bx, by, h->GetBinContent(bx));
    }
    delete h;
    h2->GetYaxis()->SetBinLabel(by, part);
    by++;
  }
  return h2;
}


void plotThermalTrees(const char* thermalTreeDataGlob = "data/makeThermalTree_*.root", const char* thermalTreeMcGlob = "mc/makeThermalTree_*.root"){

  gROOT->ProcessLine("#include \"FFTtools.h\"");
  
  TChain* c = new TChain("thermalTree");
  c->Add(thermalTreeDataGlob);
  TChain* c2 = new TChain("thermalTree");
  c2->Add(thermalTreeMcGlob);
  SummarySet::startProof();
  c->SetProof(1);



  bool plottingTheta = true;
  if(plottingTheta){
    RootTools::canvas(2);
    
    c->Draw("peak_theta>>hTheta_above(1024, -90, 90)", Acclaim::ThermalTree::isNotTaggedAsPulser + Acclaim::ThermalTree::passAllQualityCuts + TCut("peak_theta > 0"), "goff");
    TH1F* hTheta_above = (TH1F*) gProof->GetOutputList()->FindObject("hTheta_above");

    c->Draw("peak_theta>>hTheta_below(1024, -90, 90)", Acclaim::ThermalTree::isNotTaggedAsPulser + Acclaim::ThermalTree::passAllQualityCuts + TCut("peak_theta < 0"), "goff");
    TH1F* hTheta_below = (TH1F*) gProof->GetOutputList()->FindObject("hTheta_below");

    hTheta_above->SetTitle("Distribution of highest map peak #theta;Peak #theta (Degrees); Events per bin");
    hTheta_above->SetLineColor(kBlue);
    hTheta_above->SetFillColor(kBlue);
    hTheta_above->SetFillStyle(3345);
    hTheta_above->Draw();
    hTheta_below->Draw("same");
    hTheta_below->SetLineColor(kRed);
    hTheta_below->SetFillColorAlpha(kRed, 0.5);

    auto l = new TLegend(0.8, 0.8, 1, 1);
    double above = hTheta_above->Integral();
    double below = hTheta_below->Integral();
    // std::cout << above << "\t" << below << std::endl;
    l->AddEntry(hTheta_above, TString::Format("Thermal side band - %d events", int(above)), "fl");
    l->AddEntry(hTheta_below, TString::Format("Search region - %d events", int(below)), "fl");
    l->Draw();
    return;
  }




  
  bool plottingDedispersion = false;
  
  if(plottingDedispersion){
  
    RootTools::canvas(4);
    c->Draw("deconvolved_filtered_fracPowerWindowGradient:coherent_filtered_fracPowerWindowGradient>>hDedispersionPower1(1024, -10, 150, 1024, -10, 150)", ThermalTree::weight(Acclaim::ThermalTree::thermalSideBand + Acclaim::ThermalTree::passAllQualityCuts + Acclaim::ThermalTree::anita3QuietTime));
    TH2F* hDedispersionPower1 = (TH2F*) gProof->GetOutputList()->FindObject("hDedispersionPower1");
    hDedispersionPower1->Draw("colz");
    hDedispersionPower1->SetTitle("Above horizontal RF triggers from a quieter portion of the flight; Power Window Gradient (ns); Dedispersed Power Window Gradient (ns); Events per bin");

    c->Draw("deconvolved_filtered_fracPowerWindowGradient:coherent_filtered_fracPowerWindowGradient>>hDedispersionPower3(1024, -10, 150, 1024, -10, 150)", ThermalTree::weight(Acclaim::ThermalTree::isTaggedAsWaisPulser + Acclaim::ThermalTree::closeToWais + Acclaim::ThermalTree::passAllQualityCuts), "goff");
    TH2F* hDedispersionPower3 = (TH2F*) gProof->GetOutputList()->FindObject("hDedispersionPower3");
    // hDedispersionPower1->Draw("colz");
    // hDedispersionPower1->SetTitle("The power of dedispersion; Power Window Gradient (ns); Dedispersed Power Window Gradient (ns); Events per bin;");
  

    c2->Draw("deconvolved_filtered_fracPowerWindowGradient:coherent_filtered_fracPowerWindowGradient>>hDedispersionPower2(1024, -10, 150, 1024, -10, 150)", ThermalTree::weight(Acclaim::ThermalTree::closeToMC + Acclaim::ThermalTree::passAllQualityCuts), "goff");
    TH2F* hDedispersionPower2 = (TH2F*) gROOT->FindObject("hDedispersionPower2");
    hDedispersionPower2->SetLineColor(kRed);
    hDedispersionPower2->SetLineWidth(2);
    // hDedispersionPower2->SetFillColor(kRed);
    hDedispersionPower2->Scale(1e5);
    hDedispersionPower2->Draw("boxsame");
    hDedispersionPower3->SetLineColor(kBlack);
    hDedispersionPower3->SetLineWidth(2);
    hDedispersionPower3->Scale(1e5);
    hDedispersionPower3->Draw("boxsame");
    // hDedispersionPower3->SetFillColor(kBlack);

    // hDedispersionPower1->SetLineColor(kBlue);
    // hDedispersionPower1->SetLineWidth(2);
    // hDedispersionPower1->SetFillColor(kBlue);
    auto l0 = new TLegend(0.8, 0.8, 1, 1);
    // l0->AddEntry(hDedispersionPower1, "Quiet period, above horizontal RF triggers", "fl");
    l0->AddEntry(hDedispersionPower2, "MC neutrinos", "fl");
    l0->AddEntry(hDedispersionPower3, "WAIS pulses", "fl");  
    l0->Draw();
    return;
  }

  
  const TString fisher = ThermalTree::fisherDiscriminant;
  // const TString fisher = "0.898497+(1.929594*coherent_filtered_fracPowerWindowGradient/deconvolved_filtered_fracPowerWindowGradient)+(-0.195909*deconvolved_filtered_fracPowerWindowGradient)+(5.943355*coherent_filtered_impulsivityMeasure)+(0.826114*deconvolved_filtered_impulsivityMeasure)+(0.021763*coherent_filtered_peakHilbert)+(-0.012670*deconvolved_filtered_peakHilbert)+(-0.394201*peak_value)";
  Acclaim::RootTools::canvas(4);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.5);    
  TH2D* h2 = plotParts(c2, false, ThermalTree::fisherDiscriminant, ThermalTree::weight(Acclaim::ThermalTree::closeToMC + Acclaim::ThermalTree::passAllQualityCuts), "hSignal", 1024, -30, 30);
  h2->Draw("colz");
  gStyle->SetPalette(kRainBow);
  h2->SetTitle("Signal sample; Contribution to Fisher Score; weight #times variable;Events per bin");
  h2->GetYaxis()->SetTitleOffset(5);
  // return;  

  Acclaim::RootTools::canvas(4);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.5);  
  TH2D* h = plotParts(c, true, ThermalTree::fisherDiscriminant, ThermalTree::weight(Acclaim::ThermalTree::thermalSideBand + Acclaim::ThermalTree::passAllQualityCuts), "hBackground", 1024, -30, 30);
  h->Draw("colz");
  h->SetTitle("Background sample; Contribution to Fisher Score; weight #times variable;Events per bin");
  h->GetYaxis()->SetTitleOffset(5);  
  return;
  
  
  RootTools::canvas()->SetLogy(1);
  c->Draw(fisher + ">>hFu(1024, -30, 30)", ThermalTree::weight(Acclaim::ThermalTree::thermalSideBand + Acclaim::ThermalTree::passAllQualityCuts));
  TH1F* hFu = (TH1F*) gProof->GetOutputList()->FindObject("hFu");

  auto hFuI = RootTools::makeIntegralHist(hFu, false);
  auto c1_n2 = RootTools::canvas();
  c1_n2->SetLogy(1);
  hFu->SetLineColor(kBlue);
  hFuI->SetLineColor(kBlue);
  hFuI->Draw();




  RootTools::canvas()->SetLogy(1);
  c->Draw(fisher + ">>hFd(1024, -30, 30)", ThermalTree::weight(Acclaim::ThermalTree::analysisSample + Acclaim::ThermalTree::passAllQualityCuts));
  TH1F* hFd = (TH1F*) gProof->GetOutputList()->FindObject("hFd");

  auto hFdI = RootTools::makeIntegralHist(hFd, false);
  c1_n2->cd();
  hFd->SetLineColor(kBlack);
  hFdI->SetLineColor(kBlack);  
  // hFdI->Draw("same");

  
  RootTools::canvas()->SetLogy(1);
  c->Draw(fisher + ">>hFw(1024, -30, 30)", ThermalTree::weight(Acclaim::ThermalTree::isTaggedAsWaisPulser + Acclaim::ThermalTree::closeToWais + Acclaim::ThermalTree::passAllQualityCuts));
  TH1F* hFw = (TH1F*) gProof->GetOutputList()->FindObject("hFw");
  // return;
  auto hFwI = RootTools::makeIntegralHist(hFw, false);
  hFw->SetLineColor(kMagenta);
  hFwI->SetLineColor(kMagenta);  
  c1_n2->cd();
  hFwI->Draw("same");


  
  

  RootTools::canvas()->SetLogy(1);
  c2->Draw(fisher + ">>hFm(1024, -30, 30)", ThermalTree::weight(Acclaim::ThermalTree::closeToMC + Acclaim::ThermalTree::passAllQualityCuts));
  TH1F* hFm = (TH1F*) gROOT->FindObject("hFm");

  auto hFmI = RootTools::makeIntegralHist(hFm, false);
  c1_n2->cd();
  hFmI->Draw("histsame");
  hFmI->SetLineColor(kRed);
  hFm->SetLineColor(kRed);  

  hFuI->SetTitle("Fisher discriminant performance;Fisher discriminant score, F (no units); Fraction of pre-thermal cut sample above F (no units)");
  auto l = new TLegend(0.8, 0.8, 1, 1);
  l->AddEntry(hFuI, "Above horizontal RF triggers", "l");
  // l->AddEntry(hFdI, "Below horizontal RF triggers", "l");
  l->AddEntry(hFwI, "WAIS pulses", "l");
  l->AddEntry(hFmI, "MC neutrinos", "l");
  l->Draw();
  
  RootTools::canvas(2);

  c->Draw(fisher + ">>hFp(1024, -30, 30)", ThermalTree::weight(!Acclaim::ThermalTree::notPayloadBlast));
  TH1F* hFp = (TH1F*) gProof->GetOutputList()->FindObject("hFp");
  hFp->SetLineColor(kOrange);

  hFu->SetTitle("Fisher discriminant performance;Fisher discriminant score, F (no units); Events per bin");  
  hFu->Draw();
  hFd->Draw("histsame");
  hFm->Draw("histsame");
  hFw->Draw("histsame");
  hFp->Draw("histsame");
  
  TLegend* l2 = new TLegend(0.8, 0.6, 1, 1);
  l2->AddEntry(hFu, "Above horizontal RF triggers", "l");
  l2->AddEntry(hFd, "Below horizontal RF triggers", "l");
  l2->AddEntry(hFp, "Payload Blast", "l");  
  l2->AddEntry(hFw, "WAIS pulses", "l");
  l2->AddEntry(hFm, "MC neutrinos", "l");
  l2->Draw();
  
  // return;
  
}
