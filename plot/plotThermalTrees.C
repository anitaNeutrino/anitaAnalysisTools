#include "DrawStrings.h"
#include "RootTools.h"

const double tau_thermal = 2.8727096;
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
  gStyle->SetOptStat("mre");
  gStyle->SetCanvasPreferGL(true);
  
  TChain* c = new TChain("thermalTree");
  c->Add(thermalTreeDataGlob);
  TChain* c2 = new TChain("thermalTree");
  c2->Add(thermalTreeMcGlob);
  SummarySet::startProof();
  c->SetProof(1);


  bool handScanHighThermalSideband = false; //true;
  if(handScanHighThermalSideband){
    // TString scan = TString("run:eventNumber:realTime-1418939696:pol:") + ThermalTree::fisherDiscriminant;
    TString scan = TString("run:eventNumber:realTime-1418939696:pol:") + ThermalTree::fisherDiscriminant;
    TCut fisherCut = ThermalTree::fisherCut;
    c->Scan(scan, fisherCut + ThermalTree::thermalSideBand + ThermalTree::passAllQualityCuts);
    return;
  }  

  bool handScanEdge = false;//true; //false; //true;
  if(handScanEdge){
    // TString scan = TString("run:eventNumber:realTime-1418939696:pol:") + ThermalTree::fisherDiscriminant;
    TString scan = TString("run:eventNumber:pol:") + ThermalTree::fisherDiscriminant;
    TCut fisherCut = TCut(TString(ThermalTree::fisherDiscriminant + " >7 && " + ThermalTree::fisherDiscriminant + " > 4 && run > 200").Data());
    c->Scan(scan, fisherCut + ThermalTree::analysisSample + ThermalTree::passAllQualityCuts);
    return;
  }  
  
  bool plottingTheta = false;
  if(plottingTheta){
    RootTools::canvas(2);
    
    c->Draw("peak_theta>>hTheta_above(300, -90, 90)", ThermalTree::isNotTaggedAsPulser + ThermalTree::passAllQualityCuts + TCut("peak_theta > 0"), "goff");
    TH1F* hTheta_above = (TH1F*) gProof->GetOutputList()->FindObject("hTheta_above");

    c->Draw("peak_theta>>hTheta_below(300, -90, 90)", ThermalTree::isNotTaggedAsPulser + ThermalTree::passAllQualityCuts + TCut("peak_theta < 0"), "goff");
    TH1F* hTheta_below = (TH1F*) gProof->GetOutputList()->FindObject("hTheta_below");

    hTheta_above->SetTitle("Distribution of highest map peak #theta;Peak #theta (Degrees); Events per bin");
    hTheta_above->SetLineColor(kBlue);
    hTheta_above->SetFillColorAlpha(kBlue, 0.2);
    hTheta_above->SetFillStyle(1);
    hTheta_above->Draw();
    hTheta_below->Draw("same");
    hTheta_below->SetLineColor(kRed);
    hTheta_below->SetFillColorAlpha(kRed, 0.2);

    std::cout << "theta above = " << hTheta_above->Integral() << std::endl;
    std::cout << "theta below = " << hTheta_below->Integral() << std::endl;    
    std::cout << "theta above = " << int(hTheta_above->Integral()) << std::endl;
    std::cout << "theta below = " << int(hTheta_below->Integral()) << std::endl;    

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
    c->Draw("deconvolved_filtered_fracPowerWindowGradient:coherent_filtered_fracPowerWindowGradient>>hDedispersionPower1(300, -10, 150, 300, -10, 150)", ThermalTree::weight(ThermalTree::thermalSideBand + ThermalTree::passAllQualityCuts + ThermalTree::anita3QuietTime));
    TH2F* hDedispersionPower1 = (TH2F*) gProof->GetOutputList()->FindObject("hDedispersionPower1");
    hDedispersionPower1->Draw("colz");
    hDedispersionPower1->SetTitle("Above horizontal RF triggers from a quieter portion of the flight; Power Window Gradient (ns); Dedispersed Power Window Gradient (ns); Events per bin");

    c->Draw("deconvolved_filtered_fracPowerWindowGradient:coherent_filtered_fracPowerWindowGradient>>hDedispersionPower3(300, -10, 150, 300, -10, 150)", ThermalTree::weight(ThermalTree::isTaggedAsWaisPulser + ThermalTree::closeToWais + ThermalTree::passAllQualityCuts), "goff");
    TH2F* hDedispersionPower3 = (TH2F*) gProof->GetOutputList()->FindObject("hDedispersionPower3");
    // hDedispersionPower1->Draw("colz");
    // hDedispersionPower1->SetTitle("The power of dedispersion; Power Window Gradient (ns); Dedispersed Power Window Gradient (ns); Events per bin;");
  

    c2->Draw("deconvolved_filtered_fracPowerWindowGradient:coherent_filtered_fracPowerWindowGradient>>hDedispersionPower2(300, -10, 150, 300, -10, 150)", ThermalTree::weight(ThermalTree::closeToMC + ThermalTree::passAllQualityCuts), "goff");
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

  
  // const TString fisher = "0.898497+(1.929594*coherent_filtered_fracPowerWindowGradient/deconvolved_filtered_fracPowerWindowGradient)+(-0.195909*deconvolved_filtered_fracPowerWindowGradient)+(5.943355*coherent_filtered_impulsivityMeasure)+(0.826114*deconvolved_filtered_impulsivityMeasure)+(0.021763*coherent_filtered_peakHilbert)+(-0.012670*deconvolved_filtered_peakHilbert)+(-0.394201*peak_value)";

  bool doPlotParts = false;
  if(doPlotParts){
    RootTools::canvas(4);
    gPad->SetRightMargin(0.15);
    gPad->SetLeftMargin(0.5);    
    TH2D* h2 = plotParts(c2, false, ThermalTree::fisherDiscriminant, ThermalTree::weight(ThermalTree::closeToMC + ThermalTree::passAllQualityCuts), "hSignal", 300, -30, 30);
    h2->Draw("colz");
    gStyle->SetPalette(kRainBow);
    h2->SetTitle("Signal sample; Contribution to Fisher Score; weight #times variable;Events per bin");
    h2->GetYaxis()->SetTitleOffset(5);
    // return;  

    RootTools::canvas(4);
    gPad->SetRightMargin(0.15);
    gPad->SetLeftMargin(0.5);  
    TH2D* h = plotParts(c, true, ThermalTree::fisherDiscriminant, ThermalTree::weight(ThermalTree::thermalSideBand + ThermalTree::passAllQualityCuts), "hBackground", 300, -30, 30);
    h->Draw("colz");
    h->SetTitle("Background sample; Contribution to Fisher Score; weight #times variable;Events per bin");
    h->GetYaxis()->SetTitleOffset(5);  
    return;
  
  }


  
  const TString fdBins = "(300, -30, 30)";
  auto c1 = RootTools::canvas(2);
  c->Draw(ThermalTree::fisherDiscriminant + ">>hFu" + fdBins, ThermalTree::weight(ThermalTree::thermalSideBand + ThermalTree::passAllQualityCuts), "e");
  TH1F* hFu = (TH1F*) gProof->GetOutputList()->FindObject("hFu");

  auto hFuI = RootTools::makeIntegralHist(hFu, false, false);
  auto c1_n2 = RootTools::canvas();
  c1_n2->SetLogy(1);
  hFu->SetLineColor(kRed);
  hFuI->SetLineColor(kRed);
  hFuI->Draw();


  bool fitBackground = false;
  if(fitBackground){

     // gStyle->SetOptFit("");
    hFu->SetTitle("Fits to the background;Fisher Discriminant, F; Events per bin");
    
    
    TH1D* hBackgroundExpo = new TH1D("hBackgroundExpo", "", 25, 0, 1);
    TH1D* hBackgroundPower = new TH1D("hBackgroundPower", "", 25, 0, 1);
    TGraph* gr1 = new TGraph();
    TGraph* gr2 = new TGraph();
    double fitStart = 1.01;
    const double fitEnd = 30;
    int iter=0;
    const double chiSquarePerNdfThreshold = 1.5;
    while(fitStart < 3.5){
      TString f1Name = TString::Format("f1_%d", iter);
      auto f1 = new TF1(f1Name,"expo", fitStart, fitEnd);
      auto f1r = hFu->Fit(f1Name, "Q0", "", fitStart, fitEnd);
      double eval1 =f1->Eval(ThermalTree::fisherThreshold);
      
      f1->SetLineColor(kMagenta);
      f1->SetLineStyle(2);
      TString f2Name = TString::Format("f2_%d", iter);
      auto f2 = new TF1(f2Name, "[0]*pow(x, -[1])", fitStart, fitEnd);
      auto f2r = hFu->Fit(f2, "Q0", "", fitStart, fitEnd);
      // hFu->Fit(f2, "Q", "", fitStart, fitEnd);
      double eval2 = f2->Eval(ThermalTree::fisherThreshold);

      // std::cout << eval1/tau_thermal << "\t" << eval2/tau_thermal << std::endl;
      f2->SetLineColor(kBlue);


      if(f2->GetChisquare()/f2->GetNDF() < chiSquarePerNdfThreshold){
      
	c1->cd();
	f1->Draw("lsame");
	f2->Draw("lsame");

	// const double integralCutOff = 30;
	// double below1 = f1->Integral(fitStart, ThermalTree::fisherThreshold);
	// double above1 = f1->Integral(ThermalTree::fisherThreshold, integralCutOff);
	// double below2 = f2->Integral(fitStart, ThermalTree::fisherThreshold);
	// double above2 = f2->Integral(ThermalTree::fisherThreshold, integralCutOff);

	// std::cout << fitStart << "\t" << below1 << "\t" << above1 << "\t" << std::endl;	
	// std::cout << fitStart << "\t" << below2 << "\t" << above2 << std::endl;
	double A2 = f2->GetParameter(0);
	double p2 = f2->GetParameter(1);
	double T = ThermalTree::fisherThreshold;
	// so I reckon the integral is
	double above2_2 = (A2/(p2-1))*pow(T, 1 - p2);

	double A1 = exp(f1->GetParameter(0));
	double b1 = f1->GetParameter(1);
	double above1_2 = (-A1/b1)*exp(b1*T);
      
	// std::cout << fitStart << "\t" << below1 << "\t" << above1 << "\t" << std::endl;	
	// std::cout << fitStart << "\t" << below2 << "\t" << above2 << "\t" << above2_2 << std::endl;
	std::cout << above1_2 << std::endl;
	std::cout << above2_2 << std::endl;
	std::cout << std::endl << std::endl;
      
      
	// std::cout << f1->GetChisquare() << "\t" << f1->GetNDF() << "\t" << f2->GetChisquare() << "\t" << f2->GetNDF() << std::endl;      

	gr1->SetPoint(gr1->GetN(), fitStart, f1->GetChisquare()/f1->GetNDF());
	gr2->SetPoint(gr2->GetN(), fitStart, f2->GetChisquare()/f2->GetNDF());
      
	const double polFactor = 0.5; //1; //0.5;
	hBackgroundExpo->Fill(polFactor*eval1/tau_thermal);
	hBackgroundPower->Fill(polFactor*eval2/tau_thermal);
      }
      
      fitStart += 0.01;
      iter++;
    }

    RootTools::canvas();
    gr1->SetMaximum(chiSquarePerNdfThreshold);
    gr1->SetMinimum(0);
    gr1->SetLineColor(kMagenta);
    gr2->SetLineColor(kBlue);
    gr1->SetTitle("How well are we fitting this tail?; Fit start position (F); #chi^{2}/NDF");
    gr1->Draw("al");
    gr2->Draw("lsame");
    

    hBackgroundPower->SetLineColor(kBlue);
    hBackgroundExpo->SetLineColor(kMagenta);
    hBackgroundExpo->SetTitle("Distribution of background estimates;Thermal background estimate;Fits per bin");
    TH1D* hs[2] = {hBackgroundExpo, hBackgroundPower};
    RootTools::drawHistsWithStatsBoxes(2, hs, "", "mre");
    auto l0 = new TLegend(0.8, 0.8, 1, 1);
    l0->AddEntry(hBackgroundExpo, "Exponential fits", "l");
    l0->AddEntry(hBackgroundPower, "Power law fits", "l");
    l0->Draw();

    
    return;
  }
  
  // return;


  bool plotWaisAndDown = false;
  TH1F* hFw = NULL;
  TH1F* hFd = NULL;
  TH1F* hFp = NULL;
  TH1D* hFwI = NULL;
  TH1D* hFdI = NULL;
  if(plotWaisAndDown){
  
    RootTools::canvas(2);
    c->Draw(ThermalTree::fisherDiscriminant + ">>hFw" + fdBins, ThermalTree::weight(ThermalTree::isTaggedAsWaisPulser + ThermalTree::closeToWais + ThermalTree::passAllQualityCuts));
    hFw = (TH1F*) gProof->GetOutputList()->FindObject("hFw");
    // return;
    hFwI = RootTools::makeIntegralHist(hFw, false);
    hFw->SetLineColor(kMagenta);
    hFwI->SetLineColor(kMagenta);  
    c1_n2->cd();
    hFwI->Draw("same");


    RootTools::canvas(2);
    c->Draw(ThermalTree::fisherDiscriminant + ">>hFd" + fdBins, ThermalTree::weight(ThermalTree::analysisSample + ThermalTree::passAllQualityCuts));
    hFd = (TH1F*) gProof->GetOutputList()->FindObject("hFd");

    hFdI = RootTools::makeIntegralHist(hFd, false);
    c1_n2->cd();
    hFd->SetLineColor(kBlack);
    hFdI->SetLineColor(kBlack);  
    // hFdI->Draw("same");

    RootTools::canvas(2);
    c->Draw(ThermalTree::fisherDiscriminant + ">>hFp" + fdBins, ThermalTree::weight(!ThermalTree::notPayloadBlast));
    hFp = (TH1F*) gProof->GetOutputList()->FindObject("hFp");
    hFp->SetLineColor(kOrange);
    
  }

  

  RootTools::canvas(2);
  c2->Draw(ThermalTree::fisherDiscriminant + ">>hFm" + fdBins, ThermalTree::weight(ThermalTree::closeToMC + ThermalTree::passAllQualityCuts));
  TH1F* hFm = (TH1F*) gROOT->FindObject("hFm");
  
  
  auto hFmI = RootTools::makeIntegralHist(hFm, false);
  c1_n2->cd();
  hFmI->Draw("histsame");
  hFmI->SetLineColor(kBlue);
  hFm->SetLineColor(kBlue);
  
  
  // hFuI->SetTitle("Analysis #bf{B} Fisher discriminant distributions;Fisher discriminant score, F (no units); Fraction of pre-thermal cut sample above F (no units)");  
  auto l = new TLegend(0.8, 0.8, 1, 1);
  l->AddEntry(hFuI, "Non-impulsive sideband", "lf");
  if(hFdI){
    l->AddEntry(hFdI, "Below horizontal RF triggers", "l");
  }
  if(hFwI){
    l->AddEntry(hFwI, "WAIS pulses", "l");
  }
  l->AddEntry(hFmI, "MC neutrinos", "lf");
  l->Draw();
  

  RootTools::canvas(2);
  // hFu->SetTitle("Fisher discriminant performance;Fisher discriminant score, F (no units); Events per bin");
  hFm->Scale(1./hFm->Integral());
  hFu->Scale(1./hFu->Integral());

  hFu->SetTitle(";Fisher Discriminant (no units);Fraction of population per bin");
  hFm->SetTitle(";Fisher Discriminant (no units);Fraction of population per bin");  

  hFu->SetLineColorAlpha(hFu->GetLineColor(), 0.9);
  hFu->SetFillColorAlpha(hFu->GetLineColor(), 0.2);  
  hFu->SetFillStyle(1001);
  hFu->SetLineWidth(2);

  hFm->SetLineColorAlpha(hFm->GetLineColor(), 0.9);  
  hFm->SetFillColorAlpha(hFm->GetLineColor(), 0.2);
  hFm->SetFillStyle(1001);
  hFm->SetLineWidth(2);


  hFm->Draw("hist");
  hFu->Draw("histsame");
  {
    const double histMax = 1;
    hFu->SetMaximum(histMax);
    hFm->SetMaximum(histMax);    
  }

  // TLatex* fuTitle = new TLatex(-30, 2*hFu->GetMaximum(), "Analysis #bf{B} Fisher discriminant distributions");
  // fuTitle->SetTextSize(0.05);  
  // fuTitle->Draw();

  

  TLegend* l2 = new TLegend(0.75, 0.79, 0.99, 0.99);
  l2->AddEntry(hFu, "Non-impulsive sideband", "fl");  
  if(hFd){
    hFd->Draw("histsame");    
    l2->AddEntry(hFd, "Below horizontal RF triggers", "l");
  }
  if(hFp){
    hFp->Draw("histsame");
    l2->AddEntry(hFp, "Payload Blast", "l");      
  }
  if(hFw){
    hFw->Draw("histsame");
    l2->AddEntry(hFw, "WAIS pulses", "l");
  } 
  l2->AddEntry(hFm, "MC neutrinos", "fl");
  l2->Draw();
  
  // return;


  // hFm->Rebin(5);
  
  
  // return;
  
}
