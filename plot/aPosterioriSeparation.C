#include "DrawStrings.h"

using namespace Acclaim;

void aPosterioriSeparation(const char* thermalTreeDataGlob = "../cutOptimization/thermalTrees/data/makeThermalTree_*.root",
			   const char* thermalTreeMcGlob = "../cutOptimization/thermalTrees/mc/makeThermalTree_*.root"){

  gStyle->SetOptStat(0);
  gROOT->ProcessLine("#include \"FFTtools.h\"");
  gStyle->SetOptStat("mre");
  gStyle->SetCanvasPreferGL(true);
  
  TChain* c = new TChain("thermalTree");
  c->Add(thermalTreeDataGlob);
  TChain* c2 = new TChain("thermalTree");
  c2->Add(thermalTreeMcGlob);


  auto c1 = new TCanvas();  
  

  // const TString fishName("F_{2}");
  // const TString naiveAttempt1 = "-3.802951+(0.376006*coherent_filtered_fracPowerWindowGradient/deconvolved_filtered_fracPowerWindowGradient)+(0.094475*deconvolved_filtered_fracPowerWindowGradient)+(8.039617*coherent_filtered_impulsivityMeasure)+(-4.146862*deconvolved_filtered_impulsivityMeasure)+(-0.010768*coherent_filtered_peakHilbert)+(0.005617*deconvolved_filtered_peakHilbert)+(-0.709091*peak_value)";
  // const TCut defaultCut = ThermalTree::passAllQualityCuts + !ThermalTree::isAboveHorizontal + ThermalTree::fisherCut;

  const TString fishName("F_{3}");
  const TString naiveAttempt1 = "-4.625069+(0.031321*coherent_filtered_fracPowerWindowGradient)+(5.786568*coherent_filtered_impulsivityMeasure)+(-0.000420*coherent_filtered_peakHilbert)+(-0.486653*peak_value)";
  const TCut defaultCut = ThermalTree::passAllQualityCuts + !ThermalTree::isAboveHorizontal + ThermalTree::fisherCut;  

  // const TString fishName("F_{3}");
  // // const TString naiveAttempt1 = "-8.297505+(0.350286*coherent_filtered_fracPowerWindowGradient/deconvolved_filtered_fracPowerWindowGradient)+(0.015905*deconvolved_filtered_fracPowerWindowGradient)+(8.323876*coherent_filtered_impulsivityMeasure)+(2.289690*deconvolved_filtered_impulsivityMeasure)+(-0.017208*coherent_filtered_peakHilbert)+(0.009563*deconvolved_filtered_peakHilbert)+(-2.568178*peak_value)";
  // const TString naiveAttempt1 = "-1.899616+(0.422179*coherent_filtered_fracPowerWindowGradient/deconvolved_filtered_fracPowerWindowGradient)+(0.077730*deconvolved_filtered_fracPowerWindowGradient)+(7.322789*coherent_filtered_impulsivityMeasure)+(-6.175697*deconvolved_filtered_impulsivityMeasure)";
  // const TCut defaultCut = ThermalTree::passAllQualityCuts + !ThermalTree::isAboveHorizontal + ThermalTree::fisherCut;


  
  
  const TCut candidateCR = "eventNumber==7613856||eventNumber==9097075||eventNumber==11116669||eventNumber==11989349||eventNumber==15717147||eventNumber==16952229||eventNumber==19459851||eventNumber==22345215||eventNumber==23695286||eventNumber==27142546||eventNumber==32907848||eventNumber==33484995||eventNumber==41529195||eventNumber==48837708||eventNumber==58592863||eventNumber==62273732||eventNumber==65187079||eventNumber==66313844||eventNumber==68298837||eventNumber==70013898||eventNumber==71171108||eventNumber==71766273||eventNumber==73726742||eventNumber==74592579||eventNumber==75277769||eventNumber==80840274";

  //const TCut smallClusterEvent = "eventNumber==19567580||eventNumber==19567583||eventNumber==19567584||eventNumber==19606349||eventNumber==41128241||eventNumber==41475569||eventNumber==66672878||eventNumber==66683426||eventNumber==66702609||eventNumber==66703250||eventNumber==66703284||eventNumber==66703286||eventNumber==66704135||eventNumber==66704530||eventNumber==66967345||eventNumber==66967699||eventNumber==66977163||eventNumber==66993765||eventNumber==66998800||eventNumber==70029023||eventNumber==70124461||eventNumber==84650299||eventNumber==84653556";
    const TCut smallClusterEvent = "eventNumber==19567580||eventNumber==19567583||eventNumber==19567584||eventNumber==19606349||eventNumber==41128241||eventNumber==41475569||eventNumber==54431371||eventNumber==54506115||eventNumber==54601697||eventNumber==54635203||eventNumber==54710583||eventNumber==54716960||eventNumber==54755395||eventNumber==54760838||eventNumber==54762764||eventNumber==54765746||eventNumber==54799065||eventNumber==66672878||eventNumber==66683426||eventNumber==66702609||eventNumber==66703250||eventNumber==66703284||eventNumber==66703286||eventNumber==66704135||eventNumber==66704530||eventNumber==66967345||eventNumber==66967699||eventNumber==66977163||eventNumber==66993765||eventNumber==66998800||eventNumber==70029023||eventNumber==70124461||eventNumber==84650299||eventNumber==84653556";
  
  const TCut candidateNu = "eventNumber==83139414";
  c->Draw(naiveAttempt1, candidateNu, "goff");

  TGraph* grNu = new TGraph();
  grNu->SetPoint(0, c->GetV1()[0], -1e50);
  grNu->SetPoint(1, c->GetV1()[0], 1e50);

  SummarySet::startProof();
  c->SetProof(1);

  c->Draw(naiveAttempt1 + " >> hFisherCR(300, -10, 10)", candidateCR);
  TH1F* hFisherCR = (TH1F*)gProof->GetOutputList()->FindObject("hFisherCR");

  c->Draw(naiveAttempt1 + " >> hFisherSC(300, -10, 10)", smallClusterEvent);
  TH1F* hFisherSC = (TH1F*)gProof->GetOutputList()->FindObject("hFisherSC");
  

  c->Draw(naiveAttempt1 + ">>hFisherDD(300, -10, 10)",  ThermalTree::weight(defaultCut + ThermalTree::isNotTaggedAsPulser));
  c2->Draw(naiveAttempt1 + ">>hFisherMC(300, -10, 10)", ThermalTree::weight(defaultCut + ThermalTree::isNotTaggedAsPulser));  

  TH1F* hFisherDD = (TH1F*)gProof->GetOutputList()->FindObject("hFisherDD");
  TH1F* hFisherMC = (TH1F*)gDirectory->FindObject("hFisherMC");
  
  c->Draw(naiveAttempt1 + ">>hFisherWP(300, -10, 10)", ThermalTree::weight(ThermalTree::isTaggedAsWaisPulser + defaultCut));
  TH1F* hFisherWP = (TH1F*)gProof->GetOutputList()->FindObject("hFisherWP");

  


  bool scale = true; //false; //true; //false;
  hFisherDD->SetTitle(fishName + " for various populations;" + fishName + " Fisher Discriminant;Events per bin");
  
  if(scale){
  hFisherDD->SetTitle(fishName + " for various populations;" + fishName + " Fisher Discriminant;Fraction of population per bin");
    hFisherDD->Scale(1./hFisherDD->Integral());
    hFisherCR->Scale(1./hFisherCR->Integral());
    hFisherSC->Scale(1./hFisherSC->Integral());
    hFisherWP->Scale(1./hFisherWP->Integral());
    hFisherMC->Scale(1./hFisherMC->Integral());
    hFisherDD->SetMaximum(0.15);
    hFisherDD->SetMinimum(0);
  }
  else{
    c1->SetLogy(1);
  }
  auto l = new TLegend(0.79, 0.79, 0.99, 0.99);
  
  l->AddEntry(grNu, "#nu candidate (83139414)", "l");
 
  
  hFisherDD->SetLineColor(kBlue);
  hFisherDD->SetLineWidth(2);
  hFisherDD->SetFillColorAlpha(kBlue, 0.2);
  hFisherDD->Draw();
  
  grNu->SetLineWidth(2);
  grNu->Draw("lsame");
  if(!scale){
    hFisherCR->SetLineColor(kGreen);
    hFisherCR->SetLineWidth(2);
    hFisherCR->SetFillColorAlpha(kGreen, 0.2);
    hFisherCR->Draw("same");
    l->AddEntry(hFisherCR, "CR candidates", "lf");


    hFisherSC->SetLineColor(kOrange);
    hFisherSC->SetLineWidth(2);
    hFisherSC->SetFillColorAlpha(kOrange, 0.2);
    hFisherSC->Draw("same");
    l->AddEntry(hFisherSC , "Small clusters", "lf");
    
  }
  hFisherWP->SetLineColor(kMagenta);
  hFisherWP->SetLineWidth(2);
  hFisherWP->SetFillColorAlpha(kMagenta, 0.2);
  hFisherWP->Draw("same");
  l->AddEntry(hFisherWP, "WAIS pulses", "lf");
  
  hFisherMC->SetLineColor(kRed);
  hFisherMC->SetLineWidth(2);
  hFisherMC->SetFillColorAlpha(kRed, 0.2);
  hFisherMC->Draw("same");
  l->AddEntry(hFisherMC, "Kotera Max MC#nus", "lf");
  

  l->AddEntry(hFisherDD, "Below horizontal data", "lf");
  l->Draw();
  
}
