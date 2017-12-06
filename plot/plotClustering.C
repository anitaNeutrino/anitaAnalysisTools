const int numSizeGroups = 7;
const TCut sizeGroupCuts[numSizeGroups] = {"numDataEvents >= 100",
					   "numDataEvents >= 10 && numDataEvents < 100",
					   "numDataEvents >= 5 && numDataEvents < 10",
					   "numDataEvents == 4",
					   "numDataEvents == 3",
					   "numDataEvents == 2",
					   "numDataEvents == 1"};

const EColor sizeGroupColors[numSizeGroups]  = {kRed, kOrange, kYellow, kGreen, kCyan, kBlue, kMagenta};
const EMarkerStyle sizeGroupMarkers[numSizeGroups] = {kFullCircle,kFullSquare,kFullTriangleUp,kFullTriangleDown,
						      kOpenCircle,kOpenSquare,kOpenTriangleUp};//,kOpenDiamond};
// kOpenCross,
// kFullStar,
// kOpenStar
const int numK = 2;
const TCut knownCuts[numK] = {"knownBase", "!knownBase"};

void printClusterMultiplicityTable(TFile* f){
  TTree* t = (TTree*) f->Get("clusterTree");

  const int numK = 2;
  TCut knownCuts[numK] = {"knownBase", "!knownBase"};
  //  100, 10-100, 5-10, 4,3,2,1

  std::cout << " | Cluster Multiplicity | " << knownCuts[0] << " | " << knownCuts[1] << " | " << std::endl;
  for(int s=0; s < numSizeGroups; s++){
    TString nice = sizeGroupCuts[s].GetTitle();
    nice.ReplaceAll("umDataEvents", "");
    std::cout << " | " << nice << " | ";
    for(int k=0; k < numK; k++){
      if(!(s==numSizeGroups-1 && k==1)){

	t->Draw(">>elist", knownCuts[k] + sizeGroupCuts[s], "entrylist");//");
	TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
	std::cout << elist->GetN();
      }
      else{
	std::cout << "???";
      }
      std::cout << " | ";
    }
    std::cout << std::endl;
  }

// for example over 100, 10-100, 5-10, 4,3,2,1

  // 149527

}


void drawClusters(TFile* f){
  gStyle->SetPalette(kRainBow);
  TTree* eventTree = (TTree*)f->Get("eventTree");
  TTree* mcEventTree = (TTree*)f->Get("mcEventTree");  
  
  int nEvents = eventTree->GetEntries();
  

  std::vector<TTree*> clusterTrees;
  std::vector<double> llEventCuts;
  Acclaim::Clustering::Cluster* cluster = NULL;
  {
    int llCutInd = 0;
    TTree* t = NULL;
    do {
      auto treeName = TString::Format("clusterTree%d", llCutInd);
      t = (TTree*)f->Get(treeName);
      if(t){
	clusterTrees.push_back(t);
	t->SetBranchAddress("cluster", &cluster);
	t->GetEntry(0);
	llEventCuts.push_back(cluster->llEventCut);
      }
      llCutInd++;
    } while(t!=NULL);
  }


  TGraph* grMcClusterEfficiency = NULL;

  if(mcEventTree){

    grMcClusterEfficiency = new TGraph(clusterTrees.size());
    grMcClusterEfficiency->SetName("grMcClusterEfficiency");
    grMcClusterEfficiency->SetTitle("Weighted MC clustering efficiency; -2 log (L); MC clustering efficiency");

    double sumMcWeights = 0;
    Acclaim::Clustering::McEvent* mcEvent = NULL;
    mcEventTree->SetBranchAddress("mcEvent", &mcEvent);
    for(int entry=0; entry < mcEventTree->GetEntries(); entry++){
      mcEventTree->GetEntry(entry);


      for(int z=0; z < mcEvent->nThresholds; z++){
	if(mcEvent->cluster[z]==-1){
	  grMcClusterEfficiency->GetY()[z] += mcEvent->weight;
	}
      }
      sumMcWeights += mcEvent->weight;
    }

    for(int z=0; z < clusterTrees.size(); z++){
      grMcClusterEfficiency->GetY()[z]/=sumMcWeights;
      grMcClusterEfficiency->GetX()[z] = llEventCuts.at(z);
    }
    auto c0 = new TCanvas();
    grMcClusterEfficiency->Draw();
  }
  // figure out some multiplet sizes
  std::vector<TGraph*> grs(numSizeGroups, NULL);
  std::vector<TGraph*> gr2s(numSizeGroups, NULL);

  for(auto& gr : grs){
    static int i = -1;
    i++;    
    gr = new TGraph();
    TString shortTitle = sizeGroupCuts[i].GetTitle();
    shortTitle.ReplaceAll("numDataEvents", "n");
    gr->SetTitle(shortTitle);
    gr2s[i] = new TGraph();
    gr2s[i]->SetTitle(shortTitle);
  }
  
  double maxN = 0;
  for(int llCutInd=0; llCutInd < clusterTrees.size(); llCutInd++){    
    auto t = clusterTrees.at(llCutInd);
    int numDataEventsInClusterTree = 0;
    for(int entry=0; entry < t->GetEntries(); entry++){
      t->GetEntry(entry);
      numDataEventsInClusterTree += cluster->numDataEvents;
    }
    
    for(int i=0; i < numSizeGroups; i++){
      t->Draw(">>elist", sizeGroupCuts[i], "entrylist");
      TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
      int nClusters = elist->GetN();

      int entry=0;
      UInt_t totalEvents = 0;
      while(entry > -1){
	entry = elist->Next();
	if(entry > -1){
	  t->GetEntry(entry);
	  totalEvents += cluster->numDataEvents;
	}
      }
      t->GetEntry(0); 
      grs[i]->SetPoint(llCutInd, cluster->llEventCut, nClusters);
      gr2s[i]->SetPoint(llCutInd, cluster->llEventCut, totalEvents);
      maxN = nClusters >  maxN ? nClusters : maxN;
    }
  }

  // TMultiGraph* grMG = new TMultiGraph();
  auto c1 = new TCanvas();
  double c1Max = 1.1*maxN;
  auto c2 = new TCanvas();
  

  for(auto& gr : grs){
    static int i = -1;
    i++;

    if(gr->GetN()> 0){
      c1->cd();
      gr->SetFillColor(0);      
      gr->SetMaximum(c1Max);
      gr->SetMinimum(0);
      TString opt = i == 0 ? "alp plc pmc" : "lp same plc pmc";
      // gr->SetLineColor(sizeGroupColors[i]);
      gr->SetMarkerStyle(sizeGroupMarkers[i]);      
      gr->Draw(opt);

    

      gr->GetXaxis()->SetTitle("Log-Likelihood cut value");
      gr->GetYaxis()->SetTitle("Cluster multiplicity");

      c2->cd();
      // gr2s[i]->SetLineColor(i+1);
      gr2s[i]->SetFillColor(0);
      gr2s[i]->SetMarkerStyle(sizeGroupMarkers[i]);
      gr2s[i]->SetMinimum(0.9);
      gr2s[i]->SetMaximum(1e6);
      gr2s[i]->Draw(opt);
      gr2s[i]->GetXaxis()->SetTitle("Log-Likelihood cut value");
      gr2s[i]->GetYaxis()->SetTitle("Number of events");
    }
    // std::cout << i << std::endl;
    // grMG->Add(gr);
  }
  auto l1 = c1->BuildLegend();

  c1->SetLogx(1);
  // c1->SetLogy(1);  
  c2->BuildLegend();
  c2->SetLogx(1);
  c2->SetLogy(1);

  grs[0]->SetTitle("Cluster multiplicity; -2 log (L); Number of clusters");
  c1->cd();
  if(grMcClusterEfficiency){
    TGraph* grMc1 = (TGraph*) grMcClusterEfficiency->Clone("grMcClusterEfficiencyScaled1");
    grMc1->Draw("lsame");
    grMc1->SetLineStyle(2);
    grMc1->SetLineWidth(2);
    grMc1->SetLineColor(kRed);
    for(int i=0; i < grMc1->GetN(); i++){
      grMc1->GetY()[i]*=c1Max;
    }
    c1->SetTicky(0);
    TGaxis *axis = new TGaxis(grs[0]->GetXaxis()->GetXmax(),0,
			      grs[0]->GetXaxis()->GetXmax(),c1Max,
			      0, 1,510,"+L");
    TLegendEntry* l1e = l1->AddEntry(grMc1, "MC clustering efficiency", "l");
    l1e->SetTextColor(kRed);
    axis->SetTitle("MC clustering efficiency");
    axis->SetTextColor(kRed);
    axis->SetLineColor(kRed);
    axis->SetLabelColor(kRed);
    axis->Draw();
  }  
  
  
  
  
  gr2s[0]->SetTitle("Cluster multiplicity; -2 log (L); Number of events in clusters");
  
}











void plotDebug(TFile* f){


  TH2D* h = (TH2D*) f->Get("hSingleEventTest");

  if(h){
    TGraphAntarctica* gr1 = (TGraphAntarctica*) f->Get("grTestEvent1");
    TGraphAntarctica* gr2 = (TGraphAntarctica*) f->Get("grTestEvent2");
    TGraphAntarctica* grWais = (TGraphAntarctica*) f->Get("grWaisTrue");
    TGraphAntarctica* grMinPos = (TGraphAntarctica*) f->Get("grTestMinimumPosition");    
    TGraphAntarctica* grWalk = (TGraphAntarctica*) f->Get("grTestMinimizerWalk");
    TGraphAntarctica* grValue = (TGraphAntarctica*) f->Get("grTestMinimizerValue");


    gr1->SetMarkerStyle(8); gr1->SetMarkerColor(kCyan);
    gr2->SetMarkerStyle(8); gr2->SetMarkerColor(kGreen);
    grMinPos->SetMarkerStyle(8); grMinPos->SetMarkerColor(kRed);
    grWais->SetMarkerStyle(8); grWais->SetMarkerColor(kBlack);

    auto c1 = new TCanvas();
    h->Draw("colz");
    gr1->Draw("psame");
    gr2->Draw("psame");
    grWais->Draw("psame");
    grMinPos->Draw("psame");
    grWalk->Draw("lsame");
    c1->SetLogz(1);

    auto l1 = new TLegend();
    l1->AddEntry(gr1, "Event a (nudged down)", "p");
    l1->AddEntry(gr2, "Event b (nudged down)", "p");
    l1->AddEntry(grMinPos, "Fitted position, x", "p");
    l1->AddEntry(grWais, "True WAIS", "p");
    l1->AddEntry(grWalk, "Minuit's journey", "l");
    l1->Draw();
    
    new TCanvas();
    grValue->Draw();
  }
}




void plotClustering(const char* fileName = ""){

  TFile* f = TFile::Open(fileName);
  if(!f) return;

  // drawAllClusterTGraphs(f);

  // drawEvents(f, 0);
  // new TCanvas();
  // for(int i=1; i < 16; i++){
  //   drawEvents(f, i);
  // }
  // plotDebug(f);
  drawClusters(f);

  // drawResolutionDistributions(f);

  return;
  // printClusterMultiplicityTable(f);

}
