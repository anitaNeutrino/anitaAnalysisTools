
void printClusterMultiplicityTable(TFile* f){
  TTree* t = (TTree*) f->Get("clusterTree");

  const int numK = 2;
  TCut knownCuts[numK] = {"knownBase", "!knownBase"};
  //  100, 10-100, 5-10, 4,3,2,1
  const int numSizeGroups = 7;
  TCut sizeGroupCuts[numSizeGroups] = {"numDataEvents >= 100",
				       "numDataEvents >= 10 && numDataEvents < 100",
				       "numDataEvents >= 5 && numDataEvents < 10",
				       "numDataEvents == 4",
				       "numDataEvents == 3",
				       "numDataEvents == 2",
				       "numDataEvents == 1"};

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


void drawAllClusterTGraphs(TFile* f){
  TList* l = f->GetListOfKeys();

  std::vector<TGraphAntarctica*> grs;
  TIter next(l);
  int i=0;
  while (TObject* obj = next()){

    TString name = obj->GetName();
    cout << name << endl;
    if(name.Contains("grClusterCenters")){
      grs.push_back((TGraphAntarctica*)f->Get(name));
    }
  }

  for(auto& gr : grs){
    gr->Draw();
  }
}

void drawEvents(TFile* f, int clusterInd){

  TTree* clusterTree = (TTree*)f->Get("clusterTree");
  Acclaim::Clustering::Cluster* cluster = NULL;
  clusterTree->SetBranchAddress("cluster", &cluster);
  clusterTree->GetEntry(clusterInd);

  const char* eventTreeName = (cluster->numDataEvents > 1 || cluster->knownBase > 0) ? "eventTree" : "nonBaseSingletTree";

  TTree* eventTree = (TTree*)f->Get(eventTreeName);
  Acclaim::Clustering::Event* event = NULL;
  eventTree->SetBranchAddress("event", &event);


  TString eventCutString = TString::Format("cluster[0]==%d", clusterInd);
  TCut eventCut = eventCutString.Data();
  TGraphAntarctica* grAnita = new TGraphAntarctica(eventTree, "anita.longitude", "anita.latitude", eventCut);

  TString clusterCutString = TString::Format("Entry$==%d", clusterInd);
  TCut clusterCut = clusterCutString.Data();
  TGraphAntarctica* grClusterPosition = new TGraphAntarctica(clusterTree, "longitude", "latitude", clusterCut);

  std::vector<TArrowAntarctica*> arrows;


  int nEvents = cluster->numDataEvents;
  if(nEvents > Acclaim::Clustering::LogLikelihoodMethod::SmallClusterSizeThreshold){

    TString name = TString::Format("hCluster%d", clusterInd);
    TH2DAntarctica* h = new TH2DAntarctica(name, name);

    Long64_t n = eventTree->GetEntries();
    Acclaim::ProgressBar p(n);
    for(Long64_t entry=0; entry < n; entry++){
      eventTree->GetEntry(entry);

      if(event->cluster[0]==clusterInd){
	h->Fill(event->longitude, event->latitude);
      }

      p.inc(entry, n);
    }

    h->Draw("colz");
    grAnita->SetMarkerSize(0.1);

  }
  else{

    Long64_t n = eventTree->GetEntries();
    Acclaim::ProgressBar p(n);
    for(Long64_t entry=0; entry < n; entry++){
      eventTree->GetEntry(entry);

      if(event->cluster[0]==clusterInd){
	arrows.push_back(event->makeArrowFromAnitaToEvent());
	arrows.back()->SetLineWidth(2);
      }

      p.inc(entry, n);
    }

    TGraphAntarctica* grEvents = new TGraphAntarctica(eventTree, "longitude", "latitude", eventCut);
    grEvents->Draw();

    for(auto& arrow : arrows){
      arrow->Draw();
    }

  }

  grAnita->SetMarkerColor(kGreen);
  grAnita->Draw();

  grClusterPosition->SetMarkerStyle(kFullStar);
  grClusterPosition->SetMarkerColor(kRed);

  grClusterPosition->Draw("psame");
}


void drawResolutionDistributions(TFile* f, int whichCluster = -1){

  TTree* clusterTree = (TTree*)f->Get("clusterTree");
  Acclaim::Clustering::Cluster* cluster = NULL;
  clusterTree->SetBranchAddress("cluster", &cluster);

  TTree* eventTree = (TTree*)f->Get("eventTree");
  Acclaim::Clustering::Event* event = NULL;
  eventTree->SetBranchAddress("event", &event);

  const Long64_t nClusters = clusterTree->GetEntries();
  const Long64_t nEvents = eventTree->GetEntries();

  const double thetaDegRange = 5;
  const double phiDegRange = 10;
  const double maxPlotLL = 300;

  TGraphAntarctica* grPlottedClusters = new TGraphAntarctica();

  int numInVeryLargeKnownClusters = 0;
  for(int clusterInd = 0; clusterInd < nClusters; clusterInd++){

    if(whichCluster >= 0 && clusterInd != whichCluster) continue;

    clusterTree->GetEntry(clusterInd);
    if(cluster->numDataEvents > 1000 && Acclaim::Clustering::isVaguelyNearMcMurdo(*cluster)){

      const int nBins = 128;
      TString baseName = TString::Format("lon %4.2lf, lat%4.2lf: n=%d ", cluster->longitude,  cluster->latitude, cluster->numDataEvents);
      numInVeryLargeKnownClusters += cluster->numDataEvents;

      grPlottedClusters->SetPoint(grPlottedClusters->GetN(), cluster->longitude, cluster->latitude);

      TString inCluster = TString::Format("cluster[0]==%d", clusterInd);

      TString thetaName = TString::Format("hDeltaTheta%d", clusterInd);
      TH1D* hDeltaTheta = new TH1D(thetaName, baseName + "#delta#theta", nBins, -thetaDegRange, thetaDegRange);

      TString phiName = TString::Format("hDeltaPhi%d", clusterInd);
      TH1D* hDeltaPhi = new TH1D(phiName, baseName + "#delta#phi", nBins, -phiDegRange, phiDegRange);

      TString llName = TString::Format("hLL%d", clusterInd);
      TH1D* hLL = new TH1D(llName, baseName + "Log Likelihood", nBins, 0, maxPlotLL);

      eventTree->Draw(TString("event.dThetaCluster>>") + thetaName, inCluster, "goff");
      eventTree->Draw(TString("event.dPhiCluster>>") + phiName, inCluster, "goff");
      eventTree->Draw(TString("event.logLikelihood>>") + llName, inCluster, "goff");

      hDeltaTheta->SetLineColor(kBlue);
      hDeltaPhi->SetLineColor(kRed);

      auto c1 = new TCanvas();
      c1->Divide(2);
      c1->cd(1);
      hDeltaPhi->Draw();
      hDeltaTheta->Draw("same");
      c1->cd(2);
      gPad->SetLogy(1);
      hLL->Draw();
    }
  }
  if(grPlottedClusters->GetN() > 0){
    new TCanvas();
    grPlottedClusters->SetMarkerStyle(8);
    grPlottedClusters->SetMarkerColor(kMagenta);
    grPlottedClusters->Draw("p");
  }
  else{
    delete grPlottedClusters;
  }
  int numTotal = eventTree->GetEntries();

  std::cout << "Total events " << numTotal << std::endl;
  std::cout << "Num in very large clusters " << numInVeryLargeKnownClusters << std::endl;
  std::cout << "Diff = " << numTotal - numInVeryLargeKnownClusters << std::endl;

}




void drawClusters(TFile* f){

  TTree* clusterTree = (TTree*)f->Get("clusterTree");
  Acclaim::Clustering::Cluster* cluster = NULL;
  clusterTree->SetBranchAddress("cluster", &cluster);

  TTree* eventTree = (TTree*)f->Get("eventTree");
  Acclaim::Clustering::Event* event = NULL;
  eventTree->SetBranchAddress("event", &event);

  // TTree* eventTree2 = (TTree*)f->Get("nonBaseSingletTree");
  // Acclaim::Clustering::Event* event2 = NULL;
  // eventTree2->SetBranchAddress("event", &event2);

  const Long64_t nClusters = clusterTree->GetEntries();
  const Long64_t nEvents = eventTree->GetEntries();

  const double thetaDegRange = 5;
  const double phiDegRange = 10;
  const double maxPlotLL = 300;

  TGraphAntarctica* grPlottedClusters = new TGraphAntarctica();

  int numInVeryLargeKnownClusters = 0;
  eventTree->GetEntry(0);
  const int n = event->nThresholds;
  
  std::vector<int> colInd(n, 2);
  
  std::vector<TGraphAntarctica*> grs[n];
  
  for(int clusterInd = 0; clusterInd < nClusters; clusterInd++){

    clusterTree->GetEntry(clusterInd);
  // }
    if(cluster->numDataEvents > 0){

      TString cut = TString::Format("cluster[%d]==%d", cluster->llEventCutInd, clusterInd);
      TGraphAntarctica* gr = new TGraphAntarctica(eventTree, "longitude", "latitude", cut.Data());
      //  if(gr->GetN() == 0){	
      // 	delete gr;
      // 	gr =  new TGraphAntarctica(eventTree2, "longitude", "latitude", cut.Data());
      // }

      gr->SetName(TString::Format("grCluster%d_%d", clusterInd, cluster->llEventCutInd));
      gr->SetMarkerColor(colInd[cluster->llEventCutInd]);
      colInd[cluster->llEventCutInd]++;
      grs[cluster->llEventCutInd].push_back(gr);
      
      // gr->Draw();
    }
  }

  for(int i=0; i < n; i++){
    new TCanvas();
    int numClusteredEvents = 0;
    for(unsigned j=0; j < grs[i].size(); j++){
      grs[i][j]->Draw();
      numClusteredEvents = grs[i][j]->GetN(); 
    }
  }
  
}

void plotDebug(TFile* f){


  TH2D* h = (TH2D*) f->Get("hSingleEventTest");

  if(h){
    TGraphAntarctica* gr1 = (TGraphAntarctica*) f->Get("grTestEvent1");
    TGraphAntarctica* gr2 = (TGraphAntarctica*) f->Get("grTestEvent2");
    TGraphAntarctica* grWais = (TGraphAntarctica*) f->Get("grWaisTrue");
    TGraphAntarctica* grWalk = (TGraphAntarctica*) f->Get("grTestMinimizerWalk");
    TGraphAntarctica* grValue = (TGraphAntarctica*) f->Get("grTestMinimizerValue");


    gr1->SetMarkerStyle(8); gr1->SetMarkerColor(kMagenta);
    gr2->SetMarkerStyle(8); gr2->SetMarkerColor(kCyan);
    grWais->SetMarkerStyle(8); grWais->SetMarkerColor(kRed);

    new TCanvas();
    h->Draw("colz");
    gr1->Draw("psame");
    gr2->Draw("psame");
    grWais->Draw("psame");
    grWalk->Draw("lsame");

    new TCanvas();
    grValue->Draw();
  }
}




void plotClustering(const char* fileName = "doClustering_2017-11-09_12-29-34.root"){

  TFile* f = TFile::Open(fileName);
  if(!f) return;

  // drawAllClusterTGraphs(f);

  // drawEvents(f, 0);
  // new TCanvas();
  // for(int i=1; i < 16; i++){
  //   drawEvents(f, i);
  // }
  plotDebug(f);
  // drawClusters(f);

  // drawResolutionDistributions(f);

  return;
  printClusterMultiplicityTable(f);

}
