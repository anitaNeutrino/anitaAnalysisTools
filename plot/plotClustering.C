
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

  int nEvents = cluster->numDataEvents;
  

  if(nEvents < )
  
  TString name = TString::Format("hCluster%d", clusterInd);
  TH2DAntarctica* h = new TH2DAntarctica(name, name);
  TGraphAntarctica* grAnita = new TGraphAntarctica();

  TTree* clusterTree = (TTree*)f->Get("clusterTree");
  TString cutString = TString::Format("Entry$==%d", clusterInd);
  TCut cut = cutString.Data();
  TGraphAntarctica* grClusterPosition = new TGraphAntarctica(clusterTree, "longitude", "latitude", cut);

  
  TTree* clusteredDataTree = (TTree*)f->Get("clusteredDataTree");


  
  Acclaim::Clustering::Event* event = NULL;
  clusteredDataTree->SetBranchAddress("event", &event);

  
  Long64_t n = clusteredDataTree->GetEntries();
  Acclaim::ProgressBar p(n);
  for(Long64_t entry=0; entry < n; entry++){
    clusteredDataTree->GetEntry(entry);

    if(event->cluster==clusterInd){
      h->Fill(event->longitude, event->latitude);
      grAnita->SetPoint(grAnita->GetN(), event->anita.longitude, event->anita.latitude);
    }
    
    p.inc(entry, n);
  }

  h->Draw("colz");
  grAnita->SetMarkerColor(kGreen);
  grAnita->SetMarkerSize(0.1);  
  grAnita->Draw("psame");
  
  grClusterPosition->SetMarkerStyle(kFullStar);
  grClusterPosition->SetMarkerColor(kRed);
  
  grClusterPosition->Draw("psame");
  
}


void plotClustering(const char* fileName = "doClustering_2017-10-18_10-35-07.root"){

  TFile* f = TFile::Open(fileName);
  if(!f) return;

  drawAllClusterTGraphs(f);

  new TCanvas();
  drawEvents(f, 0);


  return;
  printClusterMultiplicityTable(f);
  
}
