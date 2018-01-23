#include "RootTools.h"
using namespace Acclaim;

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
const int numK = 3;
const TCut knownCuts[numK] = {"knownBase", "!knownBase", "knownBase || !knownBase"};

void printClusterMultiplicityTable(TFile* f){
  TTree* t = (TTree*) f->Get("clusterTree");


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


void drawClusters2(TFile* f){

  gStyle->SetPalette(kRainBow);
  TTree* eventTree = (TTree*)f->Get("eventTree");
  TTree* mcEventTree = (TTree*)f->Get("mcEventTree");

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

  const int numKnownBases = 437;
  const int nz = llEventCuts.size();

  const int nL = 1;
  const int largeSizes[nL] = {1};
  std::vector<TGraph*> grLarges(nL, NULL);
  for(int i=0; i < nL; i++){
    grLarges[i] = new TGraph();
    grLarges[i]->SetName(TString::Format("grLarge%d", i));
    grLarges[i]->SetTitle(TString::Format("n > %d", largeSizes[i]));
  }
  TH2D* hKnown = new TH2D("hKnown", "Clusters associated with known bases; Event cut size; Number of events", nz, 0, nz, 10, 1, 11);
  TH2D* hUnknown = new TH2D("hUnknown", "Clusters not associated with known bases; Event cut size; Number of events", nz, 0, nz, 10, 1, 11);

  ProgressBar p(nz);
  for(int z=0; z <nz; z++){
    TTree* clusterTree = clusterTrees[z];
    const int nc = clusterTree->GetEntries();
    
    TString histName = TString::Format("hc_%d", z);
    TH1D* hc = new TH1D(histName, histName, nc, 0, nc);

    TString draw = TString::Format("cluster[%d]>>" + histName, z);
    // TString selfLLCut = TString::Format("selfLogLikelihood < %lf", llEventCuts[0]);
    // TString selfLLCut = TString::Format("selfLogLikelihood < 1e-3");
    // TString selfLLCut = TString::Format("thetaAdjustmentRequired==0 && latitude < -75 && latitude > -80 && longitude > 150 && longitude < 180");
    TString selfLLCut = TString::Format("thetaAdjustmentRequired==0");        
    // TString selfLLCut = "";
    eventTree->Draw(draw, selfLLCut, "goff");

    int nLarges[nL] = {0};
    
    for(int c=0; c < nc; c++){
      int nEvents = hc->GetBinContent(c+1);
      bool knownBase = c < numKnownBases;
      if(knownBase){
	hKnown->Fill(z, nEvents);
      }
      else{
	hUnknown->Fill(z, nEvents);
      }
      for(int i=0; i < nL; i++){
	if(nEvents == largeSizes[i]){
	  nLarges[i] += nEvents;
	}
      }
    }

    for(int i=0; i < nL; i++){    
      grLarges[i]->SetPoint(z, llEventCuts.at(z), nLarges[i]);
    }
    
    delete hc;
    p.inc(z);
  }


  RootTools::canvas(3);
  for(int i=0; i < nL; i++){
    
    // double z0 = grLarges[i]->GetY()[0];
    // for(int z=0; z < grLarges[i]->GetN(); z++){
    //   grLarges[i]->GetY()[z] -= z0;
    // }

    TGraph* grGrad = Acclaim::RootTools::makeDerivativeTGraph(grLarges[i]);
    for(int z=0; z < grGrad->GetN(); z++){
      grGrad->GetY()[z] *= -1;
    }
    
    const char* opt = i == 0 ? "al" : "lsame";
    // grLarges[i]->Draw(opt);
    grGrad->Draw(opt);    
  }
  
  RootTools::canvas(4);
  hKnown->Draw("colz");
  RootTools::canvas(4);
  hUnknown->Draw("colz");
}


void drawClusters(TFile* f){
  gStyle->SetPalette(kRainBow);
  TTree* eventTree = (TTree*)f->Get("eventTree");
  TTree* mcEventTree = (TTree*)f->Get("mcEventTree");
  
  
  int nEvents = eventTree->GetEntries();
  Acclaim::Clustering::Event* event = NULL;
  if(eventTree){
    eventTree->SetBranchAddress("event", &event);
  }

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


  bool drawGraph = true;//false;
  if(drawGraph){
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

    bool printed = false;
  
    double maxN = 0;
    for(int llCutInd=0; llCutInd < clusterTrees.size(); llCutInd++){
    
      auto t = clusterTrees.at(llCutInd);
      int numDataEventsInClusterTree = 0;
      for(int entry=0; entry < t->GetEntries(); entry++){
	t->GetEntry(entry);
	numDataEventsInClusterTree += cluster->numDataEvents;
      }
    
      for(int i=0; i < numSizeGroups; i++){
	const int k = 1;
	t->Draw(">>elist", sizeGroupCuts[i] + knownCuts[k], "entrylist");

	if(!printed){
	  std::cout << knownCuts[k].GetTitle() << std::endl;
	  printed = true;
	}
	TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
	int nClusters = 0;//elist->GetN();

	int entry=0;
	UInt_t totalEvents = 0;
	while(entry > -1){
	  entry = elist->Next();

	  if(entry > -1){

	    TCut inCluster("inCluster", TString::Format("cluster[%d]==%d", llCutInd, entry));
	    eventTree->Draw(">>elist2", inCluster, "entrylist");	  
	    // TCut selfLLCut("selfLLCut", TString::Format("selfLogLikelihood < %lf", llEventCuts.at(llCutInd)));
	    // eventTree->Draw(">>elist2", inCluster + selfLLCut, "entrylist");
	    TEntryList *elist2 = (TEntryList*)gDirectory->Get("elist2");
	    int matchedEvents = elist2->GetN();
	    totalEvents += matchedEvents;
	  
	    if(matchedEvents > 0){
	      nClusters++;
	    }
	    std::cout << llCutInd << "\t" << llEventCuts.at(llCutInd) << "\t" << sizeGroupCuts[i].GetTitle() << "\t" << knownCuts[k].GetTitle() << "\t" << entry << "\t" << totalEvents << std::endl;
	  }

	  // if(entry > -1){
	  //   t->GetEntry(entry);
	  //   totalEvents += cluster->numDataEvents;
	  // }
	  
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

    // for(int p=0; p < grs.size(); p++){
    //   for(int q=0; q < grs[p]->GetN(); q++){
    //     std::cout << p << "\t" << q << "\t" << grs[p]->GetY()[q] << std::endl;
    //   }
    // }
  
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
	gr2s[i]->SetMinimum(0.1);
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
  // return;

// const TCut sizeGroupCuts[numSizeGroups] = {"numDataEvents >= 100",
// 					   "numDataEvents >= 10 && numDataEvents < 100",
// 					   "numDataEvents >= 5 && numDataEvents < 10",
// 					   "numDataEvents == 4",
// 					   "numDataEvents == 3",
// 					   "numDataEvents == 2",
// 					   "numDataEvents == 1"};


  // return;
  // for(UInt_t llInd=0; llInd < clusterTrees.size(); llInd++){
  const int sizes[numSizeGroups] = {100, 10, 5, 4, 3, 2, 1};
  bool onlySinglets = true;//false;
  int nb=0;
  
  // for(UInt_t z=21; z < 22; z++){
  // for(UInt_t z=0; z < clusterTrees.size(); z++){
  for(UInt_t z=0; z < 1; z++){
  // for(UInt_t z=1; z < 2; z++){    
    std::vector<TGraphAntarctica*> grCs(numSizeGroups,  NULL);
    std::vector<TList> arrows(numSizeGroups);
    auto l1 = new TLegend(0.8, 0.8,  1, 1);
    for(int i=0; i < numSizeGroups; i++){
      if(onlySinglets && sizes[i]!=1) continue;
      grCs[i] = new TGraphAntarctica();
      TString shortTitle = sizeGroupCuts[i].GetTitle();
      shortTitle.ReplaceAll("numDataEvents", "n");
      grCs[i]->SetTitle(shortTitle);
      // grCs[i]->SetMarkerStyle(sizeGroupMarkers[i]);

      l1->AddEntry(grCs[i], shortTitle, "p pmc");
      grCs[i]->SetName(TString::Format("grC_%d_%d", z, sizes[i]));

    }


    TH1D* hTheta = onlySinglets ? new TH1D(TString::Format("hTheta%d", z), "", 1024, -90, 90) : NULL;

    bool showAll = true; //false;
    for(UInt_t entry=0; entry < eventTree->GetEntries(); entry++){
      eventTree->GetEntry(entry);

      // if(!(event->latitude > -79 && event->latitude < -77 &&
      // 	   event->longitude > 133 && event->longitude < 138)) continue;

      Int_t clusterInd = event->cluster[z];
      clusterTrees.at(z)->GetEntry(clusterInd);
      int numEvents = cluster->numDataEvents;

      if(showAll || !(numEvents==1 && cluster->knownBase==0)){
	for(int i=0; i < numSizeGroups; i++){
	  if(onlySinglets && sizes[i]!=1) continue;
	  if((i==0 && numEvents >= sizes[i]) || (i > 0 && numEvents >= sizes[i] && numEvents < sizes[i-1])){

	    if(event->selfLogLikelihood < llEventCuts.at(z)){
	      grCs.at(i)->SetPoint(grCs.at(i)->GetN(), event->longitude, event->latitude);
	      arrows.at(i).Add(event->makeArrowFromAnitaToEvent());
	      if(hTheta){
		// std::cout << event->selfLogLikelihood << std::endl;
		std::cout << event->run << "\t" << event->eventNumber << std::endl;
		hTheta->Fill(event->theta);
	      }
	    }
	    break;
	  }
	}
      }
    }

    int nd = 0;
    auto c = new TCanvas(); //Acclaim::RootTools::squareCanvas();
    c->SetLogz(1);
    auto h = (TH2DAntarctica*)f->Get("hEvents");
    h->Draw("colz");
    auto grFlight = RootTools::flightPath(3);
    grFlight->Draw("lsame");
    
    for(int i=0; i < numSizeGroups; i++){
      if(onlySinglets && sizes[i]!=1) continue;	        
      std::cout << sizes[i] << "\t" << grCs[i]->GetN() << std::endl;

      // if(grCs[i]->GetN() > 0){
      // const char* opt = onlySinglets || i == 0 ? "ap pmc" : "psame pmc";
      const char* opt = onlySinglets || i == 0 ? "psame pmc" : "psame pmc";
      grCs[i]->Draw(opt);

      for(int j=0; j < arrows[i].GetEntries(); j++){
	arrows[i].At(j)->Draw();
      }
      // nd++;
      // }
    }
    l1->Draw();
    auto prims = c->GetListOfPrimitives();
    TString bName = TString::Format("fAntarctica%d", nb);
    auto b = (AntarcticaBackground*) prims->FindObject(bName);
    b->SetGrayScale(true);
    b->SetIcemask(true);
    b->SetShowColorAxis(false);
    // b->
    nb++;

    // TString canFileName = TString::Format("~/blind_clusters_llEventCut%d_not_vaguely_near_mcm_runs_above_160.png", TMath::Nint(llEventCuts.at(z)));
    // c->SaveAs(canFileName);


    if(hTheta){
      Acclaim::RootTools::canvas();
      hTheta->Draw();
    }
  }

  
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
  // drawClusters(f);
  drawClusters2(f);  

  // drawResolutionDistributions(f);

  return;
  // printClusterMultiplicityTable(f);

}
