#include "RootTools.h"
using namespace Acclaim;


const TString isSalt = "(eventNumber==14545077||eventNumber==15202247||eventNumber==22377564||eventNumber==26008234||eventNumber==30044333||eventNumber==37387052||eventNumber==39512687||eventNumber==39619438||eventNumber==40040470||eventNumber==49572920||eventNumber==56674128||eventNumber==67024575||eventNumber==71536794||eventNumber==75011721)";

// const TString isSalt = "0";


// const int numSizeGroups = 7;
// const int sizeGroups[numSizeGroups] = {100, 10, 5, 4, 3, 2, 1};

const int numSizeGroups = 11;
const int sizeGroups[numSizeGroups] = {100, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1};

// const int numSizeGroups = 31;
// const int sizeGroups[numSizeGroups] = {100, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
const TCut sizeGroupCuts[numSizeGroups] = {"numDataEvents >= 100",
					   // "numDataEvents >= 30 && numDataEvents < 100",
					   // // "numDataEvents >= 20 && numDataEvents < 100",
					   "numDataEvents >= 10 && numDataEvents < 100",
					   // "numDataEvents >= 5 && numDataEvents < 10",
					   // "numDataEvents == 29",
					   // "numDataEvents == 28",
					   // "numDataEvents == 27",
					   // "numDataEvents == 26",
					   // "numDataEvents == 25",
					   // "numDataEvents == 24",
					   // "numDataEvents == 23",
					   // "numDataEvents == 22",
					   // "numDataEvents == 21",
					   // "numDataEvents == 20",
					   // "numDataEvents == 19",
					   // "numDataEvents == 18",
					   // "numDataEvents == 17",
					   // "numDataEvents == 16",
					   // "numDataEvents == 15",
					   // "numDataEvents == 14",
					   // "numDataEvents == 13",
					   // "numDataEvents == 12",
					   // "numDataEvents == 11",
					   // "numDataEvents == 10",
					   "numDataEvents == 9",
					   "numDataEvents == 8",
					   "numDataEvents == 7",
					   "numDataEvents == 6",
					   "numDataEvents == 5",
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

const int minMaxSmall = 4;
const int maxMaxSmall = 9;  


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

  TString nKm = "0km";
  TString fName = f->GetName();
  if(fName.Contains("km_")){
    std::vector<TString> nameBits;
    Acclaim::RootTools::tokenize(nameBits, fName, "_");
    for(const auto& nameBit : nameBits){
      if(nameBit.Contains("km")){
	nKm = nameBit;
	break;
      }
    }
  }
  std::cout << nKm << std::endl;
  // return;
  

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
	// std::cout << cluster->llEventCut << std::endl;
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
  TGraph* grNumClusters = new TGraph();
  TGraph* grNumSmallClusters = new TGraph();

  TH2D* hKnown = new TH2D("hKnown", "Clusters associated with known bases; Event cut size; Number of events", nz, 0, nz, 10, 1, 11);
  TH2D* hUnknown = new TH2D("hUnknown", "Clusters not associated with known bases; Event cut size; Number of events", nz, 0, nz, 10, 1, 11);

  ProgressBar p(nz);
  TGraphAntarctica* grSinglets0 = NULL;
  TGraphAntarctica* grSinglets1 = NULL;
  TGraphAntarctica* grSinglets2 = NULL;
  TGraphAntarctica* grSalt = NULL;  

  std::vector<TGraph*> grKnowns(numSizeGroups, NULL);
  std::vector<TGraph*> grUnknowns(numSizeGroups, NULL);
  std::vector<TGraph*> grKnownsN(numSizeGroups, NULL);
  std::vector<TGraph*> grUnknownsN(numSizeGroups, NULL);
  std::vector<TGraph*> grKnownsH(numSizeGroups, NULL);
  std::vector<TGraph*> grUnknownsH(numSizeGroups, NULL);
  std::vector<TGraph*> grKnownsV(numSizeGroups, NULL);
  std::vector<TGraph*> grUnknownsV(numSizeGroups, NULL);
  for(int s=0; s < numSizeGroups; s++){
    grKnowns[s] = new TGraph();
    grUnknowns[s] = new TGraph();
    grKnownsN[s] = new TGraph();
    grUnknownsN[s] = new TGraph();
    grKnownsH[s] = new TGraph();
    grUnknownsH[s] = new TGraph();
    grKnownsV[s] = new TGraph();
    grUnknownsV[s] = new TGraph();
    
    for(int z=0; z < nz; z++){
      grKnowns[s]->SetPoint(z, llEventCuts[z], 0);
      grUnknowns[s]->SetPoint(z, llEventCuts[z], 0);      
      grKnownsN[s]->SetPoint(z, llEventCuts[z], 0);
      grUnknownsN[s]->SetPoint(z, llEventCuts[z], 0);      
      grKnownsH[s]->SetPoint(z, llEventCuts[z], 0);
      grUnknownsH[s]->SetPoint(z, llEventCuts[z], 0);      
      grKnownsV[s]->SetPoint(z, llEventCuts[z], 0);
      grUnknownsV[s]->SetPoint(z, llEventCuts[z], 0);      
    }
  }
  
  std::vector<std::vector<int> > binnedClusters(numSizeGroups-1, std::vector<int>());
  
  for(int z=0; z <nz; z++){
    TTree* clusterTree = clusterTrees[z];
    
    const int nc = clusterTree->GetEntries();

    TString histName = TString::Format("hc_%d", z);
    TH1D* hc = new TH1D(histName, histName, nc, 0, nc);
    TH1D* hcH = new TH1D(histName + "H", histName + "H", nc, 0, nc);
    TH1D* hcV = new TH1D(histName + "V", histName + "V", nc, 0, nc);

    TString draw = TString::Format("cluster[%d]>>" + histName, z);
    TString selfLLCut = TString::Format("selfLogLikelihood < %lf", llEventCuts[z]);

    selfLLCut += " && !" + isSalt;
    
    TString drawV = TString::Format("cluster[%d]>>" + histName + "V", z);
    TString selfLLCutV = selfLLCut + " && pol==1";


    TString drawH = TString::Format("cluster[%d]>>" + histName + "H", z);
    TString selfLLCutH = selfLLCut + " && pol==0";
    
    // TString selfLLCut = TString::Format("selfLogLikelihood < %lf && pol==0", llEventCuts[z]);
    // TString selfLLCut = TString::Format("selfLogLikelihood < 0.01");
    // TString selfLLCut = TString::Format("thetaAdjustmentRequired==0 && latitude < -75 && latitude > -80 && longitude > 150 && longitude < 180");
    // TString selfLLCut = TString::Format("thetaAdjustmentRequired==0");
    // TString selfLLCut = "";
    
    eventTree->Draw(draw, selfLLCut, "goff");
    eventTree->Draw(drawV, selfLLCutH, "goff");
    eventTree->Draw(drawH, selfLLCutV, "goff");

    
    int nLarges[nL] = {0};    
    int numNonZeroClusters = 0;

    std::vector<Int_t> nonBaseSingletClusters;
    std::vector<Int_t> baseSingletClusters;
    
    
    for(int c=0; c < nc; c++){
      int nEvents = hc->GetBinContent(c+1);

      if(nEvents > 0){
	numNonZeroClusters++;
      }

      // TString drawPol = TString::Format("pol>>hPol(2, 0, 2)");
      // TString cutPol = TString::Format("cluster[%d]==%d", z, c);
      // eventTree->Draw(drawPol, cutPol,  "goff");
      // TH1F* hPol =  (TH1F*) gROOT->FindObject("hPol");      

      bool knownBase = c < numKnownBases;

      for(int s=0; s < numSizeGroups; s++){
	if(nEvents >= sizeGroups[s]){

	  if(llEventCuts[z]==4 && s != numSizeGroups - 1){
 	    binnedClusters[s].push_back(c);
	  }

	  if(llEventCuts[z]==4 && sizeGroups[s]==2){
	    std::cout << "The cluster you want is " << c << std::endl;
	  }

	  TGraph* grN = knownBase ? grKnownsN[s] : grUnknownsN[s];
	  TGraph* grH = knownBase ? grKnownsH[s] : grUnknownsH[s];
	  TGraph* grV = knownBase ? grKnownsV[s] : grUnknownsV[s];
	  grN->GetY()[z] += hc->GetBinContent(c+1);
	  grH->GetY()[z] += hcH->GetBinContent(c+1);
	  grV->GetY()[z] += hcV->GetBinContent(c+1);

	  TGraph* gr = knownBase ? grKnowns[s] :  grUnknowns[s];
	  gr->GetY()[z] += 1;
	  break;
	}
      }
      
      if(knownBase){
	hKnown->Fill(z, nEvents);
	if(nEvents==1){
	  baseSingletClusters.push_back(c);
	}
      }
      else{
	if(nEvents  > 1){
	  hUnknown->Fill(z, nEvents);

	}
	if(nEvents==1){
	  nonBaseSingletClusters.push_back(c);
	}
      }
      for(int i=0; i < nL; i++){
	if(nEvents == largeSizes[i]){
	  nLarges[i] += nEvents;
	}
      }
    }

    if(llEventCuts[z]==4){
      std::cout << "ll = " << llEventCuts.at(z) << std::endl;
      TString scanCut = "(";
      for(auto c : nonBaseSingletClusters){
    	if(scanCut != "("){
    	  scanCut += " || ";
    	}
    	scanCut += TString::Format("cluster[%d] == %d", z, c);
      }
      scanCut += ")";


      scanCut += " && !" + isSalt;

      TString scanCut2 = "(";
      for(auto c : baseSingletClusters){
    	if(scanCut2 != "("){
    	  scanCut2 += " || ";
    	}
    	scanCut2 += TString::Format("cluster[%d] == %d", z, c);
      }
      scanCut2 += ")";
      scanCut2 += " && !" + isSalt;
      

      // eventTree->Scan("run:eventNumber:pol:theta:longitude:latitude", scanCut, "goff");
      // new TCanvas();
      // eventTree->Draw("theta:TMath::Log10(nearestEventSurfaceDistanceKm)", scanCut, "colz");
      grSinglets0 = new TGraphAntarctica(eventTree, "longitude", "latitude", TCut(scanCut) + TCut("pol==0"));
      grSinglets0->SetName("grSinglets0");
      grSinglets1 = new TGraphAntarctica(eventTree, "longitude", "latitude", TCut(scanCut) + TCut("pol==1"));
      grSinglets1->SetName("grSinglets1");
      grSinglets2 = new TGraphAntarctica(eventTree, "longitude", "latitude", TCut(scanCut2));      
      grSinglets2->SetName("grSinglets2");
      // eventTree->Scan("run:eventNumber:pol:theta:longitude:latitude", scanCut2, "goff");

      grSalt = new TGraphAntarctica(eventTree, "longitude", "latitude", TCut(isSalt));
    }
    

    for(int i=0; i < nL; i++){
      grLarges[i]->SetPoint(z, llEventCuts.at(z), nLarges[i]);
    }
    for(int i=0; i < nL; i++){
      grNumClusters->SetPoint(z, llEventCuts.at(z), numNonZeroClusters);
    }

    delete hc;
    p.inc(z);
  }

  bool varyingSmall = true;
  
  if(varyingSmall){
    
    TH2D* hPolVsSmall = new TH2D("hPolVsSmall", "Polarization content of small clusters;", 27, 2, 29, 2, 0, 2);
    hPolVsSmall->GetYaxis()->SetBinLabel(1, "HPol");
    hPolVsSmall->GetYaxis()->SetBinLabel(2, "VPol");
    hPolVsSmall->GetXaxis()->SetTitle("Cluster size #leq n");

    TString scanCut2 = "";
    int s=binnedClusters.size() - 1;
    // for(const auto& v : binnedClusters){
    for(auto it = binnedClusters.rbegin(); it != binnedClusters.rend(); ++it){      
      // std::cout << "size = " << sizeGroups[s] << std::endl;
      for(auto c : *it){
	// std::cout << c << ", ";
	if(scanCut2 != ""){
	  scanCut2 += " || ";
	}
	scanCut2 += TString::Format("cluster[%d] == %d", 2, c);
      }
      // std::cout << std::endl << std::endl;
      eventTree->Draw("pol>>hPol(2, 0, 2)", scanCut2, "goff");
      TH1F* hPol = (TH1F*) gROOT->FindObject("hPol");
      int xBin = hPolVsSmall->GetXaxis()->FindBin(sizeGroups[s]);

      if(xBin==1){
	eventTree->Scan("run:eventNumber:pol:cluster[2]", scanCut2, "goff");
      }
      
      const double n = hPol->GetBinContent(1) + hPol->GetBinContent(2);
      std::cout << s << "\t" << xBin << "\t" << hPolVsSmall->GetXaxis()->GetBinLowEdge(xBin) << "\t" << hPol->GetBinContent(1) << "\t" << hPol->GetBinContent(2) << std::endl;
      hPolVsSmall->SetBinContent(xBin, 1, hPol->GetBinContent(1)/n);
      hPolVsSmall->SetBinContent(xBin, 2, hPol->GetBinContent(2)/n);
      s--;
    }
    RootTools::canvas();
    hPolVsSmall->Draw("colz");
    hPolVsSmall->SetMinimum(0);
    hPolVsSmall->SetMaximum(1);
    // return;
    TH1D* hTau = new TH1D("hTau", "Distribution of #tau", 10, 0, 2);
    TH1D* hTau2 = new TH1D("hTau2", "Distribution of #tau_{2}", 10, 0, 2);    
    TH1D* hBHat = new TH1D("hBHat", "Distribution of estimated clustering backgrounds; Anthropogenic background estimate", 20, 0, 4);
    TH1D* hBHat2 = new TH1D("hBHat2", "Distribution of estimated clustering backgrounds; Anthropogenic background estimate", 20, 0, 4);    
    
    TGraph* grTauVsSmall = new TGraph();
    grTauVsSmall->SetName("grTauVsSmall");
    grTauVsSmall->SetTitle("How does the background estimate vary as we change what we consider \"small\"?;Upper limit on small (2 #leq n #leq x); #tau_{anthro}");
    grTauVsSmall->SetMinimum(0);

    auto grBHat = new TGraph();
    grBHat->SetName("grBHat");
    grBHat->SetTitle("Background estimate (H+V);Upper limit on small;");

    for(int maxSmall = minMaxSmall; maxSmall < maxMaxSmall ; maxSmall++){

      const int z = 2;
      double knownInt = 0;
      double unknownInt = 0;

      for(int s=0; s < numSizeGroups; s++){
	// if(sizeGroups[s] >= 2 && sizeGroups[s] < 5){
	if(sizeGroups[s] >= 2 && sizeGroups[s] <= maxSmall){
	  knownInt += grKnowns[s]->GetY()[z];
	  unknownInt += grUnknowns[s]->GetY()[z];
	}
      }
      const double polFactor = 2; //1; //2;
      double tau = polFactor*knownInt/unknownInt;
      double tau2 = polFactor*(knownInt-1)/(unknownInt);
      double C = grKnowns[numSizeGroups-1]->GetY()[z];
      double C2 = C + 1;
      double bHat = C/tau;
      double bHat2 = C2/tau2;

      std::cout << "Interesting stuff here!" << maxSmall << "\t" << knownInt << "\t" << unknownInt << "\t" << tau << "\t" << bHat << "\t" << bHat2 << std::endl;
      if(TMath::Finite(tau)){
	grTauVsSmall->SetPoint(grTauVsSmall->GetN(), maxSmall, tau);
	grBHat->SetPoint(grBHat->GetN(), maxSmall, bHat);
      }
      hTau->Fill(tau);
      hTau2->Fill(tau2);
      hBHat->Fill(bHat);
      hBHat2->Fill(bHat2);      
    }
    RootTools::canvas();
    grTauVsSmall->Draw();
    grBHat->Draw("lsame");
    grBHat->SetLineColor(kRed);
    auto l1 = new TLegend(0.8, 0.8, 1, 1);
    l1->AddEntry(grTauVsSmall, "#tau_{anthro}", "l");
    l1->AddEntry(grBHat, "Anthropogenic background estimate", "l");
    l1->Draw();
    RootTools::canvas();
    hTau->Draw();
    // hBHat->Draw();
    hBHat2->SetLineColor(kBlue);
    hBHat2->SetFillStyle(3345);
    hBHat2->SetFillColor(kBlue);    

    hBHat->SetLineColor(kBlack);
    hBHat->SetFillStyle(3354);
    hBHat->SetFillColor(kBlack);    

    TH1D* hs[2] = {hBHat, hBHat2};
    RootTools::drawHistsWithStatsBoxes(2, hs, "", "mre");
    auto l2 = new TLegend();
    l2->AddEntry(hBHat, "Initial estimate", "fl");
    l2->AddEntry(hBHat2, "Removing 70124461", "fl");
    l2->Draw();

    RootTools::canvas();    
    TH1D* hBHatCombined = (TH1D*)hBHat->Clone("hBHatCombined");
    hBHatCombined->Add(hBHat2);
    hBHatCombined->Draw();
    

    return;
  }
  
  
  std::vector<Double_t> zeros(nz, 0);
  std::vector<Double_t> knownInt(llEventCuts.size(), 0);
  std::vector<Double_t> unknownInt(llEventCuts.size(), 0);
  for(int s=0; s < numSizeGroups; s++){
    // if(sizeGroups[s] >= 2 && sizeGroups[s] < 5){
    if(sizeGroups[s] >= 2 && sizeGroups[s] < 20){      
      for(int z=0; z < grKnowns[s]->GetN(); z++){
	knownInt[z] += grKnowns[s]->GetY()[z];
      }
      for(int z=0; z < grUnknowns[s]->GetN(); z++){
	unknownInt[z] += grUnknowns[s]->GetY()[z];
      }
    }
  }


  TGraph* grEstimate = new TGraph();
  TGraph* grTau = new TGraph();
  std::vector<double> taus;
  std::vector<double> bHats;
  // unknown 1 ~ known1 * (known small)/(unknown small)
  std::cout << "Size group\tKnown\tUnknown" << std::endl;
  for(int s=0; s < numSizeGroups; s++){
    for(int z=0; z < grKnowns[s]->GetN(); z++){
      if(llEventCuts[z]==4){
	std::cout << sizeGroups[s] << "\t" << grKnowns[s]->GetY()[z] << "\t" << grUnknowns[s]->GetY()[z] << std::endl;

	std::cout << grKnownsN[s]->GetY()[z] << "(" << grKnownsH[s]->GetY()[z] << ", " << grKnownsV[s]->GetY()[z] << ")\t"
		  << grUnknownsN[s]->GetY()[z] << "(" << grUnknownsH[s]->GetY()[z] << ", " << grUnknownsV[s]->GetY()[z] << ")\t" << std::endl;
	// std::cout << "H" << sizeGroups[s] << "\t" << grKnownsH[s]->GetY()[z] << "\t" << grUnknownsH[s]->GetY()[z] << std::endl;	
	// std::cout << "V" << sizeGroups[s] << "\t" << grKnownsV[s]->GetY()[z] << "\t" << grUnknownsV[s]->GetY()[z] << std::endl;	
	// std::cout << "Known " << sizeGroups[s] << "\t" << grKnowns[s]->GetY()[z] << std::endl;
	// std::cout << "Unknown " << sizeGroups[s] << "\t" << grUnknowns[s]->GetY()[z] << std::endl;
      }
    }
  }
  return;

  for(int s=0; s < numSizeGroups; s++){    
    if(sizeGroups[s] == 1){
      for(int z=0; z < grKnowns[s]->GetN(); z++){
	double m = grKnowns[s]->GetY()[z];
	// if(estimate == 0){estimate = 1;}
	double tau = knownInt[z]/unknownInt[z];
	// tau = TMath::Finite(tau) ?  tau : 0;
	taus.push_back(tau);


	if(TMath::Finite(tau)){
	  grTau->SetPoint(grTau->GetN(), llEventCuts[z], tau);
	}

	double bHat = m/tau;
	bHats.push_back(bHat);
	if(TMath::Finite(bHat)){
	  grEstimate->SetPoint(grEstimate->GetN(),  llEventCuts[z], bHat);
	}

	// double estimate = grKnowns[s]->GetY()[z];
	// // if(estimate == 0){estimate = 1;}
	// estimate *= unknownInt[z]/knownInt[z];
	// estimate = TMath::Finite(estimate) ?  estimate : 0;
	// grEstimate->SetPoint(z,  llEventCuts[z], estimate);
	
	// std::cout << z << "\t"  << llEventCuts[z] << "\t" << tau << "\t" << bHat << std::endl;
      }
    }
  }
  // return;

  
  bool plotEstimates = true;
  if(plotEstimates){
  

    TGraph* grEfficiency = new TGraph();
    grEfficiency->SetPoint(grEfficiency->GetN(), 1,	0.958353);
    grEfficiency->SetPoint(grEfficiency->GetN(), 2,	0.948679);
    grEfficiency->SetPoint(grEfficiency->GetN(), 4,	0.924903);
    grEfficiency->SetPoint(grEfficiency->GetN(), 7,	0.890618);
    grEfficiency->SetPoint(grEfficiency->GetN(), 10,	0.868258);
    grEfficiency->SetPoint(grEfficiency->GetN(), 20,	0.801455);
    grEfficiency->SetPoint(grEfficiency->GetN(), 40,	0.697575);
    grEfficiency->SetPoint(grEfficiency->GetN(), 70,	0.61884);
    grEfficiency->SetPoint(grEfficiency->GetN(), 100,	0.555845);
    grEfficiency->SetPoint(grEfficiency->GetN(), 150,	0.485778);
    grEfficiency->SetPoint(grEfficiency->GetN(), 200,	0.441485);
    grEfficiency->SetPoint(grEfficiency->GetN(), 250,	0.409528);
    grEfficiency->SetPoint(grEfficiency->GetN(), 350,	0.346636);
    grEfficiency->SetPoint(grEfficiency->GetN(), 400,	0.328963);
    grEfficiency->SetPoint(grEfficiency->GetN(), 450,	0.306265);
    grEfficiency->SetPoint(grEfficiency->GetN(), 500,	0.292111);
    grEfficiency->SetPoint(grEfficiency->GetN(), 600,	0.25923);
    grEfficiency->SetPoint(grEfficiency->GetN(), 700,	0.237326);
    grEfficiency->SetPoint(grEfficiency->GetN(), 800,	0.221766);
    grEfficiency->SetPoint(grEfficiency->GetN(), 900,	0.20589);
    grEfficiency->SetPoint(grEfficiency->GetN(), 1000,	0.189423);
    grEfficiency->SetPoint(grEfficiency->GetN(), 1200,	0.160723);
    grEfficiency->SetPoint(grEfficiency->GetN(), 1400,	0.131399);
    grEfficiency->SetPoint(grEfficiency->GetN(), 1600,	0.111451);

    TRolke tr;   //
    double alpha = 0.9; // Confidence Level
    tr.SetCL(alpha);

    TGraph* grUpperLimit = new TGraph();
    grUpperLimit->SetName("grUpperLimit");
    grUpperLimit->SetMinimum(0);

    TGraph* grSensitivity = new TGraph();

    TH1* hb_t = NULL;
    TH1* hb_a = NULL;
    for(int z=0; z < llEventCuts.size(); z++){
      // for(int z=0; z < 3; z++){
      //-----------------------------------------------
      // Model 4 assumes:
      //
      // Poisson uncertainty in the background estimate
      // known efficiency
      //
      int y = grKnowns[6]->GetY()[z];       // events observed in the background region
      // int x = 0;       // events in the signal region
      int x = bHats[z];       // events in the signal region    
      double tau = taus[z];     // ratio between size of signal/background region
      double e = grEfficiency->GetY()[z]; //0.8;    // efficiency
      tr.SetPoissonBkgKnownEff(x,y,tau,e);
      if(TMath::Finite(tau)){
	double ul, ll;
	// tr.GetSensitivity(ll,ul);
	tr.GetLimits(ll,ul);
	grUpperLimit->SetPoint(grUpperLimit->GetN(), llEventCuts[z], ul);
	std::cout << "\n\n\n\n\n\n\n\n\n";
	std::cout << "z = " << z << ", llThreshold = " << llEventCuts[z] << "\t[" << ll << "," << ul << "]" << ", y = " << y << ", x = " << x << ", tau = " << tau << ", e = " << e << std::endl;

	SensitivityCalculator sc(alpha, false);
	sc.setAnthroBackground(y, tau, 1);
	sc.setThermalBackground(0, 3.2);
	sc.setEfficiency(e, 0);

	sc.getLimit(x, &ll, &ul);

	if(!hb_a){
	  hb_a = sc.histBAnthro();
	}
	if(!hb_t){
	  hb_t = sc.histBThermal();
	}

	grSensitivity->SetPoint(grSensitivity->GetN(), llEventCuts[z], ul);
      }
      else{
	// std::cout << "z = " << z << ", skipping infinite tau" << std::endl;
      }    
    }

  
    RootTools::canvas(1);
    grUpperLimit->SetTitle(TString::Format("TRolke %d%% upper limit vs. log likelihood; Clustering log-likelihood threshold (no units); Expected upper limit", TMath::Nint(100*alpha)));
    grUpperLimit->SetLineWidth(2);
    grUpperLimit->SetName(TString("grUpperLimit") + nKm);
    grUpperLimit->Draw();
    grSensitivity->SetLineStyle(2);
    grSensitivity->SetName(TString("grSensitivity") + nKm);
    grSensitivity->SetLineWidth(2);
    grSensitivity->Draw("lsame");
    auto lExpectedLims = new TLegend(0.8, 0.8, 1, 1);
    lExpectedLims->AddEntry(grUpperLimit, "TRolke " + nKm);
    lExpectedLims->AddEntry(grSensitivity, "SensitivityCalculator " + nKm);
    lExpectedLims->Draw();

    // return;
    // RootTools::canvas();
    // hb_a->Draw();
    // RootTools::canvas();
    // hb_t->Draw();
    // return;


  
  
    Acclaim::RootTools::canvas(1);
    auto l = new TLegend(0.45, 0.7, 1, 1);

  
    for(int s=0; s < numSizeGroups; s++){

      const char* opt = s == 0 ? "alp" : "lpsame";

      grKnowns[s]->Draw(opt);
      grKnowns[s]->SetTitle(TString::Format("Cluster multiplicity and ABCD background estimate (%s surface clustering); -2 Log (L) threshold; Number of clusters", nKm.Data()));
      grKnowns[s]->SetMinimum(0);
      if(sizeGroups[s] != 1){
	grUnknowns[s]->Draw("lp same");
      }

      grUnknowns[s]->SetLineStyle(2);
      grKnowns[s]->SetLineWidth(2);
      grUnknowns[s]->SetLineWidth(2);

      grKnowns[s]->SetLineColor(sizeGroupColors[s]);
      grUnknowns[s]->SetLineColor(sizeGroupColors[s]);

      grKnowns[s]->SetMarkerColor(sizeGroupColors[s]);
      grUnknowns[s]->SetMarkerColor(sizeGroupColors[s]);
      grKnowns[s]->SetMarkerStyle(sizeGroupMarkers[s]);
      grUnknowns[s]->SetMarkerStyle(sizeGroupMarkers[s]);

      for(int isKnown=1; isKnown >= 0; isKnown--){
	TString knownWord = isKnown ? "Known" : "Unknown";
	TGraph* gr = isKnown ? grKnowns[s] : grUnknowns[s];
	if(!(isKnown==0 && sizeGroups[s]==1)){
	  if(s > 0){
	    if(sizeGroups[s] + 1 == sizeGroups[s-1]){
	      l->AddEntry(gr, knownWord + TString::Format(" n == %d", sizeGroups[s]), "lp");
	    }
	    else{
	      l->AddEntry(gr, knownWord + TString::Format(" n >= %d && n < %d", sizeGroups[s], sizeGroups[s-1]), "lp");
	    }
	  }
	  else{
	    l->AddEntry(gr, knownWord + TString::Format(" n >= %d", sizeGroups[s]), "lp");
	  }
	}
      }
    }
    grEstimate->Draw("lsame");
    auto lExtra = new TLegend(0.8, 0, 1, 0.2);
    lExtra->AddEntry(grEstimate, "Background estimate", "lp");
    grEstimate->SetLineWidth(3);
    grTau->Draw("lsame");
    lExtra->AddEntry(grTau, "tau", "lp");
    lExtra->SetNColumns(2);
    grTau->SetLineWidth(3);
    grTau->SetLineStyle(2);
  
    l->Draw();  
    l->SetNColumns(2);
    lExtra->Draw();
  }  
  return;

  if(grSinglets0 && grSinglets1 && grSinglets2){
    grSinglets0->SetMarkerColor(kRed);
    grSinglets1->SetMarkerColor(kCyan);
    grSinglets2->SetMarkerColor(kMagenta);
    
    new TCanvas();
    TH2DAntarctica* hEvents = (TH2DAntarctica*)f->Get("hEvents");    
    if(hEvents){

      auto l1 = new TLegend(0.8, 0.8, 1, 1);
      l1->AddEntry(grSinglets0, "Isolated HPol singlets", "p");
      l1->AddEntry(grSinglets1, "Isolated VPol singlets", "p");
      // l1->AddEntry(grSinglets1, "Isolated VPol singlets", "p");

      gPad->SetLogz(1);
      hEvents->Draw("colz");
      l1->Draw();
      
      gPad->Modified();
      gPad->Update();
      hEvents->ShowBackgroundColorAxis(false);
      hEvents->SetIcemask(true);
      hEvents->SetGrayScale(true);
      
      Acclaim::RootTools::flightPath()->Draw("lsame");
      grSinglets0->Draw("psame");
      grSinglets1->Draw("psame");
      // if(grSinglets2->GetN() > 0){
      // 	grSinglets2->Draw("psame");
      // }      
      
    }
    else{
      grSinglets0->Draw("p");
      grSinglets1->Draw("psame");
      // if(grSinglets2->GetN() > 0){
      // 	grSinglets2->Draw("psame");
      // }
    }

    if(hEvents){
      new TCanvas();
      hEvents->Draw("colz");
      grSalt->SetMarkerColor(kYellow);
      grSalt->Draw("psame");
    }
    
    // return;
  }
  // return;

  RootTools::canvas(1);
  for(int i=0; i < nL; i++){

    // double z0 = grLarges[i]->GetY()[0];
    // for(int z=0; z < grLarges[i]->GetN(); z++){
    //   grLarges[i]->GetY()[z] -= z0;
    // }

    TGraph* grGrad = Acclaim::RootTools::makeDerivativeTGraph(grLarges[i]);
    // for(int z=0; z < grGrad->GetN(); z++){
    //   grGrad->GetY()[z] *= -1;
    // }
    
    const char* opt = i == 0 ? "al" : "lsame";
    // grLarges[i]->Draw(opt);
    grGrad->Draw(opt);
  }
  RootTools::canvas();
  hKnown->Draw("colz");
  RootTools::canvas();
  hUnknown->Draw("colz");
  RootTools::canvas(1);
  grNumClusters->SetTitle("Number of clusters; -2 log L threshold (no units); Number of clusters");
  grNumClusters->Draw("alp");

  // TGraph* grGradNumClusters = Acclaim::RootTools::makeDerivativeTGraph(grNumClusters);
  // grGradNumClusters->SetTitle("Number of clusters derivative; -2 log L (no units); #frac{dN_{c}}{d (-2 log L)}");
  // grGradNumClusters->Draw("alp");

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
    // for(auto& a : arrows){
    //   static int i=0; 
    //   a.SetName(TString::Format("list%d", i));
    //   i++;
    // }
    
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
    TGraphAntarctica* grAnita = (TGraphAntarctica*) f->Get("grAnita");
    TGraphAntarctica* grFlight = Acclaim::RootTools::flightPath();
    
    auto grAnitaToEvent1 = new TGraph(1,  &grAnita->GetX()[0], &grAnita->GetY()[0]);
    grAnitaToEvent1->SetPoint(grAnitaToEvent1->GetN(), gr1->GetX()[0], gr1->GetY()[0]);
    auto grAnitaToEvent2 = new TGraph(1,  &grAnita->GetX()[1], &grAnita->GetY()[1]);
    grAnitaToEvent2->SetPoint(grAnitaToEvent2->GetN(), gr2->GetX()[0], gr2->GetY()[0]);


    gr1->SetMarkerStyle(8); gr1->SetMarkerColor(kCyan); gr1->SetLineColor(gr1->GetMarkerColor());
    gr2->SetMarkerStyle(8); gr2->SetMarkerColor(kMagenta); gr2->SetLineColor(gr2->GetMarkerColor());
    grMinPos->SetMarkerStyle(8); grMinPos->SetMarkerColor(kRed);
    grWalk->SetLineWidth(3); grWalk->SetLineColor(kRed); grWalk->SetMarkerColor(grWalk->GetLineColor()); grWalk->SetMarkerStyle(8);
    grWais->SetMarkerStyle(8); grWais->SetMarkerColor(kOrange);
    grAnita->SetMarkerStyle(8); grAnita->SetMarkerColor(kBlack);

    grAnitaToEvent1->SetLineColor(gr1->GetMarkerColor());
    grAnitaToEvent2->SetLineColor(gr2->GetMarkerColor());    
    
    
    auto c1 = Acclaim::RootTools::canvas();
    h->Draw("colz");
    gr1->Draw("psame");
    gr2->Draw("psame");
    grAnita->Draw("psame");
    grFlight->Draw("lsame");
    grAnitaToEvent1->Draw("lsame");
    grAnitaToEvent2->Draw("lsame");
    grMinPos->Draw("psame");
    grWalk->Draw("lsame");
    grWais->Draw("psame");
    
    
    c1->SetLogz(1);

    auto l1 = new TLegend();
    l1->AddEntry(grWais, "True WAIS pulser position", "p");
    l1->AddEntry(gr1, "WAIS pulse a - reconstructed position", "lp");
    l1->AddEntry(gr2, "WAIS pulse b - reconstructed position", "lp");
    l1->AddEntry(grWalk, "Minuit's path to minimum", "lp");
    l1->AddEntry(grAnita, "ANITA flight path", "lp");
    
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
