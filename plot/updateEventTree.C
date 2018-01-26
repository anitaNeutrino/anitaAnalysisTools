#include "AcclaimClustering.h"

using namespace Acclaim::Clustering;

TKDTreeID* fKDTree;                     /// ROOT's implementation of a KDTree, typedef'd for int/doubleTKDTree* fKDTree = NULL;
std::vector<Acclaim::Clustering::Event> events;
std::vector<double> fEventEastings;
std::vector<double> fEventNorthings;
bool fDebug = true;


void initKDTree(TTree* eventTree){

  events.clear();
  Acclaim::Clustering::Event* event = NULL;
  eventTree->SetBranchAddress("event", &event);
  events.reserve(eventTree->GetEntries());
  for(Long64_t entry=0; entry < eventTree->GetEntries(); entry++){
    eventTree->GetEntry(entry);
    events.push_back(*event);
  }
    

  if(fDebug){
    std::cout << "About to build KDTree" << std::endl;
  }
  fEventEastings.clear();
  fEventNorthings.clear();
  fEventEastings.reserve(events.size());
  fEventNorthings.reserve(events.size());
  for(UInt_t eventInd = 0; eventInd < events.size(); eventInd++){
    const Event& event = events.at(eventInd);
    if(event.eventEventClustering){
      fEventEastings.push_back(event.easting);
      fEventNorthings.push_back(event.northing);
    }
  }

  const int binSize = 100000; // meters... too small?
  if(fKDTree){
    delete fKDTree;
  }  
  fKDTree = new TKDTreeID(events.size(), 2, binSize);
  fKDTree->SetData(0, &fEventEastings[0]);
  fKDTree->SetData(1, &fEventNorthings[0]);
  fKDTree->Build();

  if(fDebug){
    std::cout << "Built!" << std::endl;

    const int nTest = 10;
    std::vector<int> nns(nTest);
    std::vector<double> dists(nTest);
    fKDTree->FindNearestNeighbors(&events.at(0).easting, nTest, &nns[0], &dists[0]);
    std::cout << "The ten nearest neighbours of events[0] at " << events[0].longitude << "," << events[0].latitude << " are:" << std::endl;
    for(int i=0; i < nTest; i++){
      std::cout << "events[" << nns[i] << "] at " << events[nns[i]].longitude << ", " << events[nns[i]].latitude << std::endl;
    }

    std::vector<Int_t> neighbours;
    const double rangeEN = 1000e3;
    double lookup[2] = {events.at(0).easting, events.at(0).northing};
    fKDTree->FindInRange(lookup, rangeEN, neighbours);
    std::cout << "There are " << neighbours.size() << " events within a sphere of radius " << rangeEN << std::endl;
  }  
}


void updateEventTree(const char* fileName){

  TFile* f = TFile::Open(fileName);

  TList* things = (TList*)f->GetListOfKeys();

  TString fName2 = fileName;
  fName2.ReplaceAll(".root", "_updated.root");
  TFile* f2 = new TFile(fName2, "recreate");
  bool doneEventTree = false;
  for(int i=0; i < things->GetEntries(); i++){

    TKey* k = (TKey*)things->At(i);
    TString thingName = k->GetName();

    TObject* thing = f->Get(thingName);

    TTree* tree = dynamic_cast<TTree*>(thing);
    if(tree){

      TString treeName = tree->GetName();
      if(treeName=="eventTree"){
	if(!doneEventTree){
	  initKDTree(tree);

	  TTree* eventTree = new TTree("eventTree", "eventTree");
	  Event* event = NULL;
	  eventTree->Branch("event", &event);

	  const int nTest = 2;
	  std::vector<int> nns(nTest);
	  std::vector<double> dists(nTest);
	
	  Acclaim::ProgressBar p(events.size());
	  for(Long64_t entry=0; entry < events.size(); entry++){
	    event = &events.at(entry);

	    fKDTree->FindNearestNeighbors(&event->easting, nTest, &nns[0], &dists[0]);
	    event->nearestEventSurfaceDistanceKm = 1e-3*dists.at(1);
	    eventTree->Fill();
	    // std::cout << event->eventNumber << "\t" << event->nearestEventSurfaceDistanceKm << std::endl;
	    p.inc(entry);
	  }
	  eventTree->Write();
	  delete eventTree;
	  doneEventTree = true;
	}
      }
      else{
	tree = tree->CloneTree();
	tree->Write();
	delete tree;
      }
    }
    else{
      thing->Write();
    }
  }

  
  f2->Write();
  f2->Close();

  if(fKDTree){
    delete fKDTree;
  }
}
