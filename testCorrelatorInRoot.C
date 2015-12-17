{


  /* Tedious loading of libraries to make this ROOT macro work... */
  gSystem->Load("libMathMore.so"); // prereq for Ryan's ROOT FFTW wrapper
  gSystem->Load("libRootFftwWrapper.so"); // Ryan's ROOT FFTW wrapper, required by EventReaderRoot
  gSystem->Load("libAnitaEvent.so"); // load event reader ROOT into CINT
  gSystem->Load("libBensAnitaTools.so"); // load my correlator and other lovely things into CINT




  TString dataPath = "~/UCL/ANITA/flight1415/root/";
  const Int_t run = 352;//158;




  // Open ROOTified data files...
  TString fNameEvent = dataPath + TString::Format("run%d/calEventFile%d.root", run, run);
  TFile* fEvent = TFile::Open(fNameEvent);
  TTree* tEvent = (TTree*) fEvent->Get("eventTree");
  CalibratedAnitaEvent* calEvent = NULL;
  tEvent->SetBranchAddress("event", &calEvent);

  TString fNameHead = dataPath + TString::Format("run%d/headFile%d.root", run, run);
  TFile* fHead = TFile::Open(fNameHead);
  TTree* tHead = (TTree*) fHead->Get("headTree");
  RawAnitaHeader* header = NULL;
  tHead->SetBranchAddress("header", &header);


  UInt_t eventIWant = 10832108; //10106161;
  Long64_t entryIWant = -1;
  Long64_t numEntries = tHead->GetEntries();
  for(Long64_t entry=0; entry<numEntries; entry++){
    tHead->GetEntry(entry);
    if(header->eventNumber == eventIWant){
      entryIWant = entry;
    }
  }

  tHead->GetEntry(entryIWant);
  tEvent->GetEntry(entryIWant);

  // Make cross correlator
  CrossCorrelator* cc = new CrossCorrelator();
  
  UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);//, WaveCalType::kDefault, header);
    
  cc->correlateEvent(usefulEvent); // Generates set of cross correlations
  
  // TH2D* hImage = cc->makeImage(AnitaPol::kVertical); // then generate and image
  
  TH2D* hImage = cc->makeGlobalImage(AnitaPol::kHorizontal); // then generate and image
  // hImage->SetTitle("A pulser event reconstructed by Cross Correlator");
  // hImage->SetName(TString::Format("hImage%u", eventIWant));
    
  hImage->Draw("colz"); // ... and enjoy it, isn't it pretty?

}
