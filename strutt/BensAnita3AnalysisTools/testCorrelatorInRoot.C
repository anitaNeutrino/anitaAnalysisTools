{


  /* Tedious loading of libraries to make this ROOT macro work... */

  gSystem->Load("libAnitaEvent.so"); // load event reader ROOT into CINT
  // gSystem->Load("libMinuit.so"); // required by libAnitaCorrelator...
  // gSystem->Load("libAnitaCorrelator.so"); // load Ryan/Matt Mottram's correlator into CINT  

  // Load a few prerequisite libraries into CINT...
  gSystem->Load("libMathMore.so"); // prereq for Ryan's ROOT FFTW wrapper
  gSystem->Load("libfftw3.so"); // prereq for Ryan's ROOT FFTW wrapper
  gSystem->Load("libRootFftwWrapper.so"); // Ryan's ROOT FFTW wrapper, which is a prereq for Cross Correlator

  gSystem->Load("libBensAnitaTools.so"); // load my correlator and other lovely things into CINT




  TString dataPath = "~/UCL/ANITA/anita3Data/flight1415/telem/";
  const Int_t run = 335;




  // Open ROOTified data files...
  TString fNameEvent = dataPath + TString::Format("run%d/eventFile%d.root", run, run);
  TFile* fEvent = TFile::Open(fNameEvent);
  TTree* tEvent = (TTree*) fEvent->Get("eventTree");
  RawAnitaEvent* rawEvent = NULL;
  tEvent->SetBranchAddress("event", &rawEvent);

  TString fNameHead = dataPath + TString::Format("run%d/eventHeadFile%d.root", run, run);
  TFile* fHead = TFile::Open(fNameHead);
  TTree* tHead = (TTree*) fHead->Get("headTree");
  RawAnitaHeader* header = NULL;
  tHead->SetBranchAddress("header", &header);

  // In telem data the matching of events in eventHeadFile and eventFile isn't necessarily perfect.
  // Need to take account of that by building another set of indices with event number.
  tEvent->BuildIndex("event->eventNumber");

  // Make cross correlator
  CrossCorrelator* cc = new CrossCorrelator();


  
  // Loop over data...

  Long64_t numEntries = tHead->GetEntries();
  for(Long64_t entry=0; entry<numEntries; entry++){

    // Get matching event numbered entries from header and event files
    tHead->GetEntry(entry);
    Int_t eventIndex = tEvent->GetEntryWithIndex(header->eventNumber);

    UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(rawEvent, WaveCalType::kDefault, header);
    
    cc->correlateEvent(usefulEvent); // Generates set of cross correlations
    
    TH2D* hImage = cc->makeImage(AnitaPol::kVertical); // then generate and image
    
    hImage->Draw("colz"); // ... and enjoy it, isn't it pretty?
    break;
  }


}
