





void drawAnalysisCutResults(const char* fileName, bool asPercentage=false){

  TFile* f = TFile::Open(fileName);

  if(!f){
    return -1;
  }

  const char* treeName = "analysisCutTree";
  TTree* t = (TTree*) f->Get(treeName);
  if(!t){
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", can't find tree named " << treeName << " in " << f->GetName() << std::endl;
    return 1;
  }


  std::vector<TString> cutNames;
  TObjArray* l = t->GetListOfBranches();
  TIter next(l);  
  while (TBranch* branch = (TBranch*)next()){
    TString branchName = branch->GetName();
    if(branchName == "eventNumber" || branchName == "run"){
      continue;
    }

    cutNames.push_back(branchName);
  }

  const Long64_t totalEvents = t->GetEntries();

  const char* col = "|";
  const char* bdr = "+";
  const char* row = "-";

  const char* pc = asPercentage ? "%" : "";
  const double toPercentage = 100./totalEvents;

  std::cout << "There are " << totalEvents << " events " << std::endl;
  std::cout << col << "cut" << col << "in seq" << col << "if only" << col << "if not" << col << std::endl;
  std::cout << col << row   << bdr << row      << bdr << row       << bdr << row      << col << std::endl;
  TCut inSeqCut = "1";

 
  for(const auto& cutName : cutNames){
    
    TString cutString = cutName + " > 0";
    TCut ifOnlyCut(cutString);

    inSeqCut += ifOnlyCut; //

    TCut ifNotCut = "1";
    for(auto& cutName2 : cutNames){
      if(cutName2 != cutName){
	TString cutString2 = cutName2 + " > 0";
	ifNotCut += cutString2.Data();
      }
    }

    t->Draw(">>inSeq", inSeqCut, "entrylist");
    TEntryList* inSeq = (TEntryList*) gDirectory->Get("inSeq");

    t->Draw(">>ifOnly", ifOnlyCut, "entrylist");
    TEntryList* ifOnly = (TEntryList*) gDirectory->Get("ifOnly");

    t->Draw(">>ifNot", ifNotCut, "entrylist");
    TEntryList* ifNot = (TEntryList*) gDirectory->Get("ifNot");
    
    double nInSeq  = inSeq->GetN();
    double nIfOnly = ifOnly->GetN();
    double nIfNot  = ifNot->GetN();        

    if(asPercentage){
      nInSeq  *= toPercentage;
      nIfOnly *= toPercentage;
      nIfNot  *= toPercentage;
      std::cout << std::fixed;
      std::cout << std::setprecision(2);
    }
    else{
      std::cout << std::setprecision(10);
    }
    
    
    std::cout << col << cutName << col << nInSeq << pc << col << nIfOnly << pc << col << nIfNot << pc << col << std::endl;
    
  }

  
  
  
  
}
