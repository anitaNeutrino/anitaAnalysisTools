{
  // two sets of nested brackets for the local scope,
  // so we don't get things left in ROOT memory
  {
    const int nLibs = 7;
    const char* libs[nLibs] = {"libRootFftwWrapper.so",
                               "libAnitaEvent.so",
                               "libAnitaCorrelator.so",
                               "libAnitaAnalysis.so",
                               "libAnitaAnalysisTools.so",
                               "libUCorrelator.so",
                               "libAnitaMagicDisplay.so"};
    for(int i=0; i < nLibs; i++){
      TString loaded = gSystem->GetLibraries();
      if(!loaded.Contains(libs[i])){
        gSystem->Load(libs[i]);
      }
    }
    gROOT->ProcessLine("#include \"FFTtools.h\"");
  }
}
