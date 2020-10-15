{
  // two sets of nested brackets for the local scope,
  // so we don't get things left in ROOT memory
  {
    const int nLibs = 9;
    const char* libs[nLibs] = {"${ANITA_UTIL_INSTALL_DIR}/lib/libRootFftwWrapper.so",
			       "${ANITA_UTIL_INSTALL_DIR}/lib/libAntarcticaRoot.so",
			       "${ANITA_UTIL_INSTALL_DIR}/lib/libAnitaIceMC.so",
                               "${ANITA_UTIL_INSTALL_DIR}/lib/libAnitaEvent.so",
                               "${ANITA_UTIL_INSTALL_DIR}/lib/libAnitaCorrelator.so",
                               "${ANITA_UTIL_INSTALL_DIR}/lib/libAnitaAnalysis.so",
                               "${ANITA_UTIL_INSTALL_DIR}/lib/libAnitaAnalysisTools.so",
                               "${ANITA_UTIL_INSTALL_DIR}/lib/libUCorrelator.so",
                               "${ANITA_UTIL_INSTALL_DIR}/lib/libAnitaMagicDisplay.so"};
    for(int i=0; i < nLibs; i++){
      TString loaded = gSystem->GetLibraries();
      if(!loaded.Contains(libs[i])){
        gSystem->Load(libs[i]);
      }
    }
    gROOT->ProcessLine("#include \"FFTtools.h\"");
  }
}
