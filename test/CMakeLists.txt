
set(BINARIES testAveragePowerSpectrum testFancyFFTs testFancyTTreeInterpolator testOutputConvention testProgressBar testCmdLineArgs)


foreach(binary ${BINARIES})
  MESSAGE(STATUS "Process file: ${binary}")
  add_executable (${binary} ${binary}.cxx)
  target_link_libraries (${binary} LINK_PUBLIC AnitaAnalysisTools AnitaEvent )
  
endforeach(binary)
set(testBinDir ${UTIL_BIN_DIR}/Acclaim/test)
install(DIRECTORY DESTINATION ${testBinDir})
install(TARGETS ${BINARIES} DESTINATION ${testBinDir})

