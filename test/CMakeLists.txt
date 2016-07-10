add_executable (testAveragePowerSpectrum testAveragePowerSpectrum.cxx)
target_link_libraries (testAveragePowerSpectrum LINK_PUBLIC AnitaAnalysisTools AnitaEvent )  
install (TARGETS testAveragePowerSpectrum DESTINATION ${UTIL_BIN_DIR})


add_executable (testCorrelator testCorrelator.cxx)
target_link_libraries (testCorrelator LINK_PUBLIC AnitaAnalysisTools AnitaEvent )  
install (TARGETS testCorrelator DESTINATION ${UTIL_BIN_DIR})



add_executable (testFancyFFTs testFancyFFTs.cxx)
target_link_libraries (testFancyFFTs LINK_PUBLIC AnitaAnalysisTools AnitaEvent )  
install (TARGETS testFancyFFTs DESTINATION ${UTIL_BIN_DIR})



add_executable (testFancyTTreeInterpolator testFancyTTreeInterpolator.cxx)
target_link_libraries (testFancyTTreeInterpolator LINK_PUBLIC AnitaAnalysisTools AnitaEvent )  
install (TARGETS testFancyTTreeInterpolator DESTINATION ${UTIL_BIN_DIR})



add_executable (testOutputConvention testOutputConvention.cxx)
target_link_libraries (testOutputConvention LINK_PUBLIC AnitaAnalysisTools AnitaEvent )  
install (TARGETS testOutputConvention DESTINATION ${UTIL_BIN_DIR})



add_executable (testProgressBar testProgressBar.cxx)
target_link_libraries (testProgressBar LINK_PUBLIC AnitaAnalysisTools AnitaEvent )  
install (TARGETS testProgressBar DESTINATION ${UTIL_BIN_DIR})
