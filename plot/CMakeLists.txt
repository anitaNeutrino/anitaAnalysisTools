
set(BINARIES initialDistributions makeHeatMaps)

foreach(binary ${BINARIES})
  MESSAGE(STATUS "Process file: ${binary}")
  add_executable (${binary} ${binary}.cxx)
  target_link_libraries (${binary} LINK_PUBLIC AnitaAnalysisTools AnitaEvent UCorrelator)
endforeach(binary)

set(driversBinDir ${UTIL_BIN_DIR}/Acclaim/plot)
install(DIRECTORY DESTINATION ${driversBinDir})
install(TARGETS ${BINARIES} DESTINATION ${driversBinDir})
