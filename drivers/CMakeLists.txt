
set(BINARIES analyzeWais filterDecimatedAcclaim filterWaisAcclaim filterDecimatedBrickWallSatellites filterWaisBrickWallSatellites filterDecimatedUniformMag filterWaisUniformMag filterWaisSineSub filterDecimatedSineSub)

foreach(binary ${BINARIES})
  MESSAGE(STATUS "Process file: ${binary}")
  add_executable (${binary} ${binary}.cxx)
  target_link_libraries (${binary} LINK_PUBLIC AnitaAnalysisTools AnitaEvent UCorrelator)
endforeach(binary)

set(driversBinDir ${UTIL_BIN_DIR}/Acclaim/drivers)
install(DIRECTORY DESTINATION ${driversBinDir})
install(TARGETS ${BINARIES} DESTINATION ${driversBinDir})
