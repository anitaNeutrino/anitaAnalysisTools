

set(BINARIES analyzeWais filterWaisBrickWallSatellites filterWaisAcclaim)


foreach(binary ${BINARIES})
  MESSAGE(STATUS "Process file: ${binary}")
  add_executable (${binary} ${binary}.cxx)
  target_link_libraries (${binary} LINK_PUBLIC AnitaAnalysisTools AnitaEvent profiler)  
  #install (TARGETS analyzeWais DESTINATION ${UTIL_BIN_DIR})


  #add_executable (analyzeWais analyzeWais.cxx)
  #target_link_libraries (analyzeWais LINK_PUBLIC AnitaAnalysisTools AnitaEvent )  
  #install (TARGETS analyzeWais DESTINATION ${UTIL_BIN_DIR})
  
  #MESSAGE(STATUS "Process file: ${binary}")
  #add_executable(${binary} ${binary}.cxx)
  #target_link_libraries(${binary} ${ZLIB_LIBRARIES} ${ANITA_LIBS} ${ROOT_LIBRARIES} ${FFTW_LIBRARIES})
endforeach(binary)

set(driversBinDir ${UTIL_BIN_DIR}/Acclaim/drivers)
install(DIRECTORY DESTINATION ${driversBinDir})
install(TARGETS ${BINARIES} DESTINATION ${driversBinDir})
