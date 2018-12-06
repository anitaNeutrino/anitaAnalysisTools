
set(BINARIES fitPhaseCentersA4 doWaisNoFilterA4 doNoFilters doSineSub doSineSubNoCache doSineSubPlusBrickWall doRayleighFilter doCutOptimization doAPosterioriSeparation mergeSummaryTreeFiles rmsCacheSineSub rmsCacheSineSubPlusBrickWall editSumTree doClustering doSubthresholdClustering  extendClustering doWaisClustering doDataReduction doDataReductionCalPulser doDataReductionPayloadBlast newPayloadBlastCut newEventFlags getMcReconstructionEfficiency overwriteEventFlags makeThermalTree makeHiCalTree makeSurfaceTree makeSubRun) 

foreach(binary ${BINARIES})
  MESSAGE(STATUS "Process file: ${binary}")
  add_executable (${binary} ${binary}.cxx)
  target_link_libraries (${binary} LINK_PUBLIC AntarcticaRoot AnitaAnalysisTools AnitaAnalysis AnitaCorrelator AnitaEvent UCorrelator)
endforeach(binary)

set(driversBinDir ${UTIL_BIN_DIR}/Acclaim/drivers)
install(DIRECTORY DESTINATION ${driversBinDir})
install(TARGETS ${BINARIES} DESTINATION ${driversBinDir})
