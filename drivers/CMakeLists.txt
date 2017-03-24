
add_executable (analyzeWais analyzeWais.cxx)
target_link_libraries (analyzeWais LINK_PUBLIC AnitaAnalysisTools AnitaEvent )  
install (TARGETS analyzeWais DESTINATION ${UTIL_BIN_DIR})
