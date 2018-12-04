
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace Acclaim;
#pragma link C++ class Acclaim::FancyTTreeInterpolator+;
#pragma link C++ namespace Acclaim::RootTools;
#pragma link C++ namespace Acclaim::FancyFFTs;
#pragma link C++ class Acclaim::ProgressBar+;
#pragma link C++ class Acclaim::CorrelationPair+;

#pragma link C++ namespace Acclaim::Clustering;
#pragma link C++ class Acclaim::Clustering::Event+;
#pragma link C++ class Acclaim::Clustering::McEvent+;
#pragma link C++ class Acclaim::Clustering::Cluster+;
#pragma link C++ class Acclaim::Clustering::LogLikelihoodMethod+;

#pragma link C++ class Acclaim::RayleighHist+;
#pragma link C++ enum Acclaim::RayleighHist::FitMethod;
#pragma link C++ class Acclaim::TGraphFB+;

#pragma link C++ class Acclaim::AnalysisFlow+;
#pragma link C++ class Acclaim::AnalysisReco+;
#pragma link C++ class Acclaim::InterferometricMap+;

#pragma link C++ class Acclaim::SummarySet+;
#pragma link C++ class Acclaim::SummarySelector+;
#pragma link C++ class Acclaim::SumTreeReductionSelector+;
#pragma link C++ class Acclaim::CutTreeSelector+;
//#pragma link C++ namespace Acclaim::AnalysisCuts;
#pragma link C++ namespace Acclaim::SumTree;
#pragma link C++ namespace Acclaim::ThermalTree;
#pragma link C++ class Acclaim::CutOptimizer::FisherResult+;
#pragma link C++ class Acclaim::QualityCut+;
#pragma link C++ class Acclaim::SurfSaturationCut+;
#pragma link C++ class Acclaim::PayloadBlastCut+;
#pragma link C++ class Acclaim::NumPointsCut+;

#pragma link C++ class std::vector<Acclaim::Clustering::Event>+;
#pragma link C++ class std::vector<Acclaim::Clustering::McEvent>+;
#pragma link C++ class std::vector<Acclaim::Clustering::Cluster>+;

#endif
