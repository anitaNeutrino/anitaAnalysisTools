
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace Acclaim;
#pragma link C++ class Acclaim::FancyTTreeInterpolator+;
#pragma link C++ namespace Acclaim::RootTools;
#pragma link C++ namespace Acclaim::FancyFFTs;
#pragma link C++ class Acclaim::ProgressBar+;
#pragma link C++ class Acclaim::AveragePowerSpectrum+;
#pragma link C++ class Acclaim::AnitaAveragePowerSpectrum+;
#pragma link C++ class Acclaim::AntarcticaMapPlotter+;

#pragma link C++ namespace Acclaim::Clustering;
#pragma link C++ class Acclaim::Clustering::Event+;
#pragma link C++ class Acclaim::Clustering::McEvent+;
#pragma link C++ class Acclaim::Clustering::Cluster+;
#pragma link C++ class Acclaim::Clustering::DistIndex+;
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
#pragma link C++ class Acclaim::AnalysisPlot+;
#pragma link C++ class Acclaim::AnalysisProf+;

// This has gotten out of hand...
#pragma link C++ namespace Acclaim::AnalysisCuts;
#pragma link C++ enum Acclaim::AnalysisCuts::Mode;
#pragma link C++ class Acclaim::AnalysisCuts::AnalysisCut+;
#pragma link C++ class Acclaim::AnalysisCuts::IsAboveHorizontal+;
#pragma link C++ class Acclaim::AnalysisCuts::IsTaggedAsPulser+;
#pragma link C++ class Acclaim::AnalysisCuts::IsTaggedAsWaisPulser+;
#pragma link C++ class Acclaim::AnalysisCuts::IsTaggedAsLDBPulser+;
#pragma link C++ class Acclaim::AnalysisCuts::HigherPol+;
#pragma link C++ class Acclaim::AnalysisCuts::HasSourceLocation+;
#pragma link C++ class Acclaim::AnalysisCuts::IsOnContinent+;
#pragma link C++ class Acclaim::AnalysisCuts::IsTaggedAsPayloadBlast+;
#pragma link C++ class Acclaim::AnalysisCuts::IsWithin20DegreesOfSunInPhi+;
#pragma link C++ class Acclaim::AnalysisCuts::IsGood+;
#pragma link C++ class Acclaim::AnalysisCuts::GoodGPS+;
#pragma link C++ class Acclaim::AnalysisCuts::NonZeroStokesI+;
#pragma link C++ class Acclaim::AnalysisCuts::RealSNR+;
#pragma link C++ class Acclaim::AnalysisCuts::Anita3QuietTime+;
#pragma link C++ class Acclaim::AnalysisCuts::CloseToMC+;
#pragma link C++ class Acclaim::AnalysisCuts::CloseToWais+;
#pragma link C++ class Acclaim::AnalysisCuts::IsRfTrigger+;
#pragma link C++ class Acclaim::AnalysisCuts::SmallDeltaRough+;
#pragma link C++ class Acclaim::AnalysisCuts::IsNotTaggedAsPulser+;
#pragma link C++ class Acclaim::AnalysisCuts::SignalLikeFirstStandardizedPeakMoments+;
#pragma link C++ class Acclaim::AnalysisCuts::PassesThesisCuts+;
#pragma link C++ class Acclaim::AnalysisCuts::IsNotNorth+;
#pragma link C++ class Acclaim::AnalysisCuts::AcceptableHardwareAngle+;
#pragma link C++ class Acclaim::AnalysisCuts::HigherPeakHilbertAfterDedispersion+;
#pragma link C++ class Acclaim::AnalysisCuts::HigherImpulsivityMeasureAfterDedispersion+;
#pragma link C++ class Acclaim::AnalysisCuts::LowerFracPowerWindowGradientAfterDedispersion+;
#pragma link C++ class Acclaim::AnalysisCuts::DedispersedFracPowerWindowGradientBelowThreshold+;
#pragma link C++ class Acclaim::AnalysisCuts::FisherScoreAboveThreshold+;
#pragma link C++ class Acclaim::AnalysisCuts::DoesNotPointToKnownMovingSource+;

#pragma link C++ class Acclaim::CutOptimizer::FisherResult+;
#pragma link C++ class Acclaim::QualityCut+;
#pragma link C++ class Acclaim::SurfSaturationCut+;
#pragma link C++ class Acclaim::PayloadBlastCut+;
#pragma link C++ class Acclaim::NumPointsCut+;

#pragma link C++ class std::vector<Acclaim::Clustering::Event>+;
#pragma link C++ class std::vector<Acclaim::Clustering::McEvent>+;
#pragma link C++ class std::vector<Acclaim::Clustering::Cluster>+;

#endif
