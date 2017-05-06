GitError! There was a problem with the editor '/usr/bin/emacsclient --socket-name=/tmp/emacs1000/server'.  [Type `$' for details]
Head:     master Add cut on number of points to prevent interpolation core dump
Merge:    origin/master Add cut on number of points to prevent interpolation core dump

Staged changes (7)
new file   drivers/filterDecimatedBrickWallSatellites.cxx
@@ -0,0 +1,25 @@
+#include "AnalysisFlow.h"
+#include "BasicFilters.h"
+#include "AcclaimFilters.h"
+
+using namespace Acclaim;
+
+int main(int argc, char* argv[]){
+
+  AnitaVersion::set(3);
+  int run = argc > 1 ? atoi(argv[1]) : 352;
+
+  std::map<TString, FilterStrategy*> filterStrats;
+  bool saveOutput = true;  
+  Filters::appendFilterStrategies(filterStrats, saveOutput);
+
+  FilterStrategy* stupidNotchStrat = Filters::findStrategy(filterStrats, "BrickWallSatellites");
+  if(!stupidNotchStrat){ 
+    std::cerr << "Well, this script is pointless... I give up." << std::endl;
+    return 1;
+  }
+  AnalysisFlow analysisSimpleNotchFilters(argv[0], run, AnalysisFlow::kDecimated, stupidNotchStrat);
+  analysisSimpleNotchFilters.doAnalysis();
+    
+  return 0;
+}
modified   include/AcclaimFilters.h
@@ -105,6 +105,17 @@ namespace Acclaim
     };
 
 
+    class UniformMagnitude : public UniformFilterOperation {
+    public:
+      explicit UniformMagnitude();
+      virtual ~UniformMagnitude() { ;}
+      virtual void processOne(AnalysisWaveform* fw);
+      // virtual unsigned nOutputs() const {return 0;}
+      virtual const char * tag () const {return "UniformMagnitude";};
+      virtual const char * description () const {return "Gives every frequency bin the same magnitude, keeping the phase constant";}
+    };
+    
+    
 
     class SpikeSuppressor : public UniformFilterOperation {
     protected:
modified   include/FourierBuffer.h
@@ -139,6 +139,7 @@ namespace Acclaim
     mutable std::vector<TGraphFB> grNDFs[AnitaPol::kNotAPol]; // for drawSummary
     mutable std::vector<TGraphFB> grSpectrumAmplitudes[AnitaPol::kNotAPol]; // for drawSummary
     mutable std::vector<TGraphFB> grAmplitudes[AnitaPol::kNotAPol]; // for drawSummary
+    mutable std::vector<TGraphFB> grLastAmps[AnitaPol::kNotAPol]; // for drawSummary    
     mutable std::vector<TGraphFB> grProbs[AnitaPol::kNotAPol]; // for drawSummary
 
     
modified   src/AcclaimFilters.cxx
@@ -52,6 +52,8 @@ void Acclaim::Filters::appendFilterStrategies(std::map<TString, FilterStrategy*>
   double reducedChiSquareThresh = 3;
   RayleighFilter* rf = new RayleighFilter(log10ProbThresh, reducedChiSquareThresh, 1500, alfaLowPassFreq);
 
+  UniformMagnitude* um = new UniformMagnitude();
+
   // then make the strategies
   
   FilterStrategy* stupidNotchStrat = new FilterStrategy();
@@ -65,6 +67,13 @@ void Acclaim::Filters::appendFilterStrategies(std::map<TString, FilterStrategy*>
   fs->addOperation(alfaFilter, saveOutput);
   fs->addOperation(rf, saveOutput);
   filterStrats["RayleighFilter"] = fs;
+
+
+  FilterStrategy* ufs = new FilterStrategy();
+  ufs->addOperation(alfaFilter, saveOutput);  
+  ufs->addOperation(um);
+  filterStrats["UniformMagnitude"] = ufs;
+  
 }
 
 
@@ -396,6 +405,18 @@ double Acclaim::Filters::SpikeSuppressor::interpolate_dB(double x, double xLow,
 
 
 
+Acclaim::Filters::UniformMagnitude::UniformMagnitude(){
+}
+
+void Acclaim::Filters::UniformMagnitude::processOne(AnalysisWaveform* wf){
+
+  FFTWComplex* fft = wf->updateFreq();  
+  for(int i=0; i< wf->Nfreq(); i++){
+    fft[i].setMagPhase(1, fft[i].getPhase());
+  }
+}
+
+
 
 
 
modified   src/AnalysisFlow.cxx
@@ -287,21 +287,22 @@ void Acclaim::AnalysisFlow::doAnalysis(){
       Filters::makeFourierBuffersLoadHistoryOnNextEvent(fFilterStrat);
     }
 
+
     Adu5Pat* pat = fData->gps();
     UsefulAdu5Pat usefulPat(pat);
 
-    FilteredAnitaEvent filteredEvent(usefulEvent, fFilterStrat, pat, header, false);
-
-    // since we now have rolling averages make sure the filter strategy is processed before deciding whether or not to reconstruct 
-    Bool_t selectedEvent = shouldIDoThisEvent(header, &usefulPat);
+    Bool_t needToReconstruct = shouldIDoThisEvent(header, &usefulPat);
 
-    if(selectedEvent){
+    if(needToReconstruct){
       eventSummary = new AnitaEventSummary(header, &usefulPat);
+      Bool_t isGoodEvent = QualityCut::applyAll(usefulEvent, eventSummary);
 
-      QualityCut::applyAll(usefulEvent, eventSummary);
-      
+      if(!isGoodEvent){
+	// since we now have rolling averages make sure the filter strategy is sees every event before deciding whether or not to reconstruct
+	FilteredAnitaEvent filteredEvent(usefulEvent, fFilterStrat, pat, header, false);
 	// fReco->reconstructEvent(&filteredEvent, usefulPat, eventSummary);
 	fReco->process(&filteredEvent, &usefulPat, eventSummary);
+      }
 
       fSumTree->Fill();
       delete eventSummary;
modified   src/FourierBuffer.cxx
@@ -99,10 +99,10 @@ void Acclaim::FourierBuffer::initVectors(int n, double df){
       grProbs[pol].push_back(TGraphFB(this, ant, pol, n));
 
       grAmplitudes[pol].push_back(TGraphFB(this, ant, pol, n));
+      grLastAmps[pol].push_back(TGraphFB(this, ant, pol, n));      
       grSpectrumAmplitudes[pol].push_back(TGraphFB(this, ant, pol, n));
 
 
-
       for(int freqBin=0; freqBin < n; freqBin++){
 	double f = df*freqBin;
 	grChiSquares[pol][ant].GetX()[freqBin] = f;
@@ -111,6 +111,7 @@ void Acclaim::FourierBuffer::initVectors(int n, double df){
 
 	grAmplitudes[pol][ant].GetX()[freqBin] = f;	
 	grSpectrumAmplitudes[pol][ant].GetX()[freqBin] = f;
+	grLastAmps[pol][ant].GetX()[freqBin] = f;	
 
 	grProbs[pol][ant].GetX()[freqBin] = f;
 
@@ -119,6 +120,8 @@ void Acclaim::FourierBuffer::initVectors(int n, double df){
 	grReducedChiSquares[pol][ant].GetY()[freqBin] = 0;
 	
 	grAmplitudes[pol][ant].GetY()[freqBin] = 0;
+	grLastAmps[pol][ant].GetY()[freqBin] = 0;	
+	
 	grSpectrumAmplitudes[pol][ant].GetY()[freqBin] = 0;
 	grProbs[pol][ant].GetY()[freqBin] = 0;
 	
@@ -133,7 +136,11 @@ void Acclaim::FourierBuffer::initVectors(int n, double df){
 
       // get them to know about eachother
       grSpectrumAmplitudes[pol][ant].fDerivedFrom = &grAmplitudes[pol][ant];
+      grLastAmps[pol][ant].fDerivedFrom = &grAmplitudes[pol][ant];
+
       grAmplitudes[pol][ant].fDerivatives.push_back(&grSpectrumAmplitudes[pol][ant]);
+      grAmplitudes[pol][ant].fDerivatives.push_back(&grLastAmps[pol][ant]);      
+
     }
   }
 
@@ -249,6 +256,8 @@ size_t Acclaim::FourierBuffer::add(const FilteredAnitaEvent* fEv){
 	  // }
 	  
 	  bool updated = hRays[pol][ant].at(freqInd)->add(amp);
+	  grLastAmps[pol][ant].GetY()[freqInd] = amp;
+	  
 	  if(updated){
 	    anyUpdated = true;
 	    hRays[pol][ant].at(freqInd)->getRayleighFitParams(fitAmplitudes[pol][ant][freqInd],
@@ -281,6 +290,7 @@ size_t Acclaim::FourierBuffer::add(const FilteredAnitaEvent* fEv){
 	    
 	  }
 	  
+
 	  // this is N, if the probability of getting this amplitude or higher, given the rayeligh distribution is 1/N.
 
 
@@ -686,10 +696,17 @@ void Acclaim::FourierBuffer::drawSummary(TPad* pad, SummaryOption_t summaryOpt)
 
       if(summaryOpt == FourierBuffer::RayleighAmplitude){
 	grSpectrumAmplitudes[pol][ant].SetLineColor(gr->GetLineColor());
-	grSpectrumAmplitudes[pol][ant].SetLineStyle(2);
+	grSpectrumAmplitudes[pol][ant].SetLineStyle(3);
 	grSpectrumAmplitudes[pol][ant].SetMaximum(yMax);
 	grSpectrumAmplitudes[pol][ant].SetMinimum(yMin);
 	grSpectrumAmplitudes[pol][ant].Draw("lsame");
+
+	grLastAmps[pol][ant].SetLineColor(gr->GetLineColor());
+	grLastAmps[pol][ant].SetLineStyle(5);
+	grLastAmps[pol][ant].SetMaximum(yMax);
+	grLastAmps[pol][ant].SetMinimum(yMin);
+	grLastAmps[pol][ant].Draw("lsame");
+	
       }           
     }
   }
modified   src/QualityCut.cxx
@@ -18,7 +18,11 @@ Bool_t Acclaim::QualityCut::applyAll(const UsefulAnitaEvent* usefulEvent, AnitaE
   NumPointsCut npc;
   npc.apply(usefulEvent, sum);
 
-  return (ssc.eventPassesCut && stbc.eventPassesCut && npc.eventPassesCut);
+  Bool_t isGood = (ssc.eventPassesCut && stbc.eventPassesCut && npc.eventPassesCut);
+  if(sum){
+    sum->flags.isGood = isGood;
+  }  
+  return isGood;
 }
 
 
@@ -196,6 +200,7 @@ Acclaim::NumPointsCut::NumPointsCut(){
   // I just want to stop core dumps with old root when there aren't enough points for interpolation
   // which is < 5 points, so this is more than sufficient.
   numPointsCutLow = 200;
+  description = "Checks there are a reasonable number of points in each waveform.";  
 }
 
 
@@ -209,12 +214,14 @@ void Acclaim::NumPointsCut::apply(const UsefulAnitaEvent* useful, AnitaEventSumm
     }
   }
 
-  // silly old flag name
   if(sum){
+    if(!eventPassesCut){
+      // silly old flag name      
       sum->flags.isVarner2 = 1;
     }
     else{
       sum->flags.isVarner2 = 0;
     }
   }
+}
 

