#include "QualityCut.h"

ClassImp(Acclaim::QualityCut);
ClassImp(Acclaim::SurfSaturationCut);
ClassImp(Acclaim::PayloadBlastCut);
ClassImp(Acclaim::NumPointsCut);



Bool_t Acclaim::QualityCut::applyAll(const UsefulAnitaEvent* usefulEvent, AnitaEventSummary* sum){
      
  SurfSaturationCut ssc;
  ssc.apply(usefulEvent, sum);

  PayloadBlastCut stbc;
  stbc.apply(usefulEvent, sum);

  NumPointsCut npc;
  npc.apply(usefulEvent, sum);

  Bool_t isGood = (ssc.eventPassesCut && stbc.eventPassesCut && npc.eventPassesCut);
  if(sum){
    sum->flags.isGood = isGood;
  }  
  return isGood;
}



Acclaim::SurfSaturationCut::SurfSaturationCut(){
  maxLimit = 2000;
  minLimit = -2000;
  asymLimit = 500;

  maxVolts = 0;
  maxVoltsAnt = -1;
  maxVoltsPol = AnitaPol::kNotAPol;

  minVolts = 0;
  minVoltsAnt = -1;
  minVoltsPol = AnitaPol::kNotAPol;
    
  asymVolts = 0;
  asymVoltsAnt = -1;
  asymVoltsPol = AnitaPol::kNotAPol;

  eventPassesCut = true;
  description = "Surf saturation cut, looks for max volts / min volts and the asymmetry between them";
}
    


// void Acclaim::SurfSaturationCut::apply(FilteredAnitaEvent* fEv){
void Acclaim::SurfSaturationCut::apply(const UsefulAnitaEvent* useful, AnitaEventSummary* sum){  

  maxVolts = 0;
  minVolts = 0;
  for(int pol=0; pol <AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      const int chanIndex = AnitaGeomTool::getChanIndexFromAntPol(ant, (AnitaPol::AnitaPol_t) pol);
      const double* volts = useful->fVolts[chanIndex];
      int n = useful->fNumPoints[chanIndex];
      
      double maxThisChan = TMath::MaxElement(n, volts);
      double minThisChan = TMath::MinElement(n, volts);

      if(maxThisChan > maxVolts){
	maxVolts = maxThisChan;
	maxVoltsAnt = ant;
	maxVoltsPol = (AnitaPol::AnitaPol_t) pol;
      }
      if(minThisChan < minVolts){
	minVolts = minThisChan;
	minVoltsAnt = ant;
	minVoltsPol = (AnitaPol::AnitaPol_t) pol;
      }

      double asymThisChan = TMath::Abs(maxVolts + minVolts);
      if(asymThisChan > asymVolts){
	asymVolts = asymThisChan;
	asymVoltsAnt = ant;
	asymVoltsPol = (AnitaPol::AnitaPol_t) pol;
      }	  
    }
  }

  if(maxVolts > maxLimit || minVolts < minLimit || asymVolts > asymLimit){
    eventPassesCut = false;
    // std::cerr << "failed!:" << std::endl;
    // std::cerr << "\t" << maxVolts << "\t" << maxVoltsAnt << "\t" << maxVoltsPol << std::endl;
    // std::cerr << "\t" << minVolts << "\t" << minVoltsAnt << "\t" << minVoltsPol << std::endl;
    // std::cerr << "\t" << asymVolts << "\t" << asymVoltsAnt << "\t" << asymVoltsPol << std::endl;            
  }
  else{
    eventPassesCut = true;
  }

  if(sum!=NULL){
    if(eventPassesCut){
      sum->flags.isPayloadBlast = 0;
    }
    else{
      sum->flags.isPayloadBlast = 1;
    }
  }  
}










Acclaim::PayloadBlastCut::PayloadBlastCut(){
  ratioCutHigh = 2.8;
  ratioCutLow = 1.14;
  maxRatio = 0;
  maxRatioPhi = -1;
  maxRatioPol = AnitaPol::kNotAPol;
  eventPassesCut = true;
  description = "Removes events with where any peak-to-peak in the top ring is significantly greater than the bottom ring.";
      
}



void Acclaim::PayloadBlastCut::apply(const UsefulAnitaEvent* useful, AnitaEventSummary* sum){
// void Acclaim::PayloadBlastCut::apply(FilteredAnitaEvent* fEv){  

  maxRatio = 0;
  int anitaVersion = AnitaVersion::get();

  for(int pol=0; pol < AnitaPol::kNotAPol; pol++){      
    for(int phi=0; phi < NUM_PHI; phi++){

      // skip ALFA channels/ALFA cross talk
      if(anitaVersion==3 && pol==AnitaPol::kVertical && phi==7){
	continue;
      }
      if(anitaVersion==3 && pol==AnitaPol::kHorizontal && phi==4){
	continue;
      }	  

      int topAnt = phi;
      int bottomAnt = 2*NUM_PHI + phi;
	
      TGraph* grTop = useful->getGraph(topAnt, (AnitaPol::AnitaPol_t) pol);
      TGraph* grBottom = useful->getGraph(bottomAnt, (AnitaPol::AnitaPol_t) pol);
      
      // const AnalysisWaveform* waveTop = fEv->getRawGraph(topAnt, (AnitaPol::AnitaPol_t) pol);
      // const TGraphAligned* grTop = waveTop->uneven();
      // const AnalysisWaveform* waveBottom = fEv->getRawGraph(bottomAnt, (AnitaPol::AnitaPol_t) pol);
      // const TGraphAligned* grBottom = waveBottom->uneven();

      Double_t maxYTop, maxXTop, minYTop, minXTop;
      RootTools::getLocalMaxToMin((const TGraph*)grTop, maxYTop, maxXTop, minYTop, minXTop);

      Double_t maxYBottom, maxXBottom, minYBottom, minXBottom;
      RootTools::getLocalMaxToMin((const TGraph*) grBottom, maxYBottom, maxXBottom, minYBottom, minXBottom);

      delete grTop;
      delete grBottom;
      
      double ratio = (maxYBottom - minYBottom)/(maxYTop - minYTop);

      if(ratio > maxRatio){
	maxRatio = ratio;
	maxRatioPhi = phi;
	maxRatioPol = (AnitaPol::AnitaPol_t) pol;
      }	
    }	
  }

  if(maxRatio > ratioCutHigh || maxRatio < ratioCutLow){
    // std::cerr << eventNumberDQ << "\t" << maxRatio << std::endl;
    // std::cerr << "failed!:" << std::endl;
    // std::cerr << "\t" << maxRatio << "\t" << maxRatioPhi << "\t" << maxRatioPol << std::endl;    
    eventPassesCut = false;
  }
  else{
    eventPassesCut = true;
  }
  if(sum!=NULL){
    if(eventPassesCut){
      sum->flags.isVarner = 0;
    }
    else{
      sum->flags.isVarner = 1;
    }
  }
}



Acclaim::NumPointsCut::NumPointsCut(){
  // This wasn't chosen particularly carefully
  // I just want to stop core dumps with old root when there aren't enough points for interpolation
  // which is < 5 points, so this is more than sufficient.
  numPointsCutLow = 200;
  description = "Checks there are a reasonable number of points in each waveform.";  
}


void Acclaim::NumPointsCut::apply(const UsefulAnitaEvent* useful, AnitaEventSummary* sum){
  eventPassesCut = true;
  for(int chanIndex=0; chanIndex < NUM_CHAN*NUM_SURF; chanIndex++){
    const int numPoints = useful->fNumPoints[chanIndex];
    if(numPoints < numPointsCutLow){
      eventPassesCut = false;
      break;
    }
  }

  if(sum){
    if(!eventPassesCut){
      // silly old flag name      
      sum->flags.isVarner2 = 1;
    }
    else{
      sum->flags.isVarner2 = 0;
    }
  }
}

