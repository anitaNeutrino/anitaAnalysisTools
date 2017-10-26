#include "QualityCut.h"
#include "UsefulAnitaEvent.h"

/** 
 * Work around dodgy MC samples with 260 points, the last few of which are nonsense
 * 
 * @param n is the length of the timeArray
 * @param x is a pointer to the timeArray
 * 
 * @return first point at which the time array is not monatonically increasing
 */
int quickFixMC(int n, const double* x){
  // quick hack to fix this version of Linda's MC
  double lastX = -DBL_MAX;
  int firstBadPoint = n;
  for(int i=0; i < n; i++){
    double thisX = x[i];
    if(thisX < lastX){
      firstBadPoint = i;
      break;
    }
    lastX = thisX;
  }
  return firstBadPoint;
}


ClassImp(Acclaim::QualityCut);
ClassImp(Acclaim::SurfSaturationCut);
ClassImp(Acclaim::PayloadBlastCut);
ClassImp(Acclaim::NumPointsCut);


/** 
 * @brief Applies all the event quality cuts in succession, this should be the primary interface
 *
 * Each cut sets a different flag inside the AnitaEventSummary,
 * however this function also sets the master flag "isGood" if the event passes all of the cuts.
 * 
 * @param usefulEvent is the event whose quality we which to check
 * @param sum is the AnitaEventSummary, the isGood flag is set if the event passes all quality cuts
 * 
 * @return returns true if the event passes all quality cuts, false otherwise.
 */
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



/** 
 * @brief Reads the flags set in AnitaEventSummary, returns true if the event passed all quality cuts
 * 
 * @param sum points to the AnitaEventSummary to evalutate
 * @param describe if true prints a message to stdout describing the reason the event failed a quality cut (default false)
 * 
 * @return returns true if the event passes all quality cuts, false otherwise.
 */
Bool_t Acclaim::QualityCut::passedAll(const AnitaEventSummary* sum, bool describe){
      
  Bool_t isGood = sum->flags.isGood;
  if(describe){
    const char* p = "passed";
    const char* f = "failed";
    std::cout << "Event Summary run " << sum->run << " eventNumber " << sum->eventNumber << " failed a quality cut: " << std::endl;
    std::cout << "SURF Saturation cut = " << (sum->flags.isVarner ? f : p) << ", ";
    std::cout << "Payload Blast cut = " << (sum->flags.isPayloadBlast ? f : p) << ", ";
    std::cout << "Num points cut = " << (sum->flags.isVarner2 ? f : p) << std::endl;
  }  
  return isGood;
}


/** 
 * @brief Constructor for the SurfSaturation cut, sets some hard coded default values
 *
 * If the values are changed please increment the ClassDef counter in the header file.
 */
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
    


/** 
 * @brief Applies the SurfSaturation cut
 *
 * If the event does not pass this quality cut the summary flag isVarner is set to true.
 * (This wasn't the original use of this flag in Abby's analysis, but I'm coopting the flag for this purpose)
 * 
 * @param useful is the event whose quality we wish to characterise
 * @param sum is the AnitaEventSummary
 */
void Acclaim::SurfSaturationCut::apply(const UsefulAnitaEvent* useful, AnitaEventSummary* sum){  

  maxVolts = 0;
  minVolts = 0;
  for(int pol=0; pol <AnitaPol::kNotAPol; pol++){
    for(int ant=0; ant < NUM_SEAVEYS; ant++){
      const int chanIndex = AnitaGeomTool::getChanIndexFromAntPol(ant, (AnitaPol::AnitaPol_t) pol);
      const double* volts = useful->fVolts[chanIndex];
      const int n = useful->fNumPoints[chanIndex] == NUM_SAMP ? quickFixMC(NUM_SAMP, useful->fTimes[chanIndex]) : useful->fNumPoints[chanIndex];

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
      sum->flags.isVarner = 0;
    }
    else{
      sum->flags.isVarner = 1;
    }
  }
}











/** 
 * Constructor for the payload blast cut, contains some hard coded numbers
 * 
 * If the values are changed please increment the ClassDef counter in the header file.
 */
Acclaim::PayloadBlastCut::PayloadBlastCut(){
  ratioCutHigh = 2.8;
  // ratioCutLow = -1; // no longer used
  maxRatio = 0;
  maxRatioPhi = -1;
  maxRatioPol = AnitaPol::kNotAPol;
  eventPassesCut = true;
  description = "Removes events with where any peak-to-peak in the top ring is significantly greater than the bottom ring.";
      
}






/** 
 * Applies the payload blast cut. 
 *
 * Sets the pre-existings flag isPayloadBlast to true if the event does not pass the cut
 * 
 * @param useful is the event whose quality we wish to characterise
 * @param sum is the AnitaEventSummary
 */
void Acclaim::PayloadBlastCut::apply(const UsefulAnitaEvent* useful, AnitaEventSummary* sum){

  maxRatio = 0;
  int anitaVersion = AnitaVersion::get();

  if(sum){
    for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
      sum->flags.maxBottomToTopRatio[polInd] = -9999;
      sum->flags.maxBottomToTopRatioSector[polInd] = -1;
      sum->flags.minBottomToTopRatio[polInd] = 9999;
      sum->flags.minBottomToTopRatioSector[polInd] = -1;
    }
  }

  AnitaPol::AnitaPol_t maxPol = AnitaPol::kNotAPol;
  AnitaRing::AnitaRing_t maxRing = AnitaRing::kNotARing;
  int maxPhi = -1;

  int nBottomOrMidGreaterThanTop = 0;

  double largestPeakToPeak = -DBL_MAX;
  double correspondingTopPeakToPeak = -DBL_MAX;

  for(int pol=0; pol < AnitaPol::kNotAPol; pol++){

    double largestPeakToPeakPol = -DBL_MAX; // stored per polarisation
    double correspondingTopPeakToPeakPol = -DBL_MAX;

    for(int phi=0; phi < NUM_PHI; phi++){

      // skip ALFA channel
      if(anitaVersion==3 && pol==AnitaPol::kHorizontal && phi==4){
	continue;
      }
      // ALFA cross talk
      if(anitaVersion==3 && pol==AnitaPol::kHorizontal && phi==12){
	continue;
      }
      // and channel with occasionally broken amplifer
      if(anitaVersion==3 && pol==AnitaPol::kVertical && phi==7){
	continue;
      }
      // ... luckily none of those are next to each other
      // and real blasts have enough power do many phi-sectors


      // first get the top ring
      TGraph* grTop = useful->getGraph(phi, (AnitaPol::AnitaPol_t) pol);
      const int nT = grTop->GetN() == NUM_SAMP ? quickFixMC(grTop->GetN(), grTop->GetX()) : grTop->GetN();
      Double_t maxY, maxX, minY, minX;
      RootTools::getLocalMaxToMinWithinLimits(grTop, maxY, maxX, minY, minX, grTop->GetX()[0], grTop->GetX()[nT-1]+1e-10);

      Double_t peakToPeakTop = (maxY - minY);
      delete grTop;


      // then get the mid/bottom rings
      for(int ring=AnitaRing::kMiddleRing; ring < AnitaRing::kNotARing; ring++){

	int ant = ring*NUM_PHI + phi;

	TGraph* gr = useful->getGraph(ant, (AnitaPol::AnitaPol_t) pol);
	const int n = gr->GetN() == NUM_SAMP ? quickFixMC(gr->GetN(), gr->GetX()) : gr->GetN();
	Double_t maxY, maxX, minY, minX;
      	RootTools::getLocalMaxToMinWithinLimits(gr, maxY, maxX, minY, minX, gr->GetX()[0], gr->GetX()[n-1]+1e-10);

	Double_t peakToPeak = (maxY - minY);
	if(peakToPeak > largestPeakToPeakPol){
	  maxPhi = phi;
	  maxRing = (AnitaRing::AnitaRing_t) ring;
	  largestPeakToPeakPol = peakToPeak;
	  correspondingTopPeakToPeakPol = peakToPeakTop;
	}

	if(peakToPeak > peakToPeakTop){
	  nBottomOrMidGreaterThanTop++;
	}

	delete gr;
      }
    }

    if(sum){
      sum->flags.maxBottomToTopRatio[pol] = largestPeakToPeakPol/correspondingTopPeakToPeakPol;
      sum->flags.maxBottomToTopRatioSector[pol] = maxPhi + maxRing*NUM_PHI;
    }

    if(largestPeakToPeakPol > largestPeakToPeak){
      maxPol = (AnitaPol::AnitaPol_t) pol;
      largestPeakToPeak = largestPeakToPeakPol;
      correspondingTopPeakToPeak = correspondingTopPeakToPeakPol;
    }
  }


  double maxRatio = largestPeakToPeak/correspondingTopPeakToPeak;
  if(maxRatio > ratioCutHigh || maxRatio < ratioCutLow){
    eventPassesCut = false;
  }
  else{
    eventPassesCut = true;
  }
  if(sum!=NULL){
    sum->flags.nSectorsWhereBottomExceedsTop = nBottomOrMidGreaterThanTop;
    if(eventPassesCut){
      sum->flags.isPayloadBlast = 0;
    }
    else{
      sum->flags.isPayloadBlast = 1;
    }
  }  
}





/** 
 * @brief Constructor for the number of points cut, contains some hard coded numbers
 *
 * If the values are changed please increment the ClassDef counter in the header file.
 */
Acclaim::NumPointsCut::NumPointsCut(){
  // This wasn't chosen particularly carefully
  // I just want to stop core dumps with old root when there aren't enough points for interpolation
  // which is < 5 points, so this is more than sufficient.
  numPointsCutLow = 200;
  description = "Checks there are a reasonable number of points in each waveform.";  
}





/** 
 * @brief Apply the number of points cut
 *
 * This cut was really only designed to prevent core dumps when the gsl akima interpolator than ROOT wraps fails,
 * The failure occurs when n < 5. 
 *  
 * @param useful is the event whose quality we wish to characterise
 * @param sum is the AnitaEventSummary
 */
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
