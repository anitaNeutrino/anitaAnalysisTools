#include "RootTools.h"
#include <algorithm>
#include <numeric>

//---------------------------------------------------------------------------------------------------------
/**
 * @brief Silly test function to generate a UsefulAnitaEvent
 *
 * All channels have gaussian distributed voltage values with a mean of 0, rms of 1, with dt = 1/2.6 ns
 * This isn't a good fake event but may be useful for testing things.
 *
 * @param eventNumber is used as the event number in the fake event and the seed for a TRandom3
 * @returns a pointer to a UsefulAnitaEvent
*/

UsefulAnitaEvent* Acclaim::RootTools::makeGaussianEvent(UInt_t eventNumber){

  TRandom3 randy(eventNumber);
  const double meanGauss = 0;
  const double sigmaGauss = 1;
  const double dt = 1./2.6;

  UsefulAnitaEvent* useful = new UsefulAnitaEvent();
  useful->eventNumber = eventNumber;
  for(int surf=0; surf < NUM_SURF; surf++){
    for(int chan = 0; chan < NUM_CHAN - 1; chan++){
      int chanIndex = surf*NUM_CHAN + chan;
      double t = 0;
      useful->fNumPoints[chanIndex] = NUM_SAMP;
      for(int samp=0; samp < NUM_SAMP; samp++){
	useful->fVolts[chanIndex][samp] = randy.Gaus(meanGauss, sigmaGauss);
	useful->fTimes[chanIndex][samp] = t;
	t += dt;
      }
    }
  }
  return useful;

}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Adds up all y-axis values of input TGraph.
 *
 * @param gr is the TGraph you want to sum the y-axis values of.
 * @returns sum is a Double_t.
*/
Double_t Acclaim::RootTools::getSumOfYVals(const TGraph* gr){
  Double_t sum = 0;
  for(int i=0; i<gr->GetN(); i++){
    sum += gr->GetY()[i];
  }
  return sum;
}



/**
 * @brief I use this function to define my min bias sample
 *
 * Used in several places, and I've changed it at least once, so it gets its own function.
 * Since the PPS2 was only switched on around MCM for ANITA-3 I'm not going to use it
 *
 * @param header is the RawAnitaHeader of the event in question.
 *
 * @return
 */
Int_t Acclaim::RootTools::isMinBiasSampleEvent(const RawAnitaHeader* header){
  // Not using PPS2 for ANITA-3 since it was only switched on around McM
  // so for the purposes of defining a thermal sample, it's essentially useless

  // PPS2 is the header->getTriggerBitG12(), which I'm not using
  Int_t isMinBias = header->getTriggerBitADU5() || header->getTriggerBitSoftExt();

  return isMinBias;
}



//---------------------------------------------------------------------------------------------------------
/**
 * @brief Print all the elements of a c-styles array to the screen for debugging.
 *
 * @param n is the length of the array.
 * @param array is the first element of the array.
 * @param delimiter is for decoration, ", " by default, change as you please.
 * @param start is for decoration, "{" by default, change as you please.
 * @param end is for decoration, "}\n" by default, change as you please.
*/
void Acclaim::RootTools::printArray(int n, double* array, TString delimiter, TString start ,TString end){
  std::cout << start.Data();
  for(int i=0; i<n; i++){
    std::cout << array[i];
    if(i<n-1) std::cout << delimiter.Data();
  }
  std::cout << end.Data();
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Print all the the y-axis values of a TGraph to the screen for debugging.
 *
 * @param gr is the TGraph you want information about.
 * @param delimiter is for decoration, ", " by default, change as you please.
 * @param start is for decoration, "{" by default, change as you please.
 * @param end is for decoration, "}\n" by default, change as you please.
 * @sa Acclaim::RootTools::printArray(int n, double* array, TString delimiter, TString start ,TString end)
*/
void Acclaim::RootTools::printYVals(const TGraph* gr, TString delimiter, TString start, TString end){
  printArray(gr->GetN(), gr->GetY(), delimiter, start, end);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Print all the the x-axis values of a TGraph to the screen for debugging.
 *
 * @param gr is the TGraph you want information about.
 * @param delimiter is for decoration, ", " by default, change as you please.
 * @param start is for decoration, "{" by default, change as you please.
 * @param end is for decoration, "}\n" by default, change as you please.
 * @sa Acclaim::RootTools::printArray(int n, double* array, TString delimiter, TString start ,TString end)
*/
void Acclaim::RootTools::printXVals(const TGraph* gr, TString delimiter, TString start, TString end){
  printArray(gr->GetN(), gr->GetX(), delimiter, start, end);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Makes a derivative TGraph, contains gr->GetN() - 1 points.
 *
 * @param gr is the TGraph from which to make the derivative
 * @returns Pointer to TGraph containing derivative of input gr.
*/
TGraph* Acclaim::RootTools::makeDerivativeTGraph(const TGraph* gr){
  TGraph* grDer = new TGraph();
  Int_t nDer = gr->GetN()-1;
  grDer->Set(nDer);
  for(int samp=0; samp<nDer; samp++){
    Double_t dy = gr->GetY()[samp+1]-gr->GetY()[samp];
    Double_t dx = gr->GetX()[samp+1]-gr->GetX()[samp];
    grDer->SetPoint(samp, gr->GetX()[samp], dy/dx);
  }
  return grDer;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Unwraps a circularly correlated graph so the -ve time lies behind the +ve time.
 *
 * @param gr is the wrapped correlation graph.
 * @returns Pointer to TGraph containing unwrapped correlation graph.
*/
TGraph* Acclaim::RootTools::makeUnwrappedCorrelationGraph(const TGraph* gr){
  /* Unwraps a circularly correlated array so the n > 2 points get -ve values */
  Int_t n = gr->GetN();
  Double_t dt = gr->GetX()[1] - gr->GetX()[0];
  TGraph* grUnwrapped = new TGraph();
  grUnwrapped->Set(n);
  Int_t samp2=0;
  for(Int_t samp=n/2; samp<n; samp++){
    grUnwrapped->SetPoint(samp2, (samp - n)*dt, gr->GetY()[samp]);
    samp2++;
  }
  for(Int_t samp=0; samp<n/2; samp++){
    grUnwrapped->SetPoint(samp2, samp*dt, gr->GetY()[samp]);
    samp2++;
  }

  return grUnwrapped;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Scan long TGraph and update the references to the maximum and minimum input values.
 *
 * @param gr is the TGraph you want information about.
 * @param maxY is a reference to a Double_t, updated with the maximum y-axis value of the TGraph.
 * @param minY is a reference to a Double_t, updated with the minimum y-axis value of the TGraph.
*/
void Acclaim::RootTools::getMaxMin(TGraph* gr, Double_t& maxY, Double_t& minY){

  Double_t maxX=0, minX=0;
  Acclaim::RootTools::getMaxMin(gr, maxY, maxX, minY, minX);
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Scan long TGraph and update the references to the maximum and minimum input values.
 *
 * @param gr is the TGraph you want information about.
 * @param maxY is a reference to a Double_t, updated with the maximum y-axis value of the TGraph.
 * @param maxX is a reference to a Double_t, updated with the x-axis value at the maxY of the TGraph.
 * @param minY is a reference to a Double_t, updated with the minimum y-axis value of the TGraph.
 * @param minX is a reference to a Double_t, updated with the x-axis value at the minY of the TGraph.
*/
void Acclaim::RootTools::getMaxMin(TGraph* gr, Double_t& maxY, Double_t& maxX,
			  Double_t& minY, Double_t& minX){

  maxY=gr->GetY()[0];
  maxX=gr->GetX()[0];
  minY=gr->GetY()[0];
  minX=gr->GetX()[0];

  for(int i=0; i<gr->GetN(); i++){
    if(gr->GetY()[i] > maxY){
      maxY = gr->GetY()[i];
      maxX = gr->GetX()[i];
    }
    if(gr->GetY()[i] < minY){
      minY = gr->GetY()[i];
      minX = gr->GetX()[i];
    }
  }
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Scan long TGraph and update the references to the maximum and minimum input values.
 *
 * @param gr is the TGraph you want information about.
 * @param maxY is a reference to a Double_t, updated with the maximum y-axis value of the TGraph.
 * @param maxX is a reference to a Double_t, updated with the x-axis value at the maxY of the TGraph.
 * @param minY is a reference to a Double_t, updated with the minimum y-axis value of the TGraph.
 * @param minX is a reference to a Double_t, updated with the x-axis value at the minY of the TGraph.
 * @param lowerLimit is a Double_t containing the lower value of the x-axis to search from.
 * @param upperLimit is a Double_t containing the upper value of the x-axis to search to.
*/
void Acclaim::RootTools::getMaxMinWithinLimits(TGraph* gr, Double_t& maxY, Double_t& maxX,
                                               Double_t& minY, Double_t& minX,
                                               Double_t lowerLimit, Double_t upperLimit){

  // Macro values from <cfloat> header
  Double_t minPoss = -DBL_MAX;
  Double_t maxPoss = DBL_MAX;
  maxY=minPoss;
  maxX=minPoss;
  minY=maxPoss;
  minX=maxPoss;

  for(int i=0; i<gr->GetN(); i++){
    if(gr->GetX()[i] >= lowerLimit && gr->GetX()[i] <=upperLimit){
      if(gr->GetY()[i] > maxY){
	maxY = gr->GetY()[i];
	maxX = gr->GetX()[i];
      }
      if(gr->GetY()[i] < minY){
	minY = gr->GetY()[i];
	minX = gr->GetX()[i];
      }
    }
  }

  if(maxY==minPoss){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << "!" << std::endl;
    std::cerr << "Unable to find value in " << gr->GetName() << " > " << minPoss << std::endl;
  }
  if(minY==maxPoss){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << "!" << std::endl;
    std::cerr << "Unable to find value in " << gr->GetName() << " < " << maxPoss << std::endl;
  }
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Subtract a constant value from the y-axis values at each point.
 *
 * @param gr is the TGraph you want to manipulate.
 * @param offset is the value to subtract from each point.
*/
void Acclaim::RootTools::subtractOffset(TGraph* gr, Double_t offset){
  for(int i=0; i<gr->GetN(); i++){
    gr->GetY()[i] -= offset;
  }
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Do angle1-angle2 (in degrees) and +- 360 such that the result lies in -180 < deltaAngle < 180.
 *
 * @param angle1 is the first angle
 * @param angle2 is the second angle
 * @returns deltaAngle is the difference in the range -180 < deltaAngle < 180
 *
 * If it looks like there's an insane angle being passed in will return -9999.
*/
Double_t Acclaim::RootTools::getDeltaAngleDeg(Double_t angle1, Double_t angle2){

  Double_t deltaAngle = angle1 - angle2;
  Int_t loopCount = 0;
  const Int_t maxLoopCount = 5;
  if(deltaAngle > 180){
    while(deltaAngle >= 180){
      deltaAngle -= 360;
      loopCount++;
      if(loopCount >= maxLoopCount){
	deltaAngle = -9999;
	break;
      }
    }
  }
  else{
    while(deltaAngle < -180){
      deltaAngle += 360;
      loopCount++;
      if(loopCount >= maxLoopCount){
	deltaAngle = -9999;
	break;
      }
    }
  }

  return deltaAngle;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Modify gr to set mean=0 and rms=1.
 *
 * @param gr is the TGraph you want to manipulate.
*/
void Acclaim::RootTools::normalize(TGraph* gr){
  double mean, rms;
  normalize(gr, mean, rms);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Modify gr to set mean=0 and rms=1.
 *
 * @param gr is the TGraph you want to manipulate.
 * @param mean is a reference to a Double_t, updated with the mean value of gr before modification.
 * @param rms is a reference to a Double_t, updated with the rms value of gr before modification.
*/
void Acclaim::RootTools::normalize(TGraph* gr, Double_t& mean, Double_t& rms){
  Acclaim::RootTools::getMeanAndRms(gr, mean, rms);
  if(rms>0){ // Don't make any NaNs
    for(int i=0; i<gr->GetN(); i++){
      gr->GetY()[i] -= mean;
      gr->GetY()[i] /= rms;
    }
  }
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Create normalized copy of gr, where mean=0 and rms=1 (input unmodified).
 *
 * @param gr is the TGraph you want to manipulate.
*/
TGraph* Acclaim::RootTools::makeNormalized(TGraph* gr){
  /* Copies the TGraph and normalizes that */
  TGraph* grCopy = (TGraph*) gr->Clone();
  normalize(grCopy);
  return grCopy;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Create normalized copy of gr, where mean=0 and rms=1 (input unmodified).
 *
 * @param gr is the TGraph you want to manipulate.
 * @param mean is a reference to a Double_t, updated with the mean value of gr.
 * @param rms is a reference to a Double_t, updated with the rms value of gr.
*/
TGraph* Acclaim::RootTools::makeNormalized(TGraph* gr, Double_t& mean, Double_t& rms){
  /* Copies the TGraph and normalizes that */
  TGraph* grCopy = (TGraph*) gr->Clone();
  normalize(grCopy, mean, rms);
  return grCopy;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Get mean and rms of gr
 *
 * @param gr is the TGraph you want information about.
 * @param mean is a reference to a Double_t, updated with the mean value of gr before modification.
 * @param rms is a reference to a Double_t, updated with the rms value of gr before modification.
*/
void Acclaim::RootTools::getMeanAndRms(TGraph* gr, Double_t& mean, Double_t& rms){
  Double_t sum = 0;
  Double_t square = 0;
  for(int i=0; i<gr->GetN(); i++){
    sum += gr->GetY()[i];
    square += gr->GetY()[i]*gr->GetY()[i];
  }
  mean = sum/gr->GetN();
  rms = TMath::Sqrt(square/gr->GetN() - mean*mean);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Find indices where input is not a number (for debugging).
 *
 * @param len is the length of arr.
 * @param arr is a pointer to an array of doubles of length len.
 * @returns std::vector<Int_t> containing the indices of the TGraph where is value is not a number.
*/
Int_t Acclaim::RootTools::getIndexOfMaximum(Int_t len, Double_t* arr){
  Double_t max=-DBL_MAX;
  Int_t maxIndex = -1;
  for(Int_t i=0; i < len; ++i){
    if(arr[i] > max){
      max = arr[i];
      maxIndex = i;
    }
  }
  if(maxIndex==-1){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << "!" << std::endl;
    std::cerr << "Could not find an index > " << max << std::endl;
  }

  return maxIndex;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Find indices where input is not a number (for debugging)
 *
 * @param gr is the TGraph you want information about.
 * @returns std::vector<Int_t> containing the indices of the TGraph where is value is not a number
*/
std::vector<Int_t> Acclaim::RootTools::getIndicesOfNans(TGraph* gr){
  std::vector<Int_t> nanIndices;
  for(Int_t i=0; i<gr->GetN(); i++){
    if(TMath::IsNaN(gr->GetY()[i])){
      nanIndices.push_back(i);
    }
  }
  return nanIndices;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Print summary information about gr: name, title, number of points, all x and y values.
 *
 * @param gr is the TGraph you want information about.
*/
void Acclaim::RootTools::printTGraphInfo(const TGraph* gr){
  std::cout << "******************************************************************************" << std::endl;
  std::cout << "Acclaim::RootTools::printTGraphValues(TGraph* gr = " << gr << "):" << std::endl;
  std::cout << "Name: " << gr->GetName() << std::endl;
  std::cout << "Title: " << gr->GetTitle() << std::endl;
  std::cout << "Num points: " << gr->GetN() << std::endl;
  std::cout << "Xvals: ";
  for(Int_t i=0; i<gr->GetN(); i++){
    std::cout << gr->GetX()[i];
    if(i<gr->GetN()-1){
      std::cout << ", ";
    }
  }
  std::cout << std::endl;
  std::cout << "Yvals: ";
  for(Int_t i=0; i<gr->GetN(); i++){
    std::cout << gr->GetY()[i];
    if(i<gr->GetN()-1){
      std::cout << ", ";
    }
  }
  std::cout << std::endl;
  std::cout << "******************************************************************************" << std::endl;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Plot a time (or any TTree branch) ordered TGraph of some TTree branch.
 *
 * @param tree is a pointer to a TTree containing the variables to be plotted.
 * @param drawText is passed to TTree->Draw() in order to plot the data points (e.g. "eventNumber:heading") for a plot of heading as a function of event number.
 * @param cutString is passed to TTree->Draw() for selection (e.g. "eventNumber > 30000000")
 * @param wrapValue e.g. with wrapValue=360, if neighbouring points go 359->1 degrees, then draw them as 359->361. i.e. remove discontinuities.
*/
TGraph* Acclaim::RootTools::makeSortedTGraph(TTree* tree, TString drawText, TString cutString, Double_t wrapValue){

  // Draw
  const Int_t nEntries = tree->Draw(drawText, cutString, "goff"); // "goff" means graphics off

  // Sort
  std::vector<Int_t> sortedIndices(nEntries);
  TMath::Sort(nEntries, tree->GetV2(), &sortedIndices.front(), kFALSE);
  std::vector<Double_t> newX(nEntries);
  std::vector<Double_t> newY(nEntries);

  for(int i=0; i<nEntries; i++){
    newX.at(i) = tree->GetV2()[sortedIndices.at(i)];
    newY.at(i) = tree->GetV1()[sortedIndices.at(i)];
  }

  // Unwrap data here so interpolation works smoothly
  // will unwrap when getting entries from graph
  if(wrapValue != 0){
    for(int i=1; i<nEntries; i++){
      // If y[i] >> y[i-1] => then we went below zero to the wrap value
      // can only modify y[i], so we should subtract wrapValue enough times
      // that the graph becomes smooth
      while (newY.at(i) - newY.at(i-1) > wrapValue/2){
	newY.at(i) -= wrapValue;
      }

      // If y[i] << y[i-1] => then we went add zero to the wrap value
      // can only modify y[i], so we should add wrapValue enough times
      // that the graph becomes smooth
      while (newY.at(i) - newY.at(i-1) < -wrapValue/2){
	newY.at(i) += wrapValue;
      }

    }
  }

  // sorted TGraph
  TGraph* gr(new TGraph(nEntries,&newX.front(), &newY.front()));
  TString title = cutString.Length() == 0 ? drawText : drawText + ", " + cutString;
  gr->SetTitle(title);
  TObjArray* tokens = drawText.Tokenize(":");
  gr->GetXaxis()->SetTitle(((TObjString*) tokens->At(1))->String());
  gr->GetYaxis()->SetTitle(((TObjString*) tokens->At(0))->String());
  delete tokens;

  return gr;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Put Z-axis values of a TH2D histogram into a TH1D histogram.
 *
 * @param h2 is a pointer to the TH2D you want to examine.
 * @param hName is a the name of the new histogram.
 * @param nBins is the number of bins for the new histogram.
 * @param xMin is the x-axis minimum for the new histogram.
 * @param xMax is the x-axis maximum for the new histogram.
*/
TH1D* Acclaim::RootTools::plotsZaxisDist(TH2* h2, TString hName, Int_t nBins, Double_t xMin, Double_t xMax){
  TH1D* h = new TH1D(hName, hName, nBins, xMin, xMax);
  for(int xBin=1; xBin<=h2->GetNbinsX(); xBin++){
    for(int yBin=1; yBin<=h2->GetNbinsY(); yBin++){
      Double_t val = h2->GetBinContent(xBin, yBin);
      h->Fill(val);
    }
  }
  return h;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Add a bunch of zeros to the end of a TGraph.
 *
 * @param gr is the TGraph you want to pad.
 * @param newLen is the desired length of gr.
 * @param dt
*/
void Acclaim::RootTools::zeroPadTGraph(TGraph* gr, Int_t newLen, Double_t dt){
  Int_t oldLen = gr->GetN();
  gr->Set(newLen);
  if(dt >= 0){
    Double_t x0 = gr->GetX()[oldLen-1];
    for(Int_t samp=oldLen; samp<newLen; samp++){
      gr->GetX()[samp] = x0 + samp*dt;
    }
  }
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Make interpolated TGraph using linear interpolation between points.
 *
 * @param grIn is the TGraph you want to interpolate.
 * @param dt is the time (in x-axis units of grIn) between points in grOut.
 * @returns grOut, a new TGraph.
*/

TGraph* Acclaim::RootTools::makeLinearlyInterpolatedGraph(TGraph* grIn, Double_t dt){

  Int_t nIn = grIn->GetN();
  Int_t newPoints = Int_t((grIn->GetX()[nIn-1] - grIn->GetX()[0])/dt) + 1;
  std::vector<Double_t> newTimes(newPoints);
  std::vector<Double_t> newVolts(newPoints);

  Double_t time = grIn->GetX()[0];
  Double_t lastTime = grIn->GetX()[nIn-1];

  Int_t sampOut = 0; // Will iterate through
  Int_t sampIn = 0;

  Double_t y1 = grIn->GetY()[sampIn];
  Double_t x1 = grIn->GetX()[sampIn];
  Double_t x2 = grIn->GetX()[sampIn+1];
  Double_t m = (grIn->GetY()[sampIn+1]-y1)/(grIn->GetX()[sampIn+1]-x1);

  while(time < lastTime){
    if(time > x2){

      Int_t countLoop = 0;
      while(time > x2){
	if(sampIn < (nIn - 1)){ // if last point is equal don't need to increase?
	  sampIn++;
	}
	x2 = grIn->GetX()[sampIn+1];
	countLoop++;
	if(countLoop==1000){
	  std::cerr << "You can't do logic very well" << std::endl;
	  std::cerr << time << "\t" << x1 << "\t" << x2 << "\t" << nIn << "\t" << sampIn << std::endl;
	  exit(0);
	}
      }
      x1 = grIn->GetX()[sampIn];
      y1 = grIn->GetY()[sampIn];

      m = (grIn->GetY()[sampIn+1]-y1)/(x2-x1);

    }

    // std::cout << time << "\t" << x1 << "\t" << x2;

    // if(time < x1 || time > x2 ){
    //   std::cout << "!!!!!!!!!!!!!!!!!!!!!!!";
    // }
    // std::cout << std::endl;

    Double_t newY = y1 + m * (time - x1);

    newTimes.at(sampOut) = time;
    newVolts.at(sampOut) = newY;
    sampOut++;
    time += dt;
  }

  TGraph* grOut = new TGraph(newTimes.size(), &newTimes[0], &newVolts[0]);
  // std::cout << grOut->GetN() - newPoints << " " << grIn->GetN() << " " << sampIn << " "
  // 	    << grIn->GetX()[0] << " " << grIn->GetX()[nIn-1] << " " << dt << std::endl;




  return grOut;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Wrapper function for ROOT's interpolator, can zero pad the front to start from a particular time.
 *
 * @param grIn points to the TGraph containing the waveform to interpolate
 * @param startTime is the start time for interpolation: zero pads if this is earlier than the TGraph start time.
 * @param dt is the time between samples
 * @param nSamp is the time number of samples
*/
TGraph* Acclaim::RootTools::interpolateWithStartTime(TGraph* grIn, Double_t startTime, Double_t dt, Int_t nSamp){

  std::vector<Double_t> newTimes = std::vector<Double_t>(nSamp, 0);
  std::vector<Double_t> newVolts = std::vector<Double_t>(nSamp, 0);
  Double_t thisStartTime = grIn->GetX()[0];
  Double_t lastTime = grIn->GetX()[grIn->GetN()-1];


  // Quantizes the start and end times so data poInt_ts lie at Int_teger multiples of nominal sampling
  // startTime = correlationDeltaT*TMath::Nint(startTime/correlationDeltaT + 0.5);
  // lastTime = correlationDeltaT*TMath::Nint(lastTime/correlationDeltaT - 0.5);
  startTime = dt*TMath::Nint(startTime/dt + 0.5);
  lastTime = dt*TMath::Nint(lastTime/dt - 0.5);

   //ROOT Int_terpolator object constructor takes std::vector objects
  std::vector<Double_t> tVec(grIn->GetX(), grIn->GetX() + grIn->GetN());
  std::vector<Double_t> vVec(grIn->GetY(), grIn->GetY() + grIn->GetN());

  // This is ROOT's Int_terpolator object
  ROOT::Math::Interpolator chanInterp(tVec,vVec,ROOT::Math::Interpolation::kAKIMA);

  // Put new data Int_to arrays
  Double_t time = startTime;
  for(Int_t samp = 0; samp < nSamp; samp++){
    newTimes.at(samp) = time;
    if(time >= thisStartTime && time <= lastTime){
      newVolts.at(samp) = chanInterp.Eval(time);
    }
    else{
      newVolts.at(samp) = 0;
    }
    time += dt;
  }
  return new TGraph(nSamp, &newTimes[0], &newVolts[0]);

}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Draws an array of histograms with a rainbow on a single TCanvas
 *
 * @param hs is a pointer to an array of pointers to the histograms.
 * @param numHists is the number of histograms.
 * @param can is the TCanvas to draw on, if one is not given a new one is created.
 * @param colWeights is a pointer to a Double_t array, which will weight to colors drawn.
 * @returns the TCanvas the histograms are drawn on.
*/

TCanvas* Acclaim::RootTools::drawArrayOfHistosPrettily(TH1D* hs[], Int_t numHists, TCanvas* can, Double_t* colWeights){

  if(can==NULL){
    can = new TCanvas();
  }
  can->cd();

  Double_t minCol = 0;
  Double_t maxCol = 255;
  if(colWeights!=NULL){
    maxCol = TMath::MaxElement(numHists, colWeights);
    minCol = TMath::MinElement(numHists, colWeights);
  }

  Int_t numDrawn = 0;
  Int_t first = 0;
  Double_t max = -1;
  for(Int_t histInd=0; histInd < numHists; histInd++){
    if(hs[histInd] != NULL){
      if(numDrawn==0){
	first = histInd;
      }
      TString opt = numDrawn == 0 ? "" : "same";
      if(colWeights!=NULL){
	hs[histInd]->SetLineColor(gStyle->GetColorPalette(Int_t(250.*(colWeights[histInd] - minCol)/(maxCol-minCol))));
	// std::cout << colWeights[histInd] << "\t" << maxCol << "\t" << minCol << "\t"
	// 	  << ((colWeights[histInd] - minCol)/(maxCol-minCol)) << "\t"
	// 	  << Int_t(250.*(colWeights[histInd] - minCol)/(maxCol-minCol))
	// 	  << std::endl;
      }
      else{
	hs[histInd]->SetLineColor(gStyle->GetColorPalette(histInd*Int_t(255./numHists)));
      }
      hs[histInd]->Draw(opt);
      if(hs[histInd]->GetMaximum() > max){
	max = hs[histInd]->GetMaximum()*1.1;
	hs[first]->SetMaximum(max);
      }
      numDrawn++;
    }
  }
  return can;
}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Adds a set of offsets to the x-axis values of each respective TGraph
 *
 * @param numGrs is the length of both arrays.
 * @param grs is a vector of pointers to TGraphs.
 * @param offsets a vector of Double_ts to add to each x-axis point of the respective TGraph.

  Assumes that both grs and offsets point to arrays of length numGrs
*/
void Acclaim::RootTools::offsetTGraphXAxes(Int_t numGrs, TGraph* grs[], Double_t offsets[]){
  for(Int_t grInd=0; grInd < numGrs; grInd++){
    for(Int_t samp=0; samp < grs[grInd]->GetN(); samp++){
      grs[grInd]->GetX()[samp] += offsets[grInd];
    }
  }
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Multiples a set of factors to the y-axis values of each respective TGraph
 *
 * @param numGrs is the length of both arrays.
 * @param grs is an array of pointers to TGraphs.
 * @param factors an array of Double_ts to add to each x-axis point of the respective TGraph.

  Assumes that both grs and factors point to arrays of length numGrs
*/
void Acclaim::RootTools::multiplyTGraphYAxes(Int_t numGrs, TGraph* grs[], Double_t factors[]){
  for(Int_t grInd=0; grInd < numGrs; grInd++){
    for(Int_t samp=0; samp < grs[grInd]->GetN(); samp++){
      grs[grInd]->GetY()[samp] *= factors[grInd];
    }
  }
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Draws an array of TGraphs with a rainbow on a single TCanvas
 *
 * @param grs is a pointer to an array of pointers to the TGraphs.
 * @param numGrs is the number of TGraphs.
 * @param can is the TCanvas to draw on, if one is not given a new one is created.
 * @param drawOpt is the draw option, select "p" for markers and "l" for lines.
 * @param colWeights is a pointer to an array of doubles to be used as weights for colors.
 * @returns the TCanvas the TGraphs are drawn on.
*/
TCanvas* Acclaim::RootTools::drawArrayOfTGraphsPrettily(TGraph* grs[], Int_t numGrs,
					       TString drawOpt, TCanvas* can, Double_t* colWeights){

  if(can==NULL){
    can = new TCanvas();
  }
  can->cd();

  Double_t colFactor = 1;
  if(colWeights!=NULL){
    Double_t maxCol = TMath::MaxElement(numGrs, colWeights);
    Double_t minCol = TMath::MinElement(numGrs, colWeights);
    colFactor = 254./(maxCol-minCol);
  }

  Double_t max = -DBL_MAX;
  Double_t min = DBL_MAX;
  Int_t numDrawn = 0;
  Int_t first = 0;
  for(Int_t grInd=0; grInd < numGrs; grInd++){
    if(grs[grInd]!=NULL){
      if(grs[grInd]->GetN() > 0){
	TString opt = numDrawn == 0 ? "a" + drawOpt : drawOpt + "same";
	if(numDrawn==0){
	  first = grInd;
	}

	if(colWeights!=NULL){
	  grs[grInd]->SetLineColor(gStyle->GetColorPalette(Color_t(colWeights[grInd]*colFactor)));
	  grs[grInd]->SetMarkerColor(gStyle->GetColorPalette(Color_t(colWeights[grInd]*colFactor)));
	  grs[grInd]->SetMarkerColor(gStyle->GetColorPalette(Color_t(colWeights[grInd]*colFactor)));
	}
	else{
	  grs[grInd]->SetLineColor(gStyle->GetColorPalette(grInd*Int_t(254./(numGrs-1))));
	  grs[grInd]->SetMarkerColor(gStyle->GetColorPalette(grInd*Int_t(254./(numGrs-1))));
	  grs[grInd]->SetMarkerColor(gStyle->GetColorPalette(grInd*Int_t(254./(numGrs-1))));
	}

	grs[grInd]->Draw(opt);

	Double_t thisMax = TMath::MaxElement(grs[grInd]->GetN(), grs[grInd]->GetY());
	if(thisMax > max){
	  Double_t fact = thisMax > 0 ? 1.1 : 0.9;
	  max = thisMax*fact;
	  grs[first]->SetMaximum(max);
	}

	Double_t thisMin = TMath::MinElement(grs[grInd]->GetN(), grs[grInd]->GetY());
	if(thisMin < min){
	  Double_t fact = thisMin < 0 ? 1.1 : 0.9;
	  min = thisMin*fact;
	  grs[first]->SetMinimum(min);
	}
	numDrawn++;
      }
    }
  }

  return can;

}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Generates a TLegend from input arrays of TGraphs and titles (TStrings)
 *
 * @param grs is a pointer to an array of pointers to the TGraphs.
 * @param numGrs is the number of TGraphs.
 * @param titles is an array of TStrings containing the titles
 * @param opt is the legend option ("l" by default)
 * @param minX is the Pad coordinates, 0.8 by default.
 * @param minY is the Pad coordinates, 0.8 by default.
 * @param maxX is the Pad coordinates, 1 by default.
 * @param maxY is the Pad coordinates, 1 by default.
 * @returns the newly generated TLegend
*/

TLegend* Acclaim::RootTools::makeLegend(TGraph* grs[], Int_t numGrs,
			       TString titles[], TString opt,
			       Double_t minX, Double_t minY,
			       Double_t maxX, Double_t maxY){

  TLegend* l = new TLegend(minX, minY, maxX, maxY);
  for(Int_t grInd=0; grInd < numGrs; grInd++){
    if(grs[grInd]!=NULL){
      l->AddEntry(grs[grInd], titles[grInd], opt);
    }
  }
  return l;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Generates a TLegend from input arrays of TGraphs and titles (TStrings).
 *
 * @param hs is a pointer to an array of pointers to the histograms (TH1Ds).
 * @param numHists is the number of histgorams in the array hs.
 * @param titles is an array of TStrings containing the titles
 * @param opt is the legend option ("l" by default)
 * @param minX is the Pad coordinates, 0.8 by default.
 * @param minY is the Pad coordinates, 0.8 by default.
 * @param maxX is the Pad coordinates, 1 by default.
 * @param maxY is the Pad coordinates, 1 by default.
 * @returns the newly generated TLegend
*/
TLegend* Acclaim::RootTools::makeLegend(TH1D* hs[], Int_t numHists,
			       TString titles[], TString opt,
			       Double_t minX, Double_t minY,
			       Double_t maxX, Double_t maxY){

  TLegend* l = new TLegend(minX, minY, maxX, maxY);
  for(Int_t histInd=0; histInd < numHists; histInd++){
    if(hs[histInd]!=NULL){
      l->AddEntry(hs[histInd], titles[histInd], opt);
    }
  }
  return l;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Sets the name of a TGraph to name and then writes it to the current file.
 *
 * @param gr is a pointer to the TGraph you want to save.
 * @param name is the name you want to give the TGraph before saving it.
*/
void Acclaim::RootTools::writeTGraph(TGraph* gr, TString name){
  // Taking laziness to a new level... replace two lines with a one line function.
  gr->SetName(name);
  gr->Write();
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Updates input based on absolute largest local maximum to local minimum difference. For pulse finding.
 *
 * @param gr is the TGraph you want information about.
 * @param maxY is a reference to a Double_t, updated with the local maximum y-axis value of the TGraph.
 * @param maxX is a reference to a Double_t, updated with the x-axis value at the maxY of the TGraph.
 * @param minY is a reference to a Double_t, updated with the local minimum y-axis value of the TGraph.
 * @param minX is a reference to a Double_t, updated with the x-axis value at the minY of the TGraph.
*/

void Acclaim::RootTools::getLocalMaxToMin(const TGraph* gr,
				 Double_t& maxY, Double_t& maxX,
				 Double_t& minY, Double_t& minX){
  getLocalMaxToMinWithinLimits(gr, maxY, maxX, minY, minX, gr->GetX()[0], gr->GetX()[gr->GetN()-1]);

}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Updates input based on absolute largest local maximum to local minimum difference. For pulse finding.
 *
 * @param gr is the TGraph you want information about.
 * @param maxY is a reference to a Double_t, updated with the local maximum y-axis value of the TGraph.
 * @param maxX is a reference to a Double_t, updated with the x-axis value at the maxY of the TGraph.
 * @param minY is a reference to a Double_t, updated with the local minimum y-axis value of the TGraph.
 * @param minX is a reference to a Double_t, updated with the x-axis value at the minY of the TGraph.
 * @param lowerLimit is a Double_t containing the lower value of the x-axis to search from.
 * @param upperLimit is a Double_t containing the upper value of the x-axis to search to.
*/

void Acclaim::RootTools::getLocalMaxToMinWithinLimits(const TGraph* gr,
					     Double_t& maxY, Double_t& maxX,
					     Double_t& minY, Double_t& minX,
					     Double_t lowerLimit, Double_t upperLimit){


  // assumes unique value at local minima and maxima, e.g. 0.1, 0.2, 0.1 for maxima.
  // if it goes 0.1, 0.2, 0.2, 0.1 then this will need to get more complicated
  std::vector<Int_t> extremaSamps;

  for(Int_t samp=1; samp < gr->GetN()-1; samp++){
    if(gr->GetX()[samp] >= lowerLimit && gr->GetX()[samp] <=upperLimit){
      Double_t y0 = gr->GetY()[samp-1];
      Double_t y1 = gr->GetY()[samp];
      Double_t y2 = gr->GetY()[samp+1];


      if(y0 > y1 && y1 < y2){
	// Is a local minimum
	extremaSamps.push_back(samp);
      }
      else if(y0 < y1 && y1 > y2){
	// Is a local maximum
	extremaSamps.push_back(samp);
      }
    }
  }

  Double_t maxOffset = -DBL_MAX;
  Int_t maxSampInd = -1;
  for(Int_t sampInd = 0; sampInd < ((Int_t)extremaSamps.size())-1; sampInd++){
    Double_t offset = TMath::Abs(gr->GetY()[extremaSamps.at(sampInd)] - gr->GetY()[extremaSamps.at(sampInd+1)]);
    if(offset > maxOffset){
      maxOffset = offset;
      maxSampInd = sampInd;
    }
  }

  if(maxSampInd > -1){
    if(gr->GetY()[extremaSamps.at(maxSampInd)] > gr->GetY()[extremaSamps.at(maxSampInd+1)]){
      maxY = gr->GetY()[extremaSamps.at(maxSampInd)];
      maxX = gr->GetX()[extremaSamps.at(maxSampInd)];

      minY = gr->GetY()[extremaSamps.at(maxSampInd+1)];
      minX = gr->GetX()[extremaSamps.at(maxSampInd+1)];
    }
    else{
      maxY = gr->GetY()[extremaSamps.at(maxSampInd+1)];
      maxX = gr->GetX()[extremaSamps.at(maxSampInd+1)];

      minY = gr->GetY()[extremaSamps.at(maxSampInd)];
      minX = gr->GetX()[extremaSamps.at(maxSampInd)];
    }
  }
}



//---------------------------------------------------------------------------------------------------------
/**
 * @brief My function to save a TCanvas for talks, and save an editable version.
 *
 * @param c is a pointer to the TCanvas
 * @param fileName is the file name. Note that suffixes are added in the function!
*/
void Acclaim::RootTools::saveCanvas(TCanvas* c, TString fileName){

  std::cout << "Saving this canvas as an .eps, .png, and .C file..." << std::endl;
  TString fName = fileName + ".eps";
  c->SaveAs(fName);
  fName = fileName + ".png";
  c->SaveAs(fName);
  fName = fileName + ".C";
  c->SaveAs(fName);
  std::cout << "...Complete!" << std::endl;
  c->Update();
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Assumes a nice distribution with a single peak, finds the full width half max.
 *
 * @param h is a histogram
 * @returns the full width half max
*/
Double_t Acclaim::RootTools::getFullWidthHalfMax(TH1D* h){

  Int_t maxBin = Acclaim::RootTools::getPeakBinOfHistogram(h);
  Int_t n = h->GetNbinsX();
  Double_t max = h->GetBinContent(maxBin);


  Double_t halfMaxPos = 0;
  for(Int_t bin=maxBin; bin<=n; bin++){
    Double_t val = h->GetBinContent(bin);
    if(val < max/2){
      halfMaxPos = h->GetBinLowEdge(bin);
      break;
    }
  }

  Double_t halfMaxNeg = 0;
  for(Int_t bin=maxBin; bin>=1; bin--){
    Double_t val = h->GetBinContent(bin);
    if(val < max/2){
      halfMaxNeg = h->GetBinLowEdge(bin);
      break;
    }
  }

  return halfMaxPos - halfMaxNeg;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Finds the bin containing the maximum value of a TH1D
 *
 * @param h is a histogram
 * @returns the peak bin (in ROOT bin counting starts at 1)
*/
Int_t Acclaim::RootTools::getPeakBinOfHistogram(TH1D* h){

  Int_t n = h->GetNbinsX();
  Double_t max = -DBL_MAX;
  Int_t maxBin = 0;
  for(Int_t bin=1; bin<n; bin++){
    Double_t val = h->GetBinContent(bin);
    if(val > max){
      max = val;
      maxBin = bin;
    }
  }
  return maxBin;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Finds the bin containing the maximum value of a TH2D
 *
 * @param hist is a histogram
 * @param binx is a reference to the x-axis bin number (in ROOT bin counting starts at 1)
 * @param biny is a reference to the y-axis bin number (in ROOT bin counting starts at 1)
 * @returns the value of the histogram peak bin
*/
Double_t Acclaim::RootTools::getPeakBinOfHistogram(TH2D* hist, Int_t& binx, Int_t& biny){

  Int_t nx = hist->GetNbinsX();
  Int_t ny = hist->GetNbinsY();
  Double_t histPeak = -DBL_MAX;
  for(Int_t by = 1; by<=ny; by++){
    for(Int_t bx = 1; bx<=nx; bx++){
      Double_t val = hist->GetBinContent(bx, by);
      if(val > histPeak){
	histPeak = val;
	binx = bx;
	biny = by;
      }
    }
  }

  return histPeak;

}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Finds the bin containing the maximum value of a TH1D
 *
 * @param h is a histogram
 * @returns the peak bin (in ROOT bin counting starts at 1)
*/
Double_t Acclaim::RootTools::getLowBinEdgeOfHistogramPeak(TH1D* h){

  Int_t peakBin = getPeakBinOfHistogram(h);
  Double_t lowBinEdgeOfPeak = h->GetXaxis()->GetBinLowEdge(peakBin);
  return lowBinEdgeOfPeak;
}






//---------------------------------------------------------------------------------------------------------
/**
 * @brief For nice plotting of 2D dists about 0, makes max = -min.
 *
 * @param h is a 2D histogram
*/
void Acclaim::RootTools::makeZaxisScaleEqualAboutZero(TH2D* h){

  Double_t min = h->GetMinimum();
  Double_t max = h->GetMaximum();

  if(TMath::Abs(min) > max){
    h->SetMaximum(-min);
  }
  else{
    h->SetMinimum(-max);
  }
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief For decoding bit masks
 *
 * @param bitIndex is the bit to get (count from 0)
 * @param bitMask is the bitMask
*/
Int_t Acclaim::RootTools::getBit(UInt_t bitIndex, UInt_t bitMask){
  return ((bitMask >> bitIndex) & 1);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief For counting how many bits are set to one in a bitmask.
 *
 * @param numBitsToCheck is the number of little endian bits to check.
 * @param bitMask is the object to check the bits of.
*/
Int_t Acclaim::RootTools::getNumBitsSet(Int_t numBitsToCheck, UInt_t bitMask){
  Int_t numBitsHigh = 0;
  for(Int_t bitInd=0; bitInd < numBitsToCheck; bitInd++){
    if(getBit(bitInd, bitMask) > 0){
      numBitsHigh++;
    }
  }
  return numBitsHigh;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Set color scale where white is in the middle.
 *
 * You need to draw things with Acclaim::RootTools::draw2D(TH2D* hist, TString opt);
*/
void Acclaim::RootTools::setWhiteZeroColorScale(){
  const int NRGBs = 3, NCont = 999;
  gStyle->SetNumberContours(NCont);
  Double_t stops[NRGBs] = { 0.00, 0.50, 1.00};
  Double_t red[NRGBs]   = { 0.00, 1.00, 1.00};
  Double_t green[NRGBs] = { 0.00, 1.00, 0.00};
  Double_t blue[NRGBs]  = { 1.00, 1.00, 0.00};
  TColor color;
  color.InitializeColors();
  color.CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
}






//---------------------------------------------------------------------------------------------------------
/**
 * @brief Draw 2D histogram and set bin limits to be symmetrical about zero.
 *
 * @param hist is the TH2D* to draw.
 * @param opt is the draw option.
*/
void Acclaim::RootTools::draw2D(TH2D* hist, TString opt){
  hist->Draw(opt);
  Double_t max = TMath::Abs(hist->GetMaximum());
  Double_t min = TMath::Abs(hist->GetMinimum());
  Double_t limit = min > max ? min : max;
  hist->SetMaximum(limit);
  hist->SetMinimum(-limit);
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Get a TChain of the ANITA-3 header data.
 *
 * @param firstRun is the first run.
 * @param lastRun is the last run.
 * @param headPtr a reference to a pointer to a RawAnitaHeader (yikes!).
*/
TChain* Acclaim::RootTools::getHeadChain(Int_t firstRun, Int_t lastRun, RawAnitaHeader*& headPtr){
  TChain* c = new TChain("headTree");
  for(Int_t run=firstRun; run <= lastRun; run++){
    c->Add(TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run));
  }
  c->SetBranchAddress("header", &headPtr);
  return c;
}





//---------------------------------------------------------------------------------------------------------
/**
 * @brief Get a TChain of the ANITA-3 adu5Pat data.
 *
 * @param firstRun is the first run.
 * @param lastRun is the last run.
 * @param pat a reference to a pointer to an Adu5Pat (yikes!).
*/
TChain* Acclaim::RootTools::getAdu5PatChain(Int_t firstRun, Int_t lastRun, Adu5Pat*& pat){
  TChain* c = new TChain("adu5PatTree");
  for(Int_t run=firstRun; run <= lastRun; run++){
    c->Add(TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsFile%d.root", run, run));
  }
  c->SetBranchAddress("pat", &pat);
  return c;
}






//---------------------------------------------------------------------------------------------------------
/**
 * @brief Use polarization and index to get the antenna name (1st phi-sector called 1, not 0).
 *
 * @param pol is the polarization
 * @param antInd is the antenna index (counting from 0)
*/

TString Acclaim::RootTools::getAntName(AnitaPol::AnitaPol_t pol, Int_t antInd){
  // Assumes antInd starts at 0. Returns the name w/ +1
  Int_t phi = (antInd % NUM_PHI) + 1;

  Int_t ring = antInd / NUM_PHI;

  TString antName = TString::Format("%d", phi);

  antName += AnitaRing::ringAsChar(AnitaRing::AnitaRing_t(ring));
  antName += AnitaPol::polAsChar(pol);

  return antName;
}






//---------------------------------------------------------------------------------------------------------
/**
 * @brief Gets the color a fraction through the current palette.
 *
 * @param index is the numerator of the fraction
 * @param maxVal is the denominator of the fraction
*/
Int_t Acclaim::RootTools::getColorFracThroughPalette(Int_t index, Int_t maxVal){
  return gStyle->GetColorPalette(index*Int_t(255./maxVal));
}






//---------------------------------------------------------------------------------------------------------
/**
 * @brief Draw histograms on the same (new) canvas with nice stats boxes
 *
 * @param numHists is, as you might guess, the number of histograms (the lengths of the array hs[])
 * @param hs[] is the array of pointers to TH1Ds to be plotted
 * @param drawOpt is the draw option in ROOT, e.g. "e" for error bars.
 * @param statsOption sets the stats to be displayed. I usually use "mre" for mean, RMS and entries.
*/
TCanvas* Acclaim::RootTools::drawHistsWithStatsBoxes(Int_t numHists, TH1D* hs[], TString drawOpt, TString statsOption){

  gStyle->SetOptStat(statsOption);
  TCanvas* theCan = new TCanvas();

  Double_t boxSize = 0.2; // In normalized coordinates
  const Int_t maxNumBoxes = 1./boxSize;

  Int_t divisor = (numHists/maxNumBoxes) + 1;

  // This is dumb and will get messy If I have more than 2*maxnumBoxes
  boxSize/=divisor;

  for(Int_t histInd=0; histInd < numHists; histInd++){

    TString thisDrawOpt = histInd == 0 ? drawOpt : drawOpt+"sames";
    hs[histInd]->Draw(thisDrawOpt);

    theCan->Update();
    TPaveStats* sbox = (TPaveStats*) hs[histInd]->FindObject("stats");
    sbox->SetTextColor(hs[histInd]->GetLineColor());
    sbox->SetX1NDC(0.8);
    sbox->SetY1NDC(1 - (histInd+1)*boxSize);
    sbox->SetX2NDC(0.98);
    sbox->SetY2NDC(0.98 - histInd*boxSize);
  }

  theCan->Update();
  return theCan;
}





/**
 * @brief Generic function to set the bin labels from the low edge values.
 *
 * @param ax is the histogram axis
 * @param format is the format specifier for TString/sprintf, should match histogram type, e.g. %2.0lf if hist is TH1D
 */
void Acclaim::RootTools::setBinLabelsFromLowEdge(TAxis* ax, const char* format){

  Int_t n = ax->GetNbins();
  for(int bin=1; bin <= n; bin++){
    ax->SetBinLabel(bin, TString::Format(format, ax->GetBinLowEdge(bin)));
  }
  ax->SetLabelOffset(0.01);
}





/** 
 * Helper function for drawSummary
 * draws a new TPad inside the parent TPad (cds into it too) with 
 * @param parentPad the pad to draw inside of
 * @param xlow is relative to parent pad
 * @param ylow is relative to parent pad
 * @param xup is relative to parent pad
 * @param yup is relative to parent pad
 * @param suffix is the suffix to append to the pad name
 * 
 * @return the new subPad
 */
TPad* Acclaim::RootTools::makeSubPad(TPad* parentPad, double xlow, double ylow, double xup, double yup, TString suffix){

  parentPad->cd(); // go into parent pad (assume it's drawn, would that matter?)
  TString subPadName = TString::Format("%s_%s", parentPad->GetName(), suffix.Data());
  TPad* subPad = new TPad(subPadName, subPadName, xlow, ylow, xup, yup);
  subPad->Draw();
  subPad->cd();
  return subPad;
}






/** 
 * Utility function to make sure getTimeIntegratedPower matches findSmallestWindowContainingFracOfPower
 * NO BOUNDS CHECKING IS DONE!
 * 
 * @param gr is the graph
 * @param samp is the sample
 * 
 * @return 
 */
double inline samplePower(const TGraph* gr, int samp){
  return gr->GetY()[samp]*gr->GetY()[samp]*(gr->GetX()[samp+1]-gr->GetX()[samp]);
}



/** 
 * @brief Does \sum_{i} V[i]*V[i]*(t[i+1] - t[i]).
 * 
 * To get the actual power, one would have to take into account units and impedance.
 * So really this is some quantity proportional to the power.
 * Ignores the last sample as there is no t[i+1]
 * 
 * @param gr is the wave to find the time integrated power
 * @param firstSamp is the sample to begin the sum from (default = 0)
 * @param lastSamp is the sample to end the sum on (default = gr->GetN()-1)
 * 
 * @return (quantity proportional to) the time integrated power
 */
double Acclaim::RootTools::getTimeIntegratedPower(const TGraph* gr, int firstSamp, int lastSamp){

  // here actually set the lastSamp default, we can't reference the gr in the function prototype
  lastSamp == -1 ? lastSamp = gr->GetN() - 1 : lastSamp;

  if(firstSamp < 0 || firstSamp >= gr->GetN() ){
    int newFirstSamp = 0;
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", firstSamp = " << firstSamp
              << ", for " << gr->GetName() << " with " << gr->GetN()
              << "samples, setting to " << newFirstSamp << std::endl;
    firstSamp = newFirstSamp;
  }

  if(lastSamp < 0 || lastSamp >= gr->GetN()){
    int newLastSamp = gr->GetN() - 1;
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", lastSamp = " << lastSamp
              << ", for " << gr->GetName() << " with " << gr->GetN()
              << "samples, setting to " << newLastSamp << std::endl;
    lastSamp = newLastSamp;
  }
  
  double totalPower = 0;
  for(int samp=firstSamp; samp < lastSamp; samp++){    
    totalPower += samplePower(gr, samp);
  }
  return totalPower;
}




/** 
 * Find the smallest window containing a desired fraction of power
 * 
 * @param grPow is a time domain, squared waveform
 * @param fracOfPowerInWindow is the desired amount of power to find in the waveform
 * @param peakTime is the time of the peak of power, if known
 * 
 * @return a pair of times corresponding to the smallest window around the peak containing the desired fraction of power
 */
std::pair<double, double> Acclaim::RootTools::findSmallestWindowContainingFracOfPower(const TGraph* gr, double fracOfPowerInWindow){

  // A very naive implementation is O(N^{2}) which is actually quite slow for large N
  // So I'm going to try to make this O(N) with a little harder-to-read code
  
  // First, find total power and the power that we want
  const double totalPower = getTimeIntegratedPower(gr);
  const double desiredPower = totalPower*fracOfPowerInWindow;

  double minWindowWidth = DBL_MAX;
  double minWinStart = kUnknownPeakTime; // this is negative
  double minWinEnd = -kUnknownPeakTime; // this is positive

  // This little bit *=MUST=* match what is done in getTimeIntegratedPower()
  // expand initial window until we find the desired fraction of power
  int lastJ=0;
  double windowPower = 0;
  while(windowPower < desiredPower && lastJ <= gr->GetN() - 1){ // Add -1 here to account for not using last sample in getTimeIntegratedPower
    windowPower += samplePower(gr, lastJ);
    lastJ++;
  }

  // If we doun
  if(windowPower >= desiredPower){
    minWinStart = gr->GetX()[0];
    minWinEnd = gr->GetX()[lastJ];
    minWindowWidth = minWinEnd - minWinStart;
  }
  else{
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", couldn't find initial window containing "
              << 100*fracOfPowerInWindow << "% of the power. Giving up." << std::endl;
    return std::pair<double, double>(minWinStart, minWinEnd);
  }

  // now we sweep through the remaining waveform,
  // removing the first sample and expanding the window
  // until we get back to the desired sample again
  for(int i=0; i < gr->GetN() - 2; i++){
    windowPower -= samplePower(gr, i);

    while(windowPower < desiredPower && lastJ < gr->GetN() - 2){ // Add -1 here to account for not using last sample in getTimeIntegratedPower
      windowPower += samplePower(gr, lastJ);
      lastJ++;
    }
    if(windowPower >= desiredPower){
      double windowWidth =  gr->GetX()[lastJ] - gr->GetX()[0];
      if(windowWidth < minWindowWidth){
        minWinStart = gr->GetX()[0];
        minWinEnd = gr->GetX()[lastJ];
        minWindowWidth = minWinEnd - minWinStart;
      }
    }
    else{
      // we've reached the end of the waveform
      return std::pair<double, double>(minWinStart, minWinEnd);
    }
  }

  // stop the compiler complaining, but we shouldn't get here...
  // I'll print a warning and move on
  std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", you shouldn't get here, something is wrong!" << std::endl;
  std::cerr << windowPower << "\t" << totalPower << "\t" << gr->GetN() << "\t" << lastJ << "\t" << minWinStart << "\t" << minWinEnd << "\t" << std::endl;
  
  return std::pair<double, double>(minWinStart, minWinEnd);

  // for(int i=0; i < gr->GetN() - 2; i++){ // Add -2 here to account for not using last sample in getTimeIntegratedPower, and have i < j
  //   Bool_t foundFracWidthThisJ = false;
  //   for(int j=lastJ; j < gr->GetN() - 1; j++){
  //     double winStart = gr->GetX()[i];
  //     double winEnd = gr->GetX()[j];

  //     double windowPower = getTimeIntegratedPower(gr, i, j);

  //     if(windowPower >= desiredPower){
  //       double windowWidth = winEnd - winStart;

  //       if(windowWidth < minWindowWidth){
  //         minWinStart = winStart;
  //         minWinEnd = winEnd;
  //         lastJ = j;
  //         foundFracWidthThisJ = true;
  //         // std::cerr << i << "\t" << j << std::endl;
  //       }
  //     }

  //     if(foundFracWidthThisJ){
  //       // no point checking higher values of j
  //       // since that window will be larger for the same i
  //       break;
  //     }
  //   }
  // }

  // if(minWinStart==kUnknownPeakTime && minWinEnd==-kUnknownPeakTime){
  //   std::cerr << "Warning in " << __PRETTY_FUNCTION__
  //             << ", couldn't find a window containing "
  //             << fracOfPowerInWindow << " of the power."
  //             << std::endl;
  // }

  // return std::pair<double, double>(minWinStart, minWinEnd);
}
