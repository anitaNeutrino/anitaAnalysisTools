#include "RootTools.h"

double RootTools::getSumOfYVals(TGraph* gr){
  double sum = 0;
  for(int i=0; i<gr->GetN(); i++){
    sum += gr->GetY()[i];
  }
  return sum;
}

void RootTools::printArray(int n, double* array, TString delimiter, TString start ,TString end){
  std::cout << start.Data();
  for(int i=0; i<n; i++){
    std::cout << array[i];
    if(i<n-1) std::cout << delimiter.Data();
  }  
  std::cout << end.Data();
}

void RootTools::printYVals(TGraph* gr, TString delimiter, TString start, TString end){
  printArray(gr->GetN(), gr->GetY(), delimiter, start, end);
}

void RootTools::printXVals(TGraph* gr, TString delimiter, TString start, TString end){
  printArray(gr->GetN(), gr->GetX(), delimiter, start, end);
}

TGraph* RootTools::makeDerivativeTGraph(TGraph* gr){
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

void RootTools::getMaxMin(TGraph* gr, Double_t& maxY, Double_t& minY){
  Double_t maxX=0;
  Double_t minX=0;
  RootTools::getMaxMin(gr, maxY, maxX, minY, minX);
}


void RootTools::getMaxMin(TGraph* gr, Double_t& maxY, Double_t& maxX, Double_t& minY, Double_t& minX){

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


void RootTools::subtractOffset(TGraph* gr, Double_t offset){
  for(int i=0; i<gr->GetN(); i++){
    gr->GetY()[i] -= offset;
  }
}


void RootTools::normalize(TGraph* gr){
  double mean, rms;
  normalize(gr, mean, rms);
}


TGraph* RootTools::makeNormalized(TGraph* gr){
  /* Copies the TGraph and normalizes that */
  TGraph* grCopy = (TGraph*) gr->Clone();
  normalize(grCopy);
  return grCopy;
}


TGraph* RootTools::makeNormalized(TGraph* gr, Double_t& mean, Double_t& rms){
  /* Copies the TGraph and normalizes that */
  TGraph* grCopy = (TGraph*) gr->Clone();
  normalize(grCopy, mean, rms);
  return grCopy;
}



void RootTools::normalize(TGraph* gr, Double_t& mean, Double_t& rms){
  RootTools::getMeanAndRms(gr, mean, rms);
  if(rms>0){ // Don't make any NaNs
    for(int i=0; i<gr->GetN(); i++){    
      gr->GetY()[i] -= mean;
      gr->GetY()[i] /= rms;
    }
  }
}

void RootTools::getMeanAndRms(TGraph* gr, Double_t& mean, Double_t& rms){
  Double_t sum = 0;
  Double_t square = 0;
  for(int i=0; i<gr->GetN(); i++){
    sum += gr->GetY()[i];
    square += gr->GetY()[i]*gr->GetY()[i];
  }
  mean = sum/gr->GetN();
  rms = TMath::Sqrt(square/gr->GetN() - mean*mean);
}


std::vector<Int_t> RootTools::getIndicesOfNans(TGraph* gr){
  std::vector<Int_t> nanIndices;
  for(Int_t i=0; i<gr->GetN(); i++){
    if(TMath::IsNaN(gr->GetY()[i])){
      nanIndices.push_back(i);
    }
  }
  return nanIndices;
}


void RootTools::printTGraphInfo(TGraph* gr){
  std::cout << "******************************************************************************" << std::endl;
  std::cout << "RootTools::printTGraphValues(TGraph* gr = " << gr << "):" << std::endl;
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

TGraph* RootTools::makeSortedTGraph(TTree* tree, TString drawText, TString cutString, Double_t wrapValue){
  /*
    This function:
       Uses TTree::Draw to select the data.
       Sorts it using TMath::Sort.
       Returns the graph.
       Note: This does not add it to the std::map of graphs that this class uses for storage.
             To create a graph and add it to the std::map use FancyTTreeInterpolator::add.
  */
  
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


TH1D* RootTools::plotsZaxisDist(TH2* h2, TString hName, Int_t nBins, Double_t xMin, Double_t xMax){
  TH1D* h = new TH1D(hName, hName, nBins, xMin, xMax);
  for(int xBin=1; xBin<=h2->GetNbinsX(); xBin++){
    for(int yBin=1; yBin<=h2->GetNbinsY(); yBin++){
      Double_t val = h2->GetBinContent(xBin, yBin);
      h->Fill(val);
    }
  }
  return h;
}



void RootTools::zeroPadTGraph(TGraph* gr, Int_t newLen, Double_t dt){
  Int_t oldLen = gr->GetN();
  gr->Set(newLen);
  if(dt >= 0){
    Double_t x0 = gr->GetX()[oldLen-1];
    for(Int_t samp=oldLen; samp<newLen; samp++){
      gr->GetX()[samp] = x0 + samp*dt;
    }
  }
}


TGraph* RootTools::makeLinearlyInterpolatedGraph(TGraph* grIn, Double_t dt){

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
