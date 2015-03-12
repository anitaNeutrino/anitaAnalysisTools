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
