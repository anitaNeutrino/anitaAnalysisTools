// -*- C++ -*-.
/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Give this class a TTree plus a branch name to serve as an 'x-axis', something like unixtime.
             FancyTTreeInterpolator generates TGraphs of any TTree branch variable as a function of your x-axis.
             That allows us to do two slightly clever things.
             Firstly, it sorts the data so it is x-axis (i.e. time) ordered.
             Secondly, you can interpolate the data to values between TTree entries.
             This interpolation is equivalent to TGraph::Eval.
*************************************************************************************************************** */

#include "FancyTTreeInterpolator.h"
// #include "TFile.h"


FancyTTreeInterpolator::FancyTTreeInterpolator(TTree* t, TString xAxisText){

  fXAxisText = xAxisText;
  fTree = t;
  std::shared_ptr<TGraph> gr = makeSortedTGraph(xAxisText + ":" + xAxisText);
  fXmin = gr->GetX()[0];
  fXmax = gr->GetX()[gr->GetN()-1];

  // TFile* fout = new TFile("ftti.root", "recreate");
  // gr->Write();
  // fout->Close();

  /* For error reporting, want precision of unsigned int ~10 digits*/
  std::cerr.precision(10);

}


FancyTTreeInterpolator::~FancyTTreeInterpolator(){
  
}



std::shared_ptr<TGraph> FancyTTreeInterpolator::makeSortedTGraph(TString drawText){
  return makeSortedTGraph(drawText, "", 0);
}

std::shared_ptr<TGraph> FancyTTreeInterpolator::makeSortedTGraph(TString drawText, TString cutString){
  return makeSortedTGraph(drawText, cutString, 0);
}

std::shared_ptr<TGraph> FancyTTreeInterpolator::makeSortedTGraph(TString drawText, Double_t wrapValue){
  return makeSortedTGraph(drawText, "", wrapValue);
}

std::shared_ptr<TGraph> FancyTTreeInterpolator::makeSortedTGraph(TString drawText, TString cutString, Double_t wrapValue){
  /*
    This function:
       Uses TTree::Draw to select the data.
       Sorts it using TMath::Sort.
       Returns the graph.
       Note: This does not add it to the std::map of graphs that this class uses for storage.
             To create a graph and add it to the std::map use FancyTTreeInterpolator::add.
  */
  
  // Draw
  const Int_t nEntries = fTree->Draw(drawText, cutString, "goff"); // "goff" means graphics off

  // Sort
  std::vector<Int_t> sortedIndices(nEntries);
  TMath::Sort(nEntries, fTree->GetV2(), sortedIndices.data(), kFALSE);
  std::vector<Double_t> newX(nEntries);
  std::vector<Double_t> newY(nEntries);
  
  for(int i=0; i<nEntries; i++){
    newX.at(i) = fTree->GetV2()[sortedIndices.at(i)];
    newY.at(i) = fTree->GetV1()[sortedIndices.at(i)];
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
  std::shared_ptr<TGraph> gr(new TGraph(nEntries,newX.data(), newY.data()));
  gr->SetTitle(drawText + ", " + cutString);
  gr->GetXaxis()->SetTitle(fXAxisText);
  gr->GetYaxis()->SetTitle(drawText);

  return gr;

}

void FancyTTreeInterpolator::add(TString yAxisText){
  add(yAxisText, "", 0);
}

void FancyTTreeInterpolator::add(TString yAxisText, TString cutString){
  add(yAxisText, cutString, 0);
}

void FancyTTreeInterpolator::add(TString yAxisText, Double_t wrapValue){
  add(yAxisText, "", wrapValue);
}

void FancyTTreeInterpolator::add(TString yAxisText, TString cutString, Double_t wrapValue){
  TString drawText = yAxisText + ":" + fXAxisText;
  std::shared_ptr<TGraph> gr = makeSortedTGraph(drawText, cutString, wrapValue);
  fStringToGraph[yAxisText] = gr;
  fStringToWrapValue[yAxisText] = wrapValue;
}




std::shared_ptr<TGraph> FancyTTreeInterpolator::get(TString yAxisText){
  /* Use this to access a graph you've made with the add function */

  if(fStringToGraph.count(yAxisText)==0){
    std::cerr << "Can't find TGraph in FancyTTreeInterpolator created with text " + yAxisText << std::endl;
    throw std::invalid_argument("Can't find TGraph in FancyTTreeInterpolator created with text " + yAxisText);
  }
  else{
    return fStringToGraph.find(yAxisText)->second;
  }
  return nullptr;
}


Double_t FancyTTreeInterpolator::interp(TString yAxisText, Double_t xAxisValue){
  /* This function handles the getting of values from the appropriate TGraph */

  std::shared_ptr<TGraph> gr = get(yAxisText);

  if(xAxisValue >= fXmin && xAxisValue <= fXmax){
    Double_t tempVal = gr->Eval(xAxisValue);

    // do the unwrapping if required
    Double_t wrap = fStringToWrapValue.find(yAxisText)->second;
    if(wrap != 0){
      while (tempVal < 0){
        tempVal += wrap;
      }
      while (tempVal >= wrap){
        tempVal -= wrap;	
      }
    }
    return tempVal;
  }
  else{
    std::cerr << "Value " << xAxisValue << " lies outside range of " << fXAxisText << std::endl;
    std::cerr << "xMin = " << fXmin << ", xMax = " << fXmax << std::endl;
    throw std::domain_error("");
  }
  return 0;

}
