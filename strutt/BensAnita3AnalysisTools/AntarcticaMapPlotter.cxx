#include "AntarcticaMapPlotter.h"


/*!
  \brief Constructor
  \param name is the histogram name
  \param title is the histogram title
  \param nBinsX is the number of bins on the x-axis
  \param nBinsY is the number of bins on the y-axis
*/
AntarcticaMapPlotter::AntarcticaMapPlotter(TString name, TString title, Int_t nBinsX, Int_t nBinsY){
  initializeInternals();
  addHistogram(name, title, nBinsX, nBinsY);
}


/*!
  \brief Destructor
*/
AntarcticaMapPlotter::~AntarcticaMapPlotter(){
  std::map<TString, TH2D*>::iterator histItr;
  for(histItr=hists.begin(); histItr != hists.end(); ++histItr){
    delete histItr->second;
  }
  std::map<TString, TGraph*>::iterator grItr;
  for(grItr=grs.begin(); grItr != grs.end(); ++grItr){
    delete grItr->second;
  }
}


/*!
  \brief Workhorse function which initalizes internal variables
  \param latitude is the latitude
  \param longitude is the longitude
  \param x is the internal histogram variable x
  \param y is the internal histogram variable y
*/
void AntarcticaMapPlotter::initializeInternals(){
  /* The numbers in this function also come from Matt Mottram or Ryan */
  
  TrueScaleLat=71;
  CentralMeridian=0;
  RadiusOfEarth=6378.1e3; //Metres
  xOffest=375;
  yOffset=312.5; 
  scale=271.5/2.19496e+06;
  xSize=750;
  ySize=625;
  

  img = TImage::Open("~/Repositories/Install/share/anitaMap/antarcticaIceMapBW.png");
  img->SetConstRatio(kFALSE);  
}

/*!
  \brief The magic function that maps latitude and longitude onto histogram x and y coordinates
  \param latitude is the latitude
  \param longitude is the longitude
  \param x is the internal histogram variable x
  \param y is the internal histogram variable y
*/
void AntarcticaMapPlotter::getRelXYFromLatLong(Double_t latitude, Double_t longitude,Double_t &x, Double_t &y){
  /* This function, which really does all the hard work, was made either by Matt Mottram or Ryan */
  
  // Negative longitude is west
  // All latitudes assumed south
  Double_t absLat=TMath::Abs(latitude);
  Double_t r=RadiusOfEarth*TMath::Cos((90.-TrueScaleLat)*TMath::DegToRad())*TMath::Tan((90-absLat)*TMath::DegToRad());
  y=r*TMath::Cos(longitude*TMath::DegToRad());
  x=r*TMath::Sin(longitude*TMath::DegToRad());

  y*=scale;
  y+=yOffset;
  y/=ySize;
  x*=scale;
  x+=xOffest;
  x/=xSize;
}


/*!
  \brief How to fill the histogram, converts latitude and longitude into x/y bins first.
  \param latitude is the latitude
  \param longitude is the longitude
*/
Int_t AntarcticaMapPlotter::Fill(Double_t latitude, Double_t longitude, Double_t weight){
  Double_t x, y;
  getRelXYFromLatLong(latitude, longitude, x, y);
  return hCurrent->Fill(x, y, weight); 
}


/*!
  \brief Draws the canvas, image and histogram
*/
void AntarcticaMapPlotter::DrawHist(TString opt){

  TString canName = TString::Format("can%s", hCurrent->GetName());
  TString canTitle = TString::Format("Canvas of %s", hCurrent->GetTitle());
  TCanvas* can = new TCanvas(canName, canTitle, (Int_t) xSize, (Int_t) ySize);
  can->Draw();
  can->SetTopMargin(0.03);
  can->SetBottomMargin(0.03);
  can->SetLeftMargin(0.03);
  can->SetRightMargin(0.15);
  img->Draw("same");  
  hCurrent->Draw(opt + "same");
}


/*!
  \brief Draws the canvas, image and histogram
*/
void AntarcticaMapPlotter::DrawTGraph(TString opt){

  TString canName = TString::Format("can%s", grCurrent->GetName());
  TString canTitle = TString::Format("Canvas of %s", grCurrent->GetTitle());
  TCanvas* can = new TCanvas(canName, canTitle, (Int_t) xSize, (Int_t) ySize);
  can->Draw();
  can->SetTopMargin(0.03);
  can->SetBottomMargin(0.03);
  can->SetLeftMargin(0.03);
  can->SetRightMargin(0.15);
  img->Draw("same");  
  grCurrent->Draw(opt + "same");
}



/*!
  \brief Sets the current histogram pointer \a hCurrent.
  \returns 0 on failure (TH2D with name not in \a hists) and 1 on success.
*/
Int_t AntarcticaMapPlotter::setCurrentHistogram(TString name){
  
  std::map<TString, TH2D*>::iterator histItr = hists.find(name);
  Int_t successState = 0;
  if(histItr == hists.end()){
    successState = 0;
  }
  else{
    hCurrent = histItr->second;
    successState = 1;
  }
  return successState;
}

/*!
  \brief Creates a new histogram to go in the internal map and sets it as the current histogram
  \param name is the histogram name
  \param title is the histogram title
  \param nBinsX is the number of bins on the x-axis
  \param nBinsY is the number of bins on the y-axis
*/
void AntarcticaMapPlotter::addHistogram(TString name, TString title, Int_t nBinsX, Int_t nBinsY){

  if(setCurrentHistogram(name)==0){
    TH2D* theHist = new TH2D(name, title, nBinsX, 0, 1, nBinsY, 0, 1);
    hists[name] = theHist;
    hCurrent = theHist;
  }
  else{
    std::cerr << "Warning in " << __FILE__ << "! Histogram with name " << name.Data() << ", already exists!" << std::endl;
  }
}


/*!
  \brief Creates a new graph to go in the internal map and sets it as the current histogram
  \param name is the histogram name
  \param title is the histogram title
  \param n is the number of points
  \param latitudes is a pointer to an array of latitude values
  \param longitudes is a pointer to an array of longitude values
*/
  
void AntarcticaMapPlotter::addTGraph(TString name, TString title, Int_t n,
				     Double_t* latitudes, Double_t* longitudes){

  if(setCurrentTGraph(name)==0){
    std::vector<Double_t> xs(n);
    std::vector<Double_t> ys(n);
    for(Int_t i=0; i<n; i++){
      getRelXYFromLatLong(latitudes[i], longitudes[i], xs[i], ys[i]);
    }
    TGraph* theGraph = new TGraph(n, &xs[0], &ys[0]);
    
    theGraph->SetName(name);
    theGraph->SetTitle(title);
    grs[name] = theGraph;
    grCurrent = theGraph;
  }
  else{
    std::cerr << "Warning in " << __FILE__ << "! TGraph with name "
	      << name.Data() << ", already exists!" << std::endl;
    if(n>0){
      std::cerr << "No points added to " << name.Data() << "!" << std::endl;      
    }
  }
}



/*!
  \brief Sets the current histrogram pointer \a hCurrent.
  \returns 0 on failure (TGraph with name not in \a grs) and 1 on success.
*/
Int_t AntarcticaMapPlotter::setCurrentTGraph(TString name){
  
  std::map<TString, TGraph*>::iterator grItr = grs.find(name);
  Int_t successState = 0;
  if(grItr == grs.end()){
    successState = 0;
  }
  else{
    grCurrent = grItr->second;
    successState = 1;
  }
  return successState;
}



/*!
  \brief Return pointer to current TGraph
  \returns pointer to current TGraph
*/
TGraph* AntarcticaMapPlotter::getCurrentTGraph(){
  return grCurrent;
}

/*!
  \brief Return pointer to current histogram
  \returns pointer to current histogram
*/
TH2D* AntarcticaMapPlotter::getCurrentHistogram(){
  return hCurrent;
}
