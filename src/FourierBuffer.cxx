#include "FourierBuffer.h"
#include "TMath.h"
#include "FFTWComplex.h"
#include "TSpectrum.h"
#include "RayleighHist.h"
#include "TF1.h"

const char* rayleighFuncText = "([0]*x/([1]*[1]))*exp(-x*x/(2*[1]*[1]))";
const char* riceFuncText = "([0]*x/([1]*[1]))*exp(-(x*x+[2]*[2])/(2*[1]*[1]))*TMath::BesselI0([2]*x/([1]*[1]))";


Acclaim::FourierBuffer::FourierBuffer(double timeScaleSeconds, Int_t theAnt, AnitaPol::AnitaPol_t thePol){
  timeScale = timeScaleSeconds;
  ant = theAnt;
  pol = thePol;
  df = -1;
  fDrawFreqBin = 100; //26;

  // will initialize this dynamically to get around this no-copy-constructor bullshit
  spectrum = NULL;

  TString funcName = TString::Format("fRay_%d_%d", ant, pol);
  fRay = new TF1(funcName, rayleighFuncText, 0, 1);

  doneVectorInit = false;
}



void Acclaim::FourierBuffer::initVectors(int n){
  sumPower.resize(n, 0);
  sumAmps.resize(n, 0);
  hRays.resize(n, NULL);

  for(int freqBin=0; freqBin < n; freqBin++){
    TString name = TString::Format("hRayleigh");
    name += ant >= 0 && pol < AnitaPol::kNotAPol ? TString::Format("_%d_%d", ant, pol) : "";
    name += TString::Format("_%d", freqBin);
    hRays.at(freqBin) = new Acclaim::RayleighHist(this, name, name);
  }
  doneVectorInit = true;
}



Acclaim::FourierBuffer::~FourierBuffer(){
  if(spectrum){
    delete spectrum;
    spectrum = NULL;
  }
  
  for(unsigned i=0; i < hRays.size(); i++){
    delete hRays.at(i);
    hRays.at(i) = NULL;
  }
}






size_t Acclaim::FourierBuffer::add(const RawAnitaHeader* header, const AnalysisWaveform* wave){

  // get the power
  const TGraphAligned* grPower = wave->power();


  // do dynamic initialization and sanity checking
  if(!doneVectorInit){
    initVectors(grPower->GetN());
  }
  if(grPower->GetN() != (int)sumPower.size()){
    std::cerr << "Warning in " << __PRETTY_FUNCTION__ << ", unexpected waveform size" << std::endl;
  }
  if(df <= 0){
    df = grPower->GetX()[1] - grPower->GetX()[0];
    for(size_t freqBin=0; freqBin < hRays.size(); freqBin++){
      hRays.at(freqBin)->freqMHz = df*1e3*freqBin;
    }
  }
  
  eventNumbers.push_back(header->eventNumber);
  runs.push_back(header->run);

  // Double_t realTime= header->realTime;
  Double_t realTime= header->triggerTime;
  realTime += 1e-9*header->triggerTimeNs;

  // if(ant==0 && pol == 0){
  //   printf("%d\t%d\t%u\t%u\t%20.20lf\n", ant, pol, header->realTime, header->triggerTimeNs, realTime);
  // }


  realTimesNs.push_back(realTime);

  // for old compilers, push back copy of empty vector for speed.
  // then get reference that vector in the list
  powerRingBuffer.push_back(std::vector<double>(0));
  std::vector<double>& freqVec = powerRingBuffer.back();
  freqVec.assign(grPower->GetY(), grPower->GetY()+grPower->GetN());

  
  // update sum of power
  for(int freqInd=0; freqInd < grPower->GetN(); freqInd++){
    sumPower.at(freqInd) += grPower->GetY()[freqInd];

    double amp = TMath::Sqrt(grPower->GetY()[freqInd]);
    
    sumAmps.at(freqInd) += amp;

    const double meanAmp = sumAmps.at(freqInd)/eventNumbers.size();

    bool needRebin = hRays.at(freqInd)->axisRangeOK(meanAmp);
    if(needRebin){
      hRays.at(freqInd)->rebinAndEmptyHist(meanAmp);
      std::list<std::vector<double> >::const_iterator powVecPtr = powerRingBuffer.begin();
      while(powVecPtr!=powerRingBuffer.end()){
	double oldAmp = TMath::Sqrt(powVecPtr->at(freqInd));
	hRays.at(freqInd)->Fill(oldAmp);
	++powVecPtr;
      }
    }
    else{
      hRays.at(freqInd)->Fill(amp);
    }
  }
  
  removeOld();


  for(int freqInd=0; freqInd < grPower->GetN(); freqInd++){
    hRays.at(freqInd)->Fit();
  }
  
  return eventNumbers.size();
}



Int_t Acclaim::FourierBuffer::removeOld(){

  Int_t nPopped = 0;
  if(realTimesNs.size() > 0){

    Double_t mostRecentTime = realTimesNs.back();

    while(mostRecentTime - realTimesNs.front() > timeScale){
      if(ant==0 && pol == 0){
	printf("removing: %lf\t%20.20lf\t%20.20lf\n", mostRecentTime - realTimesNs.front(), realTimesNs.front(), mostRecentTime);
      }
      
      realTimesNs.pop_front();
      eventNumbers.pop_front();
      runs.pop_front();
      std::vector<double>& removeThisPower = powerRingBuffer.front();
      for(unsigned int freqInd=0; freqInd < removeThisPower.size(); freqInd++){
	sumPower.at(freqInd) -= removeThisPower.at(freqInd);
	double amp = TMath::Sqrt(removeThisPower.at(freqInd));
	sumAmps.at(freqInd) -= amp;
	hRays.at(freqInd)->Fill(amp, -1);
      }
      powerRingBuffer.pop_front();
      nPopped++;
    }
    if(ant==0 && pol == 0){
      std::cout << "after removal: " << mostRecentTime - realTimesNs.front() << "\t" << nPopped << "\t" << std::endl;
    }
    
  }
  return nPopped;
}






// const Acclaim::RayleighHist* Acclaim::FourierBuffer::fillRayleighInfo(Int_t freqBin, RayleighInfo* info) const{

//   const RayleighHist* hRay = getRayleighDistribution(freqBin);
//   TString fName = TString::Format("fRay%d", freqBin);
//   fName += ant >= 0 && pol < AnitaPol::kNotAPol ? TString::Format("_%d_%d", ant, pol) : "";
//   fName += TString::Format("_%d", TMath::Nint(1e3*df*freqBin));

//   Double_t xMin = hRay->GetXaxis()->GetBinLowEdge(1);
//   Double_t xMax = hRay->GetXaxis()->GetBinLowEdge(hRay->GetNbinsX()+1);

//   // TF1* fRay = fRays.at(freqBin);
//   fRay->SetRange(xMin, xMax);
  
  
//   // TF1* fRay = new TF1(fName, rayleighFuncText, xMin, xMax);

//   // Fill event level info
//   info->eventNumber = eventNumbers.back();
//   info->run = runs.back();
//   info->realTimeNs = realTimesNs.back();
  
//   info->firstRun = runs.front();
//   info->firstEventNumber = eventNumbers.front();
//   info->firstRealTimeNs = realTimesNs.front();

//   double powerBin = powerRingBuffer.back().at(freqBin);
//   info->eventAmp = TMath::Sqrt(powerBin);
  
//   // Fill histogram related histogram info
//   info->nBins = hRay->GetNbinsX();
//   info->integralPlusOverUnderFlow = hRay->Integral(0, info->nBins+1);
//   info->binWidth = hRay->GetBinWidth(1);
//   info->rayFitNorm = info->integralPlusOverUnderFlow*info->binWidth;
//   info->rayGuessAmp = hRay->GetMean()*TMath::Sqrt(2./TMath::Pi());

//   // Set fit params
//   fRay->FixParameter(0, info->rayFitNorm);
//   fRay->SetParameter(1, info->rayGuessAmp);
//   fRay->SetParLimits(1, info->rayGuessAmp, info->rayGuessAmp); // essentially infinite

//   info->rayChiSquare = hRay->Chisquare(fRay);

//   int nonZeroBins = 0;
//   for(int bx=1; bx <= hRay->GetNbinsX(); bx++){
//     if(hRay->GetBinContent(bx) > 0){
//       nonZeroBins++;
//     }
//   }
  
//   info->rayProb = TMath::Prob(info->rayChiSquare, nonZeroBins);
//   // std::cout << ant << "\t" << pol << "\t" << freqBin << "\t" << info->rayChiSquare << "\t" <<  nonZeroBins << "\t" << info->rayProb << std::endl;
  
//   // hRay->Fit(fRay, "Q0"); //, "", "", xMin, xMax);

//   // // Do fit
//   // // std::cout << ant << "\t" << pol << "\t" << eventNumbers.back() << "\t" << freqBin << "\t" << fRay->GetChisquare() << "\t";  
//   // // fRay->SetParLimits(1, 0, 1e9); // essentially infinite  
//   // // hRay->Fit(fRay, "Q0"); //, "", "", xMin, xMax);
//   // // std::cout << fRay->GetChisquare() << std::endl;
  
//   // Fill fit result params
//   // info->rayFitAmp = fRay->GetParameter(1);
//   // info->rayFitAmpError = fRay->GetParError(1);
//   // info->rayChiSquare = fRay->GetChisquare();
//   // info->rayNdf = fRay->GetNDF();
//   // info->rayProb = fRay->GetProb();

//   return hRay;
// }



// const Acclaim::RayleighHist* Acclaim::FourierBuffer::fillRiceInfo(Int_t freqBin, RiceInfo* info) const{

//   const RayleighHist* hRay = fillRayleighInfo(freqBin, info);

//   // for now...
//   if(info->rayChiSquare/info->rayNdf > 2){

//     TString fName = TString::Format("fRice%d", freqBin);
//     Double_t xMin = hRay->GetXaxis()->GetBinLowEdge(1);
//     Double_t xMax = hRay->GetXaxis()->GetBinLowEdge(hRay->GetNbinsX()+1);
//     TF1* fRice = new TF1(fName, riceFuncText, xMin, xMax);

//     fRice->FixParameter(0, info->rayFitNorm);
//     fRice->SetParameter(1, info->rayGuessAmp);
//     fRice->SetParLimits(1, 0, 1e9); // essentially infinite
//     fRice->SetParameter(2, 0);
//     fRice->SetParLimits(2, 0, 1e9); // essentially infinite

//     // hRay->Fit(fRice, "Q0");


//     info->riceFitNorm = fRice->GetParameter(0); // fixed for the fit  
//     info->riceFitAmp = fRice->GetParameter(1);
//     info->riceFitAmpError = fRice->GetParError(1);  
//     info->riceFitSignal = fRice->GetParameter(2);
//     info->riceFitSignalError = fRice->GetParError(2);  

//     info->riceChiSquare = fRice->GetChisquare();
//     info->riceNdf = fRice->GetNDF();

//     delete fRice;
//   }
  
//   return hRay;
// }







TGraphAligned* Acclaim::FourierBuffer::getAvePowSpec_dB(double thisTimeRange) const{

  TGraphAligned* gr = getAvePowSpec(thisTimeRange);
  gr->dBize();
  return gr;

}


TGraphAligned* Acclaim::FourierBuffer::getAvePowSpec(double thisTimeRange) const{

  // set default value (which is the whole range of the fourier buffer)
  thisTimeRange = thisTimeRange < 0 ? timeScale : thisTimeRange;

  int nEvents = eventNumbers.size();

  TGraphAligned* gr = new TGraphAligned(sumPower.size());
  
  TString name;
  TString title;  
  if(ant < 0){
    name = TString::Format("grAvePowSpec_%u_%u", eventNumbers.front(), eventNumbers.back());
    title = TString::Format("Average Power Spectrum event %u - %u (unknown channel)", eventNumbers.front(), eventNumbers.back());    
  }
  else{
    const char* polName = AnitaPol::kHorizontal ? "HPol" : "VPol";
    name = TString::Format("grAvePowSpec_%d_%s_%u_%u", ant, polName, eventNumbers.front(), eventNumbers.back());
    title = TString::Format("Average Power Spectrum %d %s event %u - %u", ant, polName, eventNumbers.front(), eventNumbers.back());        
  }
  gr->SetName(name);
  gr->SetTitle(title);
  

  for(unsigned freqInd=0; freqInd < sumPower.size(); freqInd++){
    gr->SetPoint(freqInd, freqInd*df, sumPower.at(freqInd)/nEvents);
  }

  return gr;
}



TGraphAligned* Acclaim::FourierBuffer::getBackground_dB(double thisTimeRange) const{
  
  TGraphAligned* gr = getBackground(thisTimeRange);
  gr->dBize();
  return gr;
}


TGraphAligned* Acclaim::FourierBuffer::getBackground(double thisTimeRange) const{

  TGraphAligned* gr = getAvePowSpec(thisTimeRange);

  if(!spectrum){
    spectrum = new TSpectrum();
  }

  // at some point this changed from doubles to floats, not added a new method, changed...
  // git blames Lorenzo Moneta... 
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)

  spectrum->Background(gr->GetY(),gr->GetN(),
                      6,TSpectrum::kBackDecreasingWindow,
                      TSpectrum::kBackOrder2,kFALSE,
                      TSpectrum::kBackSmoothing3,kFALSE);
#else
  std::vector<float> tempFloats(gr->GetN(), 0);
  for(int i=0; i < gr->GetN(); i++){
    tempFloats[i] = gr->GetY()[i];
  }
  spectrum->Background(&tempFloats[0],gr->GetN(),
                      6,TSpectrum::kBackDecreasingWindow,
                      TSpectrum::kBackOrder2,kFALSE,
                      TSpectrum::kBackSmoothing3,kFALSE);
  for(int i=0; i < gr->GetN(); i++){
    gr->GetY()[i] = tempFloats[i];
  }
#endif
  

  return gr;
}

