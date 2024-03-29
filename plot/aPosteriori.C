#include "ThermalChain.h"
#include "DrawStrings.h"
#include "RootTools.h"

using namespace Acclaim;
const int nd = 3;


inline double dSq(const double* x, const double* y){
  double sumSq = 0;
  for(int i=0; i < nd;  i++){
    double delta = x[i] - y[i];
    sumSq += delta*delta;
  }
  return sumSq;
}

void aPosteriori(){

  TGraphAntarctica* grNu = new TGraphAntarctica();
  grNu->SetMarkerColor(kRed);
  grNu->SetName("grNu");

  TH2DAntarctica* hA = new TH2DAntarctica("hA", "asdf");
  TH2DAntarctica* hB = new TH2DAntarctica("hB", "ghjk");

  ThermalChain tc("data/makeThermalTree_*.root");
  auto c = tc.getChain();
  std::cout << c->GetEntries() << std::endl;

  auto geom = AnitaGeomTool::Instance();
  std::vector<TGraph*> grRcc;
  std::vector<TGraph*> grPwg;
  std::vector<TGraph*> grIm;
  std::vector<TGraph*> grPh;
  std::vector<TGraph*> grFd;

  std::vector<TH2F*> hRccs;
  std::vector<TH2F*> hPwgs;
  std::vector<TH2F*> hIms;
  std::vector<TH2F*> hPhs;
  std::vector<TH1F*> hFds;

  std::vector<std::vector<UInt_t> > eventGroups;
  std::vector<std::vector<Double_t> > lons, lats, alts, xs, ys, zs;

  // unknown-base singlets
  eventGroups.push_back(std::vector<UInt_t>{ 7613856});
  eventGroups.push_back(std::vector<UInt_t>{ 9097075});
  eventGroups.push_back(std::vector<UInt_t>{11116669});
  eventGroups.push_back(std::vector<UInt_t>{11989349});
  eventGroups.push_back(std::vector<UInt_t>{15717147}); // Peter's
  eventGroups.push_back(std::vector<UInt_t>{16952229});
  eventGroups.push_back(std::vector<UInt_t>{19459851});
  eventGroups.push_back(std::vector<UInt_t>{22345215});
  eventGroups.push_back(std::vector<UInt_t>{23695286});
  eventGroups.push_back(std::vector<UInt_t>{27142546});
  eventGroups.push_back(std::vector<UInt_t>{32907848});
  eventGroups.push_back(std::vector<UInt_t>{33484995});
  eventGroups.push_back(std::vector<UInt_t>{41529195});
  eventGroups.push_back(std::vector<UInt_t>{48837708});
  eventGroups.push_back(std::vector<UInt_t>{58592863});
  eventGroups.push_back(std::vector<UInt_t>{62273732});
  eventGroups.push_back(std::vector<UInt_t>{65187079});
  eventGroups.push_back(std::vector<UInt_t>{66313844});
  eventGroups.push_back(std::vector<UInt_t>{68298837});
  eventGroups.push_back(std::vector<UInt_t>{70013898});
  eventGroups.push_back(std::vector<UInt_t>{71171108});
  eventGroups.push_back(std::vector<UInt_t>{71766273});
  eventGroups.push_back(std::vector<UInt_t>{73726742});
  eventGroups.push_back(std::vector<UInt_t>{74592579});
  eventGroups.push_back(std::vector<UInt_t>{75277769});
  eventGroups.push_back(std::vector<UInt_t>{80840274});
  eventGroups.push_back(std::vector<UInt_t>{83877990}); // VPol

  eventGroups.push_back(std::vector<UInt_t> {84650299, 84653556}); // known-base doublet
  eventGroups.push_back(std::vector<UInt_t> {70029023, 70124461}); // known-base doublet
  
  eventGroups.push_back(std::vector<UInt_t> {41128241, 41475569}); // unknown Doublet
  eventGroups.push_back(std::vector<UInt_t> {69846620});           // known-base singlet

  std::vector<UInt_t> allEvents;
  TCut allEventsCut;
  for(int e=0; e < eventGroups.size(); e++){
    
    const std::vector<UInt_t>& eventsOfInterest = eventGroups[e];
    for(auto e : eventsOfInterest){
      allEvents.push_back(e);
      TString cutStr = TString::Format("eventNumber!=%u", e);
      allEventsCut += TCut(cutStr.Data());
    }


    TString titleBit = "";
    for(auto ev : eventsOfInterest){
      titleBit += TString::Format("%u, ", ev);
    }

    grRcc.push_back(new TGraph());//1, tc.peak_value, tc.coherent_filtered_peakHilbert));
    grRcc.back()->SetName(TString::Format("grRcc_%u", eventsOfInterest[0]));

    grPwg.push_back(new TGraph());//1, &tc.deconvolved_filtered_fracPowerWindowGradient, &tc.coherent_filtered_fracPowerWindowGradient));  
    grPwg.back()->SetName(TString::Format("grPwg_%u", eventsOfInterest[0]));

    grIm.push_back(new TGraph());//1, &tc.deconvolved_filtered_impulsivityMeasure, &tc.coherent_filtered_impulsivityMeasure));  
    grIm.back()->SetName(TString::Format("grIm_%u", eventsOfInterest[0]));

    grPh.push_back(new TGraph()); //1, &tc.deconvolved_filtered_peakHilbert, &tc.coherent_filtered_peakHilbert));
    grPh.back()->SetName(TString::Format("grPh_%u", eventsOfInterest[0]));

    grFd.push_back(new TGraph());
    grFd.back()->SetName(TString::Format("grFd_%u", eventsOfInterest[0]));


    hRccs.push_back( new TH2F(TString::Format("hRcc_%u", eventsOfInterest[0]),  titleBit + ";Map peak (no units); peakHilbert (mV)", 64, 0, 0.3, 64, 0, 80) );
    hPwgs.push_back( new TH2F(TString::Format("hPwg_%u", eventsOfInterest[0]),  titleBit + ";De-dispersed Power Window Gradient (ns); Power Window Gradient (ns)", 64, 0, 120, 64, 0, 120) );
    hIms.push_back( new TH2F(TString::Format("hIm_%u", eventsOfInterest[0]),   titleBit + ";De-dispersed Impulsivity;Impulsivity",  64, 0, 1, 64, 0, 1) );
    hPhs.push_back( new TH2F(TString::Format("hPh_%u", eventsOfInterest[0]),   titleBit + ";De-dispersed peakHilbert (mV); peakHilbert (mV)", 64, 0, 80, 64, 0, 80) );
    hFds.push_back( new TH1F(TString::Format("hFd_%u", eventsOfInterest[0]),   titleBit + ";Fisher disciminant (no units); Events per bin", 128, -30, 30) );

    lons.push_back(std::vector<double>());
    lats.push_back(std::vector<double>());
    alts.push_back(std::vector<double>());
    xs.push_back(std::vector<double>());
    ys.push_back(std::vector<double>());
    zs.push_back(std::vector<double>());

    // std::cout << "here " << eventsOfInterest[0] <<  std::endl;
    
    for(auto ev : eventsOfInterest){

      // std::cout << "here3 " << ev << std::endl;
      tc.getEvent(ev);
      // std::cout << "here4 " << ev << std::endl;
      
      grRcc.back()->SetPoint(grRcc.back()->GetN(), tc.peak_value, tc.coherent_filtered_peakHilbert);
      grPwg.back()->SetPoint(grPwg.back()->GetN(), tc.deconvolved_filtered_fracPowerWindowGradient, tc.coherent_filtered_fracPowerWindowGradient);
      grIm.back()->SetPoint(grIm.back()->GetN(), tc.deconvolved_filtered_impulsivityMeasure, tc.coherent_filtered_impulsivityMeasure);
      grPh.back()->SetPoint(grPh.back()->GetN(), tc.deconvolved_filtered_peakHilbert, tc.coherent_filtered_peakHilbert);
      grFd.back()->SetPoint(grFd.back()->GetN(), tc.fisherDiscriminant(), 1);

      lons.back().push_back(tc.longitude);
      lats.back().push_back(tc.latitude);
      lats.back().push_back(tc.altitude);

      double nuCart[nd];
      geom->getCartesianCoords(tc.latitude, tc.longitude, tc.altitude, nuCart);

      xs.back().push_back(nuCart[0]);
      ys.back().push_back(nuCart[1]);
      zs.back().push_back(nuCart[2]);

      // std::cout << "here5 " << ev << std::endl;      

      grNu->SetPoint(grNu->GetN(), tc.longitude, tc.latitude);

      // std::cout << "here6 " << ev << std::endl;      
    }
  }
  tc.addCut(ThermalTree::passAllQualityCuts + allEventsCut);
  
  const double dThresh = 40e3; //100e3; // meters
  const double dSqThresh = dThresh*dThresh;

  // std::cout  << "here2 " << std::endl;
  
  ProgressBar p(tc.N());  
  
  for(Long64_t entry=0;  entry < tc.N(); entry++){
    tc.getEntry(entry);
    hA->Fill(tc.longitude, tc.latitude);

    double cart[nd];
    geom->getCartesianCoords(tc.latitude, tc.longitude, tc.altitude, cart);

    for(int i=0; i < eventGroups.size(); i++){
      double minDsq = DBL_MAX;

      for(int j=0; j < eventGroups[i].size(); j++){
	const double nuCart[nd] = {xs[i][j], ys[i][j], zs[i][j]};
	double dSquared = dSq(cart, nuCart);
	minDsq = TMath::Min(dSquared, minDsq);
      }

      if(minDsq < dSqThresh){
	hB->Fill(tc.longitude,
		 tc.latitude);
	hRccs[i]->Fill(tc.peak_value,
		       tc.coherent_filtered_peakHilbert);
	hPwgs[i]->Fill(tc.deconvolved_filtered_fracPowerWindowGradient,
		       tc.coherent_filtered_fracPowerWindowGradient);
	hIms[i]->Fill(tc.deconvolved_filtered_impulsivityMeasure,
		      tc.coherent_filtered_impulsivityMeasure);
	hPhs[i]->Fill(tc.deconvolved_filtered_peakHilbert,
		      tc.coherent_filtered_peakHilbert);
	hFds[i]->Fill(tc.fisherDiscriminant());
      }
    }
    p.inc(entry);
  }

  for(UInt_t e=0; e< eventGroups.size(); e++){
    auto mainCan = new TCanvas();
    mainCan->Divide(1, 2);
    mainCan->cd(1);
    TPad* subPad = (TPad*)gPad;
    subPad->Divide(4);
    subPad->cd(1);
    hRccs[e]->Draw("colz");
    grRcc[e]->SetMarkerStyle(8);
    grRcc[e]->SetMarkerColor(kRed);
    grRcc[e]->Draw("psame");
    subPad->cd(2);
    hPwgs[e]->Draw("colz");
    grPwg[e]->SetMarkerStyle(8);
    grPwg[e]->SetMarkerColor(kRed);
    grPwg[e]->Draw("psame");
  
    subPad->cd(3);
    hIms[e]->Draw("colz");
    grIm[e]->SetMarkerStyle(8);
    grIm[e]->SetMarkerColor(kRed);
    grIm[e]->Draw("psame");
  
    subPad->cd(4);
    hPhs[e]->Draw("colz");
    grPh[e]->SetMarkerStyle(8);
    grPh[e]->SetMarkerColor(kRed);
    grPh[e]->Draw("psame");

    mainCan->cd(2);
    gPad->SetLogy(1);
    hFds[e]->Draw();
    grFd[e]->SetMarkerStyle(8);
    grFd[e]->SetMarkerColor(kRed);
    grFd[e]->Draw("psame");
  }    

  new TCanvas();
  hA->Draw("colz");
  hA->SetIcemask(true);
  gPad->SetLogz(1);
  grNu->Draw("psame");

  new TCanvas();
  hB->Draw("colz");
  hB->SetIcemask(true);
  gPad->SetLogz(1);
  grNu->Draw("psame");
  
  

}
