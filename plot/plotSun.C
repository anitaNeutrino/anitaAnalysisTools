#include <cmath>
using namespace Acclaim;

void plotSun(const char* glob = "~/ANITA/anita3Analysis/sineSub100Percent/data/doSineSub_all_*.root"){

  SummarySet ss(glob);

  std::cout << ss.N() << std::endl;
  if(ss.N()==0){
    std::cerr << "No AnitaEventSummary files?" << std::endl;
    return;
  }
  ss.SetUseProof(1);

  bool justMinBias = false; true;
  const double maxSunTheta = 10;

  for(int pol=AnitaPol::kHorizontal; pol < AnitaPol::kNotAPol; pol++){
  
    TString phi   = TString::Format("sum.peak.phi[%d][0]", pol);
    TString theta = TString::Format("sum.peak.theta[%d][0]", pol);
  
    TString dPhi = "";
    dPhi += "((sum.sun.phi - " + phi + " > 180)*(sum.sun.phi - " + phi + " - 360) +";
    dPhi +=  "(sum.sun.phi - " + phi + " < -180)*(sum.sun.phi - " + phi + " + 360) +";
    dPhi += "(fabs(sum.sun.phi - " + phi + ") < 180)*(sum.sun.phi - " + phi + "))";


    TString dTheta = theta + "+sum.sun.theta";
  
    TString closeToSun = justMinBias ? "sum.flags.isRF==0"  : "sum.flags.isRF>=0";
    closeToSun += TString::Format(" && -sum.sun.theta < %lf", maxSunTheta);
    closeToSun += " && TMath::Abs(" +  dTheta + ") < 5";
    closeToSun += " && TMath::Abs("  + dPhi + ") < 5";
    
    const int nBinsTheta = justMinBias ? 16 : 128;
    const int nBinsPhi = justMinBias ? 256 : 2048;
    
    auto n_dTheta_theta = TString::Format("h_dTheta_theta_%d(%d,0,%lf,1024,-5,5)", pol,nBinsTheta, maxSunTheta);
    // auto h_dTheta_theta = new TH2D(n_dTheta_theta, "h_dTheta_theta", nBinsTheta, 0, 30,  1024, -5, 5);
    
    auto c1 = new TCanvas();
    ss.Draw(dTheta + ":-sum.sun.theta>>" + n_dTheta_theta, closeToSun, "colz");
    auto h_dTheta_theta = (TH2D*) ss.getDrawOutput();

    auto n_dTheta_phi = TString::Format("h_dTheta_phi_%d", pol);
    auto h_dTheta_phi = new TH2D(n_dTheta_phi,   "h_dTheta_phi",   nBinsPhi, -180, 180,  1024, -5, 5);
    auto c2 = new TCanvas();
    ss.Draw(dTheta + ":sum.sun.phi>>" + n_dTheta_phi,      closeToSun, "colz");

    auto n_dPhi_theta = TString::Format("h_dPhi_theta_%d", pol);
    auto h_dPhi_theta = new TH2D(n_dPhi_theta,   "h_dPhi_theta",   nBinsTheta, 0, maxSunTheta,  1024, -5, 5);    
    auto c3 = new TCanvas();
    ss.Draw(dPhi + " :-sum.sun.theta>>" + n_dPhi_theta,     closeToSun, "colz");

    auto n_dPhi_phi = TString::Format("h_dPhi_phi_%d", pol);    
    auto h_dPhi_phi  = new TH2D(n_dPhi_phi,     "h_dPhi_phi",     nBinsPhi, -180, 180,  1024, -5, 5);  
    auto c4 = new TCanvas();
    ss.Draw(dPhi + ":sum.sun.phi>>" + n_dPhi_phi,          closeToSun, "colz");

    delete c1;
    delete c2;
    delete c3;
    delete c4;
  
    TCanvas* c5 = new TCanvas();
    c5->Divide(2, 2);
    c5->cd(1);
    h_dTheta_theta->SetTitle(";Sun #theta (degrees); #delta#theta sun (Degrees)");
    h_dTheta_theta->Draw("colz");

    c5->cd(2);
    h_dTheta_phi->SetTitle(";Sun #phi (degrees); #delta#theta sun (Degrees)");
    h_dTheta_phi->Draw("colz");
  
    c5->cd(3);
    h_dPhi_theta->SetTitle(";Sun #theta (degrees); #delta#phi sun (Degrees)");
    h_dPhi_theta->Draw("colz");

    c5->cd(4);
    h_dPhi_phi->SetTitle(";Sun #phi (degrees); #delta#phi sun (Degrees)");
    h_dPhi_phi->Draw("colz");
  
  
    auto h_dTheta_theta_px = h_dTheta_theta->ProfileX();
    auto h_dTheta_phi_px   = h_dTheta_phi->ProfileX(); 
    auto h_dPhi_theta_px = h_dPhi_theta->ProfileX();
    auto h_dPhi_phi_px   = h_dPhi_phi->ProfileX();

    static auto c6 = new TCanvas();
    if(pol==0){
      c6->Divide(2, 2);
    }
    TString opt = pol == 0 ? "colz" : "same";
    EColor lc = pol == 0 ? kBlack : kBlue;
    c6->cd(1);
    h_dTheta_theta_px->SetTitle(";Sun #theta (degrees); #delta#theta sun (Degrees)");
    h_dTheta_theta_px->Draw(opt);
    h_dTheta_theta_px->SetLineColor(lc);

    c6->cd(2);
    h_dTheta_phi_px->SetTitle(";Sun #phi (degrees); #delta#theta sun (Degrees)");
    h_dTheta_phi_px->Draw(opt);
    h_dTheta_phi_px->SetLineColor(lc);
  
    c6->cd(3);
    h_dPhi_theta_px->SetTitle(";Sun #theta (degrees); #delta#phi sun (Degrees)");
    h_dPhi_theta_px->Draw(opt);
    h_dPhi_theta_px->SetLineColor(lc);
  
    c6->cd(4);
    h_dPhi_phi_px->SetTitle(";Sun #phi (degrees); #delta#phi sun (Degrees)");
    h_dPhi_phi_px->Draw(opt);
    h_dPhi_phi_px->SetLineColor(lc);
    
  }
}
