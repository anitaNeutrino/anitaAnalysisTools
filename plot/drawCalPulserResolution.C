void drawCalPulserResolution(const char* fileName = "doDataReductionCalPulser_*.root"){

  TCut goodWaisEvent = "run >= 331 && run <= 354 && flags.pulser==1";
  TCut goodLdbEvent = "run < 200 && flags.pulser==2";
  TCut goodLdbVPol = goodLdbEvent + "(run < 150 || run >= 154)";
  TCut goodLdbHPol = goodLdbEvent + "(run >= 150 || run < 154)";  
  
  Acclaim::SummarySet ss(fileName);
  // cout << ss.N() << endl;
  const int nBinsAngle = 32;
  const int maxSnr = 14;

  auto c1 = new TCanvas();  
  TH2D* hWaisDPhi = new TH2D("hWaisDPhi", "", maxSnr, 0, maxSnr, nBinsAngle, -3, 3);
  ss.Draw("sum.trainingPeak().dPhiWais():sum.trainingDeconvolvedFiltered().snr>>hWaisDPhi", goodWaisEvent, "colz");
  hWaisDPhi->FitSlicesY();
  TH1D* hSigmaWaisPhi = (TH1D*) gDirectory->FindObject("hWaisDPhi_2");
  
  // c1->SetLogy(1);

  auto c2 = new TCanvas();
  TH2D* hWaisDTheta = new TH2D("hWaisDTheta", "", maxSnr, 0, maxSnr, nBinsAngle, -1.5, 1.5);
  ss.Draw("sum.trainingPeak().dThetaWais():sum.trainingDeconvolvedFiltered().snr>>hWaisDTheta", goodWaisEvent, "colz");
  hWaisDTheta->FitSlicesY();
  TH1D* hSigmaWaisTheta = (TH1D*) gDirectory->FindObject("hWaisDTheta_2");




  // TF1* f1 = new TF1("f1", "exp([0]*x + [1])", 0, maxSnr);
  TF1* f1 = new TF1("f1", "exp([0]*x + [1]) + [2]", 0, maxSnr);
  f1->SetParameter(0, -1);
  // f1->SetParLimits(0, 0, 2);  
  f1->SetParameter(1, 0);
  f1->SetParameter(2, 0.1);
  // f1->SetParameter(2, 1);
  {
    auto cTemp = new TCanvas();
    hSigmaWaisPhi->Fit(f1);
    hSigmaWaisTheta->Fit(f1);
    delete cTemp;
    
  }
  
  // auto c3 = new TCanvas();
  // TH2D* hLdbDPhi = new TH2D("hLdbDPhi", "", maxSnr, 0, maxSnr, nBinsAngle, -3, 3);
  // ss.Draw("peak[1][0].dPhiLDB():deconvolved_filtered[1][0].snr>>hLdbDPhi", goodLdbVPol, "colz");
  // hLdbDPhi->FitSlicesY();
  // TH1D* hSigmaLdbPhi = (TH1D*) gDirectory->FindObject("hLdbDPhi_2");
  
  // // c1->SetLogy(1);

  // auto c4 = new TCanvas();
  // TH2D* hLdbDTheta = new TH2D("hLdbDTheta", "", maxSnr, 0, maxSnr, nBinsAngle, -1.5, 1.5);
  // ss.Draw("sum.peak[1][0].dThetaLDB():sum.deconvolved_filtered[1][0].snr>>hLdbDTheta", goodLdbVPol, "colz");
  // hLdbDTheta->FitSlicesY();
  // TH1D* hSigmaLdbTheta = (TH1D*) gDirectory->FindObject("hLdbDTheta_2");



  

  // hSigmaWaisPhi->Draw();
  // hSigmaWaisTheta->Draw("same");
  hSigmaWaisPhi->SetLineColor(kRed);
  hSigmaWaisTheta->SetLineColor(kBlue);
  ((TF1*)(hSigmaWaisPhi->GetListOfFunctions()->At(0)))->SetLineColor(hSigmaWaisPhi->GetLineColor());
  ((TF1*)(hSigmaWaisTheta->GetListOfFunctions()->At(0)))->SetLineColor(hSigmaWaisTheta->GetLineColor());
  hSigmaWaisTheta->SetLineColor(kBlue);
  
  TH1D* hs[2] = {hSigmaWaisPhi, hSigmaWaisTheta};  
  auto c5 = Acclaim::RootTools::drawHistsWithStatsBoxes(2, hs, "", "");
  hSigmaWaisPhi->SetTitle("Fitted angular resolution vs. SNR_{dedispersed}; SNR_{dedispersed}; Angular resolution (degrees)");
  
  // hSigmaLdbPhi->Draw("same");
  // hSigmaLdbTheta->Draw("same");
  // hSigmaLdbPhi->SetLineColor(kMagenta);
  // hSigmaLdbTheta->SetLineColor(kCyan);

  auto l5 = new TLegend();
  l5->AddEntry(hSigmaWaisPhi,   "WAIS #sigma_{#phi}"  , "l");
  l5->AddEntry(hSigmaWaisTheta, "WAIS #sigma_{#theta}", "l");  
  // l5->AddEntry(hSigmaLdbPhi,   "LDB #sigma_{#phi}"  , "l");
  // l5->AddEntry(hSigmaLdbTheta, "LDB #sigma_{#theta}", "l");  
  l5->Draw();
  
  // auto c= ss.getChain();
  // c->Scan("run:eventNumber:sum.trainingPeak().dPhiWais():sum.trainingPeak().dThetaWais()", goodWaisEvent + "TMath::Abs(sum.trainingPeak().dPhiWais())>6");
  
}
