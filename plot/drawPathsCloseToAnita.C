void drawPathsCloseToAnita(const char* fileName = "doDataReduction_2017-10-18_08-41-54.root"){


  Acclaim::SummarySet ss(fileName);
  Long64_t n = ss.N();
  Acclaim::ProgressBar p(n);

  const double maxDistKm = 800;

  std::vector<TGraphAntarctica*> grs(BaseList::getNumPaths(), NULL);
  TGraphAntarctica* grFlightPath = new TGraphAntarctica();
  grFlightPath->SetName("grFlightPath");
  grFlightPath->SetMarkerColor(kRed);
  grFlightPath->SetMarkerStyle(2);
  // for(int i=0; i <  BaseList::getNumPaths(); i++){auto gr = new TGraphAntarctica(BaseList::getPath(i), 600); gr->SetLineColor((i%10)+1);  gr->SetMarkerColor((i%10)+1); cout << gr->GetN() << endl; gr->Draw("lp");

  for(Long64_t entry=0; entry < n; entry++){

    ss.getEntry(entry);
    AnitaEventSummary* sum = ss.summary();

    Adu5Pat pat = sum->anitaLocation.pat();
    UsefulAdu5Pat usefulPat(&pat);

    for(int i=0; i <  BaseList::getNumPaths(); i++){
      AntarcticCoord c = BaseList::getPath(i).getPosition(sum->realTime).as(AntarcticCoord::WGS84);
  // Double_t getDistanceFromSource(Double_t sourceLat, Double_t sourceLong, Double_t sourceAlt); ///< Gets distance from any source in meters

      if(c.x < -60){

	Double_t distKm = 1e-3*usefulPat.getDistanceFromSource(c.x, c.y, c.z);


	if(distKm < maxDistKm){

	  if(!grs[i]){
	    grs[i] = new TGraphAntarctica();
	    grs[i]->SetName(BaseList::getPath(i).getName());
	    grs[i]->SetTitle(BaseList::getPath(i).getName());
	    grs[i]->SetMarkerColor((i%10)+1);
	  }
	  grs[i]->SetPoint(grs[i]->GetN(), BaseList::getPath(i).getPosition(sum->realTime));
	}
      }
    }

    // if((entry%1000)==0){
    grFlightPath->SetPoint(grFlightPath->GetN(), pat.longitude, pat.latitude);
    // }

    p.inc(entry, n);
  }


  grFlightPath->Draw("p");
  for(int i=0; i < grs.size(); i++){
    if(grs[i]){
      grs[i]->Draw("p");
    }
  }

}

