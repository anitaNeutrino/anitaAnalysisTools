#include "UsefulAdu5Pat.h"

void drawHorizonVsAltitude(){

  Adu5Pat anita;
  anita.longitude = 90;
  anita.latitude = -70;
  anita.altitude = 37e3;
  // anita.longitude = 0;
  // anita.latitude = -90;
  // anita.altitude = 37e3;

  const double feetToKm = 12*2.54/100/1000; // feet to inches to cm to m to km

  // c130 hercules ceiling: http://www.af.mil/About-Us/Fact-Sheets/Display/Article/104517/c-130-hercules/
  Double_t c130HercMaxAlt = 28000*feetToKm;
  // Twin otter ceiling: https://en.wikipedia.org/wiki/De_Havilland_Canada_DHC-6_Twin_Otter
  Double_t twinOtterMaxAlt = 25000*feetToKm;

  Double_t anitaEasting, anitaNorthing;
  RampdemReader::LonLatToEastingNorthing((double)anita.longitude, (double)anita.latitude, anitaEasting, anitaNorthing);

  const int nBins = 128;
  const double maxFlightAlt = 30;
  const double maxDeltaNorthing = 3000;

  TGraph* grHerc = new TGraph();
  grHerc->SetPoint(0, c130HercMaxAlt, 0);
  grHerc->SetPoint(1, c130HercMaxAlt, maxDeltaNorthing);
  grHerc->SetLineColor(kRed);
  TGraph* grTwin = new TGraph();
  grTwin->SetPoint(0, twinOtterMaxAlt, 0);
  grTwin->SetPoint(1, twinOtterMaxAlt, maxDeltaNorthing);
  grTwin->SetLineColor(kBlue);
  auto l = new TLegend(0.67, 0.67, 0.88, 0.88);//0.8, 0.8, 1, 1);
  l->AddEntry(grHerc, "C130 Hercules ceiling", "l");
  l->AddEntry(grTwin, "Twin otter ceiling", "l");

  UsefulAdu5Pat usefulPat(&anita);

  TProfile2D* hIntersects = new TProfile2D("hIntersects",
					   "Does theta/phi from the flight intersect the ground?; Flight altitude (km); Flight #delta(northing) (km); Ray intersects ground?",
					   nBins, 0, maxFlightAlt,
					   nBins, 0, maxDeltaNorthing);
  TProfile2D* hFlightRecoSeparation = new TProfile2D("hFlightRecoSeparation",
						    "If it does, how far is the intersection from the point directly below the flight?; Flight altitude (km); Flight #delta(northing) (km); Reco-flight-separation (km)",
						    nBins, 0, maxFlightAlt,
						    nBins, 0, maxDeltaNorthing);
  TProfile2D* hThetaDeg = new TProfile2D("hThetaDeg",
					 "Elevation angle (Degrees); Flight altitude (km); Flight #delta(northing) (km); #theta (Degrees)",
					 nBins, 0, maxFlightAlt,
					 nBins, 0, maxDeltaNorthing);

  for(int by=1; by<=hIntersects->GetNbinsY(); by++){
    for(int bx=1; bx<=hIntersects->GetNbinsX(); bx++){
      Double_t deltaNorthKm = hIntersects->GetYaxis()->GetBinCenter(by);

      Double_t flightNorthing = anitaNorthing + 1000*deltaNorthKm; // m?
      Double_t flightEasting = anitaEasting;

      Double_t flightLon = 0;
      Double_t flightLat = 0;
      RampdemReader::EastingNorthingToLonLat(flightEasting, flightNorthing, flightLon, flightLat);
      Double_t flightAlt = 1000*hIntersects->GetXaxis()->GetBinCenter(bx);

      double separation = deltaNorthKm;//1e-3*usefulPat.getDistanceFromSource(flightLat, flightLon, flightAlt);

      Double_t theta, phi;
      usefulPat.getThetaAndPhiWave(flightLon, flightLat, flightAlt, theta, phi);

      hThetaDeg->Fill(flightAlt/1000, separation, -1*theta*TMath::RadToDeg());

      Double_t recoLon = 0;
      Double_t recoLat = 0;
      Double_t recoAlt = 0;
      Double_t recoThetaAdjust = 0;
      const double maxThetaAdjustment = 0;
      const Int_t maxIter = 5;
      int retVal = usefulPat.traceBackToContinent(phi, theta, &recoLon, &recoLat, &recoAlt, &recoThetaAdjust, maxThetaAdjustment, maxIter);
      const char* traceReturnMeaning[] = {"never hits the ground, even with maximum adjustment",
					  "hits the ground with no adjustment",
					  "hits the ground with adjustment"};
      // std::cout << flightAlt << "\t" << flightLon << "\t" << flightLat << std::endl;
      // std::cout << traceReturnMeaning[retVal] << "\t" << recoLon << "\t" << recoLat << std::endl;

      Int_t intersectsWithContinent = retVal == 1 ? 1 : 0;
      hIntersects->Fill(flightAlt/1000, separation, retVal);//intersectsWithContinent);


      if(intersectsWithContinent){

	Adu5Pat flightShadow;
	flightShadow.latitude = flightLat;
	flightShadow.longitude = flightLon;
	flightShadow.altitude = 0;
	UsefulAdu5Pat usefulFlightShadow(&flightShadow);
	Double_t flightRecoSeparation = 1e-3*usefulFlightShadow.getDistanceFromSource(recoLat, recoLon, recoAlt);

	hFlightRecoSeparation->Fill(flightAlt/1000, separation, flightRecoSeparation);
      }
    }
  }
  auto c1 = new TCanvas();
  hIntersects->Draw("colz");
  grHerc->Draw("lsame");
  grTwin->Draw("lsame");
  l->Draw();
  auto c2 = new TCanvas();
  hFlightRecoSeparation->Draw("colz");
  grHerc->Draw("lsame");
  grTwin->Draw("lsame");
  l->Draw();
  c2->SetLogz(1);
  auto c3 = new TCanvas();
  hThetaDeg->Draw("colz");
  grHerc->Draw("lsame");
  grTwin->Draw("lsame");
  l->Draw();



}
