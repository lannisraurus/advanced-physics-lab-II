// C++ Libs
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>

// ROOT- CERN Includes
#include "TApplication.h"
#include "TSystem.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TPaveText.h"

#include "TH1F.h"
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"

int main(){

	//---------------------------------------------------------------
	//                      INPUT PARAMETERS
	//---------------------------------------------------------------
	
	// Interferometer distances and aquired interference patterns
	std::string resultFile = "bin/data_out/Interferometry_Contrast_vs_Distance_PEER_REVIEW.png";
	std::vector<std::string> dataFiles = {
		"bin/data_in/deslocamento_11.5mm.csv",
		"bin/data_in/deslocamento_10.75cm.csv",	
		"bin/data_in/deslocamento_35.65cm.csv",
		"bin/data_in/deslocamento_37.95cm.csv",
		"bin/data_in/deslocamento_73.5mm.csv"
	};
	std::vector<double> distances = { // in mm
		11.5,
		107.5,
		356.5,
		379.5,
		73.50
	};
	std::vector<double> distancesErr(9,0.5); // in mm
	double L1 = 121.0; // mm
	double L1Err = 1.0; // mm
	double L2a = 73.0; // mm
	double L2aErr = 1.0; // mm
	double D = 7.0; // mm
	double DErr = 1.0; // mm
	double L2i = 30.0; // mm
	double L2iErr = 2.5; // mm
	
	// Max/Min extraction algorithm
	double minThreshold = 0.05;
	double maxThreshold = 0.05;

	// Display Settings (ROOT CERN)
	double maxContrast = 1.2;										// Y Axis maximum value
	double deltaDistances[2] = {-200.0, 1200.0};				// X Axis maximum value
	int xNdiv = -514;											// Number of divisions in x axis, root notation
	int yNdiv = -512;											// Number of divisions in y axis, root notation
	int fNpoints = 1000;										// Number of points used in the display of the fitting function.
	float textLocation[4] = {325,0.7,700,1.15};					// Relative coordinates of the text, {x1rel,y1rel,x2rel,y2rel}
	int imageScaling = 140;

	// Fit Parameters
	double fit_a = 0.0034;
	double fit_c = 0.872;
	double fit_x0 = 0;
	double dfit_a = 0.0020;
	double dfit_c = 0.2;
	double dfit_x0 = 500;

	//---------------------------------------------------------------
	//                           ALGORITHM
	//---------------------------------------------------------------

	// Data Vars
	std::ifstream input;
		
	// Probing Vars
	double pixel, value;
	std::stringstream ss;
	double v[3] = {0.0, 0.0, 0.0};
	std::string line;

	// Store Vars
	std::vector<double> pixels, values, maximums, minimums;
	double avgMax = 0, avgMin = 0,
		   maxErr = 0, minErr = 0,
		   maxErrTemp = 0, minErrTemp = 0,
		   k = 0, kErr = 0,
		   averageValue = 0;
	double background = 0;
	bool take_background = true;

	// Useful vars
	std::vector<double> kVec;
	std::vector<double> kErrVec;
	
	// Get all K's and Kerrs
	for (const auto &dataFile : dataFiles) {
		// Clear vars
		pixels.clear(); values.clear(); maximums.clear(); minimums.clear();
		k = 0; kErr = 0; avgMax = 0; avgMin = 0; maxErr = 0; minErr = 0; maxErrTemp = 0; minErrTemp = 0;
		line = ""; pixel = 0; value = 0, averageValue = 0;
		// Open file
		input.open(dataFile);
		// Error Checking
		if(!input.is_open()){
			printf("Error: Couldn't open file.\n");
			return -1;
		}
		// Get data inside file
		std::getline(input, line); // Get labels (useless)
		while(!input.eof()){
			std::getline(input, line);
			std::replace(line.begin(), line.end(), ',', ' ');
			ss << line;
			ss >> pixel >> value;
			pixels.push_back(pixel);
			values.push_back(value);
			ss.clear();
		}
		// Find average of values
		for (int i = 0; i < pixels.size(); i++){
			averageValue += values[i];
		}
		averageValue /= pixels.size();
		// Maximum and minimum finder (ignores values at boundaries)
		for (int i = 1; i < pixels.size() - 1; i++){
			v[0] = values[i-1];
			v[1] = values[i];
			v[2] = values[i+1];
			if (v[1] > v[0] && v[1] > v[2] && v[1] > (1+maxThreshold)*averageValue ) maximums.push_back(v[1]);
			if (v[1] < v[0] && v[1] < v[2] && v[1] < (1-minThreshold)*averageValue ) minimums.push_back(v[1]);
		}
		// Find average max
		for (int i = 0; i < maximums.size(); i++){
			avgMax += maximums[i];
		}		avgMax /= maximums.size();
		// Find average min
		for (int i = 0; i < minimums.size(); i++){
			avgMin += minimums[i];
		}	avgMin /= minimums.size();
		// Find average max error
		for (int i = 0; i < maximums.size(); i++){
			maxErrTemp = std::abs(avgMax - maximums[i]);
			if (maxErrTemp > maxErr) maxErr = maxErrTemp;
		}
		// Find average min error
		for (int i = 0; i < minimums.size(); i++){
			minErrTemp = std::abs(avgMin - minimums[i]);
			if (minErrTemp > minErr) minErr = minErrTemp;
		}
		if(take_background){
			background = avgMin;
			take_background = false;
		}
		// Print results of averages max and min and their errors.
		printf("Average Maximum: %.5f +/- %.5f [light value]   |   ", avgMax, maxErr);
		printf("Average Minimum: %.5f +/- %.5f [light value]   |   ", avgMin, minErr);
		// Find contrast
		k = (avgMax - avgMin) / (avgMax + avgMin - 2*background);
		kErr = 2 * ( 1 / pow(avgMax + avgMin, 2) ) * (avgMin*maxErr + avgMax*minErr);
		printf("Contrast: %.5f +/- %.5f\n", k, kErr);
		// Close File
		input.close();
		// Add K and KERR
		kVec.push_back(k);
		kErrVec.push_back(kErr);
	}

	// Transform distances
	printf("\nDistances: ");
	for (int i = 0; i < distances.size(); i++){
		distances[i] = L2a + L2i + distances[i] + D - L1;
		distancesErr[i] = L2aErr + L2iErr + distancesErr[i] + DErr + L1Err;
		distances[i] *= 2;
		distancesErr[i] *= 2;
		printf("(%.5f +/- %.5f)   ", distances[i], distancesErr[i]);
	}printf("\n");


	// Plot and Fit
	TCanvas* C = new TCanvas("C", "Canvas", 16*imageScaling, 9*imageScaling);
    TGraphErrors* g = new TGraphErrors(distances.size(), distances.data(), kVec.data(), distancesErr.data(), kErrVec.data());
    TF1* f = new TF1("f" , "[1]*abs(cos([0]*(x-[2])))" , deltaDistances[0] , deltaDistances[1]);
	TF1* f2 = new TF1("f2" , "[1]*abs(cos([0]*(x-[2])))" , deltaDistances[0] , deltaDistances[1]);

	f->SetParameters(fit_a,fit_c,fit_x0);
	f2->SetParameters(0.0034,0.872,4);
	f->SetParLimits(0, fit_a - dfit_a, fit_a + dfit_a);
	f->SetParLimits(1, fit_c - dfit_c, fit_c + dfit_c);
	f->SetParLimits(2, fit_x0 - dfit_x0, fit_x0 + dfit_x0);
	//f->FixParameter(1,1);
	
    g->SetTitle( "Contrast vs. #DeltaL" );
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kAzure+2);
    g->SetLineColor(kBlue+2);
    g->SetMarkerSize(1.4);

    g->GetXaxis()->SetTitle("Difference in Distance [mm]");
    g->GetXaxis()->SetLimits(deltaDistances[0], deltaDistances[1]);
    g->GetXaxis()->SetNdivisions(xNdiv);
    g->GetXaxis()->SetLabelSize(0.028);
    g->GetXaxis()->SetTickLength(-0.04);
    g->GetXaxis()->SetLabelOffset(0.033);
    g->GetXaxis()->SetTitleOffset(1.3);

    g->GetYaxis()->SetTitle("Corrected Contrast [unitless]");
    g->GetYaxis()->SetLabelSize(0.028);
    g->GetYaxis()->SetRangeUser(0., maxContrast);
    g->GetYaxis()->SetNdivisions(yNdiv);

	f->SetLineColor(kGreen);
	f2->SetLineColor(kRed);
    f->SetLineWidth(2);
	f->SetNpx(fNpoints);

    g->Fit("f");
	
    TPaveText* pt = new TPaveText(textLocation[0],textLocation[1],textLocation[2],textLocation[3],"user");

    pt->SetTextSize(0.035);
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);

    pt->AddText(Form("K = A|cos(a(x-x_{0}))|"));
    pt->AddText(Form("a = %.5f %c %.5f [mm^{-1}]", f->GetParameter(0), 0xB1, f->GetParError(0)));
    pt->AddText(Form("A = (%.3f %c %.3f) [unitless]", f->GetParameter(1), 0xB1, f->GetParError(1)));
	pt->AddText(Form("x_{0} = %.0f %c %.0f [mm]", f->GetParameter(2), 0xB1, f->GetParError(2)));
    pt->AddText(Form("#chi^{2}/ndf = %.1f", float(f->GetChisquare()/f->GetNDF()) ));
    

    g->Draw("AP");
	f2->Draw("same");
    pt->Draw("same");
    
    C->Update();
    C->SaveAs(resultFile.c_str());
    C->Clear();

	return 0;

}
