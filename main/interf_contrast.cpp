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
	std::string resultFile = "bin/data_out/Interferometry_Contrast_vs_Distance.png";
	std::vector<std::string> dataFiles = {
		"bin/data_in/Interf2_calha_d3p0mm.csv",
		"bin/data_in/Interf2_calha_d51p5mm.csv",
		"bin/data_in/Interf2_calha_d74p0mm.csv",
		"bin/data_in/Interf2_calha_d125p0mm.csv",
		"bin/data_in/Interf2_calha_d196p5mm.csv",
		"bin/data_in/Interf2_calha_d293p0mm.csv",
		"bin/data_in/Interf2_calha_d427p0mm.csv",
		"bin/data_in/Interf2_calha_d501p5mm.csv",
		"bin/data_in/Interf2_calha_d567p0mm.csv"
	};
	std::vector<double> distances = { // in mm
		3.0,
		51.5,
		74.0,
		125.0,
		196.5,
		293.0,
		427.0,
		501.5,
		567.0
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
	double maxContrast = 1;										// Y Axis maximum value
	double deltaDistances[2] = {-100.0, 600.0};					// X Axis maximum value
	int xNdiv = -514;											// Number of divisions in x axis, root notation
	int yNdiv = -510;											// Number of divisions in y axis, root notation
	int fNpoints = 1000;										// Number of points used in the display of the fitting function.
	float textLocation[4] = {100,0.7,400,0.97};					// Relative coordinates of the text, {x1rel,y1rel,x2rel,y2rel}
	int imageScaling = 140;

	// Fit Parameters
	double fit_a = 0.0075;
	double fit_c = -0.2;
	double fit_x0 = 0;
	double dfit_a = 0.01;
	double dfit_c = 0.5;
	double dfit_x0 = 50;

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
		// Print results of averages max and min and their errors.
		printf("Average Maximum: %.5f +/- %.5f [light value]   |   ", avgMax, maxErr);
		printf("Average Minimum: %.5f +/- %.5f [light value]   |   ", avgMin, minErr);
		// Find contrast
		k = (avgMax - avgMin) / (avgMax + avgMin);
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
		printf("(%.5f +/- %.5f)   ", distances[i], distancesErr[i]);
	}printf("\n");


	// Plot and Fit
	TCanvas* C = new TCanvas("C", "Canvas", 16*imageScaling, 9*imageScaling);
    TGraphErrors* g = new TGraphErrors(distances.size(), distances.data(), kVec.data(), distancesErr.data(), kErrVec.data());
    TF1* f = new TF1("f" , "abs(cos([0]*(x-[2]))) + [1]" , deltaDistances[0] , deltaDistances[1]);

	f->SetParameters(fit_a,fit_c,fit_x0);
	f->SetParLimits(0, fit_a - dfit_a, fit_a + dfit_a);
	f->SetParLimits(1, fit_c - dfit_c, fit_c + dfit_c);
	f->SetParLimits(2, fit_x0 - dfit_x0, fit_x0 + dfit_x0);

	
    g->SetTitle( "Contrast vs. #DeltaL" );
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kAzure+2);
    g->SetLineColor(kBlue+2);
    g->SetMarkerSize(1.0);

    g->GetXaxis()->SetTitle("Difference in Distance [mm]");
    g->GetXaxis()->SetLimits(deltaDistances[0], deltaDistances[1]);
    g->GetXaxis()->SetNdivisions(xNdiv);
    g->GetXaxis()->SetLabelSize(0.028);
    g->GetXaxis()->SetTickLength(-0.04);
    g->GetXaxis()->SetLabelOffset(0.033);
    g->GetXaxis()->SetTitleOffset(1.3);

    g->GetYaxis()->SetTitle("Contrast [unitless]");
    g->GetYaxis()->SetLabelSize(0.028);
    g->GetYaxis()->SetRangeUser(0., maxContrast);
    g->GetYaxis()->SetNdivisions(yNdiv);

	f->SetLineColor(kGreen);
    f->SetLineWidth(2);
	f->SetNpx(fNpoints);

    g->Fit("f");
	
    TPaveText* pt = new TPaveText(textLocation[0],textLocation[1],textLocation[2],textLocation[3],"user");

    pt->SetTextSize(0.035);
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);

    pt->AddText(Form("K = |cos(a(x-x_{0}))|^{2} + C"));
    pt->AddText(Form("a = %.5f %c %.5f [mm^{-1}]", f->GetParameter(0), 0xB1, f->GetParError(0)));
    pt->AddText(Form("C = (%.3f %c %.3f) [unitless]", f->GetParameter(1), 0xB1, f->GetParError(1)));
	pt->AddText(Form("x_{0} = %.0f %c %.0f [mm]", f->GetParameter(2), 0xB1, f->GetParError(2)));
    pt->AddText(Form("#chi^{2}/ndf = %.2f", float(f->GetChisquare()/f->GetNDF()) ));
    

    g->Draw("AP");
    pt->Draw("same");
    
    C->Update();
    C->SaveAs(resultFile.c_str());
    C->Clear();

	return 0;

}
