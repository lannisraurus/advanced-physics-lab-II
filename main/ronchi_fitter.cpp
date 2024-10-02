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
	// Data file info
	std::string dataFile = "data_in/RONCHI1_d8p15mm_s0p061ms.csv";			// Input data file, csv
	std::string resultFile = "data_out/RONCHI1_d8p15mm_s0p061ms.png";	// Output fata file, png
	// Physical System
	double pixelToDist = ((3.69e-6)/0.5055);					// Conversion constant from pixels to real life distance, in metres 
	double pixelToDist_err = (3.69e-6)*0.0073/(0.5055*0.5055);	// Error of the conversion factor, in metres
	double pixelMeasureOffset = 0.007;							// Where is the center of the measure, to propagate error? in metres
	double laser_wavelength = 633e-9;							// Wavelength of the light, in metres
	double focal_length = 250e-3;								// Focal length of the fourier transform lens, in metres
	double lambdaf = laser_wavelength*focal_length;
	// Data filtering
	int noiseThreshold = 1010;
	// Fit initial parameters	
	double fit_E_0 = 6000*pixelToDist/lambdaf;
	double fit_x0 = 960*pixelToDist;
	double fit_C = 900;
	double fit_a = 0.0402/pixelToDist*lambdaf;
	double fit_d = 0.0115/pixelToDist*lambdaf;
	double fit_m = 10;
	// Fitting Freedom
	double dfit_E_0 = 2000*pixelToDist/lambdaf;
	double dfit_x0 = 200*pixelToDist;
	double dfit_C = 160;
	double dfit_a = 0.0005/pixelToDist*lambdaf;
	double dfit_d = 0.0005/pixelToDist*lambdaf;
	double dfit_m = 5;
	// Display Settings
	int maxLightValue = 85000;									// Y Axis maximum value
	double maxDistance = 1.5e-2;									// X Axis maximum value
	int xNdiv = -512;											// Number of divisions in x axis, root notation
	int yNdiv = -610;											// Number of divisions in y axis, root notation
	int fNpoints = 1000;										// Number of points used in the display of the fitting function.
	float textLocation[4] = {0.6,0.6,0.95,0.95};					// Relative coordinates of the text, {x1rel,y1rel,x2rel,y2rel}
	int imageScaling = 140;

	//---------------------------------------------------------------
	//                           ALGORITHM
	//---------------------------------------------------------------

	// Data Vars
	std::ifstream input;
	std::stringstream ss;
	input.open(dataFile);
	std::vector<double> X;
	std::vector<double> Y;
	std::vector<double> eX;
	std::vector<double> eY;
	
	// Probing Vars
	double pixel;
	double value;
	std::string line;
	
	// Error Checking
	if(!input.is_open()){
		printf("Error: Couldn't open file.\n");
		return -1;
	}

	// Retrieve data, calibration and errors
	std::getline(input, line);
	for (line; std::getline(input, line); ) {
		std::replace( line.begin(), line.end(), ',', ' ');
		ss << line;
		ss >> pixel >> value;
		ss.clear();
		if(value>noiseThreshold){
			X.push_back(pixel*pixelToDist);
			Y.push_back(value);
			eX.push_back( std::abs((pixel-pixelMeasureOffset/pixelToDist)*pixelToDist_err) );
			eY.push_back(0);	
		}	
	}
    
	// Plot and Fit
	TCanvas* C = new TCanvas("C", "Canvas", 16*imageScaling, 9*imageScaling);
    TGraphErrors* g = new TGraphErrors(X.size(), X.data(), Y.data(), eX.data(), eY.data());

    g->SetTitle("Ronchi Grid Fit");
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kAzure+2);
    g->SetLineColor(kBlue+2);
    g->SetMarkerSize(0.5);

    g->GetXaxis()->SetTitle("Distance [m]");
    g->GetXaxis()->SetLimits(0., maxDistance);
    g->GetXaxis()->SetNdivisions(xNdiv);
    g->GetXaxis()->SetLabelSize(0.028);
    g->GetXaxis()->SetTickLength(-0.04);
    g->GetXaxis()->SetLabelOffset(0.033);
    g->GetXaxis()->SetTitleOffset(1.3);

    g->GetYaxis()->SetTitle("Intensity [light value]");
    g->GetYaxis()->SetLabelSize(0.028);
    g->GetYaxis()->SetRangeUser(0., maxLightValue);
    g->GetYaxis()->SetNdivisions(yNdiv);

    TF1* f = new TF1("f" , "pow( [0]*sin(TMath::Pi()*[1]*(x-[2])/[6]) / (TMath::Pi()*(x-[2])/[6]) * sin(TMath::Pi()*[4]*[5]*(x-[2])/[6]) / sin(TMath::Pi()*(x-[2])*[4]/[6]) , 2) + [3]" , 0 , maxDistance); 

    f->SetLineColor(kBlack);
    f->SetLineWidth(2);
	
	f->SetParameters(fit_E_0, fit_a, fit_x0, fit_C, fit_d, fit_m, lambdaf);	
	f->SetParLimits(0, fit_E_0 - dfit_E_0, fit_E_0 + dfit_E_0);
	f->SetParLimits(1, fit_a - dfit_a, fit_a + dfit_a);
	f->SetParLimits(2, fit_x0 - dfit_x0, fit_x0 + dfit_x0);
	f->SetParLimits(3, fit_C - dfit_C, fit_C + dfit_C);
	f->SetParLimits(4, fit_d - dfit_d, fit_d + dfit_d);
	f->SetParLimits(5, fit_m - dfit_m, fit_m + dfit_m);	
	f->FixParameter(6,lambdaf);
    f->SetNpx(fNpoints);

    g->Fit("f");

    TPaveText* pt = new TPaveText(textLocation[0]*maxDistance, textLocation[1]*maxLightValue, textLocation[2]*maxDistance, textLocation[3]*maxLightValue, "user");

    pt->SetTextSize(0.025);
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);

    pt->AddText(Form("I = |E_{0}|^{2} (sin(#pia#nu)/(#pi#nu))^{2} (sin(#pimd#nu)/sin(#pid#nu))^{2}"));
    pt->AddText(Form("E_{0} = %.0f %c %.0f [light value]", f->GetParameter(0), 0xB1, f->GetParError(0)));
    pt->AddText(Form("a = %.6f %c %.6f [mm]", f->GetParameter(1)*1e3, 0xB1, f->GetParError(1)*1e3));
	pt->AddText(Form("x_{0} = %.5f %c %.5f [mm]", f->GetParameter(2)*1e3, 0xB1, f->GetParError(2)*1e3));
	pt->AddText(Form("C = %.0f %c %.0f [light value]", f->GetParameter(3), 0xB1, f->GetParError(3)));
	pt->AddText(Form("d = %.6f %c %.6f [mm]", f->GetParameter(4)*1e3, 0xB1, f->GetParError(4)*1e3));
	pt->AddText(Form("m = %.5f %c %.5f [unitless]", f->GetParameter(5), 0xB1, f->GetParError(5)));

    pt->AddText(Form("#chi^{2}/ndf = %.2f", float(f->GetChisquare()/f->GetNDF()) ));
     
    g->Draw("AP");
    pt->Draw("same");
    
    C->Update();
    C->SaveAs(resultFile.c_str());
    C->Clear();
    
	return 0;
}
