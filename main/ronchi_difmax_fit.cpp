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
	std::string dataFile = "bin/data_in/RONCHI1_d8p15mm_s32p8ms_diffraction_maximums.csv";		// Input data file, csv
	std::string resultFile = "bin/data_out/RONCHI1_d8p15mm_s32p8ms_diffraction_maximums.png";	// Output fata file, png
	// Physical System
	double pixelToDist = ((3.69e-6)/0.5055);					// Conversion constant from pixels to real life distance, in metres 
	double pixelToDist_err = (3.69e-6)*0.0073/(0.5055*0.5055);	// Error of the conversion factor, in metres
	double laser_wavelength = 633e-9;							// Wavelength of the light, in metres
	double focal_length = 250e-3;								// Focal length of the fourier transform lens, in metres
	double lambdaf = laser_wavelength*focal_length;
	// Display Settings
	double maxY = 0.006;								// Y Axis maximum value
	double maxX = 10;								// X Axis maximum value
	int xNdiv = 10;									// Number of divisions in x axis, root notation
	int yNdiv = -512;									// Number of divisions in y axis, root notation
	int fNpoints = 1000;								// Number of points used in the display of the fitting function.
	float textLocation[4] = {0.65,0.1,0.95,0.5};	// Relative coordinates of the text, {x1rel,y1rel,x2rel,y2rel}
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
	double n;
	double d;
	std::string line;
	
	// Error Checking
	if(!input.is_open()){
		printf("Error: Couldn't open file.\n");
		return -1;
	}

	// Retrieve data, calibration and errors
	for (line; std::getline(input, line); ) {

		ss << line;
		ss >> n >> d;
		ss.clear();

		X.push_back(n);
		Y.push_back(d*pixelToDist/2.);
		eX.push_back(0);
		eY.push_back(d*pixelToDist_err/2.);

		printf("(%f,%f)\n",X.back(),Y.back());

	}
    
	// Plot and Fit
	TCanvas* C = new TCanvas("C", "Canvas", 16*imageScaling, 9*imageScaling);
    TGraphErrors* g = new TGraphErrors(X.size(), X.data(), Y.data(), eX.data(), eY.data());

    g->SetTitle( (dataFile+" - Diffraction Maxs. Linear Fit").c_str() );
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kAzure+2);
    g->SetLineColor(kBlue+2);
    g->SetMarkerSize(0.5);

    g->GetXaxis()->SetTitle("Diffraction Index [unitless]");
    g->GetXaxis()->SetLimits(0., maxX);
    g->GetXaxis()->SetNdivisions(xNdiv);
    g->GetXaxis()->SetLabelSize(0.028);
    g->GetXaxis()->SetTickLength(-0.04);
    g->GetXaxis()->SetLabelOffset(0.033);
    g->GetXaxis()->SetTitleOffset(1.3);

    g->GetYaxis()->SetTitle("Distance [m]");
    g->GetYaxis()->SetLabelSize(0.028);
    g->GetYaxis()->SetRangeUser(0., maxY);
    g->GetYaxis()->SetNdivisions(yNdiv);

    TF1* f = new TF1("f" , "([0]/[1])*x + [2]" , 0 , maxX); 

    f->SetLineColor(kBlack);
    f->SetLineWidth(2);
	
	f->SetParameters(lambdaf,0,0);	
	f->FixParameter(0,lambdaf);
    f->SetNpx(fNpoints);

    g->Fit("f");

    TPaveText* pt = new TPaveText(textLocation[0]*maxX, textLocation[1]*maxY, textLocation[2]*maxX, textLocation[3]*maxY, "user");

    pt->SetTextSize(0.032);
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);

    pt->AddText(Form("x = (#lambdaf/d)n + b"));
    pt->AddText(Form("d = %.4f %c %.4f [mm]", f->GetParameter(1)*1e3, 0xB1, f->GetParError(1)*1e3));
    pt->AddText(Form("b = %.4f %c %.4f [mm]", f->GetParameter(2)*1e3, 0xB1, f->GetParError(2)*1e3));

    pt->AddText(Form("#chi^{2}/ndf = %.2f", float(f->GetChisquare()/f->GetNDF()) ));
     
    g->Draw("AP");
    pt->Draw("same");
    
    C->Update();
    C->SaveAs(resultFile.c_str());
    C->Clear();
    
	return 0;
}
