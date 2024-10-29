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
	std::string dataFile = "bin/data_in/AIRY_LINEAR.csv";			// Input data file, csv
	std::string resultFile = "bin/data_out/AIRY_LINEAR.png";	// Output fata file, png
	// Physical System
	double laser_wavelength = 633e-9;							// Wavelength of the light, in metres
	double focal_length = 250e-3;								// Focal length of the fourier transform lens, in metres
	double lambdaf = laser_wavelength*focal_length;
	// Display Settings
	double maxY = 0.12;							// Y Axis maximum value
	double minY = 0.00;
	double maxX = 0.7;								// X Axis maximum value
	double minX = 0;
	int xNdiv = -514;									// Number of divisions in x axis, root notation
	int yNdiv = -512;									// Number of divisions in y axis, root notation
	int fNpoints = 1000;								// Number of points used in the display of the fitting function.
	float textLocation[4] = {0.55,0.1,0.9,0.5};	// Relative coordinates of the text, {x1rel,y1rel,x2rel,y2rel}
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
	double D, D_err, R, R_err;
	std::string line;
	
	// Error Checking
	if(!input.is_open()){
		printf("Error: Couldn't open file.\n");
		return -1;
	}

	// Retrieve data, calibration and errors
	for (line; std::getline(input, line); ) {

		ss << line;
		ss >> D >> D_err >> R >> R_err;
		ss.clear();

		X.push_back(1./D);
		Y.push_back(R);
		eX.push_back(D_err/(D*D));
		eY.push_back(R_err);

		printf("(%f,%f) +- (%f,%f)\n",X.back(),Y.back(), eX.back(),eY.back());

	}
    
	// Plot and Fit
	TCanvas* C = new TCanvas("C", "Canvas", 16*imageScaling, 9*imageScaling);
    TGraphErrors* g = new TGraphErrors(X.size(), X.data(), Y.data(), eX.data(), eY.data());

    g->SetTitle( (dataFile+" - Airy Radius vs. Aperture Diameter").c_str() );
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kAzure+2);
	g->SetLineWidth(2);
    g->SetLineColor(kBlack);
    g->SetMarkerSize(1.6);

    g->GetYaxis()->SetTitle("Airy Radius [mm]");
    g->GetXaxis()->SetLimits(minX, maxX);
    g->GetXaxis()->SetNdivisions(xNdiv);
    g->GetXaxis()->SetLabelSize(0.028);
    g->GetXaxis()->SetTickLength(-0.04);
    g->GetXaxis()->SetLabelOffset(0.033);
    g->GetXaxis()->SetTitleOffset(1.3);

    g->GetXaxis()->SetTitle("(Aperture Diameter)^{-1} [mm^{-1}]");
    g->GetYaxis()->SetLabelSize(0.028);
    g->GetYaxis()->SetRangeUser(minY, maxY);
    g->GetYaxis()->SetNdivisions(yNdiv);

    TF1* f = new TF1("f" , "([0]*[1])*x + [2]" , minX , maxX); 

    f->SetLineColor(kBlack);
    f->SetLineWidth(2.5);
	
	f->SetParameters(lambdaf,0.00015,0.005);	
	f->FixParameter(0,lambdaf);
    f->SetNpx(fNpoints);

    g->Fit("f");

    TPaveText* pt = new TPaveText(textLocation[0]*maxX, textLocation[1]*maxY, textLocation[2]*maxX, textLocation[3]*maxY, "user");

    pt->SetTextSize(0.039);
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);

    pt->AddText(Form("R_{airy} = A #lambdaf (1/D) + C"));
    pt->AddText(Form("A = %.3f %c %.3f [unitless]", f->GetParameter(1)*1e-6, 0xB1, f->GetParError(1)*1e-6));
    pt->AddText(Form("C = %.5f %c %.5f [mm]", f->GetParameter(2), 0xB1, f->GetParError(2)));

    pt->AddText(Form("#chi^{2}/ndf = %.2f", float(f->GetChisquare()/f->GetNDF()) ));
    printf("CHI2: %f",float(f->GetChisquare()/f->GetNDF()) );

	g->Draw("AP*");
	//f->Draw("same");
    pt->Draw("same");
    
    C->Update();
    C->SaveAs(resultFile.c_str());
    C->Clear();
    
	return 0;
}
