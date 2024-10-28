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
	std::string dataFile = "bin/data_in/RONCHI1_d8p15mm_s0p061ms.csv";			// Input data file, csv
	std::string resultFile = "bin/data_out/RONCHI1_d8p15mm_s0p061ms_PLOT.png";	// Output fata file, png
	// Physical System
	double pixelToDist = ((3.69e-6)/0.5055);					// Conversion constant from pixels to real life distance, in metres 
	double pixelToDist_err = (3.69e-6)*0.0073/(0.5055*0.5055);	// Error of the conversion factor, in metres
	double pixelMeasureOffset = 0.007;							// Where is the center of the measure, to propagate error? in metres
	double laser_wavelength = 633e-9;							// Wavelength of the light, in metres
	double focal_length = 250e-3;								// Focal length of the fourier transform lens, in metres
	double lambdaf = laser_wavelength*focal_length;
	// Display Settings
	int maxLightValue = 55000;									// Y Axis maximum value
	double minDistance = 0.0035;
	double maxDistance = 0.0105;									// X Axis maximum value
	int xNdiv = -512;											// Number of divisions in x axis, root notation
	int yNdiv = -610;											// Number of divisions in y axis, root notation
	int fNpoints = 1000;										// Number of points used in the display of the fitting function.
	float textLocation[4] = {0.525,0.525,0.95,0.95};			// Relative coordinates of the text, {x1rel,y1rel,x2rel,y2rel}
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
		X.push_back(pixel*pixelToDist);
		Y.push_back(value);
		eX.push_back( std::abs((pixel-pixelMeasureOffset/pixelToDist)*pixelToDist_err) );
		eY.push_back(0);		
	}
    
	// Plot
	TCanvas* C = new TCanvas("C", "Canvas", 16*imageScaling, 9*imageScaling);
    TGraphErrors* g = new TGraphErrors(X.size(), X.data(), Y.data(), eX.data(), eY.data());

    g->SetTitle( (dataFile+" - Ronchi Grid Plot").c_str() );
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kAzure+2);
    g->SetLineColor(kBlue+2);
	g->SetLineWidth(2);
    g->SetMarkerSize(2);

    g->GetXaxis()->SetTitle("Distance [m]");
    g->GetXaxis()->SetLimits(minDistance, maxDistance);
    g->GetXaxis()->SetNdivisions(xNdiv);
    g->GetXaxis()->SetLabelSize(0.028);
    g->GetXaxis()->SetTickLength(-0.04);
    g->GetXaxis()->SetLabelOffset(0.033);
    g->GetXaxis()->SetTitleOffset(1.3);

    g->GetYaxis()->SetTitle("Intensity [light value]");
    g->GetYaxis()->SetLabelSize(0.028);
    g->GetYaxis()->SetRangeUser(0., maxLightValue);
    g->GetYaxis()->SetNdivisions(yNdiv);
    
    g->Draw("AL");
    
    C->Update();
    C->SaveAs(resultFile.c_str());
    C->Clear();
    
	return 0;
}
