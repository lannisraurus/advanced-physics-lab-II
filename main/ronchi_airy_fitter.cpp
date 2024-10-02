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
	std::string dataFile = "data_in/AIRY_d1p5mm.pgm";					// Input data file, pgm
	std::string resultFile = "data_out/RONCHI1_airy_d1p5mm.png";		// Output fata file, png
	double xAiry = 876, yAiry = 730;							// position of the airy pattern's center, in pixels
	double drBin = 1;											// Size of the radial bin step, in pixels
	double radius = 40;											// Radius of the pattern to gather, in pixels
	// Physical System
	double pixelToDist = ((3.69e-6)/0.5055);					// Conversion constant from pixels to real life distance, in metres 
	double pixelToDist_err = (3.69e-6)*0.0073/(0.5055*0.5055);	// Error of the conversion factor, in metres
	double laser_wavelength = 633e-9;							// Wavelength of the light, in metres
	double focal_length = 250e-3;								// Focal length of the fourier transform lens, in metres
	// Fit initial parameters
	double fit_E_0 = 88267900;									// Amplitude of the airy pattern, in light value
	double fit_a = 0.940e-3;									// Aperture radius, in metres
	double fit_C = 1845;										// Background light, in light value
	// Display Settings
	int maxLightValue = 60000;									// Y Axis maximum value
	double maxDistance = 0.3e-3;								// X Axis maximum value
	int xNdiv = -512;											// Number of divisions in x axis, root notation
	int yNdiv = -610;											// Number of divisions in y axis, root notation
	int fNpoints = 1000;										// Number of points used in the display of the fitting function.
	float textLocation[4] = {0.5,0.5,0.9,0.9};					// Relative coordinates of the text, {x1rel,y1rel,x2rel,y2rel}
	int imageScaling = 140;


	//---------------------------------------------------------------
	//                           ALGORITHM
	//---------------------------------------------------------------

	// SIMULATION DEDUCED PARAMETERS
	double lambdaf = laser_wavelength*focal_length;
	int binsVecSize = std::floor(radius/drBin);
	std::vector<double> values(binsVecSize,0);
	std::vector<double> counter(binsVecSize,0);
	
	// Open Data File
	printf("Opening file: %s\n",dataFile.c_str());
	std::ifstream inputFile;
	inputFile.open(dataFile);
	if(!inputFile.is_open()){
		printf("ERROR: Couldn't open specified file. Did you misspell the name?\n");
		return -1;
	}
	
	// Turn file into matrix
	int matrixWid = 0, matrixHei = 0;
	std::string line = "";
	std::stringstream ss;
	std::getline(inputFile,line);
	std::getline(inputFile,line);
	ss << line;
	ss >> matrixWid >> matrixHei;
	ss.clear();
	printf("Width: %d, Height: %d\n",matrixWid,matrixHei);
	std::getline(inputFile,line);
	std::vector<std::vector<int>> matrix(matrixWid,std::vector<int>(matrixHei,0));
	for(int j = 0; j < matrixHei; j++){
		for(int i = 0; i < matrixWid; i++){
			inputFile >> matrix[i][j]; 
		}
	}
	printf("Image Pixels Preview (Debugging): value(0,0)=%d, value(1,0)=%d, value(2,0)=%d\n",matrix[0][0],matrix[1][0],matrix[2][0]);

	// Define ROI
	int roiXmin = int(std::round(xAiry)) - int(std::ceil(radius));
	int roiXmax = int(std::round(xAiry)) + int(std::floor(radius));
	int roiYmin = int(std::round(yAiry)) - int(std::ceil(radius));
	int roiYmax = int(std::round(yAiry)) + int(std::floor(radius));
	if(roiXmin < 0) roiXmin = 0;
	if(roiXmax > matrixWid) roiXmax = matrixWid-1;
	if(roiYmin < 0) roiYmin = 0;
	if(roiYmax > matrixHei) roiYmax = matrixHei-1;

	// Classification Algorithm
	for(double j = roiYmin; j <= roiYmax; j++){
		for(double i = roiXmin; i <= roiXmax; i++){
			double r = double(sqrt( (xAiry-i)*(xAiry-i) + (yAiry-j)*(yAiry-j) ));
			if (r < radius){
				int bin = int(std::floor(r/drBin));
				counter[bin] += 1;
				values[bin] += matrix[i][j];
			}
		}
	}

	// Make average value vector
	std::vector<double> avgRadialValues;
	std::vector<double> avgRadialValues_err;
	std::vector<double> radialPos;
	std::vector<double> radialPos_err;
	printf("Radial Averages: [ ");
	for(int i = 0; i < binsVecSize; i++){
		if(counter[i] != 0){
			avgRadialValues.push_back(values[i]/counter[i]);
			radialPos.push_back( pixelToDist*(i+0.5)*drBin );
			avgRadialValues_err.push_back(0);
			radialPos_err.push_back( pixelToDist_err*(i+0.5)*drBin );
			printf("%.1f ", avgRadialValues.back());
		}
	}
	printf(" ]\n");

	// Plot and Fit
	TCanvas* C = new TCanvas("C", "Canvas", 16*imageScaling, 9*imageScaling);
    TGraphErrors* g = new TGraphErrors(avgRadialValues.size(), radialPos.data(), avgRadialValues.data(), radialPos_err.data(), avgRadialValues_err.data());
    TF1* f = new TF1("AIRY" , "pow([0]*[1]*[3]/x,2)*pow(TMath::BesselJ1(2*TMath::Pi()*[1]*x/[3]),2) + [2]" , 0 , pixelToDist*radius);
	
	f->SetParameters(fit_E_0,fit_a,fit_C,lambdaf); // E_0, APERTURE, VERTICAL DISPLACEMENT, lambdaf (fixed param)
	
	f->FixParameter(3,lambdaf);

    g->SetTitle( (dataFile+" - Airy Pattern Radial Average Fit").c_str() );
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kAzure+2);
    g->SetLineColor(kBlue+2);
    g->SetMarkerSize(0.6);

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

	f->SetLineColor(kBlack);
    f->SetLineWidth(2);
	f->SetNpx(fNpoints);

    g->Fit("AIRY");

    TPaveText* pt = new TPaveText(textLocation[0]*maxDistance, textLocation[1]*maxLightValue, textLocation[2]*maxDistance, textLocation[3]*maxLightValue, "user");

    pt->SetTextSize(0.04);
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);

    pt->AddText(Form("I = |(E_{0}a/#nu)J_{1}(2#pia#nu)|^{2} + C"));
    pt->AddText(Form("E_{0} = %.0f %c %.0f [light value]", f->GetParameter(0), 0xB1, f->GetParError(0)));
    pt->AddText(Form("a = (%.5f %c %.5f) [mm]", f->GetParameter(1)*1e3, 0xB1, f->GetParError(1)*1e3));
	pt->AddText(Form("C = %.0f %c %.0f [light value]", f->GetParameter(2), 0xB1, f->GetParError(2)));
    pt->AddText(Form("#chi^{2}/ndf = %.2f", float(f->GetChisquare()/f->GetNDF()) ));
    
    g->Draw("AP");
	f->Draw("same");
    pt->Draw("same");
    
    C->Update();
    C->SaveAs(resultFile.c_str());
    C->Clear();

	return 0;
}
