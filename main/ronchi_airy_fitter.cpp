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

	// SIMULATION INPUT PARAMETERS
	std::string dataFile = "data_in/AIRY_d1p5mm.pgm";
	std::string resultFile = "data_out/AIRY_d1p5mm.png";
	double xAiry = 876, yAiry = 730;		// pixels
	double drBin = 1;						// pixels
	double radius = 40;						// pixels
	double pixelToDist = (3.69e-6)/0.5055;						// m 
	double pixelToDist_err = (3.69e-6)*0.0073/(0.5055*0.5055);	// m 

	// SIMULATION DEDUCED PARAMETERS
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
	TCanvas* C = new TCanvas("C", "Canvas", 16*70, 9*70);
    TGraphErrors* g = new TGraphErrors(avgRadialValues.size(), radialPos.data(), avgRadialValues.data(), radialPos_err.data(), avgRadialValues_err.data());

    g->SetTitle("Airy Pattern Radial Average Plot");
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kAzure+2);
    g->SetLineColor(kBlue+2);
    g->SetMarkerSize(1);

    g->GetXaxis()->SetTitle("Distance [m]");
    g->GetXaxis()->SetLimits(0., pixelToDist*radius);
    g->GetXaxis()->SetNdivisions(-910);
    g->GetXaxis()->SetLabelSize(0.028);
    g->GetXaxis()->SetTickLength(-0.04);
    g->GetXaxis()->SetLabelOffset(0.033);
    g->GetXaxis()->SetTitleOffset(1.3);

    g->GetYaxis()->SetTitle("Light Value [unitless]");
    g->GetYaxis()->SetLabelSize(0.028);
    g->GetYaxis()->SetRangeUser(0., 65000.);
    g->GetYaxis()->SetNdivisions(-510);

    g->Draw("AP*");
    
    C->Update();
    C->SaveAs(resultFile.c_str());
    C->Clear();

	return 0;
}
