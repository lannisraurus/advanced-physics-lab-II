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
	// Physical System
	double pixelToDist = ((3.69e-6)/0.5055);					// Conversion constant from pixels to real life distance, in metres 
	double pixelToDist_err = (3.69e-6)*0.0073/(0.5055*0.5055);	// Error of the conversion factor, in metres
	double laser_wavelength = 633e-9;							// Wavelength of the light, in metres
	double focal_length = 250e-3;								// Focal length of the fourier transform lens, in metres
	double lambdaf = laser_wavelength*focal_length;

	// Data Vars
	std::string fileName = "data_in/RONCHI1_d8p15mm_s0p061ms.csv";
	std::string displayName = "Ronchi Grid 1 (d_{iris}= 8.15mm)";
	std::ifstream input;
	std::stringstream ss;
	input.open(fileName);
	std::vector<double> X;
	std::vector<double> Y;
	int noiseThreshhold = 1050;
	
	// Probing Vars
	double pixel;
	double value;
	std::string line;
	
	// Error Checking
	if(!input.is_open()){
		printf("Error: Couldn't open file.\n");
		return -1;
	}

	// Retrieve data
	std::getline(input, line);
	for (line; std::getline(input, line); ) {
		std::replace( line.begin(), line.end(), ',', ' ');
		ss << line;
		ss >> pixel >> value;
		ss.clear();
		if (value > noiseThreshhold){
			X.push_back(pixel);
			Y.push_back(value);
			printf("%f : %f\n",X.back(),Y.back());
		}
	}

    // Calibration and errors
    std::vector<double> eX(X.size(),0);
    std::vector<double> eY(Y.size(),0);
	for(int i = 0; i < X.size(); i++){
		eX[i] = X[i]*pixelToDist_err;
		X[i] = X[i]*pixelToDist;
	}
    
	// Plot and Fit
	TCanvas* C = new TCanvas("C", "Canvas", 16*140, 9*140);
    TGraphErrors* g = new TGraphErrors(X.size(), X.data(), Y.data(), eX.data(), eY.data());

    g->SetTitle(displayName.c_str());
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kAzure+2);
    g->SetLineColor(kBlue+2);
    g->SetMarkerSize(1);
	g->SetLineWidth(1);

    g->GetXaxis()->SetTitle("Distance [m]");
    g->GetXaxis()->SetLimits(0., 1900.*pixelToDist);
    g->GetXaxis()->SetNdivisions(-910);
    g->GetXaxis()->SetLabelSize(0.028);
    g->GetXaxis()->SetTickLength(-0.04);
    g->GetXaxis()->SetLabelOffset(0.033);
    g->GetXaxis()->SetTitleOffset(1.3);

    g->GetYaxis()->SetTitle("Intensity [light value]");
    g->GetYaxis()->SetLabelSize(0.028);
    g->GetYaxis()->SetRangeUser(0., 65000.);
    g->GetYaxis()->SetNdivisions(-510);

    TF1* f = new TF1("f" , "pow( [0]*sin(TMath::Pi()*[1]*(x-[2])*[6]) / (TMath::Pi()*(x-[2])*[6]) * sin(TMath::Pi()*[4]*[5]*(x-[2])*[6]) / sin(TMath::Pi()*(x-[2])*[6]*[4]) , 2) + [3]" , 0 , 1900*pixelToDist); 

    f->SetLineColor(kBlack);
    f->SetLineWidth(1.5);

    //f->SetParameters(6000*pixelToDist/lambdaf,		// A
	//				0.0402/pixelToDist*lambdaf,		// a
	//				965.*pixelToDist/lambdaf,		// x0
	//				900.,							// B
	//				0.0115/pixelToDist*lambdaf,		// d
	//				10.								// m
	//				);


	f->SetParameters(6000*lambdaf/pixelToDist,		// E_0
					0.0402*lambdaf/pixelToDist,		// a
					965.*pixelToDist,				// x0
					900.,							// B
					0.0115*lambdaf/pixelToDist,		// d
					10.,							// m
					lambdaf
					);
														  //    a=0.14mm, d=0.18mm, m=10
                                                          //    amped: a=42mm, d=54  
    f->SetNpx(5000);

    //g->Fit("f");

    //TPaveText* pt = new TPaveText(245, 5.35, 335, 5.65, "user");

    //pt->SetTextSize(0.03);
    //pt->SetFillColor(0);
    //pt->SetTextAlign(12);
    //pt->SetTextFont(42);

    //pt->AddText(Form("y = ax + b"));
    //pt->AddText(Form("a = %.10f %c %.10f", f->GetParameter(0), 0xB1, f->GetParError(0)));
    //pt->AddText(Form("b = %.10f %c %.10f", f->GetParameter(1), 0xB1, f->GetParError(1)));
    //pt->AddText(Form("#chi^{2} = %.10f", f->GetChisquare()));
    
    
    
    g->Draw("AP");
    f->Draw("same");
    //pt->Draw("same");
    
    C->Update();
    C->SaveAs("data_out/RONCHI1_d8p15mm_s0p061ms_FIT.png");
    C->Clear();
    


	return 0;
}
