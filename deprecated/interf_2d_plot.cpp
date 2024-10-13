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
#include "TGraph2D.h"
#include "TStyle.h"
#include "TPaveText.h"

#include "TH2F.h"
#include "TH2.h"
#include "TF2.h"
#include "TROOT.h"

int main(){

	//---------------------------------------------------------------
	//                      INPUT PARAMETERS
	//---------------------------------------------------------------
	// Data file info
	std::string dataFile = "bin/data_in/Interf2_calha_d427p0mm_P2.pgm";		// Input data file, pgm
	std::string resultFile = "bin/data_out/Interf2_calha_d427p0mm_P2_2d_plot.png";	// Output fata file, png
	// Physical System	
	double laser_wavelength = 633e-9;							// Wavelength of the light, in metres
	double focal_length = 250e-3;								// Focal length of the fourier transform lens, in metres
	// Fit initial parameters
	double fit_I_0 = 50;
	double fit_alpha = 0.1;
	double fit_C = 500;
	double fit_angle = M_PI/7.;
	// Display Settings
	int maxLightValue = 50000;									// Z Axis maximum value
	int xNdiv = -512;											// Number of divisions in x axis, root notation
	int yNdiv = -512;											// Number of divisions in y axis, root notation
	int zNdiv = -610;											// Number of divisions in z axis, root notation
	int fNpoints = 10*10;										// Number of points used in the display of the fitting function.
	float textLocation[4] = {0.525,0.525,0.95,0.95};			// Relative coordinates of the text, {x1rel,y1rel,x2rel,y2rel}
	int imageScaling = 140;

	//---------------------------------------------------------------
	//                           ALGORITHM
	//---------------------------------------------------------------

	// Data Vars
	std::ifstream input;
	std::stringstream ss;
	input.open(dataFile);
	int xMax, yMax, zMax;
	std::vector<double> X;
	std::vector<double> Y;
	std::vector<double> Z;
	
	// Probing Vars
	double z;
	std::string line;
	
	// Error Checking
	if(!input.is_open()){
		printf("Error: Couldn't open file.\n");
		return -1;
	}

	// Retrieve data
	std::getline(input, line); // Get format (useless)
	std::getline(input, line); // Get dimensions
	ss << line;
	ss >> xMax >> yMax;
	ss.clear();
	std::getline(input, line); // Get max pixel
	ss << line;
	ss >> zMax;
	ss.clear();
	printf("XMAX=%d, YMAX=%d, VALUEMAX=%d\n",xMax,yMax,zMax);
	for(int j = 0; j < yMax; j++){
		for(int i = 0; i < xMax; i++){
			X.push_back(i);
			Y.push_back(j);
			input >> z;
			Z.push_back(z);
		}
	}
	printf("DONE READING IMAGE.\n");
		 
	// Plot and Fit
	TCanvas* C = new TCanvas("C", "Canvas", 16*imageScaling, 9*imageScaling);
    TGraph2D* g = new TGraph2D(X.size(), X.data(), Y.data(), Z.data() );
		
    g->SetTitle( (dataFile+" - Interferometry 2D Plot;Pixel [pixels];Pixel [pixels];Intensity [light value]").c_str());
	g->SetFillColorAlpha(kRed, 1);
	g->SetLineWidth(1);
	g->SetNpx(fNpoints);
    // g->SetMarkerStyle(20);
    // g->SetMarkerColor(kAzure+2);
    // g->SetMarkerSize(1.5);

	//g->GetXaxis()->SetLabelOffset(9000);
	//g->GetYaxis()->SetLabelOffset(9000);
	//g->GetZaxis()->SetLabelOffset(9000);
	
	/*	
    g->GetXaxis()->SetTitle("Pixel [pixels]");
    g->GetXaxis()->SetLimits(0., xMax);
    g->GetXaxis()->SetNdivisions(xNdiv);
    g->GetXaxis()->SetLabelSize(0.028);
    g->GetXaxis()->SetTickLength(-0.04);
    g->GetXaxis()->SetLabelOffset(0.033);
    g->GetXaxis()->SetTitleOffset(1.3);
	
	g->GetYaxis()->SetTitle("Pixel [pixels]");
    g->GetYaxis()->SetLimits(0., yMax);
    g->GetYaxis()->SetNdivisions(yNdiv);
    g->GetYaxis()->SetLabelSize(0.028);
    g->GetYaxis()->SetTickLength(-0.04);
    g->GetYaxis()->SetLabelOffset(0.033);
    g->GetYaxis()->SetTitleOffset(1.3);

    g->GetZaxis()->SetTitle("Intensity [light value]");
    g->GetZaxis()->SetLabelSize(0.028);
    g->GetZaxis()->SetRangeUser(0., maxLightValue);
    g->GetZaxis()->SetNdivisions(zNdiv);
	*/

	TF2* f = new TF2("f" , "2*[0]*sin( [1] * (x*cos([3]) + y*sin([3])) + [4])  + [2]" , 0 , xMax , 0 , yMax); 
	// f->SetParameters(fit_I_0,fit_alpha,fit_C,fit_angle); // E_0, APERTURE, VERTICAL DISPLACEMENT, lambdaf (fixed param)
	f->SetFillColorAlpha(kBlack, 0.5);
	f->SetLineWidth(2);
	f->SetNpx(fNpoints);
	g->Fit("f");

	/*
    TPaveText* pt = new TPaveText(textLocation[0]*maxDistance, textLocation[1]*maxLightValue, textLocation[2]*maxDistance, textLocation[3]*maxLightValue, "user");

    pt->SetTextSize(0.032);
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
    */
    
	gStyle->SetCanvasPreferGL(true);
	//g->Draw("surf");
	f->Draw("surf1");
	g->Draw("surf3");
	
	//pt->Draw("same");
    
    C->Update();
    C->SaveAs(resultFile.c_str());
    C->Clear();
    
	return 0;
}
