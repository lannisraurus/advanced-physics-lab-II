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
	// Data Vars
	std::ifstream input;
	std::stringstream ss;
	input.open("data_in/RONCHI1_d8p15mm_s0p061ms.csv");
	std::vector<double> X;
	std::vector<double> Y;
	
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
		X.push_back(pixel);
		Y.push_back(value);
		printf("%f : %f\n",X.back(),Y.back());
	}

    // Calibration and errors
    std::vector<double> eX(X.size(),0);
    std::vector<double> eY(Y.size(),0);
    
	// Plot and Fit
	TCanvas* C = new TCanvas("C", "Canvas", 16*70, 9*70);
    TGraphErrors* g = new TGraphErrors(X.size(), X.data(), Y.data(), eX.data(), eY.data());

    g->SetTitle("Corrected linear fit");
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kAzure+2);
    g->SetLineColor(kBlue+2);
    g->SetMarkerSize(0.25);

    g->GetXaxis()->SetTitle("Distance [Ainda e Pixel]]");
    g->GetXaxis()->SetLimits(0., 1900.);
    g->GetXaxis()->SetNdivisions(-910);
    g->GetXaxis()->SetLabelSize(0.028);
    g->GetXaxis()->SetTickLength(-0.04);
    g->GetXaxis()->SetLabelOffset(0.033);
    g->GetXaxis()->SetTitleOffset(1.3);

    g->GetYaxis()->SetTitle("Energy [MeV]");
    g->GetYaxis()->SetLabelSize(0.028);
    g->GetYaxis()->SetRangeUser(0., 65000.);
    g->GetYaxis()->SetNdivisions(-510);

    TF1* f = new TF1("f" , "pow( [0]*sin(TMath::Pi()*[1]*(x-[2])) / (TMath::Pi()*(x-[2])) * sin(TMath::Pi()*[4]*[5]*(x-[2])) / sin(TMath::Pi()*(x-[2])*[4]) , 2) + [3]" , 0 , 1900); 

    f->SetLineColor(kBlack);
    f->SetLineWidth(2);

    f->SetParameters(6000.,0.0402,965.,900.,0.0115, 10.); // A, a, x0, B, d, m
                                                      //    a=0.14mm, d=0.18mm, m=10
                                                      //    amped: a=42mm, d=54
    f->FitParameters();    
    f->SetNpx(1000);

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
