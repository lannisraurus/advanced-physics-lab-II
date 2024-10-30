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

#include "TF1.h"
#include "TROOT.h"

#include "TMath.h"
#include "TGraph2D.h"
#include "TRandom.h"
#include "TF2.h"
#include "TError.h"

#include "TH1F.h"
#include "TH1.h"


//740,569  794,890 (y,x)
//1072,582   1041x790 (y,x)
bool in_solder(int i, int j) {
    return j>741 && j<1085 && i<1132 && j> (794-740)/(890-569)*(i-569)+ 740 && j< (1041-1072)/(790-582)*(i-582)+ 1072;
}



struct PGM {
    std::vector<std::vector<int>> M;
    std::string type;
    int N, nrows, ncols;
    PGM(std::vector<std::vector<int>> m, std::string c, int n, int nr, int nc): M(m), type(c), N(n), nrows(nr), ncols(nc) {}
};

PGM readPGM(std::string file) {
    std::ifstream FI(file);
    std::string line, type;
    int N, nrows, ncols;
    FI >> type;
    FI >> ncols >> nrows;
    FI >> N;
    std::vector<std::vector<int>> M(nrows, std::vector<int>(ncols));;
    int i = 0; int j = 0;
    while (FI >> M[i][j]) {
        j++;
        if (j==ncols) {j=0; i++;}
    }
    FI.close();
    PGM r(M,type,N,nrows,ncols);
    return r;
}

void writePGM(const PGM& pgm, std::string file) {
    std::ofstream FO(file);
    FO << pgm.type << "\n" << pgm.ncols << " " << pgm.nrows << "\n" << pgm.N << "\n";
    int i = 0; int j = 0;
    std::string line;
    while(true) {
        line = "";
        int charcount = 0;
        while (charcount + std::to_string(pgm.M[i][j]).size() < 69) {
            line += std::to_string(pgm.M[i][j]) + " ";
            charcount += std::to_string(pgm.M[i][j]).size()+1;
            j++;
            if (j==pgm.ncols) {j=0; i++;}
            if(i==pgm.nrows) break;
        }
        if(i==pgm.nrows) break;
        line += "\n";
        FO << line;
    }
    FO << line << "\n";
    FO.close();
}

void writeFitToPGM(const PGM& pgm, TF2 * fit, std::string file) {
    std::ofstream FO(file);
    FO << pgm.type << "\n" << pgm.ncols << " " << pgm.nrows << "\n" << pgm.N << "\n";
    int i = 0; int j = 0;
    std::string line;
    while(true) {
        line = "";
        int charcount = 0;
        while (charcount + std::to_string(fit->Eval(i,j)).size() < 69) {
            line += std::to_string(fit->Eval(i,j)) + " ";
            charcount += std::to_string(fit->Eval(i,j)).size()+1;
            j++;
            if (j==pgm.ncols) {j=0; i++;}
            if(i==pgm.nrows) break;
        }
        if(i==pgm.nrows) break;
        line += "\n";
        FO << line;
    }
    FO << line << "\n";
    FO.close();
}


void graphPGM(const PGM& pgm, TGraph2D* g, TRandom* r,double p) {
    int count = 0;
    for (int i = 0; i < pgm.nrows; i++) {
        for (int j = 0; j < pgm.ncols; j++) {
            //&& (i<400 || (i==1000 && j==900) || j<400 || j>1400)
            if(r->Rndm(j+i*pgm.ncols) < p && !in_solder(i,j)){
                count++;
                g->SetPoint(count,i,j,pgm.M[i][j]);
            } 
        }
    }
}

TGraph* graphLinePGM(const PGM& pgm, int line,int avg) {
    std::vector<double> y;
    std::vector<double> val;
    for (int j = 0; j < pgm.ncols; ++j) {
        double avgVal = 0;
        for (int i = line-avg; i < line+avg+1; ++i) {
            avgVal += pgm.M[i][j];
        }
        if(true || !in_solder(line,j)) {
            y.push_back(j);
            val.push_back(avgVal/(2*avg+1));
        }
    }
    TGraph* g = new TGraph(y.size(), y.data(), val.data());
    return g;
}

TGraph* graphMinPGM(const PGM& pgm, int start_x,int sieve_size) {
    std::vector<double> y;
    std::vector<double> x;
    double min_x = start_x;
    for (int j = 0; j < pgm.ncols; ++j) {
        if(!in_solder(min_x,j)) {
            std::vector<double> vals;
            for (int i = min_x-sieve_size; i < min_x+sieve_size+1; ++i) {
                vals.push_back(pgm.M[i][j]);
            }
            min_x = min_x-sieve_size + std::distance(std::begin(vals), std::min_element(std::begin(vals), std::end(vals)));
            y.push_back(j);
            x.push_back(min_x);
        }
    }
    TGraph* g = new TGraph(y.size(), y.data(), x.data());
    return g;
}

void drawGraph(TGraph2D* g, TCanvas* c, const char* style= "", TF2* f = nullptr, const char* stylef = "") {
    if(f) f->Draw(stylef);
    if(g) g->Draw(style);
    c->Update();
    gSystem->ProcessEvents();
    c->WaitPrimitive();
}

void drawGraph1D(TGraph* g1d, TCanvas* c, const char* style= "", TF1* f = nullptr, const char* stylef = "") {
    if(f) f->Draw(stylef);
    if(g1d) {g1d->Draw(style);} 
    c->Update();
    gSystem->ProcessEvents();
    c->WaitPrimitive();
}

double square(double *x1, double *par){
    double z = x1[0]; 
    double y = par[0];
    double phiFactor = par[1];
    double R = par[2];
    double r = sqrt(z*z+y*y);
    return r < R ? phiFactor : 0;
}

double triangle(double *x1, double *par){
    double z = x1[0]; 
    double y = par[0];
    double phiFactor = par[1];
    double R = par[2];
    double r = sqrt(z*z+y*y);
    return r < R ? phiFactor * (R-r)/R : 0;
}

double gauss(double *x1, double *par){
    double z = x1[0]; 
    double y = par[0];
    double phiFactor = par[1];
    double R = par[2];
    double r = sqrt(z*z+y*y);
    double sigma = R/2/sqrt(2);
    return r < R ? phiFactor * exp(-r*r/2/sigma/sigma) : 0;
}

double negExp(double *x1, double *par){
    double z = x1[0]; 
    double y = par[0];
    double phiFactor = par[1];
    double R = par[2];
    double r = sqrt(z*z+y*y);
    double alpha = R/4;
    return r < R ? phiFactor * exp(-r/alpha) : 0;
}

double fitFunc1D(double *x1, double *par){
    double y = x1[0];
    double A = par[0];
    double omega_x = par[1];
    double omega_y = par[2];
    double phi_0 = par[3];
    double BG = par[4];    
    double omega_y_secondary = par[5];    
    double phi_0_secondary = par[6];
    double x = par[7];
    double y_solder = par[8];
    double phiFactor = par[9];
    double R = par[10];
    TF1 f("radialDistribuition", triangle, 0, 1900, 3); //min, max, no.parameters
    f.SetParameter(0,y-y_solder);
    f.SetParameter(1,phiFactor);
    f.SetParameter(2,R);
    return A*(1+sin(omega_x*x + omega_y*y + phi_0 + f.Integral( -sqrt(std::max(R*R-(y-y_solder)*(y-y_solder),0.)),sqrt(std::max(R*R-(y-y_solder)*(y-y_solder),0.)) ) ))*sin(omega_y_secondary*y+phi_0_secondary) + BG;  //Integrate 
}

double fitFuncMins(double *x1, double *par){
    double y = x1[0];
    double a = par[0];
    double b = par[1];
    double y_solder = par[2];
    double phiFactor = par[3];
    double R = par[4];
    TF1 f("radialDistribuition", triangle, 0, 1900, 3); //min, max, no.parameters
    f.SetParameter(0,y-y_solder);
    f.SetParameter(1,phiFactor);
    f.SetParameter(2,R);
    return a*y+b + f.Integral( -sqrt(std::max(R*R-(y-y_solder)*(y-y_solder),0.)),sqrt(std::max(R*R-(y-y_solder)*(y-y_solder),0.)) );  //Integrate 
}
//fit paramters
// p0 = kinda contraste/ fator do sin 
// p1 = termo do x dentro do sin / quanto oscila ao longo de x
// p2 = termo do y dentro do sin / quanto oscila ao longo de y
// p3 = extra fase do sin ~entre 0 e 2pi
// p4 = avg fora da solda
// p5 = y_B
// p6 = y_C
// p7 = x limite da solda
// p8 = y_A
// p9 = x_C
// p10 = avg dentro da solda
TF2* fitToGraph(TGraph2D* g, double* pars, bool solder = false, int phi = 0, int parsToFix = 0, int* parsf = nullptr) {
    TF2* f; 
    if(solder) {
        if(phi==1) {
            //TF2* dist = new TF2("dist","y-([5]+[8])/2",0,1500,0,2000);
            //TF2* phi = new TF2("phi","[11]*sqrt(TMath::Max([12]*[12]*(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))-pow(dist,2), 0.0))",0,1500,0,2000); (TMath::Max(exp([13]*(x-[14]))-1
            f= new TF2("f","x>[7] && y<([5]-[6])/(1447.-[9])*(x-[9])+[6] && y>([8]-[6])/(1447.-[9])*(x-[9])+[6] ? [10] : [0]*sin([1]*x+[2]*y+[3] + [11]*sqrt(TMath::Max([12]*[12]*(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))-(y-([5]+[8])/2)*(y-([5]+[8])/2), 0.0)) - [11]*(1./2./[12]/sqrt(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))*(y-([5]+[8])/2)*(y-([5]+[8])/2)*log(([12]*sqrt(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))+sqrt(TMath::Max([12]*[12]*(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))-(y-([5]+[8])/2)*(y-([5]+[8])/2), 0.0))  )/([12]*sqrt(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))-sqrt(TMath::Max([12]*[12]*(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))-(y-([5]+[8])/2)*(y-([5]+[8])/2), 0.0))  )) )   )+[4]",0,1500,0,2000);
            f->SetParLimits(12, 0, 1000);
        }else if(phi==2){
            f= new TF2("f","x>[7] && y<([5]-[6])/(1447.-[9])*(x-[9])+[6] && y>([8]-[6])/(1447.-[9])*(x-[9])+[6] ? [10] : [0]*sin([1]*x+[2]*y+[3] + 2*[11]*sqrt(TMath::Max([12]*[12]*(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))-(y-([5]+[8])/2)*(y-([5]+[8])/2), 0.0)))+[4]",0,1500,0,2000);
        } else {
            f= new TF2("f","x>[7] && y<([5]-[6])/(1447.-[9])*(x-[9])+[6] && y>([8]-[6])/(1447.-[9])*(x-[9])+[6] ? [10] : [0]*sin([1]*x+[2]*y+[3])+[4]",0,1500,0,2000);
        }
    } else {
        f= new TF2("f","[0]*(1+sin([1]*x+[2]*y+[3]))*sin([5]*y+[6])+[4]",0,1500,0,2000);
    }
    f->SetParameters(pars);
    for(int i=0; i<parsToFix;++i) {
        f->FixParameter(parsf[i],pars[parsf[i]]);
    }
    g->Fit(f);

    return f;
}

TF1* fitToGraph1D(TGraph* g1d, double* pars, int parsToFix = 0, int* parsf = nullptr) {
    TF1* f; 
    //f= new TF2("f","x>[7] && y<([5]-[6])/(1447.-[9])*(x-[9])+[6] && y>([8]-[6])/(1447.-[9])*(x-[9])+[6] ? [10] : [0]*sin([1]*x+[2]*y+[3] + 2*[11]*sqrt(TMath::Max([12]*[12]*(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))-(y-([5]+[8])/2)*(y-([5]+[8])/2), 0.0)) - [11]*(sqrt(TMath::Max([12]*[12]*(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))-(y-([5]+[8])/2)*(y-([5]+[8])/2), 0.0)) + 1./2./[12]/sqrt(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))*(y-([5]+[8])/2)*(y-([5]+[8])/2)*log(([12]*sqrt(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))+sqrt(TMath::Max([12]*[12]*(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))-(y-([5]+[8])/2)*(y-([5]+[8])/2), 0.0))  )/([12]*sqrt(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))-sqrt(TMath::Max([12]*[12]*(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))-(y-([5]+[8])/2)*(y-([5]+[8])/2), 0.0))  )) )   )+[4]",0,1500,0,2000);
    f= new TF1("f",fitFunc1D,0,2000,11);
    //f= new TF1("f","[0]*sin([1]*[5]+[2]*x+[3]+ 2*[6]*sqrt(TMath::Max([7]*[7]-(x-[8])*(x-[8]), 0.0))   )+[4]",0,2000);
    f->SetParameters(pars);
    for(int i=0; i<parsToFix;++i) {
        f->FixParameter(parsf[i],pars[parsf[i]]);
    }
    g1d->Fit(f);

    return f;
}

TF1* fitToGraphMins(TGraph* g1d, double* pars, int parsToFix = 0, int* parsf = nullptr) {
    TF1* f; 
    //f= new TF2("f","x>[7] && y<([5]-[6])/(1447.-[9])*(x-[9])+[6] && y>([8]-[6])/(1447.-[9])*(x-[9])+[6] ? [10] : [0]*sin([1]*x+[2]*y+[3] + 2*[11]*sqrt(TMath::Max([12]*[12]*(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))-(y-([5]+[8])/2)*(y-([5]+[8])/2), 0.0)) - [11]*(sqrt(TMath::Max([12]*[12]*(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))-(y-([5]+[8])/2)*(y-([5]+[8])/2), 0.0)) + 1./2./[12]/sqrt(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))*(y-([5]+[8])/2)*(y-([5]+[8])/2)*log(([12]*sqrt(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))+sqrt(TMath::Max([12]*[12]*(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))-(y-([5]+[8])/2)*(y-([5]+[8])/2), 0.0))  )/([12]*sqrt(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))-sqrt(TMath::Max([12]*[12]*(TMath::Max(exp([13]*(x-[14]))-1,1.e-10))-(y-([5]+[8])/2)*(y-([5]+[8])/2), 0.0))  )) )   )+[4]",0,1500,0,2000);
    f= new TF1("f",fitFuncMins,0,2000,5);
    //f= new TF1("f","[0]*sin([1]*[5]+[2]*x+[3]+ 2*[6]*sqrt(TMath::Max([7]*[7]-(x-[8])*(x-[8]), 0.0))   )+[4]",0,2000);
    f->SetParameters(pars);
    for(int i=0; i<parsToFix;++i) {
        f->FixParameter(parsf[i],pars[parsf[i]]);
    }
    g1d->Fit(f);

    return f;
}

std::pair<double,double> CalcTemp(double PhiFactorPx, double ePhiFactorPx, double R, double eR, double Magnification, double eMagnification, double SolderR, double eSolderR, int fit=1) {
    double AirN = 1.000267; //https://www.wolframalpha.com/input?i=refractive+index+of+moist+air&assumption=%7B%22F%22%2C+%22MoistAirRefractiveIndex%22%2C+%22Tair%22%7D+-%3E%2225+%C2%B0C%22&assumption=%7B%22F%22%2C+%22MoistAirRefractiveIndex%22%2C+%22CCO2%22%7D+-%3E%22450+ppm%22&assumption=%7B%22F%22%2C+%22MoistAirRefractiveIndex%22%2C+%22rh%22%7D+-%3E%2250+%25%22&assumption=%7B%22F%22%2C+%22MoistAirRefractiveIndex%22%2C+%22lambda%22%7D+-%3E%22633+nm%22&assumption=%7B%22F%22%2C+%22MoistAirRefractiveIndex%22%2C+%22P%22%7D+-%3E%221.01325+bar%22 https://emtoolbox.nist.gov/wavelength/documentation.asp
    double eAirN = 1e-6;
    double TempBg = 25; // ºC
    double eTempBg = 1;
    double PixelSize = 3.69e-6; //m/pixel
    double ePixelSize =0.005e-6; 
    double Lambda = 6.33e-7; //m
    double eLambda = 0.005e-7;

    double PhiFactor = PhiFactorPx/PixelSize*Magnification; //1/m
    double ePhiFactor = ePhiFactorPx/PixelSize*Magnification + PhiFactorPx/PixelSize/PixelSize*Magnification*ePixelSize + PhiFactorPx/PixelSize*eMagnification;
    double DeltaNCenter = PhiFactor/2/M_PI*Lambda; //adim
    double eDeltaNCenter = ePhiFactor/2/M_PI*Lambda + PhiFactor/2/M_PI*eLambda;
    
    double DeltaNSolder;
    double eDeltaNSolder;
    if(fit==1) {
        DeltaNSolder = DeltaNCenter*(1-SolderR/R);
        eDeltaNSolder = eDeltaNCenter*(1-SolderR/R) + DeltaNCenter*eSolderR/R + DeltaNCenter*SolderR/R/R*eR;
    } else if(fit==2) {
        DeltaNSolder = DeltaNCenter*exp(-SolderR*SolderR*4/R/R);
        eDeltaNSolder = eDeltaNCenter*exp(-SolderR*SolderR*4/R/R) + DeltaNCenter*exp(-SolderR*SolderR*4/R/R)*SolderR*eSolderR*8/R/R + DeltaNCenter*exp(-SolderR*SolderR*4/R/R)*SolderR*SolderR*8/R/R/R*eR;
    } else if(fit==3) {
        DeltaNSolder = DeltaNCenter*exp(-SolderR*4/R);
        eDeltaNSolder = eDeltaNCenter*exp(-SolderR*4/R) + DeltaNCenter*exp(-SolderR*4/R)*eSolderR*4/R + DeltaNCenter*exp(-SolderR*4/R)*SolderR*4/R/R*eR;
    }

    double TempRatio = (AirN-1)/(AirN-1-DeltaNSolder); //adim
    double eTempRatio = eAirN/(AirN-1-DeltaNSolder)+(AirN-1)/(AirN-1-DeltaNSolder)/(AirN-1-DeltaNSolder)*eAirN + (AirN-1)/(AirN-1-DeltaNSolder)/(AirN-1-DeltaNSolder)*eDeltaNSolder;
    double Temp = TempRatio*(TempBg + 273.15) - 273.15; // ºC
    double eTemp = eTempRatio*(TempBg + 273.15) + TempRatio*eTempBg;
    std::cout << "Solder temp is " << Temp << "+-" << eTemp << ".ºC" << std::endl;
    return {Temp,eTemp};
}



int main(){
    // Root initialization
    TApplication *A = new TApplication("a",0,0);
    TCanvas *c = new TCanvas("Canvas name", "c", 400, 100, 1000, 800);
    TRandom *r = new TRandom();
    TGraph2D *g;
    TGraph *g1d;
    TF2* f;
    TF1* f1d;
    gStyle->SetPalette(1);
    gErrorIgnoreLevel = kPrint;

    // Data files temperature
    PGM Tbg1 = readPGM("bin/data_in/SOLDA/BG_L1_T0.pgm");
    PGM Tbg2 = readPGM("bin/data_in/SOLDA/BG_L1_T0_2.pgm");
    PGM T150C1 = readPGM("bin/data_in/SOLDA/L1_T1.pgm");
    PGM T150C2 = readPGM("bin/data_in/SOLDA/L1_T1_2.pgm");
    PGM T150C3 = readPGM("bin/data_in/SOLDA/L1_T1_3.pgm");
    PGM T300C1 = readPGM("bin/data_in/SOLDA/L1_T2.pgm");
    PGM T300C2 = readPGM("bin/data_in/SOLDA/L1_T2_2.pgm");
    PGM T300C3 = readPGM("bin/data_in/SOLDA/L1_T2_3.pgm");
    PGM T450C1 = readPGM("bin/data_in/SOLDA/L1_T3.pgm");
    PGM T450C2 = readPGM("bin/data_in/SOLDA/L1_T3_2.pgm");
    PGM T450C3 = readPGM("bin/data_in/SOLDA/L1_T3_3.pgm");
    PGM T450C4 = readPGM("bin/data_in/SOLDA/L1_T3_4.pgm");
    PGM T450C5 = readPGM("bin/data_in/SOLDA/L1_T3_5.pgm");
    PGM Vela = readPGM("bin/data_in/VELA/AMPLIACAO1/VELA_BAIXO_LADO.pgm");

    // soldering iron background analysis
    g = new TGraph2D();
    graphPGM(Tbg1,g,r,0.01);
    //drawGraph(g, c, "surf1");
    double pars2[8] = {20000,0.0625089,0.000388446,-0.5,5043.99, -0.000729725, 1.5};
    f = fitToGraph(g, pars2);
    //drawGraph(nullptr, c, "same p0", f, "surf1");
    writeFitToPGM(Tbg1, f,"bin/data_out/Tbg1Fit.pgm");


    // soldering iron 300C1 analysis 2D
    /*
    g = new TGraph2D();
    graphPGM(T150C1,g,r,0.1);
    drawGraph(g, c, "surf1");
    double pars3[15] = {13622.1,0.0658905,0.00127228,6.64916,19655.2,solder5,solder6,solder7,solder8,solder9, 10000,-0.02,400,0.002, 400};
    int pars3fix[15] = {5,6,7,8,9};
    f = fitToGraph(g, pars3, true, 1,5,pars3fix);
    drawGraph(nullptr, c, "same p0", f, "surf1");
    writeFitToPGM(Tbg, f,"bin/data_out/T150C1Fit.pgm");
    

    // soldering iron 450C5 analysis 1D
 */ double pars1D[12] = {20860.6,0.0625089,0.000388446,-0.475559,5043.99,-0.000729725, 1.54302, 0 ,905,0.018,700};
    int pars1Dfix[12] = {0,1,2,3,4,5,6,7,8};
    double parsMins[12] = {-0.01, 0 ,905,-0.05,500};
    int parsMinsfix[12] = {2};
    double line = 0;
    double start_x = 0;



    start_x = 150;
    g1d = graphMinPGM(T150C3,start_x,5);

    //drawGraph1D(g1d, c,"");
    parsMins[1] = start_x;
    f1d=fitToGraphMins(g1d,parsMins,0,parsMinsfix);
    //drawGraph1D(g1d, c,"",f1d);



  double pars1Dv[12] = {20860.6,0,0,0,5043.99,0, M_PI/2., -1 ,0,0.047,1100};
    int pars1Dvfix[12] = {1,2,5,6,7};
    line = 100;
    g1d = graphLinePGM(Vela,line,5);
    drawGraph1D(g1d, c, "");
    pars1D[7] = line;
    f1d=fitToGraph1D(g1d,pars1Dv,5,pars1Dvfix);
    g1d->SetTitle( "1d sine of triangular fit - Air above a lit candle");
    g1d -> GetXaxis() -> SetTitle(" y [px] ");
    g1d -> GetYaxis() -> SetTitle(" intensity [Celsius] ");
    g1d -> GetXaxis() -> SetRangeUser(0,1928);
    g1d -> GetYaxis() -> SetRangeUser (0,20000);
    drawGraph1D(g1d, c,"",f1d);

    c->SaveAs("vela.png");

    // soldering iron 450C1 analysis olhometro
    std::vector<double> z1 = {555,567,579,603,615,643,672,703,740,808,1056,1113,1167,1211,1254,1282,1311,1336,1364,776,780,795,815,835,1171,1160,1137,1119,1070,427,435,456,470,502,515,544,571,609,646,697,820,1024,1121,1203,1245,1294,1339,1370,1398,1426,1448,1469,1491};
    std::vector<double> y1 = {1393,1293,1195,1096,997,901,807,711,615,523,517,610,702,797,892,985,1082,1185,1286,1389,1287,1189,1089,992,1391,1283,1182,1081,988,1400,1300,1201,1101,1002,905,810,716,621,529,435,339,337,425,516,609,703,796,894,988,1086,1188,1293,1399};
    std::vector<double> phi1 = {M_PI,M_PI,M_PI,M_PI,M_PI,M_PI,M_PI,M_PI,M_PI,M_PI,M_PI,M_PI,M_PI,M_PI,M_PI,M_PI,M_PI,M_PI,M_PI,2*M_PI,2*M_PI,2*M_PI,2*M_PI,2*M_PI,2*M_PI,2*M_PI,2*M_PI,2*M_PI,2*M_PI,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.,M_PI/2.};
    auto ga = new TGraph2D(z1.size(), z1.data(), y1.data(), phi1.data());
    drawGraph(ga, c, "surf1");

    //slide 1- teorico formulas do paint inv abel , distribuições
    //slide 2- img fit do background + calibracao
    //slide 3 fit linha dos min, fit na linha com sin
    //slide 4- t vs x para 150 450 no mesmo grafico trig e gauss - expneg =1000+
    //slide 5- linhas infinitas 450C + conclusoes e comparacao - assimetria
    //slide 6- chama
    //slide final conclusoes
    double Magnification = 0.548; //adim
    double eMagnification = 0.032;
    double SolderR = 150./2.;  // pixel
    double eSolderR = 1;  // pixel
    double Omega_x = 0.0625089;
    double eOmega_x = 1.27672e-05;
    double PhiFactorMins = 0.123708;
    double ePhiFactorMins = 0.002654;
    double PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    double ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;
    CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,490.997, 1.41032,Magnification,eMagnification, SolderR, eSolderR,1);

    CalcTemp(0.0219296, 0.000228066,564.947, 2.65542,Magnification,eMagnification, SolderR, eSolderR);
    CalcTemp(0.0151286, 0.000291513,731.608, 7.7559,Magnification,eMagnification, SolderR, eSolderR);
    //calc of temp 900px
    /*

        g1d->SetTitle( "Light Intensity vs Y at X = 900px for 300C" );
    g1d->SetMarkerStyle(20);
    //g1d->SetMarkerColor(kAzure+2);
    g1d->SetLineColor(kBlue+2);
    g1d->SetMarkerSize(0.1);

    g1d->GetXaxis()->SetTitle("y [pixels]");
    g1d->GetXaxis()->SetLabelSize(0.028);


    g1d->GetYaxis()->SetTitle("Intensity [light value]");
    g1d->GetYaxis()->SetLabelSize(0.028);

    PhiFactorPx = 0.0125748; // units of 1/pixel  max ~0.33
    ePhiFactorPx = 0.000223054;
    PixelSize = 3.69e-6; //m/pixel
    ePixelSize =0.005e-6; //?
    Magnification = 0.434; //adim
    eMagnification = 0.002;
    PhiFactor = PhiFactorPx/PixelSize*Magnification; //1/m
    ePhiFactor = ePhiFactorPx/PixelSize*Magnification + PhiFactorPx/PixelSize/PixelSize*Magnification*ePixelSize + PhiFactorPx/PixelSize*eMagnification;
    Avogadro =  6.02214076e23; //1/mol
    Polarizability = 2.133e-29; //m^3
    ePolarizability = 0.032e-29;
    Lambda = 6.33e-7; //m
    eLambda = 0.005e-7;
    MolarMass = 28.96e-3; //kg/mol
    eMolarMass = 0.005e-3;
    DeltaRo = PhiFactor/4/M_PI/M_PI/Avogadro/Polarizability*Lambda*MolarMass; //kg/m^3
    eDeltaRo = ePhiFactor/4/M_PI/M_PI/Avogadro/Polarizability*Lambda*MolarMass + PhiFactor/4/M_PI/M_PI/Avogadro/Polarizability/Polarizability*Lambda*MolarMass*ePolarizability + PhiFactor/4/M_PI/M_PI/Avogadro/Polarizability*eLambda*MolarMass + PhiFactor/4/M_PI/M_PI/Avogadro/Polarizability*Lambda*eMolarMass;
    Ro = 1.225; // kg/m^3
    eRo = 0.0001;
    TempRatio = Ro/(Ro-DeltaRo); //adim
    eTempRatio = eRo/(Ro-DeltaRo) + Ro/(Ro-DeltaRo)/(Ro-DeltaRo)*eRo + Ro/(Ro-DeltaRo)/(Ro-DeltaRo)*eDeltaRo;
    TempBg = 25+273.15; // K
    eTempBg = 1;
    TempCentral = TempRatio*TempBg - 273.15; // ºC
    eTempCentral = eTempRatio*TempBg + TempRatio*eTempBg;
    std::cout << "Central temp is " << TempCentral << "+-" << eTempCentral << ".ºC" << std::endl;
*/

    //f = new TF2("f","x>[7] && y<([5]-[6])/(1447.-[9])*(x-[9])+[6] && y>([8]-[6])/(1447.-[9])*(x-[9])+[6] ? [10] : [0]*sin([1]*x+[2]*y+[3])+[4]",0,1500,0,2000);
    //(abs((y-([5]+[8])/2)/([9]-1447)*x-([5]+[8])/2+1447*([5]+[8])/2/([9]-1447))/sqrt(1+([5]+[8])/2/([9]-1447)*([5]+[8])/2/([9]-1447)))
    //double pars[15] = {13622.1,0.0658905,0.00127228,6.64916,19655.2,1150,985,600,825.703,300, 10000, -0.02,200,0.001, 200};
    //f = new TF2("f","x>[7] && y<([5]-[6])/(1447.-[9])*(x-[9])+[6] && y>([8]-[6])/(1447.-[9])*(x-[9])+[6] ? [10] : [0]*sin([1]*x+[2]*y+[3])+[4]",0,1500,0,2000);
    //d= dist from central line = 1447,([5]+[8])/2 to [9],[6] => y-([6]-([5]+[8])/2)/([9]-1447)*x-([5]+[8])/2+1447*([5]+[8])/2)/([9]-1447)=0
    //(abs((y-([5]+[8])/2)/([9]-1447)*x-([5]+[8])/2+1447*([5]+[8])/2/([9]-1447))/sqrt(1+([5]+[8])/2/([9]-1447)*([5]+[8])/2/([9]-1447)))
}


