// C++ Libs
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>

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
    }

    double TempRatio = (AirN-1)/(AirN-1-DeltaNSolder); //adim
    double eTempRatio = eAirN/(AirN-1-DeltaNSolder)+(AirN-1)/(AirN-1-DeltaNSolder)/(AirN-1-DeltaNSolder)*eAirN + (AirN-1)/(AirN-1-DeltaNSolder)/(AirN-1-DeltaNSolder)*eDeltaNSolder;
    double Temp = TempRatio*(TempBg + 273.15) - 273.15; // ºC
    double eTemp = eTempRatio*(TempBg + 273.15) + TempRatio*eTempBg;
    std::cout << "Solder temp is " << Temp << "+-" << eTemp << ".ºC" << std::endl;
    return {Temp,eTemp};
}

int main(){

    std::pair<double,double> T;

    std::vector<double> x = {278,386,484,584,684,782,880,974,1072,1170,1262,1362};
    std::vector<double> error_x = {5,5,5,5,5,5,5,5,5,5,5,5};
    std::vector<double> y;
    std::vector<double> error_y;
    std::vector<double> y2;
    std::vector<double> error_y2;

    //------------line1
    double Magnification = 0.548; //adim
    double eMagnification = 0.032;
    double SolderR = 290/2;  // pixel
    double eSolderR = 3;  // pixel
    double Omega_x = 0.0625089;
    double eOmega_x = 1.27672e-05;
    double PhiFactorMins = 0.384681/2;
    double ePhiFactorMins = 0.00180885/2;
    double PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    double ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,634.689, 1.54514,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //------------line2
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 288/2;  // pixel
    eSolderR = 10;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.374803/2;
    ePhiFactorMins = 0.00179866/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,619.39, 1.54251,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //------------line3
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 298/2;  // pixel
    eSolderR = 10;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.37746/2;
    ePhiFactorMins = 0.00185257/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,595.433, 1.45801,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //------------line4
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 284/2;  // pixel
    eSolderR = 10;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.378921/2;
    ePhiFactorMins = 0.00201571/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,572.631, 1.49295,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //------------line5
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 240/2;  // pixel
    eSolderR = 7;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.380034/2;
    ePhiFactorMins = 0.00223806/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,546.075, 1.49233,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //------------line6
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 204/2;  // pixel
    eSolderR = 7;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.369891/2;
    ePhiFactorMins = 0.00241007/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,520.641, 1.52401,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //------------line7
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 180/2;  // pixel
    eSolderR = 5;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.355938/2;
    ePhiFactorMins = 0.00246803/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,494.429, 1.47395,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //------------line8
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 147/2;  // pixel
    eSolderR = 5;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.337995/2;
    ePhiFactorMins = 0.00304792/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,463.552, 1.7202,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //------------line9
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 136/2;  // pixel
    eSolderR = 5;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.292207/2;
    ePhiFactorMins = 0.00311116/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,436.294, 1.8292,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //------------line10
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 80/2;  // pixel
    eSolderR = 3;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.256263/2;
    ePhiFactorMins = 0.00400291/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,390.072, 2.13863,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //-------------linePontinha
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 0;  // pixel
    eSolderR = 0;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.320412/2;
    ePhiFactorMins = 0.00529753/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,312.153, 2.30073,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //-------------lineFORA1
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 0;  // pixel
    eSolderR = 0;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.0617671/2;
    ePhiFactorMins = 0.000993929/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,424.596, 4.6243,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);


    //gauss
    //
    //
//------------line1
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 290/2;  // pixel
    eSolderR = 3;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.382917/2;
    ePhiFactorMins = 0.00316907/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,766.686, 3.72221,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //------------line2
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 288/2;  // pixel
    eSolderR = 10;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.374148/2;
    ePhiFactorMins = 0.0022899/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,747.099, 2.42344,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //------------line3
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 298/2;  // pixel
    eSolderR = 10;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.379193/2;
    ePhiFactorMins = 0.00209599/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,716.241, 1.9462,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //------------line4
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 284/2;  // pixel
    eSolderR = 10;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.386557/2;
    ePhiFactorMins = 0.00196252/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected, 683.743, 1.6146,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //------------line5
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 240/2;  // pixel
    eSolderR = 7;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.393124/2;
    ePhiFactorMins = 0.00296677/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,648.229, 2.35362,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //------------line6
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 204/2;  // pixel
    eSolderR = 7;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.386135/2;
    ePhiFactorMins = 0.00121139/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,616.087, 0.053623,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //------------line7
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 180/2;  // pixel
    eSolderR = 5;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.381208/2;
    ePhiFactorMins = 0.00221991/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,579.362, 1.34264,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //------------line8
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 147/2;  // pixel
    eSolderR = 5;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.366369/2;
    ePhiFactorMins = 0.0031516/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,542.427, 1.88866,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //------------line9
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 136/2;  // pixel
    eSolderR = 5;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.31954/2;
    ePhiFactorMins = 0.00263604/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected, 510.107, 1.52665,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //------------line10
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 80/2;  // pixel
    eSolderR = 3;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.315451/2;
    ePhiFactorMins = 0.00303194/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,439.212, 1.17479,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //-------------linePontinha
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 0;  // pixel
    eSolderR = 0;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.257713/2;
    ePhiFactorMins = 0.00182489/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,378.816, 1.71379,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //-------------lineFORA1
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 0;  // pixel
    eSolderR = 0;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.0558898/2;
    ePhiFactorMins = 0.000606522/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T = CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,537.622, 3.27015,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);

    //
    //
    //

    TApplication *A = new TApplication("a",0,0);
    TCanvas *c = new TCanvas("Canvas name", "c", 400, 100, 1000, 800);
    TGraphErrors *g = new TGraphErrors(x.size(),x.data(),y.data(),error_x.data(),error_y.data());
    TGraphErrors *g2 = new TGraphErrors(x.size(),x.data(),y2.data(),error_x.data(),error_y2.data());

        gStyle->SetOptFit(111);

    g->SetTitle( "Temperature along the soldering iron at 450C");
    g -> GetXaxis() -> SetTitle(" x [px] ");
    g -> GetYaxis() -> SetTitle(" T [Celsius] ");
    g -> GetXaxis() -> SetRangeUser(100,1400);
    g -> GetYaxis() -> SetRangeUser (0,600);
    g -> SetMarkerStyle(8);
    g -> SetMarkerSize(2);
    g -> SetMarkerColor(kRed);
    g -> SetLineColor(kRed);
    
    g2 -> SetMarkerStyle(8);
    g2 -> SetMarkerSize(2);
    g2 -> SetMarkerColor(kBlue);
    g2 -> SetLineColor(kBlue);
   g->SetFillColor(1);
   g->SetFillStyle(3001);
      g2->SetFillColor(2);
   g2->SetFillStyle(3001);
    g->Draw("AP");
    g2->Draw("P");
    c->Update();
    gSystem->ProcessEvents();
    c->WaitPrimitive();
    c->SaveAs("temp450.png");


}