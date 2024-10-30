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

    std::vector<double> x = {144,278,386,484,584,684,782,880,974,1072,1170,1262,1362};
    std::vector<double> error_x = {5,5,5,5,5,5,5,5,5,5,5,5,5};
    std::vector<double> y;
    std::vector<double> error_y;
    std::vector<double> y2;
    std::vector<double> error_y2;


    //------------line1
    double Magnification = 0.548; //adim
    double eMagnification = 0.032;
    double SolderR = 286/2;  // pixel
    double eSolderR = 10;  // pixel
    double Omega_x = 0.0625089;
    double eOmega_x = 1.27672e-05;
    double PhiFactorMins = 0.179572/2;
    double ePhiFactorMins = 0.00257488/2;
    double PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    double ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,510.969, 3.22017,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //------------line2
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 292/2;  // pixel
    eSolderR = 10;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.181654/2;
    ePhiFactorMins = 0.00274667/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,499.78, 3.24311,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //------------line3
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 298/2;  // pixel
    eSolderR = 10;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.180568/2;
    ePhiFactorMins = 0.00318638/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,482.623, 3.64431,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //------------line4
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 282/2;  // pixel
    eSolderR = 10;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.181516/2;
    ePhiFactorMins = 0.00278004/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,466.844, 2.90996,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //------------line5
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 270/2;  // pixel
    eSolderR = 7;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.172466/2;
    ePhiFactorMins = 0.00274871/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,449.95, 2.79788,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //------------line6
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 232/2;  // pixel
    eSolderR = 7;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.172466/2;
    ePhiFactorMins = 0.00274871/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,449.95, 2.79788,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //------------line7
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 222/2;  // pixel
    eSolderR = 5;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.155736/2;
    ePhiFactorMins = 0.00262469/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,443.584, 2.89887,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //------------line8
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 163/2;  // pixel
    eSolderR = 5;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.133506/2;
    ePhiFactorMins = 0.00287105/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,443.916, 3.70191,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //------------line9
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 144/2;  // pixel
    eSolderR = 5;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.124238/2;
    ePhiFactorMins = 0.00262729/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,434.634, 3.45774,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //------------line10
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 99/2;  // pixel
    eSolderR = 3;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.116578/2;
    ePhiFactorMins = 0.00330836/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,402.853, 4.11226,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //-------------linePontinha
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 0;  // pixel
    eSolderR = 0;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.0796616;
    ePhiFactorMins = 0.00434;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,344.281, 7.00448,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //-------------lineFORA1
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 0;  // pixel
    eSolderR = 0;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.049636;
    ePhiFactorMins = 0.00116529;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,373.47, 5.66407,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);
    //-------------lineFORA2
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 0;  // pixel
    eSolderR = 0;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.0165562;
    ePhiFactorMins = 0.000370531;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,692.695, 11.0706,Magnification,eMagnification, SolderR, eSolderR,1);
    y.push_back(T.first);
    error_y.push_back(T.second);

    
    //gauss
    //
    //

    //------------line1
     Magnification = 0.548; //adim
     eMagnification = 0.032;
     SolderR = 286/2;  // pixel
     eSolderR = 10;  // pixel
     Omega_x = 0.0625089;
     eOmega_x = 1.27672e-05;
     PhiFactorMins = 0.195473/2;
     ePhiFactorMins = 0.00145831/2;
     PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
     ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,592.456, 1.09208,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //------------line2
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 292/2;  // pixel
    eSolderR = 10;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.199526/2;
    ePhiFactorMins = 0.00134974/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,577.971, 0.134368,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //------------line3
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 298/2;  // pixel
    eSolderR = 10;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.200589/2;
    ePhiFactorMins = 0.00489528/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,557.016, 5.81116,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //------------line4
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 282/2;  // pixel
    eSolderR = 10;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.206894/2;
    ePhiFactorMins = 0.00153103/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,533.809, 0.186316,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //------------line5
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 270/2;  // pixel
    eSolderR = 7;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.194012/2;
    ePhiFactorMins = 0.00157317/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,518.333, 0.382636,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //------------line6
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 232/2;  // pixel
    eSolderR = 7;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.17918/2;
    ePhiFactorMins = 0.00304859/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,506.477, 3.16138,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //------------line7
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 222/2;  // pixel
    eSolderR = 5;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.159107/2;
    ePhiFactorMins = 0.00178607/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,499.884, 0.638899,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //------------line8
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 163/2;  // pixel
    eSolderR = 5;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.148258/2;
    ePhiFactorMins = 0.00644615/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,490.228, 8.34491,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //------------line9
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 144/2;  // pixel
    eSolderR = 5;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.146099/2;
    ePhiFactorMins = 0.0017021/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,448.659, 0.22238,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //------------line10
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 99/2;  // pixel
    eSolderR = 3;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.0948546;
    ePhiFactorMins = 0.000715138/2;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,475.629, 2.53778,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //-------------linePontinha
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 0;  // pixel
    eSolderR = 0;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.0960447;
    ePhiFactorMins = 0.00365568;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,387.126, 4.92896,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //-------------lineFORA1
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 0;  // pixel
    eSolderR = 0;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.0470579;
    ePhiFactorMins = 0.000668326;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,453.894,2.38037,Magnification,eMagnification, SolderR, eSolderR,2);
    y2.push_back(T.first);
    error_y2.push_back(T.second);
    //---------------lineFORA2
    Magnification = 0.548; //adim
    eMagnification = 0.032;
    SolderR = 0;  // pixel
    eSolderR = 0;  // pixel
    Omega_x = 0.0625089;
    eOmega_x = 1.27672e-05;
    PhiFactorMins = 0.016619;
    ePhiFactorMins = 0.000285149;
    PhiFactorMinsCorrected = PhiFactorMins*Omega_x;
    ePhiFactorMinsCorrected = ePhiFactorMins*Omega_x+PhiFactorMins*eOmega_x;

    T=CalcTemp(PhiFactorMinsCorrected, ePhiFactorMinsCorrected,801.976,3.50188,Magnification,eMagnification, SolderR, eSolderR,2);
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

    g->SetTitle( "Temperature along the soldering iron at 150C");
    g -> GetXaxis() -> SetTitle(" x [px] ");
    g -> GetYaxis() -> SetTitle(" T [Celsius] ");
    g -> GetXaxis() -> SetRangeUser(100,1400);
    g -> GetYaxis() -> SetRangeUser (0,200);
    g -> SetMarkerStyle(8);
    g -> SetMarkerSize(2);
    g -> SetMarkerColor(kRed);
    g -> SetLineColor(kRed);
    //gStyle->SetTitleFontSize(0.1);
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
    c->SaveAs("temp150.png");


}