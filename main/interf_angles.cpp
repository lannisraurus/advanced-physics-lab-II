#include <iostream>
#include <vector>
#include <math.h>

void calculate(double lambda, double pixelToDist, std::vector<double> ds, int angleNumber, std::vector<double> Ds, double L, double L_err){
	// Using line distances
	double d = 0;
	double d_err = 0;
	for (int i = 0; i < ds.size(); i++){
		d += ds[i];
	}	d /= ds.size();
	for (int i = 0; i < ds.size(); i++){
		double d_err_temp = d - ds[i];
		if (d_err_temp > d_err) d_err = d_err_temp;
	}
	d *= pixelToDist;
	d_err *= pixelToDist;
	double a = lambda/d;
	double a_err = lambda*d_err/(d*d);
	printf("Angle %d (line spacing): %f +- %f (rad) = %f +- %f (degrees).\n",angleNumber,a,a_err,a*180/M_PI,a_err*180/M_PI);
	// Using system lengths
	double D = 0;
	double D_err = 0;
	for (int i = 0; i < Ds.size(); i++){
		D += Ds[i];
	}	D /= Ds.size();
	for (int i = 0; i < Ds.size(); i++){
		double D_err_temp = D - Ds[i];
		if (D_err_temp > D_err) D_err = D_err_temp;
	}
	D *= pixelToDist;
	D_err *= pixelToDist;
	double a2 = D/L;
	double a2_err = D_err/L + D*L_err/(L*L);
	printf("Angle %d (system info): %f +- %f (rad) = %f +- %f (degrees).\n",angleNumber,a2,a2_err,a2*180/M_PI,a2_err*180/M_PI);


}

int main(){

	// INPUT PARAMETERS
	double pixelToDist = 3.69e-6;	// Convert pixels to distances (m) (amplification factor, CCD pixels)
	double lambda = 633e-9;     	// Wavelength of the laser.
	double L = 40.1e-2;					// Length of the interferometer arm (m).
	double L_err = 0.1e-2;				// Error in the measure of the interferometer arm (m).
	// CALCULATE FOR MULTIPLE ANGLES
	// dXs: distance between lines; DXs: distance foci;
	
	// Angle 1 
	std::vector<double> d1s = {198.28, 198.82, 204.35, 177.91, 184.88, 194.83};
	std::vector<double> D1s = {82.22,  86.21,  85.15,  87.14,  86.02 };
	calculate(lambda,pixelToDist,d1s,1,D1s,L,L_err);

	// Angle 2
	std::vector<double> d2s = {101.41, 99.05,  95.15,  91.32,  88.84,  99.36};
	std::vector<double> D2s = {197.41, 198.21, 197.06, 194.31, 200.66 };
	calculate(lambda,pixelToDist,d2s,1,D2s,L,L_err);

	// Angle 3
	std::vector<double> d3s = {243.40, 246.26, 248.07, 244.97, 223.79, 241.40};
	std::vector<double> D3s = {71.59,  75.17,  74.33,  76.28,  74.20};
	calculate(lambda,pixelToDist,d3s,1,D3s,L,L_err);

	// Angle 3
	std::vector<double> d4s = {52.72,  50.20,  53.90,  54.94,  54.30, 48.67};
	std::vector<double> D4s = {397.68, 395.36, 396.76, 396.43, 401.42};
	calculate(lambda,pixelToDist,d4s,1,D4s,L,L_err);

	return 0;
}
