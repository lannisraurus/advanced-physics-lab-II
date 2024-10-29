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
	double L = 12.1e-2;					// Length of the interferometer arm (m).
	double L_err = 0.1e-2;				// Error in the measure of the interferometer arm (m).
	// CALCULATE FOR MULTIPLE ANGLES
	// dXs: distance between lines; DXs: distance foci;
	std::vector<double> d1s = {198.28, 198.82, 204.35, 177.91, 184.88, 194.83};
	std::vector<double> D1s = {90.35,  82.22,  86.21,  85.15,  87.14,  86.02 };
	calculate(lambda,pixelToDist,d1s,1,D1s,L,L_err);

	return 0;
}
