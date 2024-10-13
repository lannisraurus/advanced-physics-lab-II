// C++ Libs
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>

int main(){

	//---------------------------------------------------------------
	//                      INPUT PARAMETERS
	//---------------------------------------------------------------
	
	// Data file info
	std::string dataFile = "bin/data_in/Interf2_calha_d3p0mm.csv";	// Input data file, pgm

	//---------------------------------------------------------------
	//                           ALGORITHM
	//---------------------------------------------------------------

	// Data Vars
	std::ifstream input;
	input.open(dataFile);
	
	// Probing Vars
	double pixel, value;
	double v[3] = {0.0, 0.0, 0.0};
	std::string line;

	// Store Vars
	std::vector<double> pixels;
	std::vector<double> values;
	std::vector<double> maximums;
	std::vector<double> minimums;
	double avgMax, avgMin;
	double maxErr = 0, minErr = 0,
		   maxErrTemp, minErrTemp;
	double k, kErr;
	
	// Error Checking
	if(!input.is_open()){
		printf("Error: Couldn't open file.\n");
		return -1;
	}

	// Read File - Load information into pixels and values
	std::getline(input, line); // Get labels (useless)
	while(!input.eof()){
		input >> pixel >> value;
		pixels.push_back(pixel);
		values.push_back(value);
	}

	// Maximum and minimum finder (ignores values at boundaries)
	for (int i = 1; i < pixels.size() - 1; i++){
		v[0] = values[i-1];
		v[1] = values[i];
		v[2] = values[i+1];
		if (v[1] > v[0] && v[1] > v[2]) maximums.push_back(v[1]);
		if (v[1] < v[0] && v[0] < v[2]) minimums.push_back(v[1]);
	}
	
	// Find average max
	for (int i = 0; i < maximums.size(); i++){
		avgMax += maximums[i];
	}	avgMax /= maximums.size();
	
	// Find average min
	for (int i = 0; i < minimums.size(); i++){
		avgMin += minimums[i];
	}	avgMin /= minimums.size();

	// Find average max error
	for (int i = 0; i < maximums.size(); i++){
		maxErrTemp = std::abs(avgMax - maximums.size());
		if (maxErrTemp > maxErr) maxErr = maxErrTemp;
	}

	// Find average min error
	for (int i = 0; i < minimums.size(); i++){
		minErrTemp = std::abs(avgMin - minimums.size());
		if (minErrTemp > minErr) minErr = minErrTemp;
	}

	// Print results of averages max and min and their errors.
	printf("Average Maximum: %f.5 +/- %f.5 [light value]\n", avgMax, avgMaxErr);
	printf("Average Minimum: %f.5 +/- %f.5 [light value]\n", avgMin, avgMinErr);

	// Find contrast
	k = (avgMax - avgMin) / (avgMax + avgMin);


	return 0;

}
