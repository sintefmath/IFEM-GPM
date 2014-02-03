/**********************************************************************************//**
 * \file main_enforceRHS.cpp
 *
 * \author Kjetil A. Johannessen
 * \date March 2012
 *
 *************************************************************************************/

#include <iostream>
#include <fstream>
#include <string.h>

#include "SplineModel.h"
#include "TopologySet.h"
#include "primitives.h"

#include <GoTools/trivariate/SplineVolume.h>

using namespace std;
using namespace Go;

bool storeFile         = false;
const char *outputFile = new char[1024];
string fileUsage       = "\
File usage: enforceRHS [-out=<outfile>] <inputFile> \n\
  \n\
  Arguments\n\
    <inputFile> : one or more .g2-files describing the spline volumes \n\
  FLAGS\n\
    -out        : stores the resulting spline model in the file <outfile>\n\
    -help       : display this help screen";
SplineModel model;

void processParameters(int argc, char** argv) {
	bool fileRead = false;
	for(int argi=1; argi<argc; argi++) {
		const char *arg = argv[argi];
		if(strcmp(arg, "-help") == 0) {
			cout << fileUsage << endl;
			exit(0);
		} else if(strncmp(arg, "-out=", 5) == 0) {
			outputFile = arg+5;
			storeFile  = true;
		} else {
			ifstream inFile;
			inFile.open(arg);
			if(!inFile.good()) {
				cerr << "Error reading input file \"" << arg << "\"" << endl;
				exit(1);
			}
			model.readSplines(inFile, false);
			inFile.close();
			fileRead = true;
		}
	}
	if(!fileRead) {
		cout << fileUsage << endl;
		exit(1);
	}
}

/**********************************************************************************//**
 * \brief Main program for enforcing right-hand side parametrization
 *
 * This evalutates the determinant of the jacobian at the midpoint of all spline patches.
 * For all patches with negative determinant, the last parameter direction is reversed
 * (i.e. for surfaces v-direction is reversed, and for volumes the w-direction). 
 *
 *************************************************************************************/

int main(int argc, char **argv) {
	processParameters(argc, argv);
	if( model.enforceRightHandSystem() ) {
		// so some information print?
	}
	if(storeFile) {
		ofstream outFile;
		outFile.open(outputFile);
		model.writeSplines(outFile);
		outFile.close();
	} else {
		model.writeSplines(cout);
	}

	return 0;
}

