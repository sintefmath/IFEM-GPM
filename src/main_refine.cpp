/**********************************************************************************//**
 * \file main_getGNO.cpp
 *
 * \author Kjetil A. Johannessen
 * \date August 2010
 *
 *************************************************************************************/

#include <iostream>
#include <fstream>
#include <string.h>

#include "SplineModel.h"
#include "TopologySet.h"
#include "primitives.h"

#include <GoTools/trivariate/SplineVolume.h>

typedef unsigned int uint;

using namespace std;
using namespace Go;
using boost::shared_ptr;

vector<double> uKnot       ;
vector<double> vKnot       ;
vector<double> wKnot       ;
vector<int>    eff_volumes ;
bool   uniform             = false;
string fileUsage           = "\
File usage: gpm [-u <arg>] [-v <arg>] [-w <arg>] [-all] <inputFile> \n\
  \n\
  Arguments\n\
    <inputFile> : one or more .g2-files describing the spline volumes \n\
  FLAGS\n\
    -help        : display file usage (this screen) \n\
    -all         : uniform refinement (not working yet) \n\
    -u <knot>    : insert <knot> into the first parametric knot vector  \n\
    -v <knot>    : insert <knot> into the second parametric knot vector  \n\
    -w <knot>    : insert <knot> into the third parametric knot vector  \n\
    -vol <index> : only apply refinement to the volume <index> \n";
SplineModel model;

void processParameters(int argc, char** argv) {
	bool fileRead = false;
	for(int argi=1; argi<argc; argi++) {
		const char *flag = argv[argi];
		char *flag_arg;
		if(strcmp(flag, "-u") == 0) {
			flag_arg = argv[++argi];
			char *oneArg;
			oneArg = strtok(flag_arg, ",");
			while(oneArg != NULL) {
				uKnot.push_back(atof(oneArg));
				oneArg = strtok(NULL, ",");
			}
		} else if(strcmp(flag, "-v") == 0) {
			flag_arg = argv[++argi];
			char *oneArg;
			oneArg = strtok(flag_arg, ",");
			while(oneArg != NULL) {
				vKnot.push_back(atof(oneArg));
				oneArg = strtok(NULL, ",");
			}
		} else if(strcmp(flag, "-w") == 0) {
			flag_arg = argv[++argi];
			char *oneArg;
			oneArg = strtok(flag_arg, ",");
			while(oneArg != NULL) {
				wKnot.push_back(atof(oneArg));
				oneArg = strtok(NULL, ",");
			}
		} else if(strcmp(flag, "-vol") == 0) {
			flag_arg = argv[++argi];
			char *oneArg;
			oneArg = strtok(flag_arg, ",");
			while(oneArg != NULL) {
				eff_volumes.push_back(atof(oneArg));
				oneArg = strtok(NULL, ",");
			}
		} else if(strcmp(flag, "-all") == 0) {
			uniform = true;
		} else if(strcmp(flag, "-help") == 0) {
			cout << fileUsage << endl;
			exit(0);
		} else {
			ifstream inFile;
			inFile.open(flag);
			if(!inFile.good()) {
				cerr << "Error reading input file \"" << flag << "\"" << endl;
				exit(1);
			}
			model.readSplines(inFile);
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
 * \brief Main program for command-line patch refinement
 *************************************************************************************/

int main(int argc, char **argv) {
	processParameters(argc, argv);
	
	vector<shared_ptr<SplineVolume> > volumes = model.getSplines();
	for(uint i=0; i<volumes.size(); i++) {
		if(eff_volumes.size() > 0) {
			bool found = false;
			for(uint j=0; j<eff_volumes.size(); j++)
				if(eff_volumes[j] == (int) i)
					found = true;

			if(!found)
				continue;
		}
		if(uniform) {
			cerr << "Uniform refinement not implemented yet. Come back later\n";
			exit(0);
		}
		for(uint j=0; j<uKnot.size(); j++)
			volumes[i]->insertKnot(0, uKnot[j]);
		for(uint j=0; j<vKnot.size(); j++)
			volumes[i]->insertKnot(1, vKnot[j]);
		for(uint j=0; j<wKnot.size(); j++)
			volumes[i]->insertKnot(2, wKnot[j]);
	}

	model.writeSplines(cout);
	return 0;
}

