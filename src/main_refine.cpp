/**********************************************************************************//**
 * \file main_refine.cpp
 *
 * \author Kjetil A. Johannessen
 * \date November 2010
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

int    eff_patch          = -1;
int    prefDir            = -1;
int    prefAmount         = -1;
int    hrefDir            = -1;
int    edge               = -1;
double r                  = -1;
double nBoundary          = -1;
vector<double> hrefKnot   ;
bool   uniform_u          = false;
bool   uniform_v          = false;
bool   uniform_w          = false;
bool   regularize         = false;
bool   glob_refine        = false;
bool   local_refine       = false;
string fileUsage          = "\
File usage: refine [-href <dir> <knot>] [-pref <dir> <n>] [-patch <index>] \n\
                   [-boundary <edge> <r>] \n\
                   [-unif <dir>] [-regul] <inputFile> \n\
  \n\
  Arguments\n\
    <inputFile> : one or more .g2-files describing the spline volumes \n\
  FLAGS\n\
    -help                : display file usage (this screen) \n\
  GLOBAL REFINEMENT FLAGS\n\
    -unif <dir>          : uniform refinement (h-refinement) in the <dir> parametric\n\
	                       direction (-1 for all).\n\
    -regul               : regularize biased patches by ensuring that no element \n\
                           is larger than twice the size of the smallest element \n\
                           as well as making all patches of equal degree (highest) \n\
  LOCAL REFINEMENT FLAGS\n\
    -patch <index>       : apply refinement to the patch <index> \n\
    -href <dir> <knot>   : insert <knot(s)> into the <dir> parametric knot vector  \n\
    -pref <dir> <n>      : raise the order in parametric direction <dir> by <n> \n\
    -bndry <edge> <r> <n>: resolve boundary layer by splitting the elements nearest \n\
                           local edge number <edge> in <n> new knot lines an aspect \n\
						   ratio of <r>. With <r> being between 0 and 1\n";
	/*
    -u <knot>    : insert <knot> into the first parametric knot vector  \n\
    -v <knot>    : insert <knot> into the second parametric knot vector  \n\
    -w <knot>    : insert <knot> into the third parametric knot vector  \n\
	*/
SplineModel model;

void processParameters(int argc, char** argv) {
	bool fileRead = false;
	for(int argi=1; argi<argc; argi++) {
		const char *flag = argv[argi];
		char *flag_arg;
		if(strcmp(flag, "-href") == 0) {
			flag_arg = argv[++argi];
			hrefDir = atoi(flag_arg);
			char *oneArg;
			flag_arg = argv[++argi];
			oneArg = strtok(flag_arg, ",");
			while(oneArg != NULL) {
				hrefKnot.push_back(atof(oneArg));
				oneArg = strtok(NULL, ",");
			}
			local_refine = true;

		} else if(strcmp(flag, "-pref") == 0) {
			flag_arg = argv[++argi];
			prefDir = atoi(flag_arg);
			flag_arg = argv[++argi];
			prefAmount = atoi(flag_arg);
			local_refine = true;

		} else if(strcmp(flag, "-patch") == 0) {
			flag_arg = argv[++argi];
			eff_patch = atof(flag_arg);

		} else if(strcmp(flag, "-bndry") == 0) {
			edge      = atoi(argv[++argi]);
			r         = atof(argv[++argi]);
			nBoundary = atoi(argv[++argi]);
			local_refine = true;

		} else if(strcmp(flag, "-regul") == 0) {
			regularize  = true;
			glob_refine = true;
		} else if(strcmp(flag, "-unif") == 0) {
			flag_arg = argv[++argi];
			int dir = atoi(flag_arg);
			if(dir==-1) {
				uniform_u = true;
				uniform_v = true;
				uniform_w = true;
			} else if(dir==0) {
				uniform_u = true;
			} else if(dir==1) {
				uniform_v = true;
			} else if(dir==2) {
				uniform_w = true;
			}
				
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
	
	if(glob_refine && local_refine) {
		cerr << "Can't both perform local and global refinements\n";
		exit(2);
	}
	if((uniform_u || uniform_v || uniform_w)  && regularize) {
		cerr << "Can't do both uniform refinement and regularization at once\n";
		exit(1);
	}
	if(model.isVolumetricModel()) {
		cerr << "Volumetric model refinement needs some love. Avoid using it for the moment\n";
		exit(1154152);
		vector<shared_ptr<SplineVolume> > volumes = model.getSplineVolumes();
		if(eff_patch == -1) {
			for(uint i=0; i<volumes.size(); i++) {
				if(uniform_u || uniform_v || uniform_w) {
					vector<double> uniqueKnots;
					for(int dir=0; dir<3; dir++) {
						uniqueKnots.clear();
						volumes[i]->basis(dir).knotsSimple(uniqueKnots);
						for(uint j=0; j<uniqueKnots.size()-1; j++)
							volumes[i]->insertKnot(dir, (uniqueKnots[j]+uniqueKnots[j+1])/2);
					}
				} else if(regularize) {
					cerr << "Regularization not implemented yet. Come back later\n";
					exit(0);
				} else {
					if(hrefDir > -1)
						for(vector<double>::iterator it=hrefKnot.begin(); it!=hrefKnot.end(); it++)
							volumes[i]->insertKnot(hrefDir, *it);
					if(prefDir == 0 )
						volumes[i]->raiseOrder(prefAmount,0,0);
					else if(prefDir == 1 )
						volumes[i]->raiseOrder(0,prefAmount, 0);
					else if(prefDir == 2 )
						volumes[i]->raiseOrder(0,0,prefAmount);
				}
			}
		}
	} else {
		vector<shared_ptr<SplineSurface> > surfaces = model.getSplineSurfaces();
		if(eff_patch==-1) {
			for(uint i=0; i<surfaces.size(); i++) {
				if(uniform_u || uniform_v) {
					if(!(uniform_u && uniform_v)) // specified global uniform reffinement in one parametric direction
						cerr << "WARNING: global uniform refinement in one direction is not guaranteed to give a consistent model back, due to potential local orientations. Uniform refinement in all parametric directions was performed\n";
					vector<double> uniqueKnots_u;
					surfaces[i]->basis_u().knotsSimple(uniqueKnots_u);
					for(uint j=0; j<uniqueKnots_u.size()-1; j++)
						surfaces[i]->insertKnot_u((uniqueKnots_u[j]+uniqueKnots_u[j+1])/2);
					vector<double> uniqueKnots_v;
					surfaces[i]->basis_v().knotsSimple(uniqueKnots_v);
					for(uint j=0; j<uniqueKnots_v.size()-1; j++)
						surfaces[i]->insertKnot_v((uniqueKnots_v[j]+uniqueKnots_v[j+1])/2);
				}
				if(regularize) {
					cerr << "Regularization not implemented yet. Come back later\n";
					exit(0);
				} else {
					/*
					if(hrefDir == 0)
						for(vector<double>::iterator it=hrefKnot.begin(); it!=hrefKnot.end(); it++)
							surfaces[i]->insertKnot_u(*it);
					else if(hrefDir == 1)
						for(vector<double>::iterator it=hrefKnot.begin(); it!=hrefKnot.end(); it++)
							surfaces[i]->insertKnot_v(*it);
					*/
					if(prefDir == 0 )
						surfaces[i]->raiseOrder(prefAmount,0);
					else if(prefDir == 1 )
						surfaces[i]->raiseOrder(0,prefAmount);
				}
			}
		} else {
			if(hrefDir > -1) {
				model.knot_insert(eff_patch, hrefDir, hrefKnot[0]);
			}
			if(edge > -1) {
				model.boundary_layer_refinement(eff_patch, 1-edge/2, 1-edge%2, r, nBoundary);
			}
			if(uniform_u) {
				vector<double> uniqueKnots_u;
				surfaces[eff_patch]->basis_u().knotsSimple(uniqueKnots_u);
				for(uint j=0; j<uniqueKnots_u.size()-1; j++)
					model.knot_insert(eff_patch, 0, (uniqueKnots_u[j]+uniqueKnots_u[j+1])/2);
			}
			if(uniform_v) {
				vector<double> uniqueKnots_v;
				surfaces[eff_patch]->basis_v().knotsSimple(uniqueKnots_v);
				for(uint j=0; j<uniqueKnots_v.size()-1; j++)
					model.knot_insert(eff_patch, 1, (uniqueKnots_v[j]+uniqueKnots_v[j+1])/2);
			}
		}
	}

	cout << setprecision(17);
	model.writeSplines(cout);
	return 0;
}

