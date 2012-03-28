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
int    hrefN              = -1;
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
                   [-unif <dir> <n>] [-regul] <inputFile> \n\
  \n\
  Arguments\n\
    <inputFile> : one or more .g2-files describing the spline volumes \n\
  FLAGS\n\
    -help                : display file usage (this screen) \n\
  GLOBAL REFINEMENT FLAGS\n\
    -unif <dir> <n>      : uniform refinement (h-refinement) in the <dir> parametric\n\
	                       direction (-1 for all). Inserts <n> knots in each knot span\n\
    -regul               : regularize biased patches by ensuring that no element \n\
                           is larger than twice the size of the smallest element \n\
                           as well as making all patches of equal degree (highest) \n\
  LOCAL REFINEMENT FLAGS\n\
    -patch <index>       : apply refinement only to the patch <index> \n\
    -href <dir> <knot>   : insert <knot(s)> into the <dir> parametric knot vector  \n\
    -pref <dir> <n>      : raise the order in parametric direction <dir> (-1 for all\n\
	                       parameteric directions) by <n> polynomial degrees \n\
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
			hrefDir  = atoi(flag_arg);
			char *oneArg;
			flag_arg = argv[++argi];
			oneArg = strtok(flag_arg, ",");
			while(oneArg != NULL) {
				hrefKnot.push_back(atof(oneArg));
				oneArg = strtok(NULL, ",");
			}
			local_refine = true;

		} else if(strcmp(flag, "-pref") == 0) {
			prefDir    = atoi(argv[++argi]);
			prefAmount = atoi(argv[++argi]);
			local_refine = true;

		} else if(strcmp(flag, "-patch") == 0) {
			eff_patch = atof(argv[++argi]);
		} else if(strcmp(flag, "-bndry") == 0) {
			edge      = atoi(argv[++argi]);
			r         = atof(argv[++argi]);
			nBoundary = atoi(argv[++argi]);
			local_refine = true;

		} else if(strcmp(flag, "-regul") == 0) {
			regularize  = true;
			glob_refine = true;
		} else if(strcmp(flag, "-unif") == 0) {
			int dir   = atoi(argv[++argi]);
			hrefN     = atoi(argv[++argi]);
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
		vector<shared_ptr<SplineVolume> > volumes = model.getSplineVolumes();
		if(eff_patch == -1) {
			if(prefDir != -1) {// specified global uniform reffinement in one parametric direction
				cerr << "WARNING: global order elevation in one direction is not guaranteed to give a consistent model back, due to potential local orientations.";
				cerr << "         Uniform order elevations in all parametric directions was performed\n";
				prefDir = -1;
			}
			if(uniform_u || uniform_v || uniform_w) {
				if(!(uniform_u && uniform_v && uniform_w)) {// specified global uniform reffinement in one parametric direction
					cerr << "WARNING: global uniform refinement in one direction is not guaranteed to give a consistent model back, due to potential local orientations.";
					cerr << "         Uniform refinement in all parametric directions was performed\n";
					uniform_u = uniform_v = uniform_w = true;
				}
			}

			for(uint i=0; i<volumes.size(); i++) {
				if(uniform_u || uniform_v || uniform_w) {
					vector<double> uniqueKnots;
					for(int dir=0; dir<3; dir++) {
						uniqueKnots.clear();
						volumes[i]->basis(dir).knotsSimple(uniqueKnots);
						for(uint j=0; j<uniqueKnots.size()-1; j++)
							for(int k=0; k<hrefN; k++)
								volumes[i]->insertKnot(dir, uniqueKnots[j]*(k+1)/(hrefN+1) + uniqueKnots[j+1]*(hrefN-k)/(hrefN+1));
					}
				} else if(regularize) {
					cerr << "Regularization not implemented yet. Come back later\n";
					exit(0);
				} else {
					// if(hrefDir > -1)
						// for(vector<double>::iterator it=hrefKnot.begin(); it!=hrefKnot.end(); it++)
							// volumes[i]->insertKnot(hrefDir, *it);
					if(prefAmount>0 && (prefDir==-1 || prefDir == 0 ))
						volumes[i]->raiseOrder(prefAmount,0,0);
					if(prefAmount>0 && (prefDir==-1 || prefDir == 1 ))
						volumes[i]->raiseOrder(0,prefAmount, 0);
					if(prefAmount>0 && (prefDir==-1 || prefDir == 2 ))
						volumes[i]->raiseOrder(0,0,prefAmount);
				}
			}
		} else {
			if(hrefDir > -1) {
				model.knot_insert(eff_patch, hrefDir, hrefKnot[0]);
			}
			if(edge > -1) {
				model.boundary_layer_refinement(eff_patch, edge/2, 1-edge%2, r, nBoundary);
			}
			if(uniform_u) {
				vector<double> uniqueKnots_u;
				volumes[eff_patch]->basis(0).knotsSimple(uniqueKnots_u);
				for(uint j=0; j<uniqueKnots_u.size()-1; j++)
					for(int k=0; k<hrefN; k++)
						model.knot_insert(eff_patch, 0, uniqueKnots_u[j]*(k+1)/(hrefN+1) + uniqueKnots_u[j+1]*(hrefN-k)/(hrefN+1));
							
			}
			if(uniform_v) {
				vector<double> uniqueKnots_v;
				volumes[eff_patch]->basis(1).knotsSimple(uniqueKnots_v);
				for(uint j=0; j<uniqueKnots_v.size()-1; j++)
					for(int k=0; k<hrefN; k++)
						model.knot_insert(eff_patch, 1, uniqueKnots_v[j]*(k+1)/(hrefN+1) + uniqueKnots_v[j+1]*(hrefN-k)/(hrefN+1));
			}
			if(uniform_w) {
				vector<double> uniqueKnots_w;
				volumes[eff_patch]->basis(2).knotsSimple(uniqueKnots_w);
				for(uint j=0; j<uniqueKnots_w.size()-1; j++)
					for(int k=0; k<hrefN; k++)
						model.knot_insert(eff_patch, 2, uniqueKnots_w[j]*(k+1)/(hrefN+1) + uniqueKnots_w[j+1]*(hrefN-k)/(hrefN+1));
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
						for(int k=0; k<hrefN; k++)
							surfaces[i]->insertKnot_u(uniqueKnots_u[j]*(k+1)/(hrefN+1) + uniqueKnots_u[j+1]*(hrefN-k)/(hrefN+1));
					vector<double> uniqueKnots_v;
					surfaces[i]->basis_v().knotsSimple(uniqueKnots_v);
					for(uint j=0; j<uniqueKnots_v.size()-1; j++)
						for(int k=0; k<hrefN; k++)
							surfaces[i]->insertKnot_v(uniqueKnots_v[j]*(k+1)/(hrefN+1) + uniqueKnots_v[j+1]*(hrefN-k)/(hrefN+1));
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
					if(prefAmount>0 && (prefDir == -1 || prefDir == 0 ))
						surfaces[i]->raiseOrder(prefAmount,0);
					if(prefAmount>0 && (prefDir == -1 || prefDir == 1 ))
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
					for(int k=0; k<hrefN; k++)
						model.knot_insert(eff_patch, 0, uniqueKnots_u[j]*(k+1)/(hrefN+1) + uniqueKnots_u[j+1]*(hrefN-k)/(hrefN+1));
							
			}
			if(uniform_v) {
				vector<double> uniqueKnots_v;
				surfaces[eff_patch]->basis_v().knotsSimple(uniqueKnots_v);
				for(uint j=0; j<uniqueKnots_v.size()-1; j++)
					for(int k=0; k<hrefN; k++)
						model.knot_insert(eff_patch, 1, uniqueKnots_v[j]*(k+1)/(hrefN+1) + uniqueKnots_v[j+1]*(hrefN-k)/(hrefN+1));
			}
		}
	}

	cout << setprecision(17);
	model.writeSplines(cout);
	return 0;
}

