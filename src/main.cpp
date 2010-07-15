#include <iostream>
#include <fstream>
#include <string.h>

#include "TopologySet.h"
#include "primitives.h"

#include <GoTools/trivariate/SplineVolume.h>
#include <GoTools/geometry/ObjectHeader.h>

typedef unsigned int uint;

using namespace std;
using namespace Go;
using boost::shared_ptr;

bool   verbose         = false;
string fileUsage       = "\
File usage: gpm [-v] <inputFile> \n\
  \n\
  Arguments\n\
    <inputFile> : one or more .g2-files describing the spline volumes \n\
  FLAGS\n\
    -v          : verbose output";

typedef struct globNumber {
	int vertex[8];
	int edge[12];
	int edge_incr[12];
	int surface[6];
	int surface_incr_i[6];
	int surface_incr_j[6];
	int volume;
} globNumber;

bool enforceRightHandSystem(vector<shared_ptr<SplineVolume> > &volumes) ;
vector<shared_ptr<SplineVolume> > readInput(const char *filename);
void writeReparameterization(vector<shared_ptr<SplineVolume> > &volumes, const char* filename);
vector<shared_ptr<SplineVolume> > processParameters(int argc, char** argv);

vector<shared_ptr<SplineVolume> > processParameters(int argc, char** argv) {
	bool fileRead = false;
	vector<shared_ptr<SplineVolume> > results;
	for(int argi=1; argi<argc; argi++) {
		const char *arg = argv[argi];
		if(strcmp(arg, "-v") == 0) {
			verbose = true;
		} else {
			results = readInput(arg);
			fileRead = true;
		}
	}
	if(!fileRead) {
		cout << fileUsage << endl;
		exit(1);
	}
	return results;
}

bool enforceRightHandSystem(vector<shared_ptr<SplineVolume> > &volumes) {
	bool anything_switched = false;
	for(uint i=0; i<volumes.size(); i++) {
		// do test evaluation in the middle of the parametric domain
		double u = (volumes[i]->startparam(0) + volumes[i]->endparam(0)) / 2;
		double v = (volumes[i]->startparam(1) + volumes[i]->endparam(1)) / 2;
		double w = (volumes[i]->startparam(2) + volumes[i]->endparam(2)) / 2;
		vector<Point> results(4); // one position and three derivatives
		volumes[i]->point(results, u, v, w, 1);
		double jacobian = (results[1] % results[2])*results[3];
		if(jacobian < 0) {
			anything_switched = true;
			volumes[i]->reverseParameterDirection(2);
		}
	}
	return anything_switched;
}

void writeReparameterization(vector<shared_ptr<SplineVolume> > &volumes, const char* filename) {
	ofstream outFile;
	outFile.open(filename);
	for(uint i=0; i<volumes.size(); i++) {
		volumes[i]->writeStandardHeader(outFile);
		volumes[i]->write(outFile);
	}
	outFile.close();
}

vector<shared_ptr<SplineVolume> > readInput(const char *filename) {
	vector<shared_ptr<SplineVolume> > results;
	ifstream inFile;
	ObjectHeader head;
	inFile.open(filename);
	if(!inFile.good()) {
		cerr << "Error reading input file \"" << filename << "\"" << endl;
		exit(1);
	}
	while(!inFile.eof()) {
		head.read(inFile);
		if(head.classType() == Class_SplineVolume) {
			shared_ptr<SplineVolume> v(new SplineVolume());
			v->read(inFile);
			results.push_back(v);
			
			if(verbose) {
				cout << "SplineVolume succesfully read" << endl;
			}
		// } else if(head.classType == Class_SplineSurface) {  // should add support for this eventually
		} else {
			fprintf(stderr, "File \"%s\" contains unknown or unsupported class object\n", filename);
			exit(1);
		}
		ws(inFile); // eats up as many whitespaces as it can
	}
	inFile.close();
	return results;
}

int main(int argc, char **argv) {
	vector<shared_ptr<SplineVolume> > volumes = processParameters(argc, argv);
	if( enforceRightHandSystem(volumes) ) {
		cerr << "WARNING: system reparameterized to strict right-hand-system. \n";
		cerr << "         stored in \"reparameterized.g2\"\n";
		writeReparameterization(volumes, "reparameterized.g2");

	}
	TopologySet topology(volumes); // epsilon 1e-4
	topology.buildTopology();
	if(verbose) cout << "Model build complete\n";

	// debug output:
	/*
	set<Vertex*>::iterator vit;
	set<Line*>::iterator lit;
	set<Face*>::iterator fit;
	// cout << setiosflags(ios::fixed) << setprecision(3);
	cout << "All topology vertices:\n";
	for(vit=topology.vertex_begin(); vit!=topology.vertex_end(); vit++) {
		cout << " (" << (*vit)->cp << ") \t Bordering volumes: ";
		set<Volume*>::iterator volit;
		for(volit=(*vit)->volume.begin(); volit != (*vit)->volume.end(); volit++)
			cout << (*volit)->id << ", ";
		cout << endl;
	}
	for(lit=topology.line_begin(); lit!=topology.line_end(); lit++) {
		cout << "Line " << **lit << "\t Bordering volumes: ";
		set<Volume*>::iterator volit;
		// cout << "Distance at this line: " << (*lit)->v1->cp.dist((*lit)->v2->cp) << endl;
		for(volit=(*lit)->volume.begin(); volit != (*lit)->volume.end(); volit++)
			cout << (*volit)->id << ", ";
		cout << endl;
	}
	for(fit=topology.face_begin(); fit != topology.face_end(); fit++) {
		cout << "Face: " << **fit << "\t Bordering volumes: ";
		if((*fit)->v1) cout << (*fit)->v1->id << ", ";
		if((*fit)->v2) cout << (*fit)->v2->id;
		cout << endl;
	}
	*/
	if(verbose) {
		cout << "Total number of vertices: " << topology.numbVertices() << endl;
		cout << "Total number of lines   : " << topology.numbLines() << endl;
		cout << "Total number of faces   : " << topology.numbFaces() << endl;
		cout << "Total number of volumes : " << topology.numbVolumes() << endl;
	}



	/*************************************************************

	************  START COPY-PASTE FROM PREV ED. HERE *************

	**************************************************************/

	globNumber l2g[volumes.size()];
	// initalize all l2g-variables 
	for(uint i=0; i<volumes.size(); i++) {
		for(int j=0; j<8; j++)
			l2g[i].vertex[j] = -1;
		for(int j=0; j<12; j++) {
			l2g[i].edge[j] = -1;
			l2g[i].edge_incr[j] = 0;
		}
		for(int j=0; j<6; j++) {
			l2g[i].surface[j] = -1;
			l2g[i].surface_incr_i[j] = 0;
			l2g[i].surface_incr_j[j] = 0;
		}
		l2g[i].volume = -1;
	}

	set<Vertex*>::iterator   v_it;
	set<Line*>::iterator     l_it ;
	set<Face*>::iterator     f_it ;

	vector<shared_ptr<SplineVolume> >::iterator vol_it;

	// Enumerate vertices, lines, surfaces and finally volumes, in that order
	int glob_i = 0;

	vector<int>::iterator pos;
	set<Volume*>::iterator v;
	
	// for all VERTICES, assign number
	for(v_it=topology.vertex_begin(); v_it != topology.vertex_end(); v_it++) {
		for(v=(*v_it)->volume.begin(); v != (*v_it)->volume.end(); v++) {
			vector<int> corners = (*v)->getVertexEnumeration(*v_it);
			for(pos=corners.begin(); pos!=corners.end(); pos++)
				l2g[(*v)->id].vertex[*pos] = glob_i;
		}
		glob_i++;
	}

	// initialize all (possible) degenerate lines
	for(uint i=0; i<volumes.size(); i++) {
		int lineCount = 0;
		for(int parDir=0; parDir<3; parDir++) {
			for(int u2=0; u2<2; u2++) {
				for(int u1=0; u1<2; u1++) {
					int v_start, v_stop;
					if(parDir==0)      v_start =      2*u1 + 4*u2;
					else if(parDir==1) v_start = u1        + 4*u2;
					else if(parDir==2) v_start = u1 + 2*u2       ;

					if(parDir==0)      v_stop  =  1 + 2*u1 + 4*u2;
					else if(parDir==1) v_stop  = u1 +  2   + 4*u2;
					else if(parDir==2) v_stop  = u1 + 2*u2 +  4  ;

 					/* note that this test is NOT sufficient to verify degenerate lines as they may circle and stop at the same 
					 * place as they started. However if this is the case, then the values here will be overwritten at a later point.
					 */
					if(l2g[i].vertex[v_start] == l2g[i].vertex[v_stop] ) {
						// cout << "Detected line with same start- and end-point\n";
						l2g[i].edge[lineCount]      = l2g[i].vertex[v_start];
						l2g[i].edge_incr[lineCount] = 0;
					}
					lineCount++;
				}
			}
		}
	}

	// for all EDGES, assign startnumber and increment
	for(l_it=topology.line_begin(); l_it != topology.line_end(); l_it++) {
		vector<int> numb;
		vector<int> parDir;
		vector<int> parStep;
		Volume *first_vol = (*(*l_it)->volume.begin());
		shared_ptr<SplineVolume> sv = volumes[first_vol->id]; // this edge SHOULD have at least one volume-connection
		first_vol->getEdgeEnumeration(*l_it, numb, parDir, parStep);
		int coefs_here = sv->numCoefs(parDir[0])-2; // coefs_here SHOULD be identical for all lines in all volumes for this line
		if(coefs_here > 0) {
			for(v=(*l_it)->volume.begin(); v != (*l_it)->volume.end(); v++) {
				(*v)->getEdgeEnumeration(*l_it, numb, parDir, parStep); // no worries: numb, parDir & parStep all get cleaned inside function call before returned
				for(uint i=0; i<numb.size(); i++) {
					l2g[(*v)->id].edge[numb[i]]      = (parStep[i]==1) ? glob_i : glob_i+coefs_here-1;
					l2g[(*v)->id].edge_incr[numb[i]] =  parStep[i] ;
				}
			}
			glob_i += coefs_here;
		}
	}

	// initialize all (possible) degenerate surfaces
	for(uint i=0; i<volumes.size(); i++) {
		int surfCount = 0;
		int v_step   = 1; // starting vertex for face corresponding to parametric max value
		for(int parDir=0; parDir<3; parDir++) {
			if(parDir > 0) v_step *= 2; // v_step = 1,2,4
			for(int parVal=0; parVal<2; parVal++) {  // parVal = max or min
				int incr1 = (parDir==0) ? 2 : 1; 
				int incr2 = (parDir==2) ? 2 : 4;
				int v_start = parVal * v_step;

				if(l2g[i].vertex[v_start] == l2g[i].vertex[v_start+incr1]) { 
					// degenerate in local i-direction, store the line in the local j-direction
					// cout << "Detected surface degenerated in the i-direction\n";
					int line_i = Line::getLineEnumeration(v_start, v_start+incr2);
					l2g[i].surface[surfCount]        = l2g[i].edge[line_i];
					l2g[i].surface_incr_j[surfCount] = l2g[i].edge_incr[line_i];
				}
				if(l2g[i].vertex[v_start] == l2g[i].vertex[v_start+incr2]) {
					// vica versa
					//cout << "Detected surface degenerated in the j-direction\n";
					int line_i = Line::getLineEnumeration(v_start, v_start+incr1);
					l2g[i].surface[surfCount]        = l2g[i].edge[line_i];
					l2g[i].surface_incr_i[surfCount] = l2g[i].edge_incr[line_i];
				}

				surfCount++;
			}
		}
	}

	// for all SURFACES, assign startnumber and increments
	for(f_it=topology.face_begin(); f_it != topology.face_end(); f_it++) {
		/*  Face Index (logic behind numCoef_u/v)
		 *  0-1 : surface umin/umax - numcCefs(1) X numCoefs(2)
		 *  2-3 : surface vmin/vmax - numcCefs(0) X numCoefs(2)
		 *  4-5 : surface wmin/wmax - numcCefs(0) X numCoefs(1)
		 */
		if((*f_it)->v2) { // inner face w/ two adjacent volumes
			int id1   = (*f_it)->v1->id;
			int face1 = (*f_it)->face1;
			int id2   = (*f_it)->v2->id;
			int face2 = (*f_it)->face2;
			shared_ptr<SplineVolume> s_volume = volumes[id1];
			int numCoef_u = s_volume->numCoefs(  (face1 < 2) ) - 2;
			int numCoef_v = s_volume->numCoefs(2-(face1 > 3) ) - 2;
			int coefs_here  = numCoef_u * numCoef_v;
			if(coefs_here > 0) {
				l2g[id1].surface[face1]        = glob_i;
				l2g[id1].surface_incr_i[face1] = 1;
				l2g[id1].surface_incr_j[face1] = numCoef_u;

				bool uv_flip   = (*f_it)->uv_flip;
				bool u_reverse = (*f_it)->u_reverse;
				bool v_reverse = (*f_it)->v_reverse;

				if(uv_flip) {
					l2g[id2].surface_incr_i[face2] = numCoef_u;
					l2g[id2].surface_incr_j[face2] = 1;
				} else {
					l2g[id2].surface_incr_i[face2] = 1;
					l2g[id2].surface_incr_j[face2] = numCoef_u;
				}
				if(u_reverse)
					l2g[id2].surface_incr_i[face2] *= -1;
				if(v_reverse)
					l2g[id2].surface_incr_j[face2] *= -1;

				if(!u_reverse && !v_reverse)
					l2g[id2].surface[face2] = glob_i ;
				else if(u_reverse && v_reverse)
					l2g[id2].surface[face2] = glob_i + coefs_here - 1;
				else if(v_reverse == uv_flip)
					l2g[id2].surface[face2] = glob_i + numCoef_u - 1;
				else // u_reverse == uv_flip
					l2g[id2].surface[face2] = glob_i + coefs_here - numCoef_u;

				glob_i += coefs_here;
			}
		} else { // boundary face
			int id    = (*f_it)->v1->id;
			int face1 = (*f_it)->face1;
			shared_ptr<SplineVolume> s_volume = volumes[id];
			int numCoef_u = s_volume->numCoefs(  (face1 < 2) ) - 2;
			int numCoef_v = s_volume->numCoefs(2-(face1 > 3) ) - 2;
			int coefs_here  = numCoef_u * numCoef_v;
			if(coefs_here > 0) {
				l2g[id].surface[face1]        = glob_i;
				l2g[id].surface_incr_i[face1] = 1;
				l2g[id].surface_incr_j[face1] = numCoef_u;
				glob_i += coefs_here;
			}
		}
	}

	// for all VOLUMES assign startnumber
	for(uint i=0; i<volumes.size(); i++) {
		shared_ptr<SplineVolume> s_volume = volumes[i];
		int numCoef_u = s_volume->numCoefs(0)-2;
		int numCoef_v = s_volume->numCoefs(1)-2;
		int numCoef_w = s_volume->numCoefs(2)-2;
		int coefs_here  = numCoef_u * numCoef_v * numCoef_w;
		if(coefs_here > 0) {
			l2g[i].volume = glob_i;
			glob_i += coefs_here;
		}
	}

	if(verbose) { // human friendly output
		for(uint i=0; i<volumes.size(); i++) {
			cout << "IBLOCK " << i << endl;
			cout << "  Global vertex numbers:\n    ";
			for(int j=0; j<8; j++)
				cout << l2g[i].vertex[j] << " ";
			cout << endl;
			cout << "  Global line numbers and iterator:\n";
			for(int j=0; j<12; j++)
				cout << "    " << l2g[i].edge[j] << " " << l2g[i].edge_incr[j] << endl;
			cout << "  Global face numbers and iterators:\n";
			for(int j=0; j<6; j++)
				cout << "    " << l2g[i].surface[j] << " " << l2g[i].surface_incr_i[j] << " " << l2g[i].surface_incr_j[j] << endl;
			cout << "  Internal volume start:\n";
			cout << "    " << l2g[i].volume << endl;
		}
	} else { // production output
		for(uint i=0; i<volumes.size(); i++) {
			cout << i << endl;
			for(int j=0; j<8; j++)
				cout << l2g[i].vertex[j] << " ";
			cout << endl;
			for(int j=0; j<12; j++)
				cout << l2g[i].edge[j] << " " << l2g[i].edge_incr[j] << endl;
			for(int j=0; j<6; j++)
				cout << l2g[i].surface[j] << " " << l2g[i].surface_incr_i[j] << " " << l2g[i].surface_incr_j[j] << endl;
			cout << l2g[i].volume << endl;
		}
	}

	return 0;
}
