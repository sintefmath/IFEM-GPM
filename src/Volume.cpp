
#include "primitives.h"
#include <vector>

using namespace std;

Volume::Volume(int id) {
	this->id = id;
}

vector<int> Volume::getVertexEnumeration(Vertex *v) {
	vector<int> results;
	for(int i=0; i<8; i++)
		if(v == corner[i])
			results.push_back(i);
	return results;
}

void Volume::getEdgeEnumeration(Line *l, vector<int> &numb, vector<int> &parDir, vector<int> &parStep) {
	numb.clear();
	parDir.clear();
	parStep.clear();
	int lineCount = 0;
	for(int p=0; p<3; p++) { // parametric direction
		for(int u2=0; u2<2; u2++) {
			for(int u1=0; u1<2; u1++) {
				if(l == line[lineCount]) {
					int v_start = -1; // index of line starting vertex
					if(p==0)      v_start  =      2*u1 + 4*u2;
					else if(p==1) v_start  = u1        + 4*u2;
					else if(p==2) v_start  = u1 + 2*u2       ;

					// is the line oriented the same way globaly as localy?
					if(l->v1 == corner[v_start])
						parStep.push_back(1);
					else
						parStep.push_back(-1);
					numb.push_back(lineCount);
					parDir.push_back(p);
				}
				lineCount++;
			}
		}
	}
}

vector<int> Volume::getSurfaceEnumeration(Face *f) {
	vector<int> results;
	if(f->v1 == this)
		results.push_back(f->face1);
	if(f->v2 != NULL && f->v2 == this)
		results.push_back(f->face2);
	return results;
}

