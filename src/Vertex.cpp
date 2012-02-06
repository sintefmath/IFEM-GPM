#include "primitives.h"

using namespace std;

void Vertex::getVertexEnumerationOnFace(int line, int &v1, int &v2) {
	if(line == 0) {
		v1 = 0;
		v2 = 2;
	} else if(line == 1) {
		v1 = 1;
		v2 = 3;
	} else if(line == 2) {
		v1 = 0;
		v2 = 1;
	} else if(line == 3) {
		v1 = 2;
		v2 = 3;
	}
}

void Vertex::getVertexEnumerationOnVolume(int line, int &v1, int &v2) {
	int parDir = line/4;
	int n      = line%4;
	if(parDir==0)      v1 = 2*n;
	else if(parDir==1) v1 = n + 2*(n>1);
	else if(parDir==2) v1 = n;

	if(parDir==0)      v2 = v1+1;
	else if(parDir==1) v2 = v1+2;
	else if(parDir==2) v2 = v1+4;
}

std::vector<int> Vertex::getVertexEnumeration(int face) {
	vector<int> results;
	results.push_back( (face%2) * (1 << (face/2) ) );
	results.push_back(results[0] + (face<2) + 1);
	results.push_back(results[1] + 2*(face<4) - (face<2) + 1);
	results.push_back(results[2] + (face<2) + 1);
	return results;
}

