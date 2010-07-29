#include "primitives.h"

using namespace std;


std::vector<int> Vertex::getVertexEnumeration(int face) {
	vector<int> results;
	results.push_back( (face%2) * (1 << (face/2) ) );
	results.push_back(results[0] + (face<2) + 1);
	results.push_back(results[1] + 2*(face<4) - (face<2) + 1);
	results.push_back(results[2] + (face<2) + 1);
	return results;
}

