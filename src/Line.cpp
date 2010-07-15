#include "primitives.h"
#include <vector>
#include <GoTools/utils/Point.h>

using namespace std;
using namespace Go;

bool Line::equals(Line *l, double tol) {
	vector<Point>::iterator cp_new  = l->cp.begin();
	if(cp_new->dist( cp.front() ) < tol ) {
		vector<Point>::iterator this_cp ;
		for(this_cp=cp.begin(); this_cp!=cp.end() && cp_new!=l->cp.end(); this_cp++, cp_new++)
			if( this_cp->dist(*cp_new) > tol )
				return false;
		return true;
	} else if(cp_new->dist( cp.back() ) < tol ) {
		vector<Point>::reverse_iterator this_cp ;
		for(this_cp=cp.rbegin(); this_cp!=cp.rend() && cp_new!=l->cp.end(); this_cp++, cp_new++)
			if( this_cp->dist(*cp_new) > tol )
				return false;
		return true;
	} else {
		return false;
	}
}

void Line::write(ostream &os) const {
	os << "(" << v1->cp << ") -> (" << v2->cp << ")";
}

void Line::read(std::istream &is) {
	// do nothing
}

int Line::getLineEnumeration(int vert1, int vert2) {
	if(vert2<vert1) {
		int tmp = vert2;
		vert2   = vert1;
		vert1   = tmp;
	}
	if(vert2-vert1 == 1) // u-line
		return vert1/2;
	else if(vert2-vert1==2) // v-line
		return (vert2+vert1)/4 + 4;
	else if(vert2-vert1==4) // w-line
		return vert1 + 8;
	else // no valid vertex combination
		return -1;
}
