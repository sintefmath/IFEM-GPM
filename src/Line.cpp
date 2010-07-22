/**********************************************************************************//**
 * \file Line.cpp
 *
 * \author Kjetil A. Johannessen
 *
 * \date July 2010
 *
 *************************************************************************************/

#include "primitives.h"
#include <vector>
#include <GoTools/utils/Point.h>

using namespace std;
using namespace Go;

Line::Line() {
	v1      = NULL;
	v2      = NULL;
	bc_code = 0;
}

/**********************************************************************************//**
 * \brief Checks for line equality by comparing all controlpoints
 * \param l The line which is to be compared
 * \param tol The tolerance given to the euclidean  distance between the controlpoints
 * \return True if all control points are equal (within tolerance)
 *
 * This function will try both orientations of the control points since a line primitive may contain the 
 * start at either v1 or v2
 *************************************************************************************/
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

/**********************************************************************************//**
 * \brief User friendly output format
 * \param os output stream to be written to
 *************************************************************************************/
void Line::write(ostream &os) const {
	os << "(" << v1->cp << ") -> (" << v2->cp << ")";
}

void Line::read(std::istream &is) {
	// do nothing
}

/**********************************************************************************//**
 * \brief get the enumeration of the four edge lines corresponding to one face
 * \param face local face enumeration
 * \return vector containing the number to the four edges
 *************************************************************************************/
std::vector<int> Line::getLineEnumeration(int face) {
	vector<int> results;

	results.push_back( (face<2)*4 + face%2 + (face>4) );
	results.push_back( 1+(face<4)+4*(face<2) + face%2 + (face>4) );
	results.push_back( 4+(face<4)*4 + face%2 + (face>2)*(face%2) );
	results.push_back( 5+(face<4)*4+(face<2) + face%2 + (face>2)*(face%2) );
	return results;
}

/**********************************************************************************//**
 * \brief get the line enumeration corresponding to two corner numbers
 * \param vert1 First corner
 * \param vert2 Second corner
 * \return Local line enumeration as defined in Volume
 *************************************************************************************/
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
