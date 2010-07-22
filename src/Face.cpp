/**********************************************************************************//**
 * \file Face.cpp
 *
 * \author Kjetil A. Johannessen
 *
 * \date July 2010
 *
 *************************************************************************************/
#include "primitives.h"
#include <iostream>

using namespace std;

/**********************************************************************************//**
 * \brief Basic constructor for the Face class
 * \param n1 Number of internal controlpoints in the first parametric direction
 * \param n2 Number of internal controlpoints in the second parametric direction
 *************************************************************************************/
Face::Face(int n1, int n2) {
	cp.resize(n1);
	for(int u=0; u<n1; u++)
		cp[u].resize(n2);

	v1    = NULL;
	v2    = NULL;
	face1 = -1;
	face2 = -1;
	uv_flip    = false;
	u_reverse  = false;
	v_reverse  = false;
	bc_code    = 0;
}

/**********************************************************************************//**
 * \brief Test face equality
 * \param f face to be compared with this
 * \param tol Equality-tolerance for the euclidean distance between the control points
 * \return True if control points are equal (up to tolerance)
 * 
 * In the case of equality being detected, the function will set the uv_flip, u_reverse
 * and v_reverse variables to their correct state. These are being set with respect to 
 * *this' relation to *f.
 *************************************************************************************/
bool Face::equals(Face *f, double tol) {
	int n1  = cp.size();
	int n2  = cp[0].size();
	int fn1 = f->cp.size();
	int fn2 = f->cp[0].size();
	for(int uv_flip=0; uv_flip<2; uv_flip++) {
		if( !( ( uv_flip && n1==fn2 && n2==fn1) ||
	       (!uv_flip && n1==fn1 && n2==fn2) ) )
			continue;
		for(int u_reverse=0; u_reverse<2; u_reverse++) {
			for(int v_reverse=0; v_reverse<2; v_reverse++) {
				int fu_start, fv_start, du, dv;
				if((u_reverse && !uv_flip) ||
				   (v_reverse && uv_flip) ) {
					fu_start = n1-1;
					du = -1;
				} else {
					fu_start = 0;
					du = 1;
				}
				if((v_reverse && !uv_flip) ||
				   (u_reverse && uv_flip) ) {
					fv_start = n2-1;
					dv = -1;
				} else {
					fv_start = 0;
					dv = 1;
				}
				bool match = true;
				int fv = fv_start;
				for(int v=0; v<n2 && match; v++, fv+=dv) {
					int fu = fu_start;
					for(int u=0; u<n1; u++, fu+=du) {
						if( ( uv_flip && cp[u][v].dist(f->cp[fv][fu]) > tol) ||
						    (!uv_flip && cp[u][v].dist(f->cp[fu][fv]) > tol) ) {
							match = false;
							break;
						}
					}
				}
				if(match) {
					this->uv_flip   = uv_flip;
					this->u_reverse = u_reverse;
					this->v_reverse = v_reverse;
					return true;
				}
			}
		}
	}
	return false;
}

/*! \brief Is face degenerated to point or line */
bool Face::isDegen() {
	return degen1 || degen2;
}


void Face::write(std::ostream &os) const {
	// do nothing
}

void Face::read(std::istream &is) {
	// do nothing
}
