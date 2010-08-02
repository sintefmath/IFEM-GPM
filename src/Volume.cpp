/**********************************************************************************//**
 * \file Volume.cpp
 *
 * \author Kjetil A. Johannessen
 *
 * \date July 2010
 *
 *************************************************************************************/

#include "primitives.h"
#include <vector>

using namespace std;

/**********************************************************************************//**
 * Constructor for the Volume class
 *
 * \param id Unique identification number
 *************************************************************************************/
Volume::Volume(int id) {
	this->id      = id;
	material_code = 0;
}

/**********************************************************************************//**
 * \brief get local enumeration of a vertex
 *
 * \param v Vertex to be tested
 * \return A vector containing all corner positions where v appears
 *************************************************************************************/
vector<int> Volume::getVertexEnumeration(Vertex *v) {
	vector<int> results;
	for(int i=0; i<8; i++)
		if(v == corner[i])
			results.push_back(i);
	return results;
}

/**********************************************************************************//**
 * \brief get local enumeration of a Line
 *
 * \param [in]  l       line to be tested
 * \param [out] numb    local line number
 * \param [out] parDir  0,1 or 2, giving the parametric direction over which the line is running
 * \param [out] parStep +1 or -1 depending on whether the line is oriented the same way as l
 * 
 * This function searches for all occurrences of *l among *this̈́' edge lines. For each occurrence,
 * the relation between *l and the detected line will be stored and returned as three vectors.
 *************************************************************************************/
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

/**********************************************************************************//**
 * \brief get local enumeration of a Face
 *
 * \param f face to be tested
 * \return all occurrences of the face *f among *this' stored faces, returned as their local
 *         face enumeration
 *************************************************************************************/
vector<int> Volume::getSurfaceEnumeration(Face *f) {
	vector<int> results;
	for(uint i=0; i<f->volume.size(); i++)
		if(f->volume[i] == this)
			results.push_back(f->face[i]);
	return results;
}

