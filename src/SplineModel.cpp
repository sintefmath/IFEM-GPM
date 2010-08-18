/**********************************************************************************//**
 * \file TopologySet.cpp
 *
 * \author Kjetil A. Johannessen
 *
 * \date July 2010
 *
 *************************************************************************************/
#include "SplineModel.h"
#include "TopologySet.h"
#include "primitives.h"

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include <GoTools/geometry/ObjectHeader.h>

using namespace std;
using namespace Go;
using boost::shared_ptr;

SplineModel::SplineModel() {
	l2g = new globNumber[0];
	topology = new TopologySet();
}

/**********************************************************************************//**
 * \brief Constructor
 * \param spline_volumes All spline volumes to be considered part of this model
 *************************************************************************************/
SplineModel::SplineModel(std::vector<boost::shared_ptr<Go::SplineVolume> > &spline_volumes) {
	l2g = new globNumber[0];
	topology = new TopologySet(spline_volumes);
	topology->buildTopology();
	spline_volumes_ = spline_volumes;
}

//! \brief Destructor
SplineModel::~SplineModel() {
	if(topology)
		delete topology;
	delete[] l2g;
}

TopologySet* SplineModel::getTopology() {
	return topology;
}

std::vector<boost::shared_ptr<Go::SplineVolume> >& SplineModel::getSplines() {
	return spline_volumes_;
}

/**********************************************************************************//**
 * \brief Make the parametrization of all spline-volumes to be right-hand-system
 * \returns true if any volumes was reparameterized
 *
 * This ensures that all jacobians are evaluated positive, assuming of course that you
 * have a model descirption which is not wrapping itself inside out at any points (i.e.
 * jacboian must be non-negative or non-positive in the entire domain prior to
 * reparametrization).
 *************************************************************************************/
bool SplineModel::enforceRightHandSystem() {
	bool anything_switched = false;
	for(uint i=0; i<spline_volumes_.size(); i++) {
		// do test evaluation in the middle of the parametric domain
		double u = (spline_volumes_[i]->startparam(0) + spline_volumes_[i]->endparam(0)) / 2;
		double v = (spline_volumes_[i]->startparam(1) + spline_volumes_[i]->endparam(1)) / 2;
		double w = (spline_volumes_[i]->startparam(2) + spline_volumes_[i]->endparam(2)) / 2;
		vector<Point> results(4); // one position and three derivatives
		spline_volumes_[i]->point(results, u, v, w, 1);
		double jacobian = (results[1] % results[2])*results[3];
		if(jacobian < 0) {
			anything_switched = true;
			spline_volumes_[i]->reverseParameterDirection(2);
		}
	}
	if(anything_switched)
		topology->buildTopology();
	return anything_switched;
}

/**********************************************************************************//**
 * \brief sets a (material) property code to a volume
 * \param volId     volume identifier (natural numbering if read from a .g2-file)
 * \param propCode  property code
 * \param inclusive Whether the corners, lines and faces associated with the volume
 *                  should get the same code
 * \returns true if any corners, lines or faces had previously been tagged with a
 *          different code and has now been overwritten (only applicable if
 *          inclusive=true)
 *************************************************************************************/
bool SplineModel::addVolumePropertyCode(int volId, int propCode, bool inclusive) {
	Volume* v = topology->getVolume(volId);
	v->material_code = propCode;
	if(inclusive) {
		bool changes = false;
		for(int i=0; i<6; i++) {
			Face *f = v->face[i];
			changes = changes || (f->bc_code!=propCode);
			f->bc_code = propCode;
		}
		for(int i=0; i<12; i++) {
			Line *l = v->line[i];
			changes = changes || (l->bc_code!=propCode);
			l->bc_code = propCode;
		}
		for(int i=0; i<8; i++) {
			Vertex *c = v->corner[i];
			changes = changes || (c->bc_code!=propCode);
			c->bc_code = propCode;
		}
		return changes;
	}
	return false;
}
/**********************************************************************************//**
 * \brief sets a (boundary condition) property code to a face
 * \param volId     volume identifier (natural numbering if read from a .g2-file)
 * \param faceId    local face number corresponding to the volume volId, as defined in
 *                  the class Volume (values 0-5)
 * \param propCode  property code
 * \param inclusive Whether the corners or lines associated with the face should get 
 *                  the same code
 * \returns true if any corners or lines had previously been tagged with a different 
 *          code and has now been overwritten (only applicable if inclusive=true)
 *************************************************************************************/
bool SplineModel::addFacePropertyCode(int volId, int faceId, int propCode, bool inclusive) {
	Volume* v = topology->getVolume(volId);
	Face*   f = v->face[faceId];
	f->bc_code = propCode;
	if(inclusive) {
		bool changes = false;
		
		// tracking edge lines & points through the first volume in *f - not neccessary equal to
		// *v (but the lines & points should be equal)
		vector<int> lineNumb = Line::getLineEnumeration( f->face[0] );
		for(int i=0; i<4; i++) {
			Line* l = f->volume[0]->line[lineNumb[i]];
			changes = changes || (l->bc_code!=propCode);
			l->bc_code = propCode;
			changes = changes || (l->v1->bc_code!=propCode);
			l->v1->bc_code = propCode;
			changes = changes || (l->v2->bc_code!=propCode);
			l->v2->bc_code = propCode;
		}
		return changes;
	}
	return false;
}
/**********************************************************************************//**
 * \brief sets a (boundary condition) property code to a line
 * \param volId     volume identifier (natural numbering if read from a .g2-file)
 * \param lineId    local line number corresponding to the volume volId, as defined in
 *                  the class Volume (values 0-11)
 * \param propCode  property code
 * \param inclusive Whether the corners associated with the face should get the same code
 * \returns true if any corners had previously been tagged with a different code and has 
 *          now been overwritten (only applicable if inclusive=true)
 *************************************************************************************/
bool SplineModel::addLinePropertyCode(int volId, int lineId, int propCode, bool inclusive) {
	Volume* v = topology->getVolume(volId);
	Line*   l = v->line[lineId];
	l->bc_code = propCode;
	if(inclusive) {
		bool changes = false;
		l->bc_code = propCode;
		changes = changes || (l->v1->bc_code!=propCode);
		l->v1->bc_code = propCode;
		changes = changes || (l->v2->bc_code!=propCode);
		l->v2->bc_code = propCode;
		return changes;
	}
	return false;
}
/**********************************************************************************//**
 * \brief sets a (boundary condition) property code to a vertex
 * \param volId     volume identifier (natural numbering if read from a .g2-file)
 * \param vertId    local vertex number corresponding to the volume volId, as defined in
 *                  the class Volume (values 0-7)
 * \param propCode  property code
 *************************************************************************************/
void SplineModel::addVertexPropertyCode(int volId, int vertId, int propCode) {
	Volume* v = topology->getVolume(volId);
	v->corner[vertId]->bc_code = propCode;
}

int SplineModel::getVolumePropertyCode(int volId) {
	Volume* v = topology->getVolume(volId);
	return v->material_code;
}

int SplineModel::getFacePropertyCode(int volId, int faceId) {
	Volume* v = topology->getVolume(volId);
	return v->face[faceId]->bc_code;
}

int SplineModel::getLinePropertyCode(int volId, int lineId) {
	Volume* v = topology->getVolume(volId);
	return v->line[lineId]->bc_code;
}

int SplineModel::getVertexPropertyCode(int volId, int vertId) {
	Volume* v = topology->getVolume(volId);
	return v->corner[vertId]->bc_code;
}

//! \brief generate the local to global enumerations
void SplineModel::generateGlobalNumbers() {

	delete[] l2g;
	l2g = new globNumber[spline_volumes_.size()];
	// initalize all l2g-variables 
	for(uint i=0; i<spline_volumes_.size(); i++) {
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

	set<Vertex*>::iterator v_it;
	set<Line*>::iterator   l_it ;
	set<Face*>::iterator   f_it ;

	vector<shared_ptr<SplineVolume> >::iterator vol_it;

	// Enumerate vertices, lines, surfaces and finally volumes, in that order
	int glob_i = 0;

	vector<int>::iterator pos;
	set<Volume*>::iterator v;
	
	// for all VERTICES, assign number
	for(v_it=topology->vertex_begin(); v_it != topology->vertex_end(); v_it++) {
		for(v=(*v_it)->volume.begin(); v != (*v_it)->volume.end(); v++) {
			vector<int> corners = (*v)->getVertexEnumeration(*v_it);
			for(pos=corners.begin(); pos!=corners.end(); pos++)
				l2g[(*v)->id].vertex[*pos] = glob_i;
		}
		glob_i++;
	}

	// for all EDGES, assign startnumber and increment
	for(l_it=topology->line_begin(); l_it != topology->line_end(); l_it++) {
		vector<int> numb;
		vector<int> parDir;
		vector<int> parStep;
		Volume *first_vol = (*(*l_it)->volume.begin());
		shared_ptr<SplineVolume> sv = spline_volumes_[first_vol->id]; // this edge SHOULD have at least one volume-connection
		first_vol->getEdgeEnumeration(*l_it, numb, parDir, parStep);
		int coefs_here = sv->numCoefs(parDir[0])-2; // coefs_here SHOULD be identical for all lines in all volumes for this line
		if(coefs_here > 0) {
			for(v=(*l_it)->volume.begin(); v != (*l_it)->volume.end(); v++) {
				(*v)->getEdgeEnumeration(*l_it, numb, parDir, parStep); // no worries: numb, parDir & parStep all get cleaned inside function call before returned
				for(uint i=0; i<numb.size(); i++) {
					if((*l_it)->degen) {
						vector<int> corners = (*v)->getVertexEnumeration((*l_it)->v1);
						int corner_pos = corners[0]; // should be the same global number on all elements of "corners"
						l2g[(*v)->id].edge[numb[i]]      = l2g[(*v)->id].vertex[corner_pos];
						l2g[(*v)->id].edge_incr[numb[i]] = 0;
					} else {
						l2g[(*v)->id].edge[numb[i]]      = (parStep[i]==1) ? glob_i : glob_i+coefs_here-1;
						l2g[(*v)->id].edge_incr[numb[i]] =  parStep[i] ;
					}
				}
			}
			if( !(*l_it)->degen )
				glob_i += coefs_here;
		}
	}

	// for all SURFACES, assign startnumber and increments
	for(f_it=topology->face_begin(); f_it != topology->face_end(); f_it++) {
		Face *f = *f_it;
		/*  Face Index (logic behind numCoef_u/v)
		 *  0-1 : surface umin/umax - numcCefs(1) X numCoefs(2)
		 *  2-3 : surface vmin/vmax - numcCefs(0) X numCoefs(2)
		 *  4-5 : surface wmin/wmax - numcCefs(0) X numCoefs(1)
		 */
		shared_ptr<SplineVolume> s_volume = spline_volumes_[f->volume[0]->id];
		int numCoef_u  = s_volume->numCoefs(  (f->face[0] < 2) ) - 2;
		int numCoef_v  = s_volume->numCoefs(2-(f->face[0] > 3) ) - 2;
		int coefs_here = numCoef_u*numCoef_v;

		for(uint i=0; i<f->volume.size(); i++) {
			Volume *vol = f->volume[i];
			int faceId  = f->face[i];
			bool uv_flip   = f->uv_flip[i];
			bool u_reverse = f->u_reverse[i];
			bool v_reverse = f->v_reverse[i];


			if(f->degen1 && f->degen2) {
				vector<int> corners = Vertex::getVertexEnumeration(faceId);
				l2g[vol->id].surface[faceId]        = l2g[vol->id].vertex[corners[0]];
				l2g[vol->id].surface_incr_i[faceId] = 0;
				l2g[vol->id].surface_incr_j[faceId] = 0;
			} else if(f->degen1) {
				vector<int> lines = Line::getLineEnumeration(faceId);
				if(uv_flip) {
					l2g[vol->id].surface[faceId]        = l2g[vol->id].edge[lines[0]];
					l2g[vol->id].surface_incr_i[faceId] = l2g[vol->id].edge_incr[lines[0]];
					l2g[vol->id].surface_incr_j[faceId] = 0;
				} else {
					l2g[vol->id].surface[faceId]        = l2g[vol->id].edge[lines[2]];
					l2g[vol->id].surface_incr_i[faceId] = 0;
					l2g[vol->id].surface_incr_j[faceId] = l2g[vol->id].edge_incr[lines[2]];
				}
			} else if(f->degen2) {
				vector<int> lines = Line::getLineEnumeration(faceId);
				if(uv_flip) {
					l2g[vol->id].surface[faceId]        = l2g[vol->id].edge[lines[2]];
					l2g[vol->id].surface_incr_i[faceId] = 0;
					l2g[vol->id].surface_incr_j[faceId] = l2g[vol->id].edge_incr[lines[2]];
				} else {
					l2g[vol->id].surface[faceId]        = l2g[vol->id].edge[lines[0]];
					l2g[vol->id].surface_incr_i[faceId] = l2g[vol->id].edge_incr[lines[0]];
					l2g[vol->id].surface_incr_j[faceId] = 0;
				}
			} else {
				if(!u_reverse && !v_reverse)
					l2g[vol->id].surface[faceId] = glob_i ;
				else if(u_reverse && v_reverse)
					l2g[vol->id].surface[faceId] = glob_i + coefs_here - 1;
				else if(v_reverse == uv_flip)
					l2g[vol->id].surface[faceId] = glob_i + numCoef_u - 1;
				else // u_reverse == uv_flip
					l2g[vol->id].surface[faceId] = glob_i + coefs_here - numCoef_u;

				if(uv_flip) {
					l2g[vol->id].surface_incr_i[faceId] = numCoef_u;
					l2g[vol->id].surface_incr_j[faceId] = 1;
				} else {
					l2g[vol->id].surface_incr_i[faceId] = 1;
					l2g[vol->id].surface_incr_j[faceId] = numCoef_u;
				}
				if(u_reverse)
					l2g[vol->id].surface_incr_i[faceId] *= -1;
				if(v_reverse)
					l2g[vol->id].surface_incr_j[faceId] *= -1;
			}
		}
		if(!f->isDegen()) 
			glob_i += coefs_here;
	}

	// for all VOLUMES assign startnumber
	for(uint i=0; i<spline_volumes_.size(); i++) {
		shared_ptr<SplineVolume> s_volume = spline_volumes_[i];
		int numCoef_u = s_volume->numCoefs(0)-2;
		int numCoef_v = s_volume->numCoefs(1)-2;
		int numCoef_w = s_volume->numCoefs(2)-2;
		int coefs_here  = numCoef_u * numCoef_v * numCoef_w;
		if(coefs_here > 0) {
			l2g[i].volume = glob_i;
			glob_i += coefs_here;
		}
	}
}

void SplineModel::writeGlobalNumberOrdering(std::ostream &os) const {
	for(uint i=0; i<spline_volumes_.size(); i++) {
		os << i << endl;
		for(int j=0; j<8; j++)
			os << l2g[i].vertex[j] << " ";
		os << endl;
		for(int j=0; j<12; j++)
			os << l2g[i].edge[j] << " " << l2g[i].edge_incr[j] << endl;
		for(int j=0; j<6; j++)
			os << l2g[i].surface[j] << " " << l2g[i].surface_incr_i[j] << " " << l2g[i].surface_incr_j[j] << endl;
		os << l2g[i].volume << endl;
	}
}

void SplineModel::writeModelProperties(std::ostream &os) const {
	set<Volume*>::iterator v_it;
	set<Face*>::iterator   f_it;
	set<Line*>::iterator   l_it;
	set<Vertex*>::iterator c_it; //corner iterator

	/***************     WRITE IT IN THE INPUT-FORM, NOT REALLY USEFULL FOR PRODUCTION    ***********************/

	for(v_it=topology->volume_begin(); v_it != topology->volume_end(); v_it++) {
		Volume* v = *v_it;
		if(v->material_code != 0)
			// os << "Volume " << v->id << " " << v->material_code << endl;
			os << v->material_code << " " << v->id << " 3 " << endl;
	}
	for(f_it=topology->face_begin(); f_it != topology->face_end(); f_it++) {
		Face *f = *f_it;
		if(f->bc_code != 0) 
			// os << "Face " << f->volume[0]->id << " " << f->face[0] << " " << f->bc_code << " 0" << endl;
			os << f->bc_code << " " << f->volume[0]->id << " 2 " << f->face[0] << endl;
	}
	for(l_it=topology->line_begin(); l_it != topology->line_end(); l_it++) {
		Line *l = *l_it;
		vector<int> numb, parDir, parStep;
		(*l->volume.begin())->getEdgeEnumeration(l, numb, parDir, parStep);
		if(l->bc_code != 0)
			// os << "Line " << (*l->volume.begin())->id << " " << numb[0] << " " << l->bc_code << " 0" << endl;
			os << l->bc_code << " " << (*l->volume.begin())->id << " 1 " << numb[0] << endl;
	}
	for(c_it=topology->vertex_begin(); c_it != topology->vertex_end(); c_it++) {
		Vertex *c = *c_it;
		vector<int> corner_id = (*c->volume.begin())->getVertexEnumeration(c);
		if(c->bc_code != 0)
			// os << "Vertex " << (*c->volume.begin())->id << " " << corner_id[0] << " " << c->bc_code << endl;
			os << c->bc_code << " " << (*c->volume.begin())->id << " 0 " << corner_id[0] << endl;
	}
	/************************************************************************************************************/

/*
	globNumber propCodes[spline_volumes_.size()];
	for(v_it=topology->volume_begin(); v_it != topology->volume_end(); v_it++)
		propCodes[(*v_it)->id].volume = (*v_it)->material_code;

	for(f_it=topology->face_begin(); f_it != topology->face_end(); f_it++) {
		for(uint i=0; i<(*f_it)->volume.size(); i++) {
			Volume *vol = (*f_it)->volume[i];
			vector<int> fIds = vol->getSurfaceEnumeration(*f_it);
			for(uint j=0; j<fIds.size(); j++)
				propCodes[vol->id].surface[fIds[j]] = (*f_it)->bc_code;
		}
	}
	
	for(l_it=topology->line_begin(); l_it != topology->line_end(); l_it++) {
		for(v_it=(*l_it)->volume.begin(); v_it != (*l_it)->volume.end(); v_it++) {
			vector<int> numb, parDir, parStep;
			(*v_it)->getEdgeEnumeration(*l_it, numb, parDir, parStep);
			for(uint i=0; i<numb.size(); i++)
				propCodes[(*v_it)->id].edge[numb[i]] = (*l_it)->bc_code;
		}
	}

	for(c_it=topology->vertex_begin(); c_it != topology->vertex_end(); c_it++) {
		for(v_it=(*c_it)->volume.begin(); v_it != (*c_it)->volume.end(); v_it++) {
			vector<int> cIds = (*v_it)->getVertexEnumeration(*c_it);
			for(uint i=0; i<cIds.size(); i++)
				propCodes[(*v_it)->id].vertex[cIds[i]] = (*c_it)->bc_code;
		}
	}

	for(uint i=0; i<spline_volumes_.size(); i++) {
		os << propCodes[i].volume << endl;
		for(uint j=0; j<8; j++)
			os << propCodes[i].vertex[j] << " ";
		os << endl;
		for(uint j=0; j<12; j++)
			os << propCodes[i].edge[j] << " ";
		os << endl;
		for(uint j=0; j<6; j++)
			os << propCodes[i].surface[j] << " ";
		os << endl;
	}
*/


}

void SplineModel::writeSplines(std::ostream &os) const {
	for(uint i=0; i<spline_volumes_.size(); i++) {
		spline_volumes_[i]->writeStandardHeader(os);
		spline_volumes_[i]->write(os);
	}
}

void SplineModel::readSplines(std::istream &is) {
	ObjectHeader head;
	while(!is.eof()) {
		head.read(is);
		if(head.classType() == Class_SplineVolume) {
			shared_ptr<SplineVolume> v(new SplineVolume());
			v->read(is);
			spline_volumes_.push_back(v);
			topology->addVolume(v);
		// } else if(head.classType == Class_SplineSurface) {  // should add support for this eventually
		} else {
			fprintf(stderr, "unknown or unsupported class object\n");
			exit(1);
		}
		ws(is); // eats up as many whitespaces as it can
	}
	topology->buildTopology();
}

void SplineModel::readModelProperties(std::istream &is) {
	if(!is.good())
		// some error code 
		;
	string str_line, primitive;
	int volId, locId, code, inclusive;
	while( getline(is, str_line) ) {
		istringstream ss(str_line);
		ss >> primitive;
		if(primitive.compare("Volume") == 0 ||
		   primitive.compare("volume") == 0 )  {
			ss >> volId;
			ss >> code;
			if( !(ss >> inclusive) )
				inclusive = 0;
			addVolumePropertyCode(volId, code, inclusive);
		} else if(primitive.compare("Face") == 0 ||
		          primitive.compare("face") == 0 )  {
			ss >> volId;
			ss >> locId;
			ss >> code;
			if( !(ss >> inclusive) )
				inclusive = 1;
			addFacePropertyCode(volId, locId, code, inclusive);
		} else if(primitive.compare("Line") == 0 ||
		          primitive.compare("line") == 0 )  {
			ss >> volId;
			ss >> locId;
			ss >> code;
			if( !(ss >> inclusive) )
				inclusive = 1;
			addLinePropertyCode(volId, locId, code, inclusive);
		} else if(primitive.compare("Vertex") == 0 ||
		          primitive.compare("vertex") == 0 )  {
			ss >> volId;
			ss >> locId;
			ss >> code;
			addVertexPropertyCode(volId, locId, code);
		} else if(primitive[0] == '/' && primitive[1] == '/') { // comment - continue without doing anything
			;
		} else {
			cerr << "Syntax error near: \"" << str_line << "\"\n";
			exit(1);
		}
	}
}

