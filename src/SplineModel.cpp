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
#include <stack>

#include <GoTools/geometry/ObjectHeader.h>

using namespace std;
using namespace Go;
using boost::shared_ptr;

SplineModel::SplineModel() {
	sl2g = new surfGlobNumber[0];
	vl2g = new volGlobNumber[0];
	sl2g = new surfGlobNumber[0];
	topology = new TopologySet();
	volumetric_model = false;
}

/**********************************************************************************//**
 * \brief Constructor
 * \param spline_surfaces All spline volumes to be considered part of this model
 *************************************************************************************/
SplineModel::SplineModel(std::vector<boost::shared_ptr<Go::SplineSurface> > &spline_surfaces) {
	sl2g = new surfGlobNumber[0];
	vl2g = new volGlobNumber[0];
	topology = new TopologySet(spline_surfaces);
	topology->buildTopology();
	spline_surfaces_ = spline_surfaces;
	volumetric_model = false;
}

/**********************************************************************************//**
 * \brief Constructor
 * \param spline_volumes All spline volumes to be considered part of this model
 *************************************************************************************/
SplineModel::SplineModel(std::vector<boost::shared_ptr<Go::SplineVolume> > &spline_volumes) {
	sl2g = new surfGlobNumber[0];
	vl2g = new volGlobNumber[0];
	topology = new TopologySet(spline_volumes);
	topology->buildTopology();
	spline_volumes_ = spline_volumes;
	volumetric_model = true;
}

//! \brief Destructor
SplineModel::~SplineModel() {
	if(topology)
		delete topology;
	delete[] vl2g;
	delete[] sl2g;
}

//! \brief Checks to see if the model is built up from trivariate volume splines
bool SplineModel::isVolumetricModel() const {
	return volumetric_model;
}

TopologySet* SplineModel::getTopology() {
	return topology;
}

std::vector<boost::shared_ptr<Go::SplineSurface> >& SplineModel::getSplineSurfaces() {
	return spline_surfaces_;
}

std::vector<boost::shared_ptr<Go::SplineVolume> >& SplineModel::getSplineVolumes() {
	return spline_volumes_;
}

/**********************************************************************************//**
 * \brief Make the parametrization of all spline-volumes to be right-hand-system
 * \returns true if any patches was reparameterized
 *
 * This ensures that all jacobians are evaluated positive, assuming of course that you
 * have a model descirption which is not wrapping itself inside out at any points (i.e.
 * jacboian must be non-negative or non-positive in the entire domain prior to
 * reparametrization).
 *************************************************************************************/
bool SplineModel::enforceRightHandSystem() {
	bool anything_switched = false;
	if(volumetric_model) {
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
	} else {
		for(uint i=0; i<spline_surfaces_.size(); i++) {
			// do test evaluation in the middle of the parametric domain
			double u = (spline_surfaces_[i]->startparam_u() + spline_surfaces_[i]->endparam_u()) / 2;
			double v = (spline_surfaces_[i]->startparam_v() + spline_surfaces_[i]->endparam_v()) / 2;
			vector<Point> results(3); // one position and two derivatives
			spline_surfaces_[i]->point(results, u, v, 1);
			Point normal = results[1] % results[2];
			double jacobian = normal[2];
			if(jacobian < 0) {
				anything_switched = true;
				spline_surfaces_[i]->reverseParameterDirection(true);
			}
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

/**********************************************************************************//**
 * \brief performs the minimum number of refinements required to contain an analysis
 *        suitable model
 * \param patchId  the patch requested for the initial refinement
 * \param parDir   parametric direction of the refinement
 * \param knot     knot vector inserted
 *
 * This method will insert the requested knot in the patch, as well as propagte this
 * refinement to all neighbouring patches and their neigbours and so on to ensure that
 * the model remains "analysis suitable" in a "corner-to-corner" fashion.
 *************************************************************************************/
void SplineModel::knot_insert(int patchId, int parDir, double knot) {
	if(volumetric_model) {
		cerr << "Refinement for volumes is not yet implemented\n";
		exit(142990);
	/*
		vector<bool> patchIds(topology->numbVolumes(), false);

		vector<bool> minParDirU(topology->numbVolumes(), false);
		vector<bool> minParDirV(topology->numbVolumes(), false);
		vector<bool> minParDirW(topology->numbVolumes(), false);
		vector<bool> maxParDirU(topology->numbVolumes(), false);
		vector<bool> maxParDirV(topology->numbVolumes(), false);
		vector<bool> maxParDirW(topology->numbVolumes(), false);

		vector<int>  parameterDirection(topology->numbVolumes(), -1);
		vector<bool> swapDir(topology->numbVolumes(), false);
		int otherSide[] = {1,0,  3,2,  5,4}; // given face index i, the other side is index otherSide[i]

		double knotStart = spline_volumes_[patchId]->startparam(parDir);
		double knotEnd   = spline_volumes_[patchId]->endparam(parDir);
		newKnotRatio     = (knot - knotStart) / (knotEnd - knotStart);

		stack<Volume*> volStack ;
		Volume *vol = topology->getVolume(patchId);
		Face   *f   = vol->face[parDir*2];
		vector<Volume*> neighb =  f->volume;
		for(uint i=0; i<neighb.size(); i++) {
		}
	*/
	} else {
		Line *l;
		set<Face*>::iterator it;
		vector<bool> minorRefine[2];
		vector<bool> majorRefine[2];
		minorRefine[0].insert(minorRefine[0].begin(), topology->numbFaces(), false);
		minorRefine[1].insert(minorRefine[1].begin(), topology->numbFaces(), false);
		majorRefine[0].insert(majorRefine[0].begin(), topology->numbFaces(), false);
		majorRefine[1].insert(majorRefine[1].begin(), topology->numbFaces(), false);

		int otherSide[] = {1,0,  3,2}; // given face index i, the other side is index otherSide[i]
		double knotStart = (parDir==0) ? spline_surfaces_[patchId]->startparam_u() : spline_surfaces_[patchId]->startparam_v();
		double knotEnd   = (parDir==0) ? spline_surfaces_[patchId]->endparam_u()   : spline_surfaces_[patchId]->endparam_v();
		double knotRel   = (knot - knotStart) / (knotEnd - knotStart);
		Face *f = topology->getFace(patchId);
		minorRefine[parDir][f->id] = true;

		stack<Face*> faceStack;
		stack<bool> minor;
		stack<int> edgeIn;

		l = f->line[parDir*2];
		vector<int> numb, runningParDir, parStep;
		if(! l->degen) {
			// l is degenerate or has one or two corresponding faces
			bool rightWay, sameWay;
			for(it=l->face.begin(); it!=l->face.end(); it++) {
				(*it)->getEdgeEnumeration(l, numb, runningParDir, parStep);
				if(*it == f)
					rightWay = (parStep[0]==1);
			}
			for(it=l->face.begin(); it!=l->face.end(); it++) {
				if(*it == f)
					continue;
				(*it)->getEdgeEnumeration(l, numb, runningParDir, parStep);
				int refParDir = runningParDir[0];
				sameWay = rightWay == (parStep[0]==1);
				if(sameWay) {
					if(!minorRefine[refParDir][(*it)->id]) {
						minorRefine[refParDir][(*it)->id] = true;
						faceStack.push(*it);
						minor.push(true);
						edgeIn.push(numb[0]);
					}
				} else {
					if(!majorRefine[refParDir][(*it)->id]) {
						majorRefine[refParDir][(*it)->id] = true;
						faceStack.push(*it);
						minor.push(false);
						edgeIn.push(numb[0]);
					}
				}
			}
		}

		faceStack.push(f);
		minor.push(true);
		edgeIn.push(parDir*2);

		while(!faceStack.empty()) {
			f = faceStack.top();
			int eIn = edgeIn.top();
			bool isMinor = minor.top();
			faceStack.pop();
			edgeIn.pop();
			minor.pop();

			l = f->line[otherSide[eIn]];
			vector<int> numb, runningParDir, parStep;
			if(! l->degen) {
				// l is degenerate or has one or two corresponding faces
				bool rightWay, sameWay;
				for(it=l->face.begin(); it!=l->face.end(); it++) {
					(*it)->getEdgeEnumeration(l, numb, runningParDir, parStep);
					if(*it == f)
						rightWay = (parStep[0]==1);
				}
				for(it=l->face.begin(); it!=l->face.end(); it++) {
					if(*it == f)
						continue;
					(*it)->getEdgeEnumeration(l, numb, runningParDir, parStep);
					int refParDir = runningParDir[0];
					sameWay = rightWay == (parStep[0]==1);
					if(sameWay == isMinor) {
						if(!minorRefine[refParDir][(*it)->id]) {
							minorRefine[refParDir][(*it)->id] = true;
							faceStack.push(*it);
							minor.push(true);
							edgeIn.push(numb[0]);
						}
					} else {
						if(!majorRefine[refParDir][(*it)->id]) {
							majorRefine[refParDir][(*it)->id] = true;
							faceStack.push(*it);
							minor.push(false);
							edgeIn.push(numb[0]);
						}
					}
				}
			}
		} // end faceStack loop
		for(uint i=0; i<spline_surfaces_.size(); i++) {
			double start_u = spline_surfaces_[i]->startparam_u();
			double stop_u  = spline_surfaces_[i]->endparam_u();
			double start_v = spline_surfaces_[i]->startparam_v();
			double stop_v  = spline_surfaces_[i]->endparam_v();
			if(minorRefine[0][i])
				spline_surfaces_[i]->insertKnot_u( (1-knotRel)*start_u + knotRel*stop_u);
			if(minorRefine[1][i])
				spline_surfaces_[i]->insertKnot_v( (1-knotRel)*start_v + knotRel*stop_v);
			if(majorRefine[0][i])
				spline_surfaces_[i]->insertKnot_u( knotRel*start_u + (1-knotRel)*stop_u);
			if(majorRefine[1][i])
				spline_surfaces_[i]->insertKnot_v( knotRel*start_v + (1-knotRel)*stop_v);
		}
	} // end if surfaceModel
}

/**********************************************************************************//**
 * \brief resolve boundary layer refinement
 * \param patchId  the patch requested for the initial refinement
 * \param parDir   parametric direction of the refinement
 * \param start    if the boundary layer is at the parametric start or end value
 * \param scale    ratio between the new innermost element and the previous innermost
 *                 element
 * \param n        number of new knot lines inserted
 *
 * This method will refine the elements closest to the requested boundary by splitting
 * with a suitable ratio. As with knot_insert, this refinement may propagate throughout 
 * the other patches to ensure a consistent model.
 *************************************************************************************/
void SplineModel::boundary_layer_refinement(int patchId, int parDir, bool start, double scale, int n) {
	if(n<1)
		return;
	if(volumetric_model) {
		cerr << "Refinement for volumes is not yet implemented\n";
		exit(142983);
	} else {
		vector<double> simpleKnots;
		if(parDir==0)
			spline_surfaces_[patchId]->basis_u().knotsSimple(simpleKnots);
		else
			spline_surfaces_[patchId]->basis_v().knotsSimple(simpleKnots);
		double pow = 1;
		double sum = 0;
		for(int i=0; i<=n; i++) {
			sum += pow;
			pow *= scale;
		}
		double minor   = (start) ? simpleKnots[1] : simpleKnots[simpleKnots.size()-2];
		double major   = (start) ? simpleKnots[0] : simpleKnots[simpleKnots.size()-1];
		double alpha   = 1.0/sum;
		sum = 0;
		pow = 1;
		for(int i=0; i<n; i++) {
			sum += alpha*pow;
			pow *= scale;
			double newKnot = sum*major + (1-sum)*minor;
			knot_insert(patchId, parDir, newKnot);
		}
	}
}

/**********************************************************************************//**
 * \brief uniform h-refinement on all patches in all directions
 *************************************************************************************/
void SplineModel::uniform_h_refine() {
	if(volumetric_model) {
		for(uint i=0; i<spline_volumes_.size(); i++) {
			vector<double> uniqueKnots;
			for(int dir=0; dir<3; dir++) {
				uniqueKnots.clear();
				spline_volumes_[i]->basis(dir).knotsSimple(uniqueKnots);
				for(uint j=0; j<uniqueKnots.size()-1; j++)
					spline_volumes_[i]->insertKnot(dir, (uniqueKnots[j]+uniqueKnots[j+1])/2);
			}
		}
	} else {
		for(uint i=0; i<spline_surfaces_.size(); i++) {
			vector<double> uniqueKnots_u;
			vector<double> uniqueKnots_v;
			spline_surfaces_[i]->basis_u().knotsSimple(uniqueKnots_u);
			for(uint j=0; j<uniqueKnots_u.size()-1; j++)
				spline_surfaces_[i]->insertKnot_u((uniqueKnots_u[j]+uniqueKnots_u[j+1])/2);
			spline_surfaces_[i]->basis_v().knotsSimple(uniqueKnots_v);
			for(uint j=0; j<uniqueKnots_v.size()-1; j++)
				spline_surfaces_[i]->insertKnot_v((uniqueKnots_v[j]+uniqueKnots_v[j+1])/2);
		}
	}
}

/**********************************************************************************//**
 * \brief uniform p-refinement on all patches in all directions
 *************************************************************************************/
void SplineModel::uniform_p_refine() {
	if(volumetric_model)
		for(uint i=0; i<spline_volumes_.size(); i++)
			spline_volumes_[i]->raiseOrder(1,1,1);
	else 
		for(uint i=0; i<spline_surfaces_.size(); i++)
			spline_surfaces_[i]->raiseOrder(1,1);
}




//! \brief generate the local to global enumerations
void SplineModel::generateGlobalNumbersPETSc() 
{
  uint iu;
  int  i, j, e;

  vector<int>::iterator pos;
  set<Volume*>::iterator v1, v2;

  delete[] vl2g;
  vl2g = new volGlobNumber[spline_volumes_.size()];

  // initialize all vl2g-variables 
  for(iu = 0;iu < spline_volumes_.size();iu++) {
    for(j = 0;j < 8;j++)
      vl2g[iu].vertex[j] = -1;
    for(j = 0;j < 12;j++) {
      vl2g[iu].edge[j] = -1;
      vl2g[iu].edge_incr[j] = 0;
    }
    for(j = 0;j < 6;j++) {
      vl2g[iu].surface[j] = -1;
      vl2g[iu].surface_incr_i[j] = 0;
      vl2g[iu].surface_incr_j[j] = 0;
    }
    vl2g[iu].volume = -1;
  }

  // Enumerate global nodes
  int glob_i = 0;

  // Assign global node numbers for all spline volumes
  for (v1=topology->volume_begin();v1 != topology->volume_end();v1++) {
    // Assign node numbers for vertices
    for (i = 0;i < 8;i++) 
      if (vl2g[(*v1)->id].vertex[i] < 0) {
	Vertex* v = (*v1)->corner[i];
	
	bool dup = false;
	int id   = glob_i;
	for (j = 0;j < i;j++)
	  if (v == (*v1)->corner[j]) {
	    dup = true;
	    id = vl2g[(*v1)->id].vertex[j];
	    break;
	  }
	
	for (v2 = v->volume.begin();v2 != v->volume.end();v2++) {
	  vector<int> corners = (*v2)->getVertexEnumeration(v);
	  for(pos = corners.begin(); pos != corners.end(); pos++)
	    vl2g[(*v2)->id].vertex[*pos] = id;
	}
	
	if (!dup)
	  glob_i++;
      }
  
    // Assign node numbers for edges
    for (e = 0;e < 12;e++) 
      if (vl2g[(*v1)->id].edge[e] < 0) {
	vector<int> edges, dir, step;
	
	Line* l1 = (*v1)->line[e];	
	Volume *first_vol = *(l1->volume.begin());
	shared_ptr<SplineVolume> sv = spline_volumes_[first_vol->id]; // this edge SHOULD have at least one volume-connection
	first_vol->getEdgeEnumeration(l1, edges, dir, step);
	int numCoeffs = sv->numCoefs(dir[0])-2; // coefs_here SHOULD be identical for all lines in all volumes for this line
	
	if(numCoeffs > 0) {
	  for (v2 = l1->volume.begin();v2 != l1->volume.end();v2++) {
	    (*v2)->getEdgeEnumeration(l1,edges,dir,step);
	    
	    for(iu = 0;iu < edges.size();iu++) 
	      if(l1->degen) {
		vector<int> corners = (*v2)->getVertexEnumeration(l1->v1);
		int c0 = corners[0]; // should be the same global number on all elements of "corners"
		vl2g[(*v2)->id].edge[edges[iu]]      = vl2g[(*v2)->id].vertex[c0];
		vl2g[(*v2)->id].edge_incr[edges[iu]] = 0; 
	      } 
	      else {
		vl2g[(*v2)->id].edge[edges[iu]]      = (step[iu]==1) ? glob_i : glob_i + numCoeffs - 1;
		vl2g[(*v2)->id].edge_incr[edges[iu]] =  step[iu] ;
	      }
	  }
	}
	
	if( !(l1->degen) )
	  glob_i += numCoeffs;
      }
 
    // Assign node numbers for faces
    for (i = 0;i < 6;i++) 
      if (vl2g[(*v1)->id].surface[i] < 0) {
	Face *f = (*v1)->face[i];
	/*  Face Index (logic behind numCoef_u/v)
	 *  0-1 : surface umin/umax - numcCefs(1) X numCoefs(2)
	 *  2-3 : surface vmin/vmax - numcCefs(0) X numCoefs(2)
	 *  4-5 : surface wmin/wmax - numcCefs(0) X numCoefs(1)
	 */
	shared_ptr<SplineVolume> s_volume = spline_volumes_[f->volume[0]->id];
	int numCoeff_u  = s_volume->numCoefs(  (f->face[0] < 2) ) - 2;
	int numCoeff_v  = s_volume->numCoefs(2-(f->face[0] > 3) ) - 2;
	int numCoeff    = numCoeff_u*numCoeff_v;
	
	for(iu = 0;iu < f->volume.size();iu++) {
	  Volume *vol = f->volume[iu];
	  int faceId  = f->face[iu];
	  bool uv_flip   = f->uv_flip[iu];
	  bool u_reverse = f->u_reverse[iu];
	  bool v_reverse = f->v_reverse[iu];
	  
	  if(f->degen1 && f->degen2) {
	    vector<int> corners = Vertex::getVertexEnumeration(faceId);
	    vl2g[vol->id].surface[faceId]        = vl2g[vol->id].vertex[corners[0]];
	    vl2g[vol->id].surface_incr_i[faceId] = 0;
	    vl2g[vol->id].surface_incr_j[faceId] = 0;
	  } 
	  else if(f->degen1) {
	    vector<int> lines = Line::getLineEnumeration(faceId);
	    
	    if(uv_flip) {
	      vl2g[vol->id].surface[faceId]        = vl2g[vol->id].edge[lines[0]];
	      vl2g[vol->id].surface_incr_i[faceId] = vl2g[vol->id].edge_incr[lines[0]];
	      vl2g[vol->id].surface_incr_j[faceId] = 0;
	    } 
	    else {
	      vl2g[vol->id].surface[faceId]        = vl2g[vol->id].edge[lines[2]];
	      vl2g[vol->id].surface_incr_i[faceId] = 0;
	      vl2g[vol->id].surface_incr_j[faceId] = vl2g[vol->id].edge_incr[lines[2]];
	    }
	  } 
	  else if(f->degen2) {
	    vector<int> lines = Line::getLineEnumeration(faceId);
	    if(uv_flip) {
	      vl2g[vol->id].surface[faceId]        = vl2g[vol->id].edge[lines[2]];
	      vl2g[vol->id].surface_incr_i[faceId] = 0;
	      vl2g[vol->id].surface_incr_j[faceId] = vl2g[vol->id].edge_incr[lines[2]];
	    } 
	    else {
	      vl2g[vol->id].surface[faceId]        = vl2g[vol->id].edge[lines[0]];
	      vl2g[vol->id].surface_incr_i[faceId] = vl2g[vol->id].edge_incr[lines[0]];
	      vl2g[vol->id].surface_incr_j[faceId] = 0;
	    }
	  } 
	  else {
	    if(!u_reverse && !v_reverse)
	      vl2g[vol->id].surface[faceId] = glob_i ;
	    else if(u_reverse && v_reverse)
	      vl2g[vol->id].surface[faceId] = glob_i + numCoeff - 1;
	    else if(v_reverse == uv_flip)
	      vl2g[vol->id].surface[faceId] = glob_i + numCoeff_u - 1;
	    else // u_reverse == uv_flip
	      vl2g[vol->id].surface[faceId] = glob_i + numCoeff - numCoeff_u;
	    
	    if(uv_flip) {
	      vl2g[vol->id].surface_incr_i[faceId] = numCoeff_u;
	      vl2g[vol->id].surface_incr_j[faceId] = 1;
	    } 
	    else {
	      vl2g[vol->id].surface_incr_i[faceId] = 1;
	      vl2g[vol->id].surface_incr_j[faceId] = numCoeff_u;
	    }
	    if(u_reverse)
	      vl2g[vol->id].surface_incr_i[faceId] *= -1;
	    if(v_reverse)
	      vl2g[vol->id].surface_incr_j[faceId] *= -1;
	  }
	}

	if(!f->isDegen()) 
	  glob_i += numCoeff;
      }
    
    // Assign node numbers to volumes
    shared_ptr<SplineVolume> s_volume = spline_volumes_[(*v1)->id];
    int numCoeff_u = s_volume->numCoefs(0)-2;
    int numCoeff_v = s_volume->numCoefs(1)-2;
    int numCoeff_w = s_volume->numCoefs(2)-2;
    int numCoeff   = numCoeff_u * numCoeff_v * numCoeff_w;
    if(numCoeff > 0) {
      vl2g[(*v1)->id].volume = glob_i;
      glob_i += numCoeff;
    }
  }
}


//! \brief generate the local to global enumerations
void SplineModel::generateGlobalNumbers() {
	
	if(volumetric_model) {
		delete[] vl2g;
		vl2g = new volGlobNumber[spline_volumes_.size()];
		// initalize all vl2g-variables 
		for(uint i=0; i<spline_volumes_.size(); i++) {
			for(int j=0; j<8; j++)
				vl2g[i].vertex[j] = -1;
			for(int j=0; j<12; j++) {
				vl2g[i].edge[j] = -1;
				vl2g[i].edge_incr[j] = 0;
			}
			for(int j=0; j<6; j++) {
				vl2g[i].surface[j] = -1;
				vl2g[i].surface_incr_i[j] = 0;
				vl2g[i].surface_incr_j[j] = 0;
			}
			vl2g[i].volume = -1;
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
					vl2g[(*v)->id].vertex[*pos] = glob_i;
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
							vl2g[(*v)->id].edge[numb[i]]      = vl2g[(*v)->id].vertex[corner_pos];
							vl2g[(*v)->id].edge_incr[numb[i]] = 0;
						} else {
							vl2g[(*v)->id].edge[numb[i]]      = (parStep[i]==1) ? glob_i : glob_i+coefs_here-1;
							vl2g[(*v)->id].edge_incr[numb[i]] =  parStep[i] ;
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
					vl2g[vol->id].surface[faceId]        = vl2g[vol->id].vertex[corners[0]];
					vl2g[vol->id].surface_incr_i[faceId] = 0;
					vl2g[vol->id].surface_incr_j[faceId] = 0;
				} else if(f->degen1) {
					vector<int> lines = Line::getLineEnumeration(faceId);
					if(uv_flip) {
						vl2g[vol->id].surface[faceId]        = vl2g[vol->id].edge[lines[0]];
						vl2g[vol->id].surface_incr_i[faceId] = vl2g[vol->id].edge_incr[lines[0]];
						vl2g[vol->id].surface_incr_j[faceId] = 0;
					} else {
						vl2g[vol->id].surface[faceId]        = vl2g[vol->id].edge[lines[2]];
						vl2g[vol->id].surface_incr_i[faceId] = 0;
						vl2g[vol->id].surface_incr_j[faceId] = vl2g[vol->id].edge_incr[lines[2]];
					}
				} else if(f->degen2) {
					vector<int> lines = Line::getLineEnumeration(faceId);
					if(uv_flip) {
						vl2g[vol->id].surface[faceId]        = vl2g[vol->id].edge[lines[2]];
						vl2g[vol->id].surface_incr_i[faceId] = 0;
						vl2g[vol->id].surface_incr_j[faceId] = vl2g[vol->id].edge_incr[lines[2]];
					} else {
						vl2g[vol->id].surface[faceId]        = vl2g[vol->id].edge[lines[0]];
						vl2g[vol->id].surface_incr_i[faceId] = vl2g[vol->id].edge_incr[lines[0]];
						vl2g[vol->id].surface_incr_j[faceId] = 0;
					}
				} else {
					if(!u_reverse && !v_reverse)
						vl2g[vol->id].surface[faceId] = glob_i ;
					else if(u_reverse && v_reverse)
						vl2g[vol->id].surface[faceId] = glob_i + coefs_here - 1;
					else if(v_reverse == uv_flip)
						vl2g[vol->id].surface[faceId] = glob_i + numCoef_u - 1;
					else // u_reverse == uv_flip
						vl2g[vol->id].surface[faceId] = glob_i + coefs_here - numCoef_u;
					
					if(uv_flip) {
						vl2g[vol->id].surface_incr_i[faceId] = numCoef_u;
						vl2g[vol->id].surface_incr_j[faceId] = 1;
					} else {
						vl2g[vol->id].surface_incr_i[faceId] = 1;
						vl2g[vol->id].surface_incr_j[faceId] = numCoef_u;
					}
					if(u_reverse)
						vl2g[vol->id].surface_incr_i[faceId] *= -1;
					if(v_reverse)
						vl2g[vol->id].surface_incr_j[faceId] *= -1;
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
				vl2g[i].volume = glob_i;
				glob_i += coefs_here;
			}
		}
	} else {
		delete[] sl2g;
		sl2g = new surfGlobNumber[spline_surfaces_.size()];
		// initalize all sl2g-variables 
		for(uint i=0; i<spline_surfaces_.size(); i++) {
			for(int j=0; j<4; j++)
				sl2g[i].vertex[j] = -1;
			for(int j=0; j<4; j++) {
				sl2g[i].edge[j] = -1;
				sl2g[i].edge_incr[j] = 0;
			}
			sl2g[i].surface = -1;
		}
		
		set<Vertex*>::iterator v_it;
		set<Line*>::iterator   l_it ;
		set<Face*>::iterator   f_it;
		
		vector<shared_ptr<SplineSurface> >::iterator surf_it;
		
		// Enumerate vertices, lines and finally surfaces, in that order
		int glob_i = 0;
		
		vector<int>::iterator pos;
		
		// for all VERTICES, assign number
		for(v_it=topology->vertex_begin(); v_it != topology->vertex_end(); v_it++) {
			for(f_it=(*v_it)->face.begin(); f_it != (*v_it)->face.end(); f_it++) {
				vector<int> corners = (*f_it)->getVertexEnumeration(*v_it);
				for(pos=corners.begin(); pos!=corners.end(); pos++)
					sl2g[(*f_it)->id].vertex[*pos] = glob_i;
			}
			glob_i++;
		}
		
		// for all EDGES, assign startnumber and increment
		for(l_it=topology->line_begin(); l_it != topology->line_end(); l_it++) {
			vector<int> numb;
			vector<int> parDir;
			vector<int> parStep;
			Face *first_face = (*(*l_it)->face.begin());
			shared_ptr<SplineSurface> sf = spline_surfaces_[first_face->id]; // this edge should have at least one face-connection
			first_face->getEdgeEnumeration(*l_it, numb, parDir, parStep);
			int coefs_here = (parDir[0]==0) ? sf->numCoefs_u()-2 : sf->numCoefs_v()-2; // coefs_here should be identical for all lines in all faces for this line
			if(coefs_here > 0) {
				for(f_it=(*l_it)->face.begin(); f_it != (*l_it)->face.end(); f_it++) {
					(*f_it)->getEdgeEnumeration(*l_it, numb, parDir, parStep); // no worries: numb, parDir & parStep all get cleaned inside function call before returned
					for(uint i=0; i<numb.size(); i++) {
						if((*l_it)->degen) {
							vector<int> corners = (*f_it)->getVertexEnumeration((*l_it)->v1);
							int corner_pos = corners[0]; // should be the same global number on all elements of "corners"
							sl2g[(*f_it)->id].edge[numb[i]]      = sl2g[(*f_it)->id].vertex[corner_pos];
							sl2g[(*f_it)->id].edge_incr[numb[i]] = 0;
						} else {
							sl2g[(*f_it)->id].edge[numb[i]]      = (parStep[i]==1) ? glob_i : glob_i+coefs_here-1;
							sl2g[(*f_it)->id].edge_incr[numb[i]] =  parStep[i] ;
						}
					}
				}
				if( !(*l_it)->degen )
					glob_i += coefs_here;
			}
		}
		
		// for all FACES assign startnumber
		for(uint i=0; i<spline_surfaces_.size(); i++) {
			shared_ptr<SplineSurface> ss = spline_surfaces_[i];
			int numCoef_u = ss->numCoefs_u()-2;
			int numCoef_v = ss->numCoefs_v()-2;
			int coefs_here  = numCoef_u * numCoef_v;
			if(coefs_here > 0) {
				sl2g[i].surface = glob_i;
				glob_i += coefs_here;
			}
		}
	}
}

void SplineModel::writeGlobalNumberOrdering(std::ostream &os) const {
	if(volumetric_model) {
		for(uint i=0; i<spline_volumes_.size(); i++) {
			os << i << endl;
			for(int j=0; j<8; j++)
				os << vl2g[i].vertex[j] << " ";
			os << endl;
			for(int j=0; j<12; j++)
				os << vl2g[i].edge[j] << " " << vl2g[i].edge_incr[j] << endl;
			for(int j=0; j<6; j++)
				os << vl2g[i].surface[j] << " " << vl2g[i].surface_incr_i[j] << " " << vl2g[i].surface_incr_j[j] << endl;
			os << vl2g[i].volume << endl;
		}
	} else {
		for(uint i=0; i<spline_surfaces_.size(); i++) {
			os << i << endl;
			for(int j=0; j<4; j++)
				os << sl2g[i].vertex[j] << " ";
			os << endl;
			for(int j=0; j<4; j++)
				os << sl2g[i].edge[j] << " " << sl2g[i].edge_incr[j] << endl;
			os << sl2g[i].surface << endl;
		}
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
	volGlobNumber propCodes[spline_volumes_.size()];
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
	if(volumetric_model) {
		for(uint i=0; i<spline_volumes_.size(); i++) {
			spline_volumes_[i]->writeStandardHeader(os);
			spline_volumes_[i]->write(os);
		}
	} else {
		for(uint i=0; i<spline_surfaces_.size(); i++) {
			spline_surfaces_[i]->writeStandardHeader(os);
			spline_surfaces_[i]->write(os);
		}
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
			topology->addPatch(v);
			volumetric_model = true;
		} else if(head.classType() == Class_SplineSurface) {
			shared_ptr<SplineSurface> s(new SplineSurface());
			s->read(is);
			spline_surfaces_.push_back(s);
			topology->addPatch(s);
			surface_model = true;
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

