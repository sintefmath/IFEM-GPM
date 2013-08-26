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
SplineModel::SplineModel(std::vector<SurfacePointer> &spline_surfaces) {
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
SplineModel::SplineModel(std::vector<VolumePointer> &spline_volumes) {
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

void SplineModel::setTopologyTolerance(double tol) {
	topology->setTolerance(tol);
}

//! \brief Checks to see if the model is built up from trivariate volume splines
bool SplineModel::isVolumetricModel() const {
	return volumetric_model;
}


void SplineModel::buildTopology(std::vector<bool>* periodic) {
  topology->buildTopology(periodic);
}


TopologySet* SplineModel::getTopology() {
	return topology;
}

std::vector<SurfacePointer>& SplineModel::getSplineSurfaces() {
	return spline_surfaces_;
}

std::vector<VolumePointer>& SplineModel::getSplineVolumes() {
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
			vector<Go::Point> results(4); // one position and three derivatives
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
			vector<Go::Point> results(3); // one position and two derivatives
			spline_surfaces_[i]->point(results, u, v, 1);
                        Go::Point normal = results[1] % results[2];
			double jacobian = normal[2];
			if(jacobian < 0) {
				anything_switched = true;
				spline_surfaces_[i]->reverseParameterDirection(false);
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
bool SplineModel::addVolumePropertyCode(int volId, const char* propCode, bool inclusive) {
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
 * \param patchId   patch identifier (natural numbering if read from a .g2-file)
 * \param faceId    local face number corresponding to the patch patchId, as defined in
 *                  the class Volume (volumeModel: values 0-5, surfaceModel: value IGNORED)
 * \param propCode  property code
 * \param inclusive Whether the corners or lines associated with the face should get 
 *                  the same code
 * \returns true if any corners or lines had previously been tagged with a different 
 *          code and has now been overwritten (only applicable if inclusive=true)
 *************************************************************************************/
bool SplineModel::addFacePropertyCode(int patchId, int faceId, const char* propCode, bool inclusive) {
	Face* f;
	if(volumetric_model) {
		Volume* v = topology->getVolume(patchId);
		f = v->face[faceId];
	} else {
		f = topology->getFace(patchId);
	}
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
 * \param patchId   volume identifier (natural numbering if read from a .g2-file)
 * \param lineId    local line number corresponding to the patch patchId, as defined in
 *                  the class Volume (volumeModel: values 0-11, surfaceModel: 0-3)
 * \param propCode  property code
 * \param inclusive Whether the corners associated with the face should get the same code
 * \returns true if any corners had previously been tagged with a different code and has 
 *          now been overwritten (only applicable if inclusive=true)
 *************************************************************************************/
bool SplineModel::addLinePropertyCode(int patchId, int lineId, const char* propCode, bool inclusive) {
	Line* l;
	if(volumetric_model) {
		Volume* v = topology->getVolume(patchId);
		l = v->line[lineId];
	} else {
		Face* f = topology->getFace(patchId);
		l = f->line[lineId];
	}
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
 * \param patchId   volume identifier (natural numbering if read from a .g2-file)
 * \param vertId    local vertex number corresponding to the patch patchId, as defined in
 *                  the class Volume (volumeModel: values 0-7, surfaceModel: 0-3)
 * \param propCode  property code
 *************************************************************************************/
void SplineModel::addVertexPropertyCode(int patchId, int vertId, const char* propCode) {
	if(volumetric_model) {
		Volume* v = topology->getVolume(patchId);
		v->corner[vertId]->bc_code = propCode;
	} else {
		Face* f = topology->getFace(patchId);
		f->corner[vertId]->bc_code = propCode;
	}
}

const char* SplineModel::getVolumePropertyCode(int patchId) {
	if(volumetric_model) {
		Volume* v = topology->getVolume(patchId);
		return v->material_code;
	} else {
		return "";
	}
}

const char* SplineModel::getFacePropertyCode(int patchId, int faceId) {
	if(volumetric_model) {
		Volume* v = topology->getVolume(patchId);
		return v->face[faceId]->bc_code;
	} else {
		Face* f = topology->getFace(patchId);
		return f->bc_code;
	}
}

const char* SplineModel::getLinePropertyCode(int patchId, int lineId) {
	if(volumetric_model) {
		Volume* v = topology->getVolume(patchId);
		return v->line[lineId]->bc_code;
	} else {
		Face* f = topology->getFace(patchId);
		return f->line[lineId]->bc_code;
	}
}

const char* SplineModel::getVertexPropertyCode(int patchId, int vertId) {
	if(volumetric_model) {
	Volume* v = topology->getVolume(patchId);
	return v->corner[vertId]->bc_code;
	} else {
		Face* f = topology->getFace(patchId);
		return f->corner[vertId]->bc_code;
	}
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

		vector<bool> patchIds(topology->numbVolumes(),   false);

		vector<bool> minorRefine[3];
		vector<bool> majorRefine[3];
		minorRefine[0].insert(minorRefine[0].begin(), topology->numbFaces(), false);
		minorRefine[1].insert(minorRefine[1].begin(), topology->numbFaces(), false);
		minorRefine[2].insert(minorRefine[2].begin(), topology->numbFaces(), false);
		majorRefine[0].insert(majorRefine[0].begin(), topology->numbFaces(), false);
		majorRefine[1].insert(majorRefine[1].begin(), topology->numbFaces(), false);
		majorRefine[2].insert(majorRefine[2].begin(), topology->numbFaces(), false);

		// Note: "minor" in uv-,uw- and vw-pairs refers to the first parametric direction (i.e. u,u and v respectively)

		int otherSide1[6]    =   {1,0,  3,2,  5,4};   // given face index i, the opposing output face is index otherSide[i]
		int otherSide2[2][6] = { {2,2,  0,0,  0,0},   // the first side face, given input face and running minor parameter direction
		                         {4,4,  4,4,  2,2} }; // the first side face, given running direction as "major parameter"
		int otherSide3[2][6] = { {3,3,  1,1,  1,1},   // second side face, when running parameter is "minor" (refining parameter = "major")
		                         {5,5,  5,5,  3,3} }; // second side face, when running parameter is "major" (refining parameter = "minor")
		int refDir[2][6]     = { {2,2,  2,2,  1,1},   // given input face and running parameter, what parametric
		                         {1,1,  0,0,  0,0} }; // direction should be refined

		double knotStart = spline_volumes_[patchId]->startparam(parDir);
		double knotEnd   = spline_volumes_[patchId]->endparam(parDir);
		double knotRel   = (knot - knotStart) / (knotEnd - knotStart); // knotRel is the actual % from the edge which is going to be refined
		Volume *vol = topology->getVolume(patchId);

		stack<Volume*> volStack;
		stack<bool> minor;            // in v-refining, true if % is evaluated from vmin (false if % is evaluated from vmax)
		stack<bool> refiningMinor;    // in vw-face, true if refining direction is minor (i.e. v)
		stack<int> faceIn;
		Volume     *v;
		

		// start by pushing one face in, and all the rest will follow
		volStack.push(vol);
		minor.push(parDir != 0);
		refiningMinor.push(true);
		faceIn.push((parDir*2+3)%6);

		while(!volStack.empty()) {
			Volume *v       = volStack.top();
			int    fIn      = faceIn.top();
			bool   isMinor  = minor.top();
			bool   refMinor = refiningMinor.top();
			volStack.pop();
			faceIn.pop();
			minor.pop();
			refiningMinor.pop();

			// tag acutal refinement
			int parametricRefine = refDir[!isMinor][fIn];
			if(refMinor) {
				if(minorRefine[parametricRefine][v->id])
					continue; // already done with this one, keep going, no propagation
				minorRefine[parametricRefine][v->id] = true;
			} else {
				if(majorRefine[parametricRefine][v->id])
					continue; // already done with this one, keep going, no propagation
				majorRefine[parametricRefine][v->id] = true;
			}

			// propegate to 3 edge faces
			int faceOut[4] = {otherSide1[fIn],
			                  otherSide2[!isMinor][fIn],
			                  otherSide3[!isMinor][fIn],
							  fIn}; // <-- this last guy is only used for the FIRST volume to ensure 
							        //     propegation in both directions
			for(int i=0; i<4; i++) {
				Face *f = v->face[faceOut[i]];
				// face degenerate to point does not propegate refinement
				// note however that faces degenerating to lines DO cause propegation
				if(f->degen1 && f->degen2)
					continue;
				vector<Volume*>::iterator it;
				bool thisFlipUV   = false;
				bool thisReverseU = false;
				bool thisReverseV = false;

				// uv-flip and reversing is wrt volume[0] which may not be equal to *v
				int j=0;
				for(it=f->volume.begin(); it!=f->volume.end(); it++, j++) {
					if(*it == v) {
						thisFlipUV   = f->uv_flip[j];
						thisReverseU = f->u_reverse[j];
						thisReverseV = f->v_reverse[j];
					}
				}

				if( (f->degen1 && !thisFlipUV && !isMinor) ||
				    (f->degen1 &&  thisFlipUV &&  isMinor))
					continue;
				if( (f->degen2 && !thisFlipUV &&  isMinor) ||
				    (f->degen2 &&  thisFlipUV && !isMinor))
					continue;

				// for all neighbouring volumes to this (possible degenerate) face
				j=0;
				for(it=f->volume.begin(); it!=f->volume.end(); it++, j++) {
					if(*it == v && f->id == fIn)
						continue;

					// add neighbouring volume to stack
					volStack.push(*it);

					bool uv_flip   = (f->uv_flip[j]   != thisFlipUV);
					bool u_reverse = (f->u_reverse[j] != thisReverseU);
					bool v_reverse = (f->v_reverse[j] != thisReverseV);
					// is "running index" the minor one?
					bool thisMinor = (parametricRefine == 1 && faceOut[i]>3) ||
					                 (parametricRefine == 2 );
					minor.push( thisMinor != uv_flip );

					// on the neighbouting volume, should we refine major?
					bool nextRefMinor = !thisMinor;
					nextRefMinor = ( !thisFlipUV &&  nextRefMinor && ! uv_flip && !u_reverse ) ||
					               ( !thisFlipUV &&  nextRefMinor &&   uv_flip && !v_reverse ) ||
								   ( !thisFlipUV && !nextRefMinor && ! uv_flip && !v_reverse ) ||
								   ( !thisFlipUV && !nextRefMinor &&   uv_flip && !u_reverse ) ||
								   (  thisFlipUV &&  nextRefMinor && ! uv_flip && !v_reverse ) ||
								   (  thisFlipUV &&  nextRefMinor &&   uv_flip && !u_reverse ) ||
								   (  thisFlipUV && !nextRefMinor && ! uv_flip && !u_reverse ) ||
								   (  thisFlipUV && !nextRefMinor &&   uv_flip && !v_reverse ) ;
					nextRefMinor = (nextRefMinor == refMinor);
					refiningMinor.push(nextRefMinor);
					
					// get local face enumeration on the neighbouring volume
					vector<int> nextFaceIn = (*it)->getSurfaceEnumeration(f);
					faceIn.push(nextFaceIn[0]);

				}
			}
		} // end volStack loop

		// perform knot insetion 
		for(uint i=0; i<spline_volumes_.size(); i++) {
			double start_u = spline_volumes_[i]->startparam(0);
			double stop_u  = spline_volumes_[i]->endparam(0);
			double start_v = spline_volumes_[i]->startparam(1);
			double stop_v  = spline_volumes_[i]->endparam(1);
			double start_w = spline_volumes_[i]->startparam(2);
			double stop_w  = spline_volumes_[i]->endparam(2);
			if(minorRefine[0][i])
				spline_volumes_[i]->insertKnot(0, (1-knotRel)*start_u + knotRel*stop_u);
			if(minorRefine[1][i])
				spline_volumes_[i]->insertKnot(1, (1-knotRel)*start_v + knotRel*stop_v);
			if(minorRefine[2][i])
				spline_volumes_[i]->insertKnot(2, (1-knotRel)*start_w + knotRel*stop_w);
			if(majorRefine[0][i])
				spline_volumes_[i]->insertKnot(0, knotRel*start_u + (1-knotRel)*stop_u);
			if(majorRefine[1][i])
				spline_volumes_[i]->insertKnot(1, knotRel*start_v + (1-knotRel)*stop_v);
			if(majorRefine[2][i])
				spline_volumes_[i]->insertKnot(2, knotRel*start_w + (1-knotRel)*stop_w);
		}

	} else {
		Line *l;
		set<Face*>::iterator it;
		// create a list over how each patch should be refined. minorRefine[0] means refine % is given from u_min to u_max,
		// while majorRefine[0] means percentage from u_max to u_min. Likewise with [1] for the v-parametric direction
		vector<bool> minorRefine[2];
		vector<bool> majorRefine[2];
		minorRefine[0].insert(minorRefine[0].begin(), topology->numbFaces(), false);
		minorRefine[1].insert(minorRefine[1].begin(), topology->numbFaces(), false);
		majorRefine[0].insert(majorRefine[0].begin(), topology->numbFaces(), false);
		majorRefine[1].insert(majorRefine[1].begin(), topology->numbFaces(), false);

		int otherSide[] = {1,0,  3,2}; // given face index i, the other side is index otherSide[i]
		double knotStart = (parDir==0) ? spline_surfaces_[patchId]->startparam_u() : spline_surfaces_[patchId]->startparam_v();
		double knotEnd   = (parDir==0) ? spline_surfaces_[patchId]->endparam_u()   : spline_surfaces_[patchId]->endparam_v();
		double knotRel   = (knot - knotStart) / (knotEnd - knotStart); // knotRel is the actual % from the edge which is going to be refined
		Face *f = topology->getFace(patchId);
		minorRefine[parDir][f->id] = true;

		stack<Face*> faceStack;
		stack<bool> minor;
		stack<int> edgeIn;

		// pick either the bottom or left line to figure out local orientation (store in "rightWay")
		l = f->line[2-parDir*2]; 
		vector<int> numb, runningParDir, parStep;
		if(! l->degen) { // l will always be degenerate or has one or two corresponding faces 
			bool rightWay, sameWay;
			for(it=l->face.begin(); it!=l->face.end(); it++) {
				(*it)->getEdgeEnumeration(l, numb, runningParDir, parStep);
				if(*it == f)
					rightWay = (parStep[0]==1);
			}
			for(it=l->face.begin(); it!=l->face.end(); it++) {
				if(*it == f)
					continue;
				(*it)->getEdgeEnumeration(l, numb, runningParDir, parStep); // vectors returned will always be of size 1
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

		// the loop above creates search recursion down/left, while these 3 lines create up/right
		faceStack.push(f);
		minor.push(true);
		edgeIn.push(2-parDir*2);

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

	vector<double> simpleKnots;

	if(volumetric_model) {
		spline_volumes_[patchId]->basis(parDir).knotsSimple(simpleKnots);
	} else {
		if(parDir==0)
			spline_surfaces_[patchId]->basis_u().knotsSimple(simpleKnots);
		else
			spline_surfaces_[patchId]->basis_v().knotsSimple(simpleKnots);
	}
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
void SplineModel::generateGlobalNumbersPETSc(bool mixed, int start) 
{
    int mx = 0;
  if (mixed) mx = 1;

  if (volumetric_model) {
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
    int glob_i = start;
    
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
	  VolumePointer sv = spline_volumes_[first_vol->id]; // this edge SHOULD have at least one volume-connection
	  first_vol->getEdgeEnumeration(l1, edges, dir, step);
	  int numCoeffs = sv->numCoefs(dir[0])-2+mx; // coefs_here SHOULD be identical for all lines in all volumes for this line
	  
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
	  VolumePointer s_volume = spline_volumes_[f->volume[0]->id];
	  int numCoeff_u  = s_volume->numCoefs(  (f->face[0] < 2) ) - 2+mx;
	  int numCoeff_v  = s_volume->numCoefs(2-(f->face[0] > 3) ) - 2+mx;
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
      VolumePointer s_volume = spline_volumes_[(*v1)->id];
      int numCoeff_u = s_volume->numCoefs(0)-2+mx;
      int numCoeff_v = s_volume->numCoefs(1)-2+mx;
      int numCoeff_w = s_volume->numCoefs(2)-2+mx;
      int numCoeff   = numCoeff_u * numCoeff_v * numCoeff_w;
      if(numCoeff > 0) {
	vl2g[(*v1)->id].volume = glob_i;
	glob_i += numCoeff;
      }
    }
  }
  else {
    vector<int>::iterator pos;
    set<Face*>::iterator  f_it1, f_it2;

    delete[] sl2g;
    sl2g = new surfGlobNumber[spline_surfaces_.size()];
    // Initialize all sl2g-variables
    for (uint i = 0;i < spline_surfaces_.size();i++) {
      for (int j = 0;j < 4;j++) {
	sl2g[i].vertex[j] = -1;
	sl2g[i].edge[j]   = -1;
	sl2g[i].edge_incr[j] = 0;
      }
      sl2g[i].surface = -1;
    }
                
    // Enumerate vertices, lines and finally surfaces, in that order
    int glob_i = start;
    
    // Assign global node numbers for all spline volumes
    for (f_it1 = topology->face_begin();f_it1 != topology->face_end();f_it1++) {
      for (int i = 0;i < 4;i++) 
	if (sl2g[(*f_it1)->id].vertex[i] < 0) {
	  Vertex* v = (*f_it1)->corner[i];
	  
	  bool dup = false;
	  int id   = glob_i;
	  for (int j = 0;j < i;j++)
	    if (v == (*f_it1)->corner[j]) {
	      dup = true;
	      id = sl2g[(*f_it1)->id].vertex[j];
	      break;
	    }
	  
	  for (f_it2 = topology->face_begin();f_it2 != topology->face_end();f_it2++) {
	    vector<int> corners = (*f_it2)->getVertexEnumeration(v);
	    for(pos = corners.begin(); pos != corners.end(); pos++)
	      sl2g[(*f_it2)->id].vertex[*pos] = id;
	  }
	  
	  if (!dup)
	    glob_i++;
	}

      // Assign node numbers for edges
      for (int e = 0;e < 4;e++)
	if (sl2g[(*f_it1)->id].edge[e] < 0) {
	  Line* l1 = (*f_it1)->line[e];
	  
	  vector<int> numb;
	  vector<int> parDir;
	  vector<int> parStep;
	  Face *first_face = *(l1->face.begin());
	  SurfacePointer sf = spline_surfaces_[first_face->id]; // this edge should have at least one face-connection
	  first_face->getEdgeEnumeration(l1, numb, parDir, parStep);
	  int numCoeffs = (parDir[0]==0) ? sf->numCoefs_u()-2+mx : sf->numCoefs_v()-2+mx; // coefs_here should be identical for all lines in all parDir, parStep);
	  if(numCoeffs > 0) {
	    for(f_it2 = l1->face.begin(); f_it2 != l1->face.end(); f_it2++) {
	      (*f_it2)->getEdgeEnumeration(l1, numb, parDir, parStep); // no worries: numb, parDir & parStep all get cleaned inside function call before returned
	      for(uint i=0; i<numb.size(); i++) {
		if (l1->degen) {
		  vector<int> corners = (*f_it2)->getVertexEnumeration(l1->v1);
		  int corner_pos = corners[0]; // should be the same global number on all elements of "corners"
		  sl2g[(*f_it2)->id].edge[numb[i]]      = sl2g[(*f_it2)->id].vertex[corner_pos];
		  sl2g[(*f_it2)->id].edge_incr[numb[i]] = 0;
		} else {
		  sl2g[(*f_it2)->id].edge[numb[i]]      = (parStep[i]==1) ? glob_i : glob_i+numCoeffs-1;
		  sl2g[(*f_it2)->id].edge_incr[numb[i]] =  parStep[i] ;
		}
	      }
	    }
	    if( !l1->degen )
	      glob_i += numCoeffs;
	  }
	}
     
      // Assign node numbers to surfaces
      SurfacePointer s_surface = spline_surfaces_[(*f_it1)->id];
      int numCoef_u = s_surface->numCoefs_u()-2+mx;
      int numCoef_v = s_surface->numCoefs_v()-2+mx;
      int numCoeffs = numCoef_u * numCoef_v;
      if (numCoeffs > 0) {
	sl2g[(*f_it1)->id].surface = glob_i;
	glob_i += numCoeffs;
      }
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
		
		vector<VolumePointer>::iterator vol_it;
		
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
			VolumePointer sv = spline_volumes_[first_vol->id]; // this edge SHOULD have at least one volume-connection
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
			VolumePointer s_volume = spline_volumes_[f->volume[0]->id];
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
			VolumePointer s_volume = spline_volumes_[i];
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
		
		vector<SurfacePointer>::iterator surf_it;
		
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
			SurfacePointer sf = spline_surfaces_[first_face->id]; // this edge should have at least one face-connection
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
			SurfacePointer ss = spline_surfaces_[i];
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


void SplineModel::getGlobalNaturalNumbering(std::vector<std::vector<int> >& num) const
{
  this->getGlobalNumbering(num);
  this->renumberNatural(num);
}


void SplineModel::getGlobalNumbering(std::vector<std::vector<int> >& num) const
{
  if (volumetric_model)
    getGlobalNumberingVolumes(num);
  else
    getGlobalNumberingSurfaces(num);
}


void SplineModel::getGlobalNumberingSurfaces(std::vector<std::vector<int> >& num) const
{
  size_t ns = spline_surfaces_.size();
  num.resize(ns);
  
  for (size_t s = 0;s < ns;s++) {
    size_t nx = spline_surfaces_[s]->numCoefs_u();
    size_t ny = spline_surfaces_[s]->numCoefs_v();
    num[s].resize(nx*ny);
    
    // Vertex numbers
    num[s][0]         = sl2g[s].vertex[0];
    num[s][nx-1]      = sl2g[s].vertex[1];
    num[s][nx*(ny-1)] = sl2g[s].vertex[2];
    num[s][nx*ny-1]   = sl2g[s].vertex[3];

    // Edge numbers
    int gnod = sl2g[s].edge[0];
    int lnod = nx;
    for (size_t i = 1;i < ny-1;i++) {
      num[s][lnod++] = gnod;
      gnod +=  sl2g[s].edge_incr[0];
    }

    gnod = sl2g[s].edge[1];
    lnod = 2*nx-1;
    for (size_t i = 1;i < ny-1;i++) {
      num[s][lnod++] = gnod;
      gnod +=  sl2g[s].edge_incr[1];
    }

    gnod = sl2g[s].edge[2];
    lnod = 1;
    for (size_t i = 1;i < nx-1;i++) {
      num[s][lnod] = gnod;
      gnod +=  sl2g[s].edge_incr[2];
      lnod += nx;
    }

    gnod = sl2g[s].edge[3];
    lnod = nx*(ny-1)+1;
    for (size_t i = 1;i < nx-1;i++) {
      num[s][lnod] = gnod;
      gnod +=  sl2g[s].edge_incr[3];
      lnod += nx;
    }

    // Face numbers
    gnod = sl2g[s].surface;
    lnod = nx+1;
    for (size_t j = 1;j < ny-1;j++, lnod++) 
      for (size_t i = 1;i < nx-1;i++, lnod++) 
	num[s][lnod] = gnod++;
  }
}


void SplineModel::getGlobalNumberingVolumes(std::vector<std::vector<int> >& num) const
{
  size_t nv = spline_volumes_.size();
  num.resize(nv);
  
  for (size_t v = 0;v < nv;v++) {
    size_t nx = spline_volumes_[v]->numCoefs(0);
    size_t ny = spline_volumes_[v]->numCoefs(1);
    size_t nz = spline_volumes_[v]->numCoefs(2);
    num[v].resize(nx*ny*nz);

    // Vertex numbers
    num[v][0]            = vl2g[v].vertex[0];
    num[v][nx-1]         = vl2g[v].vertex[1];
    num[v][nx*(ny-1)]    = vl2g[v].vertex[2];
    num[v][nx*ny-1]      = vl2g[v].vertex[3];
    size_t offset = nx*ny*(nz-1);
    num[v][offset]           = vl2g[v].vertex[4];
    num[v][offset+nx-1]      = vl2g[v].vertex[5];
    num[v][offset+nx*(ny-1)] = vl2g[v].vertex[6];
    num[v][offset+nx*ny-1]   = vl2g[v].vertex[7];

    // Edge numbers
    int gnod = vl2g[v].edge[0];
    int lnod = 1;
    for (size_t i = 1;i < nx-1;i++) {
      num[v][lnod++] = gnod;
      gnod +=  vl2g[v].edge_incr[0];
    }

    gnod = vl2g[v].edge[1];
    lnod = nx*(ny-1) + 1;
    for (size_t i = 1;i < nx-1;i++) {
      num[v][lnod++] = gnod;
      gnod +=  vl2g[v].edge_incr[1];
    }

    gnod = vl2g[v].edge[2];
    lnod = nx*ny*(nz-1) + 1;
    for (size_t i = 1;i < nx-1;i++) {
      num[v][lnod++] = gnod;
      gnod +=  vl2g[v].edge_incr[2];
    }

    gnod = vl2g[v].edge[3];
    lnod = nx*ny*(nz-1) + nx*(ny-1) + 1;
    for (size_t i = 1;i < nx-1;i++) {
      num[v][lnod++] = gnod;
      gnod +=  vl2g[v].edge_incr[3];
    }

    gnod = vl2g[v].edge[4];
    lnod = nx;
    for (size_t i = 1;i < ny-1;i++) {
      num[v][lnod++] = gnod;
      gnod +=  vl2g[v].edge_incr[4];
    }

    gnod = vl2g[v].edge[5];
    lnod = 2*nx-1;
    for (size_t i = 1;i < ny-1;i++) {
      num[v][lnod++] = gnod;
      gnod +=  vl2g[v].edge_incr[5];
    }

    gnod = vl2g[v].edge[6];
    lnod = nx*ny*(nz-1) + nx;
    for (size_t i = 1;i < ny-1;i++) {
      num[v][lnod++] = gnod;
      gnod +=  vl2g[v].edge_incr[6];
    }

    gnod = vl2g[v].edge[7];
    lnod = nx*ny*(nz-1) + 2*nx-1;
    for (size_t i = 1;i < ny-1;i++) {
      num[v][lnod++] = gnod;
      gnod +=  vl2g[v].edge_incr[7];
    }

    gnod = vl2g[v].edge[8];
    lnod = nx*ny;
    for (size_t i = 1;i < nz-1;i++) {
      num[v][lnod++] = gnod;
      gnod +=  vl2g[v].edge_incr[8];
    }

    gnod = vl2g[v].edge[9];
    lnod = nx*ny + nx-1;
    for (size_t i = 1;i < nz-1;i++) {
      num[v][lnod++] = gnod;
      gnod +=  vl2g[v].edge_incr[9];
    }

    gnod = vl2g[v].edge[10];
    lnod = nx*(2*ny-1);
    for (size_t i = 1;i < nz-1;i++) {
      num[v][lnod++] = gnod;
      gnod +=  vl2g[v].edge_incr[10];
    }

    gnod = vl2g[v].edge[11];
    lnod = 2*nx*ny-1;
    for (size_t i = 1;i < nz-1;i++) {
      num[v][lnod++] = gnod;
      gnod +=  vl2g[v].edge_incr[11];
    }

    // Faces
    gnod = vl2g[v].surface[0];
    int lnod2 = nx*(ny+1);
    for (size_t k = 1;k < nz-1;k++, lnod2 += vl2g[v].surface_incr_j[0]) {
      lnod = lnod2;
      for (size_t j = 1;j < ny-1;j++, lnod += vl2g[v].surface_incr_i[0]) 
	num[v][lnod++] = gnod;
    }

    gnod = vl2g[v].surface[1];
    lnod2 = nx*ny + 2*nx-1;;
    for (size_t k = 1;k < nz-1;k++, lnod2 += vl2g[v].surface_incr_j[1]) {
      lnod = lnod2;
      for (size_t j = 1;j < ny-1;j++, lnod += vl2g[v].surface_incr_i[1]) 
	num[v][lnod++] = gnod;
    }

    gnod = vl2g[v].surface[2];
    lnod2 = nx*ny + 1;;
    for (size_t k = 1;k < nz-1;k++, lnod2 += vl2g[v].surface_incr_j[2]) {
      lnod = lnod2;
      for (size_t i = 1;i < nx-1;i++, lnod += vl2g[v].surface_incr_i[2]) 
	num[v][lnod++] = gnod;
    }

    gnod = vl2g[v].surface[3];
    lnod2 = nx*ny + nx*(ny-1)+1;;
    for (size_t k = 1;k < nz-1;k++, lnod2 += vl2g[v].surface_incr_j[3]) {
      lnod = lnod2;
      for (size_t i = 1;i < nx-1;i++, lnod += vl2g[v].surface_incr_i[3]) 
	num[v][lnod++] = gnod;
    }

    gnod = vl2g[v].surface[4];
    lnod2 = nx + 1;;
    for (size_t j = 1;j < ny-1;j++, lnod2 += vl2g[v].surface_incr_j[4]) {
      lnod = lnod2;
      for (size_t i = 1;i < nx-1;i++, lnod += vl2g[v].surface_incr_i[4]) 
	num[v][lnod++] = gnod;
    }

    gnod = vl2g[v].surface[5];
    lnod2 = nx*ny*(nz-1)+nx+1;;
    for (size_t j = 1;j < ny-1;j++, lnod2 += vl2g[v].surface_incr_j[5]) {
      lnod = lnod2;
      for (size_t i = 1;i < nx-1;i++, lnod += vl2g[v].surface_incr_i[5]) 
	num[v][lnod++] = gnod;
    }

    // Interior nodes
    gnod = vl2g[v].volume;
    lnod = nx*ny + nx + 1;
    for (size_t k = 1;k < nz-1;k++)
      for (size_t j = 1;j < ny-1;j++)
	for (size_t i = 1;i < nx-1;i++) 
	  num[v][lnod++] = gnod++;
  }
}



void SplineModel::renumberNatural(std::vector<std::vector<int> >& num) const
{
  size_t ns = num.size();

  size_t nnod = num[num.size()-1][num[ns-1].size()-1];
  std::vector<int> gnum;
  gnum.resize(nnod,-1);
  
  int gnod = 0;
  for (size_t s = 0;s < ns;s++) 
    for (size_t i = 0;i < num[s].size();i++) {
      size_t node = num[s][i];
      if (gnum[node] > -1)
	num[s][i] = gnum[node];
      else 
	num[s][i] = gnum[node] = gnod++;
    }
}


void SplineModel::writeModelXMLProperties(std::ostream &os) const {
	set<Volume*>::iterator v_it;
	set<Face*>::iterator   f_it;
	set<Line*>::iterator   l_it;
	set<Vertex*>::iterator c_it; //corner iterator

	string activeString = "";
	bool openSetBlock = false;

	/***************     WRITE IT IN THE INPUT-FORM, NOT REALLY USEFULL FOR PRODUCTION    ***********************/
	os << "  <boundaryconditions>" << endl;

	if(volumetric_model) {
		for(f_it=topology->face_begin(); f_it != topology->face_end(); f_it++) {
			Face *f = *f_it;
			if(f->bc_code != NULL) {
				if(activeString.compare(f->bc_code) != 0) { 
					if(openSetBlock) os << "    </set>" << endl;
					activeString = f->bc_code;
					openSetBlock = true;
					os << "    <set name=\"" << activeString << "\" type=\"face\">" << endl;
				}
				os << "      <item patch=\"" << f->volume[0]->id << "\"> " << f->face[0] << " </item>" << endl;
			}
		}
		for(l_it=topology->line_begin(); l_it != topology->line_end(); l_it++) {
			Line *l = *l_it;
			vector<int> numb, parDir, parStep;
			(*l->volume.begin())->getEdgeEnumeration(l, numb, parDir, parStep);
			if(l->bc_code != NULL) {
				if(activeString.compare(l->bc_code) != 0) { 
					if(openSetBlock) os << "    </set>" << endl;
					activeString = l->bc_code;
					openSetBlock = true;
					os << "    <set name=\"" << activeString << "\" type=\"edge\">" << endl;
				}
				os << "      <item patch=\"" << (*l->volume.begin())->id << "\"> " << numb[0] << " </item>" << endl;
			}
		}
		for(c_it=topology->vertex_begin(); c_it != topology->vertex_end(); c_it++) {
			Vertex *c = *c_it;
			vector<int> corner_id = (*c->volume.begin())->getVertexEnumeration(c);
			if(c->bc_code != NULL) {
				if(activeString.compare(c->bc_code) != 0) { 
					if(openSetBlock) os << "    </set>" << endl;
					activeString = c->bc_code;
					openSetBlock = true;
					os << "    <set name=\"" << activeString << "\" type=\"vertex\">" << endl;
				}
				os << "      <item patch=\"" << (*c->volume.begin())->id << "\"> " << corner_id[0] << " </item>" << endl;
			}
		}
	} else { // SURFACE model
		for(l_it=topology->line_begin(); l_it != topology->line_end(); l_it++) {
			Line *l = *l_it;
			vector<int> numb, parDir, parStep;
			(*l->face.begin())->getEdgeEnumeration(l, numb, parDir, parStep);
			if(l->bc_code != NULL) {
				if(activeString.compare(l->bc_code) != 0) { 
					if(openSetBlock) os << "    </set>" << endl;
					activeString = l->bc_code;
					openSetBlock = true;
					os << "    <set name=\"" << activeString << "\" type=\"edge\">" << endl;
				}
				os << "      <item patch=\"" << (*l->volume.begin())->id << "\"> " << numb[0] << " </item>" << endl;
			}
		}
		for(c_it=topology->vertex_begin(); c_it != topology->vertex_end(); c_it++) {
			Vertex *c = *c_it;
			vector<int> corner_id = (*c->face.begin())->getVertexEnumeration(c);
			if(c->bc_code != NULL) {
				if(activeString.compare(c->bc_code) != 0) { 
					if(openSetBlock) os << "    </set>" << endl;
					activeString = c->bc_code;
					openSetBlock = true;
					os << "    <set name=\"" << activeString << "\" type=\"vertex\">" << endl;
				}
				os << "      <item patch=\"" << (*c->volume.begin())->id << "\"> " << corner_id[0] << " </item>" << endl;
			}
		}
	}
	if(openSetBlock) os << "    </set>" << endl;

	os << "  </boundaryconditions>" << endl;
}

void SplineModel::writeModelProperties(std::ostream &os) const {
	set<Volume*>::iterator v_it;
	set<Face*>::iterator   f_it;
	set<Line*>::iterator   l_it;
	set<Vertex*>::iterator c_it; //corner iterator

	/***************     WRITE IT IN THE INPUT-FORM, NOT REALLY USEFULL FOR PRODUCTION    ***********************/

	if(volumetric_model) {
		for(v_it=topology->volume_begin(); v_it != topology->volume_end(); v_it++) {
			Volume* v = *v_it;
			if(v->material_code != NULL)
				os << v->material_code << "\n" << v->id << " 3 " << endl;
		}
		for(f_it=topology->face_begin(); f_it != topology->face_end(); f_it++) {
			Face *f = *f_it;
			if(f->bc_code != NULL) 
				os << f->bc_code << "\n" << f->volume[0]->id << " 2 " << f->face[0] << endl;
		}
		for(l_it=topology->line_begin(); l_it != topology->line_end(); l_it++) {
			Line *l = *l_it;
			vector<int> numb, parDir, parStep;
			(*l->volume.begin())->getEdgeEnumeration(l, numb, parDir, parStep);
			if(l->bc_code != NULL)
				os << l->bc_code << "\n" << (*l->volume.begin())->id << " 1 " << numb[0] << endl;
		}
		for(c_it=topology->vertex_begin(); c_it != topology->vertex_end(); c_it++) {
			Vertex *c = *c_it;
			vector<int> corner_id = (*c->volume.begin())->getVertexEnumeration(c);
			if(c->bc_code != NULL)
				os << c->bc_code << "\n" << (*c->volume.begin())->id << " 0 " << corner_id[0] << endl;
		}
	} else { // SURFACE model
		for(f_it=topology->face_begin(); f_it != topology->face_end(); f_it++) {
			Face *f = *f_it;
			if(f->bc_code != 0) 
				os << f->bc_code << "\n" << f->id << " 2 " << f->face[0] << endl;
		}
		for(l_it=topology->line_begin(); l_it != topology->line_end(); l_it++) {
			Line *l = *l_it;
			vector<int> numb, parDir, parStep;
			(*l->face.begin())->getEdgeEnumeration(l, numb, parDir, parStep);
			if(l->bc_code != 0)
				os << l->bc_code << "\n" << (*l->face.begin())->id << " 1 " << numb[0] << endl;
		}
		for(c_it=topology->vertex_begin(); c_it != topology->vertex_end(); c_it++) {
			Vertex *c = *c_it;
			vector<int> corner_id = (*c->face.begin())->getVertexEnumeration(c);
			if(c->bc_code != 0)
				os << c->bc_code << "\n" << (*c->face.begin())->id << " 0 " << corner_id[0] << endl;
		}
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

void SplineModel::readSplines(std::istream &is, bool buildTopology) {
	Go::ObjectHeader head;
	while(!is.eof()) {
		head.read(is);
		if(head.classType() == Go::Class_SplineVolume) {
			VolumePointer v(new Go::SplineVolume());
			v->read(is);
			spline_volumes_.push_back(v);
			topology->addPatch(v);
			volumetric_model = true;
		} else if(head.classType() == Go::Class_SplineSurface) {
			SurfacePointer s(new Go::SplineSurface());
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
	if(buildTopology)
		topology->buildTopology();
}

void SplineModel::readModelProperties(std::istream &is) {
	if(!is.good())
		// some error code 
		;
	string str_line, primitive, code;
	int volId, locId, inclusive;
	while( getline(is, code) ) {
		getline(is, str_line);
		istringstream ss(str_line);
		ss >> primitive;
		if(primitive.compare("Volume") == 0 ||
		   primitive.compare("volume") == 0 )  {
			ss >> volId;
			if( !(ss >> inclusive) )
				inclusive = 0;
			addVolumePropertyCode(volId, code.c_str(), inclusive);
		} else if(primitive.compare("Face") == 0 ||
		          primitive.compare("face") == 0 )  {
			ss >> volId;
			ss >> locId;
			if( !(ss >> inclusive) )
				inclusive = 1;
			addFacePropertyCode(volId, locId, code.c_str(), inclusive);
		} else if(primitive.compare("Line") == 0 ||
		          primitive.compare("line") == 0 )  {
			ss >> volId;
			ss >> locId;
			if( !(ss >> inclusive) )
				inclusive = 1;
			addLinePropertyCode(volId, locId, code.c_str(), inclusive);
		} else if(primitive.compare("Vertex") == 0 ||
		          primitive.compare("vertex") == 0 )  {
			ss >> volId;
			ss >> locId;
			addVertexPropertyCode(volId, locId, code.c_str());
		} else if(primitive[0] == '/' && primitive[1] == '/') { // comment - continue without doing anything
			;
		} else {
			cerr << "Syntax error near: \"" << str_line << "\"\n";
			exit(1);
		}
	}
}

