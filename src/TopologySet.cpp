/**********************************************************************************//**
 * \file TopologySet.cpp
 *
 * \author Kjetil A. Johannessen
 *
 * \date July 2010
 *
 *************************************************************************************/
#include "TopologySet.h"
#include "primitives.h"

using namespace std;

/**********************************************************************************//**
 * \brief Constructor
 * \param tol            Tolerance used when checking for control point equality (measured in euclidean distance)
 *
 *************************************************************************************/
TopologySet::TopologySet(double tol) {
	this->tol = tol;
	volumetric_model = false;
	surface_model    = false;
}

/**********************************************************************************//**
 * \brief Constructor
 * \param spline_surfaces All spline surfaces to be considered part of this model
 * \param tol             Tolerance used when checking for control point equality (measured in euclidean distance)
 *
 * Note that the constructor is NOT creating the topology. It should usually immediately be followed by a call to
 * buildTopology().
 *************************************************************************************/
TopologySet::TopologySet(std::vector<SurfacePointer> &spline_surfaces, double tol) {
	this->tol = tol;
	spline_surfaces_ = spline_surfaces;
	volumetric_model = false;
	surface_model    = true;
}

/**********************************************************************************//**
 * \brief Constructor
 * \param spline_volumes All spline volumes to be considered part of this model
 * \param tol            Tolerance used when checking for control point equality (measured in euclidean distance)
 *
 * Note that the constructor is NOT creating the topology. It should usually immediately be followed by a call to
 * buildTopology().
 *************************************************************************************/
TopologySet::TopologySet(std::vector<VolumePointer> &spline_volumes, double tol) {
	this->tol = tol;
	spline_volumes_ = spline_volumes;
	volumetric_model = true;
	surface_model    = false;
}

TopologySet::~TopologySet() {
	set<Vertex*>::iterator	v_it;
	set<Line*>::iterator    l_it;
	FaceSet::iterator    f_it;

	for(v_it = vertex_begin(); v_it != vertex_end(); v_it++)
		delete *v_it;
	for(l_it = line_begin(); l_it != line_end(); l_it++)
		delete *l_it;
	for(f_it = face_begin(); f_it != face_end(); f_it++)
		delete *f_it;
}

/**********************************************************************************//**
 * \brief sets neighbouring tolerance
 * \param tol absolute tolerance given in euclidean distance measured for each controlpoint
 *
 *************************************************************************************/
void TopologySet::setTolerance(double tol) {
	this->tol = tol;
	buildTopology();
}


/**********************************************************************************//**
 * \brief sets periodicitty in a given direction
 * \param dir periodic direction (0=X, 1=Y, 2=Z) 
 *
 *************************************************************************************/
void TopologySet::setPeriodic(int dir) {
	if (dir > 2) {
		cerr << "Topology::setPeriodic: dir > 2 not possible" << endl;    
		return;
	}
	else if ((dir > 1) && !volumetric_model) {
		cerr << "Topology::setPeriodic: dir > 1 not possible for surface model" << endl;   
		return;
	}

	// Find min/max coordinate for the diven direction
	double min = 1.0e20;
	double max = -1.0e20;

	std::vector<VolumePointer>::iterator vit;
	std::vector<SurfacePointer>::iterator sit;

	std::vector<double>::iterator start;
	std::vector<double>::iterator end;
	std::vector<double>::iterator coef; 

	// Find min and max coordinates
	if(volumetric_model)
		for (vit=spline_volumes_.begin();vit != spline_volumes_.end();vit++) {
			int dim = (*vit)->dimension();
			
			if ((*vit)->rational()) {
				dim++;

				start = (*vit)->rcoefs_begin();
				end   = (*vit)->rcoefs_end();
			}
			else {
				start = (*vit)->coefs_begin();
				end   = (*vit)->coefs_end();
			}

			for (coef = start;coef != end;coef += dim) {
				double val = *(coef+dir);
				
				if (val < min)      
					min = val;
				else if (val > max) 
					max = val;
			}
		}
	else  // Surface spline
		for (sit=spline_surfaces_.begin();sit != spline_surfaces_.end();sit++) {
			int dim = (*sit)->dimension();
			
			if ((*sit)->rational()) {
				dim++;
				
				start = (*sit)->rcoefs_begin();
				end   = (*sit)->rcoefs_end();
			}
			else {
				start = (*sit)->coefs_begin();
				end   = (*sit)->coefs_end();
			}
			
			for (coef = start;coef != end;coef += dim) {
				double val = *(coef+dir);
				
				if (val < min)      
					min = val;
				else if (val > max) 
					max = val;
			}
		}    

	// Set max coordinate equal to min coordinate
		if(volumetric_model) 
		for (vit=spline_volumes_.begin();vit != spline_volumes_.end();vit++) {
			int dim = (*vit)->dimension();

			if ((*vit)->rational()) {
				dim++;

				start = (*vit)->rcoefs_begin();
				end   = (*vit)->rcoefs_end();
			}
			else {
				start = (*vit)->coefs_begin();
				end   = (*vit)->coefs_end();
			}

			for (coef = start;coef != end;coef += dim) {
				double val = *(coef+dir);
				
				if (fabs(max-val) < tol)      
					*(coef+dir) = min;
			}
		}
		else  // Surface spline
			for (sit=spline_surfaces_.begin();sit != spline_surfaces_.end();sit++) {
				int dim = (*sit)->dimension();
				
				if ((*sit)->rational()) {
					dim++;
					
					start = (*sit)->rcoefs_begin();
					end   = (*sit)->rcoefs_end();
				}
				else {
					start = (*sit)->coefs_begin();
					end   = (*sit)->coefs_end();
				}
				
				for (coef = start;coef != end;coef += dim) {
					double val = *(coef+dir);
					
					if (fabs(max-val) < tol)      
					  *(coef+dir) = min;
				}
			}       
}


/**********************************************************************************//**
 * \brief adds another patch to the model
 * \param surf surface patch to be added to the total model
 *
 * Nice method to use if you are sequentially building up your topology model. Note that
 * this causes error if called on a volumetric model.
 *************************************************************************************/
void TopologySet::addPatch(SurfacePointer surf) {
	if(volumetric_model) {
		cerr << "Trying to add surface patch to a volumetric topology model\n";
		exit(4251);
	} else {
		spline_surfaces_.push_back(surf);
		surface_model = true;
	}
}

/**********************************************************************************//**
 * \brief adds another patch to the model
 * \param vol volume patch to be added to the total model
 *
 * Nice method to use if you are sequentially building up your topology model. Note that
 * this causes error if called on a surface model.
 *************************************************************************************/
void TopologySet::addPatch(VolumePointer vol) {
	if(surface_model) {
		cerr << "Trying to add volume patch to a surface topology model\n";
		exit(4252);
	} else {
		spline_volumes_.push_back(vol);
		volumetric_model = true;
	}
}

/*! \brief build the topology given by all vertex, line, face and volume relations */
void TopologySet::buildTopology(std::vector<bool>* periodic) {
	volume_.clear();
	face_.clear();
	line_.clear();
	vertex_.clear();
	
	// Handle periodicity
	if (periodic)
	  for (int i = 0;i < periodic->size();i++)
	    if ((*periodic)[i])
	      this->setPeriodic(i);

	if(volumetric_model) {
		for(uint i=0; i<spline_volumes_.size(); i++) {
			bool rat = spline_volumes_[i]->rational();
			int dim  = spline_volumes_[i]->dimension();
			int n1   = spline_volumes_[i]->numCoefs(0);
			int n2   = spline_volumes_[i]->numCoefs(1);
			int n3   = spline_volumes_[i]->numCoefs(2);
			vector<double>::iterator coef = (rat) ? spline_volumes_[i]->rcoefs_begin() : spline_volumes_[i]->coefs_begin();
			// cout << "Volume " << i << endl;

			// create the volume object (should never be degenerated to face, line or point)
			Volume *vol = new Volume(i);
			addVolume(vol);

			// add all corner vertices to the appropriate places
			int corner_index = 0;
			for(int w=0; w<2; w++) {
				for(int v=0; v<2; v++) {
					for(int u=0; u<2; u++) {
						int coefNmb = u*(n1-1) + v*n1*(n2-1) + w*n1*n2*(n3-1);
						Go::Point p(coef+(coefNmb*(dim+rat)), coef+((coefNmb+1)*(dim+rat)));

						if(rat)
							for(int aa=0; aa<dim+rat; aa++)
								p[aa] /= p[dim];
						// cout << "  Corner " << corner_index++ << ": " << p << endl;
						Vertex *vert = new Vertex();
						vert->cp = p;
						vert->volume.insert(vol);

						vert = addVertex(vert);
						vol->corner[corner_index++] = vert;
					}
				}
			}

			// add all lines to the right places
			int step = 1;
			int lineCount = 0;
			for(int parDir=0; parDir<3; parDir++) {
				if(parDir==1)      step *= n1;
				else if(parDir==2) step *= n2;
				for(int u2=0; u2<2; u2++) {
					for(int u1=0; u1<2; u1++) {
						Line *line = new Line();
						int start = 0;
						int v1_i  = -1;
						int v2_i  = -1;
						if(parDir==0)      start =          + u1*n1*(n2-1) + u2*n1*n2*(n3-1);
						else if(parDir==1) start = u1*(n1-1)+              + u2*n1*n2*(n3-1);
						else if(parDir==2) start = u1*(n1-1)+ u2*n1*(n2-1)                  ;

						if(parDir==0)      v1_i  =      2*u1 + 4*u2;
						else if(parDir==1) v1_i  = u1        + 4*u2;
						else if(parDir==2) v1_i  = u1 + 2*u2       ;

						if(parDir==0)      v2_i  =  1 + 2*u1 + 4*u2;
						else if(parDir==1) v2_i  = u1 +   2  + 4*u2;
						else if(parDir==2) v2_i  = u1 + 2*u2 +   4 ;

						line->v1 = vol->corner[v1_i];
						line->v2 = vol->corner[v2_i];

						bool degen = true;
						// cout << "Volume " << vol->id << " line " << lineCount << endl;
						// cout << "Line: " << *line << endl;
						for(int cpCount=0; cpCount<spline_volumes_[i]->numCoefs(parDir); cpCount++) {
							Go::Point p(coef+start*(dim+rat), coef+(start+1)*(dim+rat));
							if(rat)
								for(int aa=0; aa<dim+rat; aa++)
									p[aa] /= p[dim];
							// cout << " (" << p << ") ";
							if( line->cp.size()>0 && line->cp.back().dist(p)>tol ) {
								// cout << " DEGEN FALSE ";
								degen = false;
							}
							line->cp.push_back(p);
							start += step;
						}
						// cout << endl;
						line->degen = degen;
						line->volume.insert(vol);
						line = addLine(line);
						vol->line[lineCount] = line;
						lineCount++;
					}
				}
			}

			// add all faces to the right places
			int surfCount = 0;
			for(int parDir=0; parDir<3; parDir++) {
				for(int minmax=0; minmax<2; minmax++) { // parMin or parMax
					int coefNmb = -1;
					bool degen1  = true;
					bool degen2  = true;
					if(parDir==0)      coefNmb = minmax *(n1-1)                ;
					else if(parDir==1) coefNmb = minmax *  n1  * (n2-1)        ;
					else if(parDir==2) coefNmb = minmax *  n1  *   n2   *(n3-1);
					int du    = (parDir==0) ? n1 :1;     // step length walking in local u-dir
					int dv    = (parDir==2) ? n1 :n1*n2; // step length walking in local v-dir
					int u_siz = (parDir==0) ? n2 : n1;
					int v_siz = (parDir==2) ? n2 : n3;
					
					// cout << "Volume #" << i << " - surace " << surfCount << endl;

					// initialize Face
					Face *f  = new Face(u_siz, v_siz);
					f->volume.push_back(vol);
					f->face.push_back(surfCount);

					for(int v=0; v<v_siz; v++) {
						for(int u=0; u<u_siz; u++) {
							Go:: Point p(coef+coefNmb*(dim+rat), coef+(coefNmb+1)*(dim+rat));
							if(rat)
								for(int aa=0; aa<dim+rat; aa++)
									p[aa] /= p[dim];
							f->cp[u][v] = p;
							if( u>0 && f->cp[u-1][v].dist(p) > tol )
								degen1 = false;
							if( v>0 && f->cp[u][v-1].dist(p) > tol )
								degen2 = false;
							
							// cout << "   (" << p << ")" << endl;
							coefNmb += du;
						}
						coefNmb -= u_siz*du;
						coefNmb += dv;
					}
					f->degen1 = degen1;
					f->degen2 = degen2;
					f = addFace(f);
					vol->face[surfCount] = f;
					surfCount++;
				}
			}
		}
	} else if(surface_model) {
		for(uint i=0; i<spline_surfaces_.size(); i++) {
			bool rat = spline_surfaces_[i]->rational();
			int dim  = spline_surfaces_[i]->dimension();
			int n1   = spline_surfaces_[i]->numCoefs_u();
			int n2   = spline_surfaces_[i]->numCoefs_v();
			vector<double>::iterator coef = (rat) ? spline_surfaces_[i]->rcoefs_begin() : spline_surfaces_[i]->coefs_begin();
			// cout << "Face " << i << endl;

			// create the face object (should never be degenerated to line or point)
			Face *face = new Face(i);
			addFace(face);

			// add all corner vertices to the appropriate places
			int corner_index = 0;
			for(int v=0; v<2; v++) {
				for(int u=0; u<2; u++) {
					int coefNmb = u*(n1-1) + v*n1*(n2-1);
					Go::Point p(coef+(coefNmb*(dim+rat)), coef+((coefNmb+1)*(dim+rat)));

					if(rat)
						for(int aa=0; aa<dim+rat; aa++)
							p[aa] /= p[dim];
					// cout << "  Corner " << corner_index << ": " << p << endl;
					Vertex *vert = new Vertex();
					vert->cp = p;
					vert->face.insert(face);

					vert = addVertex(vert);
					face->corner[corner_index++] = vert;
				}
			}

			// add all lines to the right places
			int step = n1;
			int lineCount = 0;
			for(int parDir=2; parDir-->0; ) {
				if(parDir==0)      step /= n1;
				for(int u1=0; u1<2; u1++) {
					Line *line = new Line();
					int start = 0;
					int v1_i  = -1;
					int v2_i  = -1;
					if(parDir==0)      start =          + u1*n1*(n2-1);
					else if(parDir==1) start = u1*(n1-1)              ;

					if(parDir==0)      v1_i  =      2*u1;
					else if(parDir==1) v1_i  = u1       ;

					if(parDir==0)      v2_i  =  1 + 2*u1;
					else if(parDir==1) v2_i  = u1 +   2 ;

					line->v1 = face->corner[v1_i];
					line->v2 = face->corner[v2_i];

					bool degen = true;
					// cout << "Face " << face->id << " line " << lineCount << endl;
					// cout << "Line: " << *line << endl;
					for(int cpCount=0; (parDir==0 && cpCount<spline_surfaces_[i]->numCoefs_u()) ||
					                   (parDir==1 && cpCount<spline_surfaces_[i]->numCoefs_v())   ; cpCount++) {
						Go::Point p(coef+start*(dim+rat), coef+(start+1)*(dim+rat));
						if(rat)
							for(int aa=0; aa<dim+rat; aa++)
								p[aa] /= p[dim];
						// cout << " (" << p << ") ";
						if( line->cp.size()>0 && line->cp.back().dist(p)>tol ) {
							// cout << " DEGEN FALSE ";
							degen = false;
						}
						line->cp.push_back(p);
						start += step;
					}
					// cout << endl;
					line->degen = degen;
					line->face.insert(face);
					line = addLine(line);
					face->line[lineCount] = line;
					lineCount++;
				}
			}
		}
	}
}

/**********************************************************************************//**
 * \brief try to add a vertex to the set (keeping uniqueness)
 * \param v Vertex to be added
 * \return either *v or the Vertex already contained in the set and being equal to *v
 *************************************************************************************/
Vertex* TopologySet::addVertex(Vertex *v) {
	set<Vertex*>::iterator it; 
	for(it=vertex_.begin(); it != vertex_.end(); it++) {
		if( (*it)->cp.dist(v->cp) < tol) {
			if(volumetric_model)
				(*it)->volume.insert(v->volume.begin(), v->volume.end());
			else
				(*it)->face.insert(v->face.begin(), v->face.end());
			delete v;
			return *it;
		}
	}
	vertex_.insert(v);
	return v;
}

/**********************************************************************************//**
 * \brief try to add a line to the set (keeping uniqueness)
 * \param l Line to be added
 * \return either *l or the Line already contained in the set and being equal to *l
 *************************************************************************************/
Line* TopologySet::addLine(Line* l) {
	set<Line*>::iterator it; 
	for(it=line_.begin(); it != line_.end(); it++) {
		if( l->equals(*it, tol) ) {
			if(volumetric_model)
				(*it)->volume.insert(l->volume.begin(), l->volume.end());
			else
				(*it)->face.insert(l->face.begin(), l->face.end());
			delete l;
			return *it;
		}
	}
	line_.insert(l);
	return l;
}

/**********************************************************************************//**
 * \brief try to add a face to the set (keeping uniqueness)
 * \param f Face to be added
 * \return either *f or the Face already contained in the set and being equal to *f
 *
 * In the case of surface models, uniqueness or degeneracy is not tested (*f always
 * returned)
 *************************************************************************************/
Face* TopologySet::addFace(Face* f) {
	if(surface_model) {
	  face_.insert(f);
	  return f;
	}
	FaceSet::iterator it;
	for(it=face_.begin(); it != face_.end(); it++) {
		// cout << "Testing equality " << (*it)->v1->id << "(" << (*it)->face1 << ") - " << f->v1->id << "(" << f->face1 << ")\n";
		if( (*it)->equals(f, tol) ) {
			// cout << "YAY - face equals!!!\n";
			// at this point we KNOW that *f only has one neigbhouring volume
			(*it)->volume.push_back(*f->volume.begin());
			(*it)->face.push_back(*f->face.begin());
			delete f;
			return *it;
		}
	}
	face_.insert(f);
	return f;
}

/**********************************************************************************//**
 * \brief try to add a volume to the set (uniqueness not tested)
 * \param v volume to be added
 *************************************************************************************/
void TopologySet::addVolume(Volume* v) {
	volume_.insert(v);
}

/*! \brief fetch all vertices on the boundary of the model */
std::set<Vertex*> TopologySet::getBoundaryVertices() {
	set<Vertex*> vert;
	set<Line*>   line;
	FaceSet   face;
	getBoundaries(vert, line, face);
	return vert;
}

/*! \brief fetch all lines on the boundary of the model */
std::set<Line*> TopologySet::getBoundaryLines() {
	set<Line*>   line;
	if(volumetric_model) {
		set<Vertex*> vert;
		FaceSet      face;
		getBoundaries(vert, line, face);
	} else if(surface_model) {
		set<Line*>   line;
		for(Line* l : line_) {
			if(!l->degen && l->face.size() == 1) { 
				vector<int> numb,dir,step;
				Face* f = *l->face.begin();
				f->getEdgeEnumeration(l, numb,dir,step);
				if(numb.size() == 1) // don't inlcude singlepatch periodic lines
					line.insert(l);
			}
			   
		}
	}
	return line;
}

/*! \brief fetch all faces on the boundary of the model */
FaceSet TopologySet::getBoundaryFaces() {
	FaceSet results;
	FaceSet::iterator it;
	if(volumetric_model) {
		for(it=face_.begin(); it != face_.end(); it++) 
			if((*it)->volume.size() == 1) // faces neighbouring only one volume are on edges
				if(! (*it)->isDegen())   // skip completely degenerate faces
					if((*it)->volume[0]->getSurfaceEnumeration(*it).size() == 1) // periodic patch faces are not part of boundary 
						results.insert(*it);
	} else if(surface_model) {
		results = face_;
	}
	return results; // empty set
}

/**********************************************************************************//**
 * \brief fetch all boundary primitives
 * \param [out] vertices Boundary vertices
 * \param [out] lines    Boundary lines
 * \param [out] faces    Boundary faces
 *
 * An optimized version to get all boundaries. This is faster than sequentially calling
 * getBoundaryVertices(),getBoundaryLines() and getBoundaryFaces().
 *
 *************************************************************************************/
void TopologySet::getBoundaries(std::set<Vertex*>& vertices, std::set<Line*>& lines, FaceSet& faces) {
	faces.clear();
	lines.clear();
	vertices.clear();
	if(volumetric_model) {
		faces = getBoundaryFaces();
		FaceSet::iterator it;
		for(it=faces.begin(); it != faces.end(); it++) {
			vector<int> lineNumb = Line::getLineEnumeration( *(*it)->face.begin() );
			Volume *firstVol = *(*it)->volume.begin();
			for(int i=0; i<4; i++) {
				lines.insert( firstVol->line[lineNumb[i]] );
				vertices.insert( firstVol->line[lineNumb[i]]->v1 ); // yes, this will insert all vertices (at least) twice, but it is preffered
				vertices.insert( firstVol->line[lineNumb[i]]->v2 ); // to the alternative which is making more functions 
			}
		}
	} else if(surface_model) {
		lines = getBoundaryLines();
		set<Line*>::iterator it;
	}
}

VolSet::iterator TopologySet::volume_begin() {
	return volume_.begin();
}

VolSet::iterator TopologySet::volume_end() {
	return volume_.end();
}

FaceSet::iterator TopologySet::face_begin() {
	return face_.begin();
}

FaceSet::iterator TopologySet::face_end() {
	return face_.end();
}

set<Line*>::iterator TopologySet::line_begin() {
	return line_.begin();
}

set<Line*>::iterator TopologySet::line_end() {
	return line_.end();
}

set<Vertex*>::iterator TopologySet::vertex_begin() {
	return vertex_.begin();
}

set<Vertex*>::iterator TopologySet::vertex_end() {
	return vertex_.end();
}

/*! \brief Number of unique vertices in the model */
int TopologySet::numbVertices() const {
	return vertex_.size();
}

/*! \brief Number of unique lines in the model */
int TopologySet::numbLines() const {
	return line_.size();
}

/*! \brief Number of unique lines which are not degenerated to points */
int TopologySet::numbNonDegenLines() const {
	int results = 0;
	set<Line*>::iterator it;
	for(it=line_.begin(); it != line_.end(); it++)
		if( !(*it)->degen )
			results++;
	return results;
}

/*! \brief Number of unique faces in the model */
int TopologySet::numbFaces() const {
	return face_.size();
}

/*! \brief Number of unique faces which are not degenerated to lines or points */
int TopologySet::numbNonDegenFaces() const {
	int results = 0;
	FaceSet::iterator it;
	for(it=face_.begin(); it != face_.end(); it++)
		if( !(*it)->isDegen() )
			results++;
	return results;
}

/*! \brief Number of unique volumes in the model */
int TopologySet::numbVolumes() const {
	return volume_.size();
}

/*! \brief get one specific face */
Face* TopologySet::getFace(int id) {
	FaceSet::iterator it;
	for(it=face_.begin(); it != face_.end(); it++)
		if( (*it)->id == id)
			return *it;
	return NULL;
}

/*! \brief get one specific volume */
Volume* TopologySet::getVolume(int id) {
	VolSet::iterator it;
	for(it=volume_.begin(); it != volume_.end(); it++)
		if( (*it)->id == id)
			return *it;
	return NULL;
}

void TopologySet::write(ostream &out) const {
	VolSet::iterator it;
	vector<Volume*>::iterator it2;
	for(it=volume_.begin(); it != volume_.end(); it++) {
		out << "Volume #" << (*it)->id << ":  ";
		for(int i=0; i<6; i++) {
			int written=0;
			if((*it)->face[i]) {
				for(int j=0; j<(*it)->face[i]->volume.size(); j++) { 
					if((*it)->face[i]->volume[j]->id != (*it)->id) {
						if(written++ > 0)
							out << "&";
						out << (*it)->face[i]->volume[j]->id ;
					}
				}
			}
			if(written == 0)
				out << "-";
			out << " ";
		}
		out << endl;
	}
	
}

void TopologySet::writeXML(ostream &out) const {
        out << "  <topology>" << endl;
        if (volumetric_model) {
                FaceSet::iterator it;
                vector<int> faceId;
                vector<int> allFaceId;
                vector<int> allVolId;
                for(it=face_.begin(); it != face_.end(); it++) {
                        int masterVol  = (**it).volume[0]->id;
                        allFaceId      = (**it).volume[0]->getSurfaceEnumeration(*it);
                        int masterFace = allFaceId[0];
                        for(int i=1; i<allFaceId.size(); i++) { // if this happens we have a periodic geometry
                                out << "    <connection master=\"" << masterVol +1 << "\""
                                    << " mface=\"" << masterFace                +1 << "\""
                                    << " slave=\"" << masterVol                 +1 << "\""
                                    << " sface=\"" << allFaceId[i]              +1 << "\"";
                                int orient = (**it).uv_flip[i]*4 + (**it).u_reverse[i]*2 + (**it).v_reverse[i]*1; // forms a 3bit mask of local orientation
                                if(orient != 0)
                                        out << " orient=\"" << orient << "\"";
                                out << "/>\n";
                        }
                        for(int i=1; i<(**it).volume.size(); i++) {
                                if((**it).volume[i]->id == masterVol) continue;
                                allFaceId      = (**it).volume[i]->getSurfaceEnumeration(*it);
                                int slaveVol   = (**it).volume[i]->id;
                                for(int j=0; j<allFaceId.size(); j++) { // if this happens we have a periodic geometry
                                        out << "    <connection master=\"" << masterVol +1<< "\""
                                        << " mface=\"" << masterFace                +1 << "\""
                                        << " slave=\"" << slaveVol                  +1 << "\""
                                        << " sface=\"" << allFaceId[j]              +1 << "\"";
                                        int orient = (**it).uv_flip[i]*4 + (**it).u_reverse[i]*2 + (**it).v_reverse[i]*1; // forms a 3bit mask of local orientation
                                        if(orient != 0)
                                                out << " orient=\"" << orient << "\"";
                                        out << "/>\n";
                                }
                        }
                }
        } else if (surface_model) {
                for(auto& it : line_) {
                        int masterVol  = (*it->face.begin())->id;
                        std::vector<int> dummy;
                        std::vector<int> allEdgeId;
                        (*it->face.begin())->getEdgeEnumeration(const_cast<Line*>(it), allEdgeId, dummy, dummy);
                        int masterEdge = allEdgeId[0];
                        for(int i=1; i<allEdgeId.size(); i++) { // if this happens we have a periodic geometry
                                out << "    <connection master=\"" << masterVol +1 << "\""
                                    << " medge=\"" << masterEdge                +1 << "\""
                                    << " slave=\"" << masterVol                 +1 << "\""
                                    << " sedge=\"" << allEdgeId[i]              +1 << "\"";
                                int orient = 0;
                                if(orient != 0)
                                        out << " orient=\"" << orient << "\"";
                                out << "/>\n";
                        }
                        for(auto& face : it->face) {
                                int slaveVol   = face->id;
                                if(slaveVol == masterVol)
                                  continue;
                                std::vector<int> dummy, parDir;
                                face->getEdgeEnumeration(it, allEdgeId, dummy, parDir);

                                for(int j=0; j<allEdgeId.size(); j++) { // if this happens we have a periodic geometry
                                        out << "    <connection master=\"" << masterVol +1<< "\""
                                        << " medge=\"" << masterEdge                +1 << "\""
                                        << " slave=\"" << slaveVol                  +1 << "\""
                                        << " sedge=\"" << allEdgeId[j]              +1 << "\"";
                                        if(parDir[j] < 0)
                                                out << " reverse=\"true\"";
                                        out << "/>\n";
                                }
                        }
                }

        }
        out << "  </topology>" << endl;
}

