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
using namespace Go;
using boost::shared_ptr;

/**********************************************************************************//**
 * \brief Constructor
 * \param tol            Tolerance used when checking for control point equality (measured in euclidean distance)
 *
 *************************************************************************************/
TopologySet::TopologySet(double tol) {
	this->tol = tol;
}

/**********************************************************************************//**
 * \brief Constructor
 * \param spline_volumes All spline volumes to be considered part of this model
 * \param tol            Tolerance used when checking for control point equality (measured in euclidean distance)
 *
 * Note that the constructor is NOT creating the topology. It should usually immediately be followed by a call to
 * buildTopology().
 *************************************************************************************/
TopologySet::TopologySet(std::vector<boost::shared_ptr<Go::SplineVolume> > &spline_volumes, double tol) {
	this->tol = tol;
	spline_volumes_ = spline_volumes;
}

TopologySet::~TopologySet() {
	set<Vertex*>::iterator	v_it;
	set<Line*>::iterator    l_it;
	set<Face*>::iterator    f_it;
	for(v_it = vertex_begin(); v_it != vertex_end(); v_it++)
		delete *v_it;
	for(l_it = line_begin(); l_it != line_end(); l_it++)
		delete *l_it;
	for(f_it = face_begin(); f_it != face_end(); f_it++)
		delete *f_it;
}

void TopologySet::addVolume(boost::shared_ptr<Go::SplineVolume> vol) {
	spline_volumes_.push_back(vol);
}

/*! \brief build the topology given by all vertex, line, face and volume relations */
void TopologySet::buildTopology() {
	volume_.clear();
	face_.clear();
	line_.clear();
	vertex_.clear();
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
					Point p(coef+(coefNmb*(dim+rat)), coef+((coefNmb+1)*(dim+rat)));

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
						Point p(coef+start*(dim+rat), coef+(start+1)*(dim+rat));
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
					/*
					if(degen) {
						// cout << "Ignoring degenerate line\n";
						delete line;
						vol->line[lineCount] = NULL;
					} else {
						line->volume.insert(vol);
						line = addLine(line);
						vol->line[lineCount] = line;
					}
					*/
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
				f->v1    = vol;
				f->face1 = surfCount;

				for(int v=0; v<v_siz; v++) {
					for(int u=0; u<u_siz; u++) {
						Point p(coef+coefNmb*(dim+rat), coef+(coefNmb+1)*(dim+rat));
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
				/*
				if(degen1 && degen2)
					// cout << "SURFACE DEGENERATES TO POINT\n";
					;
				else if(degen1 || degen2)
					// cout << "SURFACE DEGENERATES TO LINE\n";
					;

				if(degen1 || degen2) {
					delete f;
					vol->face[surfCount] = NULL;
				} else {
					f = addFace(f);
					vol->face[surfCount] = f;
				}
				*/
				f->degen1 = degen1;
				f->degen2 = degen2;
				f = addFace(f);
				vol->face[surfCount] = f;
				surfCount++;
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
			(*it)->volume.insert(v->volume.begin(), v->volume.end());
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
			(*it)->volume.insert(l->volume.begin(), l->volume.end());
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
 *************************************************************************************/
Face* TopologySet::addFace(Face* f) {
	set<Face*>::iterator it;
	for(it=face_.begin(); it != face_.end(); it++) {
		// cout << "Testing equality " << (*it)->v1->id << "(" << (*it)->face1 << ") - " << f->v1->id << "(" << f->face1 << ")\n";
		if( (*it)->equals(f, tol) ) {
			if((*it)->v2) {
				// actually you could optimize by putting this test before
				// testing for surface equality. Just gonna leave it here
				// for now to help in verification of the code
				cerr << "Error: Face given three neighbouring volumes\n";
				exit(1);
			}
			// cout << "YAY - face equals!!!\n";
			(*it)->v2    = f->v1;
			(*it)->face2 = f->face1;
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
	set<Face*>   face;
	getBoundaries(vert, line, face);
	return vert;
}

/*! \brief fetch all lines on the boundary of the model */
std::set<Line*> TopologySet::getBoundaryLines() {
	set<Vertex*> vert;
	set<Line*>   line;
	set<Face*>   face;
	getBoundaries(vert, line, face);
	return line;
}

/*! \brief fetch all faces on the boundary of the model */
std::set<Face*> TopologySet::getBoundaryFaces() {
	set<Face*> results;
	set<Face*>::iterator it;
	for(it=face_.begin(); it != face_.end(); it++)
		if((*it)->v2 == NULL)
			results.insert(*it);
	return results;
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
 * \note that in a very special case (i.e. single-volume pawn), this method may give
 *       degenerated internal faces as "boundaries" even if they are completely contained
 *       inside the model. This is because they are only neighbouring one volume and is as
 *       such considered on the boundary. This will in turn give wrong internal lines and
 *       vertices as well.
 *************************************************************************************/
void TopologySet::getBoundaries(std::set<Vertex*> vertices, std::set<Line*> lines, std::set<Face*> faces) {
	faces.clear();
	lines.clear();
	vertices.clear();
	faces = getBoundaryFaces();
	set<Face*>::iterator it;
	for(it=faces.begin(); it != faces.end(); it++) {
		vector<int> lineNumb = Line::getLineEnumeration( (*it)->face1 );
		for(int i=0; i<4; i++) {
			lines.insert( (*it)->v1->line[lineNumb[i]] );
			vertices.insert( (*it)->v1->line[lineNumb[i]]->v1 ); // yes, this will insert all vertices (at least) twice, but it is preffered
			vertices.insert( (*it)->v1->line[lineNumb[i]]->v2 ); // to the alternative which is making more functions 
		}
	}
}

set<Volume*>::iterator TopologySet::volume_begin() {
	return volume_.begin();
}

set<Volume*>::iterator TopologySet::volume_end() {
	return volume_.end();
}

set<Face*>::iterator TopologySet::face_begin() {
	return face_.begin();
}

set<Face*>::iterator TopologySet::face_end() {
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
	set<Face*>::iterator it;
	for(it=face_.begin(); it != face_.end(); it++)
		if( !(*it)->isDegen() )
			results++;
	return results;
}

/*! \brief Number of unique volumes in the model */
int TopologySet::numbVolumes() const {
	return volume_.size();
}

/*! \brief get one specific volume */
Volume* TopologySet::getVolume(int id) {
	set<Volume*>::iterator it;
	for(it=volume_.begin(); it != volume_.end(); it++)
		if( (*it)->id == id)
			return *it;
	return NULL;
}

