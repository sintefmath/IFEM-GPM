#ifndef _TOPOLOGYSET_H
#define _TOPOLOGYSET_H

#include <set>
#include <vector>
#include "SplineModel.h"
#include "primitives.h"


template<class T> struct ltOper
{
  bool operator()(const T* f1, const T* f2) const
  {
    if (f1->id == -1)
      return true;
    return (f1->id < f2->id);
  }
};

typedef std::set<Volume*,ltOper<Volume> > VolSet;
typedef std::set<Face*,ltOper<Face> >  FaceSet;


/**********************************************************************************//**
 * \brief keeping track of the module topology and primitive (point/line/face) relations
 * 
 * This class keeps control of all topological relations needed to create a consistent
 * local-to-global enumeration mapping. It relies on comparison of control points for 
 * its construction, so it is recommended that one creates the topology first, and uses
 * this in all subsequent computations (i.e. before refinement). 
 *************************************************************************************/
class TopologySet {
	public:
		TopologySet(double tol=1e-4);
		TopologySet(std::vector<SurfacePointer> &spline_surfaces, double tol=1e-4);
		TopologySet(std::vector<VolumePointer> &spline_volumes, double tol=1e-4);
		~TopologySet();

		void addPatch(VolumePointer vol);
		void addPatch(SurfacePointer surf);
		void buildTopology(std::vector<bool>* periodic=NULL);

		int numbVertices() const;
		int numbLines() const;
		int numbNonDegenLines() const;
		int numbFaces() const;
		int numbNonDegenFaces() const;
		int numbVolumes() const;

		void setTolerance(double tol);
		void setPeriodic(int dir);

		Volume* getVolume(int id);
		Face* getFace(int id);
		
		std::set<Vertex*> getBoundaryVertices();
		std::set<Line*>   getBoundaryLines();
		FaceSet   getBoundaryFaces();
		void              getBoundaries(std::set<Vertex*> vertices, std::set<Line*> lines, FaceSet faces);

		std::set<Vertex*>::iterator vertex_begin();
		std::set<Vertex*>::iterator vertex_end();
		std::set<Line*>::iterator line_begin();
		std::set<Line*>::iterator line_end();
		FaceSet::iterator face_begin();
		FaceSet::iterator face_end();
		VolSet::iterator volume_begin();
		VolSet::iterator volume_end();

	private:
		void    addVolume(Volume* v);
		Face*   addFace(Face* f);
		Line*   addLine(Line* l);
		Vertex* addVertex(Vertex* v);

		double tol;                            //!< Control point tolerance (given in euclidean distance)
		VolSet  volume_;                       //!< All unique volumes
		FaceSet face_;                         //!< All unique (possible degenerate) faces
		std::set<Line*>   line_;               //!< All unique (possible degenerate) edge lines
		std::set<Vertex*> vertex_;             //!< All unique corner vertices
 
 		bool volumetric_model;
 		bool surface_model;
		std::vector<SurfacePointer>  spline_surfaces_; //!< Spline objects
		std::vector<VolumePointer> spline_volumes_; //!< Spline objects
};

#endif
