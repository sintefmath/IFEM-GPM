#ifndef _TOPOLOGYSET_H
#define _TOPOLOGYSET_H

#include <set>
#include <vector>
#include <GoTools/trivariate/SplineVolume.h>

class Volume;
class Face;
class Line;
class Vertex;

class TopologySet {
	public:
		TopologySet(std::vector<boost::shared_ptr<Go::SplineVolume> > &spline_volumes, double tol=1e-4);
		void buildTopology();
		Vertex* addVertex(Vertex* v);
		Line* addLine(Line* l);
		Face* addFace(Face* f);
		void addVolume(Volume* v);
		int numbVertices() const;
		int numbLines() const;
		int numbFaces() const;
		int numbVolumes() const;
		
		std::set<Vertex*>::iterator vertex_begin();
		std::set<Vertex*>::iterator vertex_end();
		std::set<Line*>::iterator line_begin();
		std::set<Line*>::iterator line_end();
		std::set<Face*>::iterator face_begin();
		std::set<Face*>::iterator face_end();

	private:
		double tol;
		std::set<Volume*> volume_;
		std::set<Face*>   face_;
		std::set<Line*>   line_;
		std::set<Vertex*> vertex_;
		std::vector<boost::shared_ptr<Go::SplineVolume> > spline_volumes_;

};

#endif
