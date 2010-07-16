#ifndef _PRIMITIVES_H
#define _PRIMITIVES_H

#include <GoTools/utils/Point.h>
#include <GoTools/geometry/Streamable.h>
#include <vector>
#include <set>
#include <iostream>


class Volume;
class Face;
class Line;
class Vertex;

/**********************************************************************************//**
 * \brief Class describing topological volume primitives
 *
 * A volume is defined by 8 corner vertices, 12 edge lines and 6 faces (collectively referred 
 * to as primitives). The primitives are stored as pointers and carefully not duplicated such that
 * equality testing can be done by address testing. Moreover, the class gives support for degenerate
 * volumes, making it possible that one unique vertex is stored in several of the 8 corner positions,
 * or likewise for the lines and faces.
 *
 *************************************************************************************/
class Volume {
	public:
		Volume(int id) ;
		std::vector<int> getVertexEnumeration(Vertex *v);
		void getEdgeEnumeration(Line *l, std::vector<int> &numb, std::vector<int> &parDir, std::vector<int> &parStep);
		std::vector<int> getSurfaceEnumeration(Face *f);

		Vertex *corner[8];
		Line   *line[12];
		Face   *face[6];
		int    id;
};

/**********************************************************************************//**
 * \brief Class describing topological face primitives
 *
 * A face is is defined by its internal controlpoints, which are also used for checking
 * uniqueness. Once two faces are found equal, the form a twin-pair with access to both
 * parent volumes. There are eight different ways that these faces can be oriented with
 * respect to each other, and these are stored in the boolean variables uv_flip, u_reverse
 * and v_reverse.
 *
 *************************************************************************************/
class Face : public Go::Streamable {
	public:
		Face(int n1, int n2);
		bool equals(Face *f, double tol);
		void write(std::ostream &os) const;
		void read(std::istream &is);

		std::vector<std::vector<Go::Point> > cp; //!< Defining control points
		Volume *v1;        //!< First neighbouring volume
		Volume *v2;        //!< Second neighbouring volume (or null in the case of edge face)
		int    face1;      //!< face index on volume 1
		int    face2;      //!< face index on volume 2
		bool   uv_flip;    //!< is the first/second parameter direction swapped on the twin face
		bool   u_reverse;  //!< is the first parameter direction reversed
		bool   v_reverse;  //!< is the second parameter direction reversed
};

/**********************************************************************************//**
 * \brief Class describing topological line primitives
 *
 * A line is defined by its internal controlpoints and contains access to all neighbouring volumes
 * as well as its start- and end-vertex
 *
 *************************************************************************************/
class Line : public Go::Streamable {
	public:
		bool equals(Line *l, double tol);
		void write(std::ostream &os) const;
		void read(std::istream &is);
		static int getLineEnumeration(int vert1, int vert2);

		std::vector<Go::Point> cp; //!< Defining control points
		std::set<Volume*> volume;  //!< Neighbouring volumes
		Vertex *v1;                //!< Start or stop vertex.
		Vertex *v2;                //!< Start or stop vertex.
};

/**********************************************************************************//**
 * \brief Class describing topological vertex primitives
 *
 * A vertex is a single topological corner point, defined by it's corresponding controlpoint
 * and has access to all neighbouring volumes
 *
 *************************************************************************************/
class Vertex {
	public:
		Go::Point cp;              //!< Defining control point
		std::set<Volume*> volume;  //!< Neighbouring volumes
};

#endif
