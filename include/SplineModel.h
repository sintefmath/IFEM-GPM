#ifndef _SPLINE_MODEL_H
#define _SPLINE_MODEL_H

#include <vector>
#include <GoTools/trivariate/SplineVolume.h>

class TopologySet;

/**********************************************************************************//**
 * \brief Global number ordering (gno) structure
 *
 * Contains all necessary information to go from any local enumeration (i,j,k) to a global enumeration.
 * One globNumber-struct should be available for each topological volume
 *************************************************************************************/
typedef struct globNumber {
	int vertex[8];          //!< global number of the 8 corner vertices
	int edge[12];           //!< start number for the 12 edge lines
	int edge_incr[12];      //!< +1 or -1 depending on whether the numbers are ascending or descending
	int surface[6];         //!< start number for the edge faces
	int surface_incr_i[6];  //!< increment by going in the first parametric direction
	int surface_incr_j[6];  //!< increment by going in the second parametric direction
	int volume;             //!< internal volume starting number
} globNumber;


/**********************************************************************************//**
 * \brief Main class for keeping track of the module topology, and properties
 *************************************************************************************/
class SplineModel {
	public:
		// constructors and destructors
		SplineModel();
		SplineModel(std::vector<boost::shared_ptr<Go::SplineVolume> > &spline_volumes);
		~SplineModel();

		// common get-functions
		TopologySet *getTopology();
		std::vector<boost::shared_ptr<Go::SplineVolume> > &getSplines();
		
		// model geometry functions
		void setTopologyTolerance(double tol);
		bool enforceRightHandSystem();

		// local-to-global mapping functions
		void generateGlobalNumbers();

		// model property functions
		bool addVolumePropertyCode(int volId, int propCode, bool inclusive=false);
		bool addFacePropertyCode(int volId, int faceId, int propCode, bool inclusive=true);
		bool addLinePropertyCode(int volId, int lineId, int propCode, bool inclusive=true);
		void addVertexPropertyCode(int volId, int vertId, int propCode);
		int getVolumePropertyCode(int volId);
		int getFacePropertyCode(int volId, int faceId);
		int getLinePropertyCode(int volId, int lineId);
		int getVertexPropertyCode(int volId, int vertId);
		
		// input-/output-functions
		void writeSplines(std::ostream &os) const;
		void writeGlobalNumberOrdering(std::ostream &os) const;
		void writeModelProperties(std::ostream &os) const;
		void readSplines(std::istream &is);
		void readGlobalNumberOrdering(std::istream &is);
		void readModelProperties(std::istream &is);

	private:
		TopologySet *topology;
		globNumber *l2g;
		std::vector<boost::shared_ptr<Go::SplineVolume> > spline_volumes_; //!< Spline objects

};

#endif

