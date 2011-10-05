#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stack>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <GoTools/geometry/ObjectHeader.h>
#include <GoTools/geometry/SplineSurface.h>
#include <GoTools/geometry/SurfaceInterpolator.h>
#include <GoTools/trivariate/SplineVolume.h>
#include <GoTools/trivariate/VolumeInterpolator.h>

using namespace Go;

/*!
  \brief Program for order and continuity elevation of NURBS finite element.

  This program increases the order and continuity of a spline finite element
  mesh in the form of a g2-file by one. The application is to mixed spline
  finite element formulations.

  The input to the program is specified through the following 
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file: Input file in g2-format
  \arg -2D:           Input file with spline surfaces
  \arg -3D:           Input file with spline volumes
*/


int main(int argc, char** argv)
{
  bool surface_model = true;
  char* infile = 0;

  for (int i = 1; i < argc; i++)
    if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;
  
  size_t i = 0;
  while (i < strlen(infile) && isspace(infile[i])) i++;
  std::ifstream isp(infile+i);

  // For spline surface models
  std::vector<SplineSurface*> in_surf;

  // For spline volume models
  std::vector<SplineVolume*> in_vol;

  ObjectHeader head;
  int n = 0;
  while (!isp.eof()) {
    head.read(isp);
    if (head.classType() == Class_SplineVolume) {
      SplineVolume* v(new SplineVolume());
      v->read(isp);
      in_vol.push_back(v);
      surface_model = false;
    }
    else if (head.classType() == Class_SplineSurface) {
      SplineSurface* s(new SplineSurface());
      s->read(isp);
      in_surf.push_back(s);
      surface_model = true;
    }
    else
      std::cerr << "Unknown spline model" << std::endl;
    
    // Ignore blanks
    ws(isp); 
  }

  if (surface_model) {
    std::vector<SplineSurface*> out_surf;

    for (size_t i = 0;i < in_surf.size();i++) {
      SplineSurface* s_it = in_surf[i];

      // basis1 should be one degree higher than basis2 and C^p-1 continuous
      int ndim = s_it->dimension();
      Go::BsplineBasis b1 = s_it->basis(0).extendedBasis(s_it->order_u()+1);
      Go::BsplineBasis b2 = s_it->basis(1).extendedBasis(s_it->order_v()+1);

      // Note: Currently this is implemented for non-rational splines only.
      // TODO: Ask the splines people how to fix this properly, that is, how
      // may be obtain the correct weights for basis1 when *surf is a NURBS?
      if (s_it->rational())
        std::cerr <<"WARNING: The geometry basis is rational (using NURBS)\n."
                  <<"         The basis for the unknown fields of one degree"
                  <<" higher will however be non-rational.\n"
                  <<"         This may affect accuracy.\n"<< std::endl;

      // Compute parameter values of the Greville points
      size_t i;
      std::vector<double> ug(b1.numCoefs()), vg(b2.numCoefs());
      for (i = 0; i < ug.size(); i++)
        ug[i] = b1.grevilleParameter(i);
      for (i = 0; i < vg.size(); i++)
        vg[i] = b2.grevilleParameter(i);

      // Evaluate the spline surface at all points
      std::vector<double> XYZ(ndim*ug.size()*vg.size());
      s_it->gridEvaluator(XYZ,ug,vg);

      // Project the coordinates onto the new basis (the 2nd XYZ is dummy here)
      SplineSurface* s;
      s = Go::SurfaceInterpolator::regularInterpolation(b1,b2,
							ug,vg,XYZ,
							ndim,false,XYZ);
      
      out_surf.push_back(s);
      s->writeStandardHeader(std::cout);
      s->write(std::cout);
    }

    for (uint i = 0;i < out_surf.size();i++) {
      delete in_surf[i];
      delete out_surf[i];
    }
  }
  else {
    std::vector<SplineVolume*> out_vol; 

    for (size_t i = 0;i < in_vol.size();i++) {
      SplineVolume* v_it = in_vol[i];

      // basis1 should be one degree higher than basis2 and C^p-1 continuous
      int ndim = v_it->dimension();
      Go::BsplineBasis b1 = v_it->basis(0).extendedBasis(v_it->order(0)+1);
      Go::BsplineBasis b2 = v_it->basis(1).extendedBasis(v_it->order(1)+1);
      Go::BsplineBasis b3 = v_it->basis(2).extendedBasis(v_it->order(2)+1);

      // Note: Currently this is implemented for non-rational splines only.
      // TODO: Ask the splines people how to fix this properly, that is, how
      // may be obtain the correct weights for basis1 when *v_it is a NURBS?
      if (v_it->rational())
        std::cerr <<"WARNING: The geometry basis is rational (using NURBS)\n."
                  <<"         The basis for the unknown fields of one degree"
                  <<" higher will however be non-rational.\n"
                  <<"         This may affect accuracy.\n"<< std::endl;

      // Compute parameter values of the Greville points
      size_t i;
      std::vector<double> ug(b1.numCoefs()), vg(b2.numCoefs()), wg(b3.numCoefs());
      for (i = 0; i < ug.size(); i++)
        ug[i] = b1.grevilleParameter(i);
      for (i = 0; i < vg.size(); i++)
        vg[i] = b2.grevilleParameter(i);
      for (i = 0; i < wg.size(); i++)
        wg[i] = b3.grevilleParameter(i);

      // Evaluate the spline surface at all points
      std::vector<double> XYZ(ndim*ug.size()*vg.size()*wg.size());
      v_it->gridEvaluator(ug,vg,wg,XYZ);

      // Project the coordinates onto the new basis (the 2nd XYZ is dummy here)
      SplineVolume* v;
      v = Go::VolumeInterpolator::regularInterpolation(b1,b2,b3,
						       ug,vg,wg,XYZ,
						       ndim,false,XYZ);
      
      out_vol.push_back(v);
      v->writeStandardHeader(std::cout);
      v->write(std::cout);
    }

    for (uint i = 0;i < out_vol.size();i++) {
      delete in_vol[i];
      delete out_vol[i];
    }
  }

  return 0;
}  
