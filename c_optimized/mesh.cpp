#define MESH_CPP
#include "mesh.hpp"
#include "Polytope.hpp"
#include "materials.hpp"
#include "particles.hpp"

#include "debug.hpp"

#include <stdio.h>
#include <iostream>
using namespace std;
Point::Point(int _mpos_id, double *_mpos, int _type)
{
	mpos_id = _mpos_id;
	mpos = _mpos;
	type = _type;
}


//Python interfaces
extern "C" void *create_mesh(double *points, int n_points, int dim,
			     Material **materials,
			     int *boundary, int nboundary,
			     int *ptype, int n_ptype, 
			     int *ntype, int n_ntype, 
			     void *kdt,
			     int gen_num,
			     double particle_weight,
			     void *outer, void *inner)
{
	if(dim == 2)
	{
		Mesh<kdtree,2> *mesh = new Mesh<kdtree,2>(points, n_points,
					  materials,
					  boundary, nboundary,
					  ntype, n_ntype,
					  ptype, n_ptype, 
					  (kdtree *)kdt,
					  gen_num,
					  particle_weight,
					  (Polytope<2> *)outer,
					  (Polytope<2> *)inner);
		return mesh;
	}
	if(dim == 3)
	{
		Mesh<kdtree3,3> *mesh = new Mesh<kdtree3,3>(
						points, n_points,
						materials,
						boundary, nboundary,
						ntype, n_ntype,
						ptype, n_ptype, 
						(kdtree3 *)kdt,
						gen_num,
						particle_weight,
						(Polytope<3> *)outer,
						(Polytope<3>*)inner);
		return mesh;

	}
	//debug::dbg<< "Invalid dimension ("<<dim<<")"<<endl;
	exit(-1);
	return 0;
}


