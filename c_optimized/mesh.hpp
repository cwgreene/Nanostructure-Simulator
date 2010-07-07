#ifndef MESH_HPP
#define MESH_HPP
#include <list>

//Types needed for other headers (pointers)
template<class KD,int _dim> class Mesh;

//Headers
#include "Polytope.hpp"
#include "particles.hpp"
#include "materials.hpp"
extern "C"{
#include "kdtree.h"
#include "kdtree3.h"
}

//Namespaces
using namespace std;
enum {P_TYPE,N_TYPE};

//Classes
class Point
{
public:
	int mpos_id;
	double *mpos;
	int type;

	//Methods
	Point(int mpos_id, double *mpos, int type);
};

template<class KD,int dim> class Mesh
{
public:
	Polytope<dim> *outer;
	Polytope<dim> *inner;
	//These are for ordered transversal, store mesh_pos_id's
	list<int> boundary;
	list<int> n_type;
	list<int> p_type;

	//Local list of particles
	//Inverse hole and electron count lookup, indexed by mesh_pos_id
	list<int> *electrons_pos;//array of lists, indexed by mesh_pos_id
	list<int> *holes_pos;//Array of lists
	KD *kdt;

	//These are indexed by mesh_pos_id	
	int *is_boundary;
	int *is_n_type;
	int *is_p_type;
	Material **materials;//There should be only be two materials,
			     //each point will be associated with one

	double *mpos;//coordinates of points
	int npoints;
	double particle_weight;
	int find_point_id(double *position);
	double current_exit(Particles *p_data,int part_id);
	double current_exit(Particles *p_data,int part_id, Face *exit_face);
	int gen_num;
	Mesh(double *points, int n_points,
		Material **materials,
		int *boundary, int nboundary,
		int *ntype, int n_ntype,
		int *ptype, int n_ptype,KD *_kdt,
		int gen_num,double particle_weight);
	Face *nearest_edge(double *point);
};

/*extern "C" Mesh *create_mesh(double *points, int n_points,
			     int *boundary, int nboundary,
			     int *ntype, int n_ntype, 
			     int *ptype, int n_ptype);*/
#ifndef MESH_CPP
#include "mesh.cpp"
#endif
#endif
