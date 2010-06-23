#ifndef MESH_HPP
#define MESH_HPP
#include <list>
#include "materials.hpp"
extern "C"{
#include "kdtree.h"
}
using namespace std;
enum {P_TYPE,N_TYPE};
class Point
{
public:
	int mpos_id;
	double *mpos;
	int type;

	//Methods
	Point(int mpos_id, double *mpos, int type);
};

template<class KD> class Mesh
{
public:
	int dim;
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
	Mesh(double *points, int n_points, int dim,
		Material **materials,
		int *boundary, int nboundary,
		int *ntype, int n_ntype,
		int *ptype, int n_ptype,KD *_kdt,double particle_weight);
};

/*extern "C" Mesh *create_mesh(double *points, int n_points,
			     int *boundary, int nboundary,
			     int *ntype, int n_ntype, 
			     int *ptype, int n_ptype);*/
#endif
