#ifndef MESH_HPP
#define MESH_HPP
#include "materials.hpp"
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

class Mesh
{
public:
	//These are for ordered transversal, store mesh_pos_id's
	list<int> boundary;
	list<int> n_type;
	list<int> p_type;

	//These are indexed by mesh_pos_id	
	int *is_boundary;
	int *is_n_type;
	int *is_p_type;
	Material **materials;//There should be only be two materials,
			     //each point will be associated with one

	double *mpos;
	int npoints;
	Mesh(double *ntype, int n_ntype, double *ptype, int p_ptype);
};

extern "C" Mesh *create_mesh(double *points, int n_points,
			     int *boundary, int nboundary,
			     int *ntype, int n_ntype, 
			     int *ptype, int n_ptype);
#endif
