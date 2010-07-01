#define MESH_CPP
#include "mesh.hpp"
#include "Polytope.hpp"
#include "materials.hpp"
#include "particles.hpp"

#include <stdio.h>
#include <iostream>
using namespace std;
Point::Point(int _mpos_id, double *_mpos, int _type)
{
	mpos_id = _mpos_id;
	mpos = _mpos;
	type = _type;
}

typedef Material *Material_Ptr;

template <class KD,int dim> 
Mesh<KD,dim>::Mesh(double *points, 
		int n_points,
		Material **materials,
		int *boundary, int nboundary,
		int *ntype, int n_ntype,
		int *ptype, int n_ptype,
		KD *_kdt,
	   	int gen_num,
	   	double particle_weight)
{
	this->npoints = n_points;
	this->mpos = new double[dim*n_points];
	this->is_boundary = new int[n_points];
	this->is_n_type = new int[n_points];
	this->is_p_type = new int[n_points];
	this->materials = new Material_Ptr[n_points];
	this->electrons_pos = new list<int>[n_points];
	this->holes_pos = new list<int>[n_points];
	this->particle_weight = particle_weight;
	this->gen_num = gen_num;

	//This is pretty hideous.
	//I should probably figure out how to nicely
	//abstract all this
	//Performance isn't all that important since this is
	//just called once, and better abstractions would be nice	
	for(int i = 0; i < dim*n_points;i+=dim)
	{
		//copy points
		for(int d = 0; d < dim;d++)
			this->mpos[i+d] = points[i+d];

		is_n_type[i/dim] = 0;
		is_p_type[i/dim] = 0;
	}
	for(int i = 0; i < nboundary;i++)
	{
		this->is_boundary[boundary[i]] = 1;
		this->boundary.push_back(boundary[i]);
	}
	for(int i = 0; i < n_ntype;i++)
	{	
		this->is_n_type[ntype[i]] = 1;
		this->n_type.push_back(ntype[i]);
	}
	for(int i = 0; i < n_ptype;i++)
	{
		this->is_p_type[ptype[i]] = 1;
		this->p_type.push_back(ptype[i]);
	}
	for(int i = 0; i < n_points;i++)
	{
		if(is_p_type[i])
		{
			this->materials[i] = materials[0];
		}
		else if(is_n_type[i])
		{
			this->materials[i] = materials[1];
		}
		else
		{
			cout << "HUH?!"<<endl;
			this->materials[i] = materials[0];
		}
	}
	//Zero electron and hole count at each point
	for(int i = 0; i< n_points;i++)
	{
		electrons_pos[i]= list<int>();
		holes_pos[i] = list<int>();
	}
	this->kdt = _kdt;
}


//find_point_id templates
template <>
int Mesh<kdtree3,3>::find_point_id(double *pos)
{
	vector3 x = {pos[0],pos[1],pos[2]};
	return kdtree_find_point3_id(kdt,&x);
}

template <>
int Mesh<kdtree,2>::find_point_id(double *pos)
{	
	vector2 x = {pos[0],pos[1]};
	return kdtree_find_point_id(kdt,&x);
}

//Python interfaces
extern "C" void *create_mesh(double *points, int n_points, int dim,
			     Material **materials,
			     int *boundary, int nboundary,
			     int *ptype, int n_ptype, 
			     int *ntype, int n_ntype, void *kdt,
			     int gen_num,
			     double particle_weight)
{
	if(dim == 2)
	{
		Mesh<kdtree,2> *mesh = new Mesh<kdtree,2>(points, n_points,
					  materials,
					  boundary, nboundary,
					  ntype, n_ntype,
					  ptype, n_ptype, (kdtree *)kdt,
					  gen_num,
					  particle_weight);
		return mesh;
	}
	if(dim == 3)
	{
		Mesh<kdtree3,3> *mesh = new Mesh<kdtree3,3>(
						points, n_points,
						materials,
						boundary, nboundary,
						ntype, n_ntype,
						ptype, n_ptype, (kdtree3 *)kdt,
						gen_num,
						particle_weight);
		return mesh;

	}
	cout << "Invalid dimension ("<<dim<<")"<<endl;
	exit(-1);
	return 0;
}

//2 Dimensional current
template <>
double Mesh<kdtree,2>::current_exit(Particles *p_data, int part_id)
{
	double *particles = p_data->pos;
	int dim = 2;
	double v_x = pknx(part_id,0)/p_data->p_mass[part_id];
	double v_y = pknx(part_id,1)/p_data->p_mass[part_id];
	return sqrt(v_x*v_x+v_y*v_y);
}

//3 Dimensional current
template <>
double Mesh<kdtree3,3>::current_exit(
				Particles *p_data, 
				int part_id,
				Face &exit_face)
{
	Vector3f velocity;
	double *particles = p_data->pos;
	int dim = 3;
	for(int i = 0; i < 3; i++)
	{
		velocity[i] = pknx(part_id,i)/p_data->p_mass[part_id];
	}
	return velocity.dot(exit_face.normal);
}

Face *Mesh::nearest_edge(double *point)
{
	double inner_dist;
	double outer_dist;
	Vector3f vpoint(point[0],point[1],point[2]);
	Face *inner_face,outer_face;
	inner_face = mesh->outer->nearest_face(vpoint,&inner_dist);
	outer_face = mesh->inner->nearest_face(vpoint,&outer_dist);
	return inner_face < outer_face ? inner_face : outer_face;
}
