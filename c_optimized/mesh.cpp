#define MESH_CPP
#include "mesh.hpp"
#include "materials.hpp"

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
template <class KD> Mesh<KD>::Mesh(double *points, int n_points, int dim,
	   Material **materials,
	   int *boundary, int nboundary,
	   int *ntype, int n_ntype,
	   int *ptype, int n_ptype,
	   KD *_kdt,double particle_weight)
{
	this->npoints = n_points;
	this->dim = dim;
	this->mpos = new double[dim*n_points];
	this->is_boundary = new int[n_points];
	this->is_n_type = new int[n_points];
	this->is_p_type = new int[n_points];
	this->materials = new Material_Ptr[n_points];
	this->electrons_pos = new list<int>[n_points];
	this->holes_pos = new list<int>[n_points];
	this->particle_weight = particle_weight;
	
	//cout << this->mpos << endl;
	cout << "c_mesh" << this << endl;
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
	printf("hi3\n");
	for(int i = 0; i < n_ntype;i++)
	{	
		this->is_n_type[ntype[i]] = 1;
		this->n_type.push_back(ntype[i]);
	}
	printf("hi4\n");
	for(int i = 0; i < n_ptype;i++)
	{
		this->is_p_type[ptype[i]] = 1;
		this->p_type.push_back(ptype[i]);
	}
	printf("hi5\n");
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
	printf("hi6\n");
}

extern "C" void *create_mesh(double *points, int n_points, int dim,
			     Material **materials,
			     int *boundary, int nboundary,
			     int *ptype, int n_ptype, 
			     int *ntype, int n_ntype, kdtree *kdt,
			     double particle_weight)
{
	Mesh<kdtree> *bob = new Mesh<kdtree>(points, n_points, dim,
				  materials,
				  boundary, nboundary,
				  ntype, n_ntype,
				  ptype, n_ptype, kdt,particle_weight);
	//cout << bob << endl;;

	return bob;
}
/*
extern "C" void *create_mesh3(double *points, int n_points, int dim, Material **materials,
			     int *boundary, int nboundary,
			     int *ptype, int n_ptype, 
			     int *ntype, int n_ntype, kdtree *kdt,
			     double particle_weight)
{
	Mesh<kdtree3> = new Mesh(points, n_points, dim,
				  materials,
				  boundary, nboundary,
				  ntype, n_ntype,
				  ptype, n_ptype, kdt,particle_weight);
	//cout << bob << endl;;

	return (void *)mesh; 
}*/
