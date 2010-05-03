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
Mesh::Mesh(double *points, int n_points,
	   Material **materials,
	   int *boundary, int nboundary,
		int *ntype, int n_ntype,
		int *ptype, int n_ptype,
	   kdtree *_kdt)
{
	this->mpos = new double[2*n_points];
	this->is_boundary = new int[n_points];
	this->is_n_type = new int[n_points];
	this->is_p_type = new int[n_points];
	this->materials = new Material_Ptr[n_points];
	this->electrons_pos = new list<int>[n_points];
	this->holes_pos = new list<int>[n_points];
	
	//cout << this->mpos << endl;
	cout << "c_mesh" << this << endl;
	for(int i = 0; i < 2*n_points;i+=2)
	{
		this->mpos[i] = points[i];
		this->mpos[i+1] = points[i+1];
		is_n_type[i/2] = 0;
		is_p_type[i/2] = 0;
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

extern "C" Mesh *create_mesh(double *points, int n_points,
			     Material **materials,
			     int *boundary, int nboundary,
			     int *ptype, int n_ptype, 
			     int *ntype, int n_ntype, kdtree *kdt)
{
	Mesh *bob = new Mesh(points, n_points,
				  materials,
				  boundary, nboundary,
				  ntype, n_ntype,
				  ptype, n_ptype, kdt);
	//cout << bob << endl;;

	return bob;
}
