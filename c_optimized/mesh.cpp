#include "mesh.hpp"

Point::Point(int _mpos_id, double *_mpos, int _type)
{
	mpos_id = _mpos_id;
	mpos = _mpos;
	type = _type;
}

Mesh::Mesh(double *points, int n_points,
	   Material **materials,
	   int *boundary, int nboundary,
		int *ntype, int n_ntype,
		int *ptype, int n_ptype)
{
	this->mpos = new double[2*n_points];
	this->boundary = new int[n_points];
	this->is_n_type = new int[n_points];
	this->is_p_type = new int[n_points];
	this->materials = new int[n_points];

	for(int i = 0; i < 2*n_points;i+=2)
	{
		this->mpos[i] = points[i];
		this->mpos[i+1] = points[i+1];
	}
	for(int i = 0; i < nboundary;i++)
	{
		this->is_boundary[boundary[i]] = 1;
		this->boundary.push_back(boundary[i]);
	}
	for(int i = 0; i < n_ntype;i++)
	{
		this->is_p_type[ptype[i]] = 1;
		this->p_type.push_back(ptype[i]);
	}
	for(int i = 0; i < p_ntype;i++)
	{
		this->is_p_type[ptype[i]] = 1;
		this->p_type.push_back(ptype[i]);
	}
	for(int i = 0; i < n_points;i++)
	{
		if(is_p_type[i])
			this->materials[i] = materials[0];
		if(is_n_type[i])
			this->materials[i] = materials[1];

	}
}

extern "C" Mesh *create_mesh(double *points, int n_points,
			     Material **materials,
			     int *boundary, int nboundary,
			     int *ntype, int n_ntype, 
			     int *ptype, int n_ptype)
{
	Mesh *new_mesh = new Mesh(double *points, int n_points,
				  Material **materials,
				  int *boundary, int nboundary,
				  int *ntype, int n_ntype,
				  int *ptype, int n_ptype);
}
