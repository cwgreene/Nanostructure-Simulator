#include "mesh.hpp"

Point::Point(int _mpos_id, double *_mpos, int _type)
{
	mpos_id = _mpos_id;
	mpos = _mpos;
	type = _type;
}

Mesh::Mesh(double *ntype, int n_ntype, double *ptype, int p_ptype)
{
	double *points = new double[2*(n_ntype+p_ptype)];
	mesh_points = new list<Point *>();
	for(int i = 0; i < 2*n_type;i+=2)
	{
		points[i] = ntype[i];
		points[i+1] = ntype[i+1];
		mesh_points.push_back(new Point(i,(points+i),N_TYPE));
	}
	for(int i = 2*n_ntype; i < 2*n_ntype;i+=2)
	{
		int j = (i/2-n_ntype);
		points[i] = ptype[j];
		points[i+1] = ptype[j+1];
		mesh_points.push_back(new Point(i,(points+i),N_TYPE));

	}
}

extern "C" Mesh *create_mesh(double *ntype, int n_ntype double *ptype, int n_ptype)
{
	Mesh *new_mesh = new Mesh(ntype,n_ntype,ptype,n_ptype);
}
