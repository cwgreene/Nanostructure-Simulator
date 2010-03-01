#ifndef MESH_HPP
#define MESH_HPP
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
	list<Point *> mesh_points;
	double *positions;
	Mesh();
};

extern "C" create_mesh(double *ntype, int n_ntype, double *ptype, int n_ptype);
#endif
