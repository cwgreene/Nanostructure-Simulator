#ifndef MESH_HPP
#define MESH_HPP
#include <list>

//Types needed for other headers (pointers)
template<class KD,int _dim> class Mesh;

//Headers
#include "Utils.hpp"
#include "Boundary.hpp"
#include "Polytope.hpp"
#include "particles.hpp"
#include "mesh.hpp"
#include "materials.hpp"

extern "C"{
#include "kdtree.h"
#include "kdtree3.h"
}

class Particles;

//Namespaces
typedef Material *Material_Ptr;
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

	Boundary<dim> *rboundary;//For reflections
	//These are for ordered transversal, store mesh_pos_id's
	vector<int> boundary;
	vector<int> n_type;
	vector<int> p_type;

	//Local list of particles
	//Inverse hole and electron count lookup, indexed by mesh_pos_id
	list<int> *electrons_pos;//array of lists, indexed by mesh_pos_id
	list<int> *holes_pos;//Array of lists

	list<int> ***local_particles;//Array for sign based indexing

	KD *kdt;

	//These are indexed by mesh_pos_id	
	int *is_boundary;
	int *is_n_type;
	int *is_p_type;
	int *is_reflect; //TODO

	//Geometric information, indexed by mesh_pos_id
	//TODO: Initialize these
	Vector2d *normals;
	Vector2d *tangents;
	
	Material **materials;//There should be only be two materials (for now),
			     //each point will be associated with one

	double *mpos;//coordinates of points
	int npoints;
	double particle_weight;
	int gen_num;

	//Methods
	int find_point_id(double *position);
	double current_exit(Particles *p_data,int part_id);
	double current_exit(Particles *p_data,int part_id, Face *exit_face);
	int has_escaped(Particles *p_data,int i);
	bool reflect(Particles *p_data,int i,double *old_pos);
	bool valid_momentum(double *dmomentum, double *dpoint);

	//Constructors
	Mesh(double *points, int n_points,
		Material **materials,
		int *boundary, int nboundary,
		int *ntype, int n_ntype,
		int *ptype, int n_ptype,KD *_kdt,
		int gen_num,double particle_weight,
		Polytope<dim> *outer, Polytope<dim> *inner);
	Face *nearest_edge(double *point);
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template <class KD,int dim>
bool Mesh<KD,dim>::valid_momentum(double *dmomentum, double *dpoint)
{
	typename Vector<dim>::Type point(dpoint);
	typename Vector<dim>::Type momentum(dmomentum);

	return true;
	return rboundary->is_inward_facing(momentum,point);
}

/*extern "C" Mesh *create_mesh(double *points, int n_points,
			     int *boundary, int nboundary,
			     int *ntype, int n_ntype, 
			     int *ptype, int n_ptype);*/
template <class KD,int dim> 
Mesh<KD,dim>::Mesh(double *points, 
		int n_points,
		Material **materials,
		int *boundary, int nboundary,
		int *ntype, int n_ntype,
		int *ptype, int n_ptype,
		KD *_kdt,
	   	int gen_num,
	   	double particle_weight,
		Polytope<dim> *outer, Polytope<dim> *inner)
{
	std::cout << "boundary location: "<<boundary<<std::endl;
	std::cout << "boundary end: "<<boundary+n_points*sizeof(int)<<std::endl;
	std::cout << "ntype location: "<<ntype<<std::endl;
	std::cout << "ntype end: "<<ntype+n_points*sizeof(int)<<std::endl;
	std::cout << "ptype location: "<<ptype<<std::endl;
	std::cout << "ptype end: "<<ptype+n_points*sizeof(int)<<std::endl;
	std::cout << "kd location: "<<_kdt<<std::endl;
	std::cout << "outer location: "<<outer<<std::endl;
	std::cout << "inner location: "<<inner<<std::endl;
	std::cout << "this location: "<<this<<std::endl;
	std::cout << "n_points: "<<n_points<<std::endl;
	this->npoints = n_points;
	this->mpos = new double[dim*n_points];
	this->is_boundary = new int[n_points];
	this->is_n_type = new int[n_points];
	this->is_p_type = new int[n_points];
	this->materials = new Material_Ptr[n_points];
	this->electrons_pos = new list<int>[n_points];
	this->holes_pos = new list<int>[n_points];
	
	this->local_particles = new list<int> **[n_points];

	this->particle_weight = particle_weight;
	this->gen_num = gen_num;

	this->inner = inner;
	this->outer = outer;

	//This is pretty hideous.
	//I should probably figure out how to nicely
	//abstract all this
	//Performance isn't all that important since this is
	//just called once, and better abstractions would be nice	
	//Note that i is incremented by dim, not 1
	cout<<"Copying is_n_type"<<std::endl;
	for(int i = 0; i < dim*n_points;i+=dim)
	{
		//copy points
		for(int d = 0; d < dim;d++){
			this->mpos[i+d] = points[i+d];
		}

		is_n_type[i/dim] = 0;
		is_p_type[i/dim] = 0;
	}
	cout<<"Copying nboundary"<<std::endl;
	for(int i = 0; i < nboundary;i++)
	{
		this->is_boundary[boundary[i]] = 1;
		this->boundary.push_back(boundary[i]);
	}
	cout<<"Copying n_ntype"<<std::endl;
	for(int i = 0; i < n_ntype;i++)
	{	
		this->is_n_type[ntype[i]] = 1;
		this->n_type.push_back(ntype[i]);
	}
	cout<<"Copying n_ptype"<<std::endl;
	for(int i = 0; i < n_ptype;i++)
	{
		this->is_p_type[ptype[i]] = 1;
		this->p_type.push_back(ptype[i]);
	}
	cout<<"Associating materials with points"<<std::endl;
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
			//debug::dbg << "HUH?!"<<endl;
			this->materials[i] = materials[0];
		}
	}
	cout<<"Associating electron/hole lists with points"<<std::endl;
	//Zero electron and hole count at each point
	for(int i = 0; i< n_points;i++)
	{
		electrons_pos[i]= list<int>();
		holes_pos[i] = list<int>();
		local_particles[i] = new list<int> *[2];

		local_particles[i][0] = electrons_pos;
		local_particles[i][1] = holes_pos;
	}
	cout<<"Attaching kdtree"<<std::endl;
	//Attach kdtree.
	this->kdt = _kdt;
	cout<<"Attaching boundaries"<<std::endl;
	double *boundaries_a = new double[nboundary*dim];
	Eigen::Matrix<double,dim,1> interior_p;
	for(int i = 0; i < nboundary;i++)
	{
		std::cout << i<<std::endl;
		int offset = i*dim;
		int boffset = boundary[i]*dim;
		for(int k = 0; k < dim;k++)
		{
			boundaries_a[offset+k] = points[boffset+k];
			interior_p[k] += points[boffset+k];
			interior_p[k] /= nboundary;
		}
	}
	std::cout << "Mesh Construction Completed" << std::endl;
	rboundary= new Boundary<dim>(boundaries_a,nboundary,interior_p);
	std::cout << "Seriously" << std::endl;
			         
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
	//debug::dbg <<"find_point_id:"<<(i++)<<" "<< pos[0] << " " << pos[1] << ": ";
	int id =  kdtree_find_point_id(kdt,&x);	
	//cout <<id<<endl;
	return id;
}

//2 Dimensional current
template <>
double Mesh<kdtree,2>::current_exit(Particles *p_data, int part_id)
{
	//double *particles = p_data->pos;
//	int dim = 2;
//	double v_x = pknx(part_id,0)/(particle_weight*p_data->p_mass[part_id]);
//	double v_y = pknx(part_id,1)/(particle_weight*p_data->p_mass[part_id]);
	return particle_weight;
}

template <>
double Mesh<kdtree3,3>::current_exit(Particles *p_data, int part_id)
{	
	//debug::dbg<< "current(Particles *p_data,int part_id) not supported"<<
	//	" for dim =3"<<endl;
	return 1;
}

//3 Dimensional current
//Another problem here, since calling function
//will always call wrong one
template <>
double Mesh<kdtree3,3>::current_exit(
				Particles *p_data, 
				int part_id,
				Face *exit_face)

{
	Vector3d velocity;
	double *particles = p_data->pos;
	int dim = 3;
	for(int i = 0; i < 3; i++)
	{
		velocity[i] = pknx(part_id,i)/p_data->p_mass[part_id];
	}
	return velocity.dot(exit_face->normal);
}

//template <class KD, int dim>
//Face *Mesh<KD,dim>nearest_face::nearest_face(double *point)

template <class KD,int dim>
Face *Mesh<KD,dim>::nearest_edge(double *point)
{
	double inner_dist;
	double outer_dist;
	Vector3d vpoint(point[0],point[1],point[2]);
	Face *inner_face,*outer_face;
	inner_face = inner->nearest_face(vpoint,&inner_dist);
	outer_face = outer->nearest_face(vpoint,&outer_dist);
	if(inner_dist < outer_dist)
		return inner_face;
	return outer_face;
}

//has_escaped
//returns id of nearest_exit
//or -1 if it has not escaped.
//This makes me unhappy, we would like 0 if not escaped so we can say
//if(has_escaped), but since the id can be zero we can't do it that way.
//the alternative would involve antoher call to nearest_exit, 
//which also makes me unhappy.
template<>
int Mesh<kdtree,2>::has_escaped(Particles *p_data, 	
				int part_id)
{
	//double *particles = p_data->pos;
	int dim = 2; //Macro

	if(outer->contains(p_data->pos+2*dim*part_id))
	{
		//If outer contains the point
		//since the inner doesn't we're still inside
		//Problematic because we might have some generalized version
		
		//Then go check to make sure we aren't in the inner thingy
		if(inner != NULL)
		{
			if(!(inner->contains(p_data->pos+2*dim*part_id)))
				return -1;
		}else
			return -1;
	}
	//Check if nearest_exit is reflecting boundary
	//int nearest_exit = find_point_id(p_data->pos+2*dim*part_id);
	int nearest_exit = p_data->p_id[part_id];
	return nearest_exit;
}

//reflect
//takes the old position, it's current positoin
//and reflects it across the boundary.
template <class KD,int dim>
bool Mesh<KD,dim>::reflect(Particles *p_data,int id, double *old_pos)
{
	typename Boundary<dim>::VectorD voldpos(old_pos);
	return false;
	//return rboundary->reflect_trajectory(p_data->pos+id*dim*2,voldpos); //should be true
}

//has_escaped
//returns id of nearest_exit
//or -1 if it has not escaped.
//This makes me unhappy, we would like 0 if not escaped so we can say
//if(has_escaped), but since the id can be zero we can't do it that way.
//the alternative would involve antoher call to nearest_exit, 
//which also makes me unhappy.
template<>
int Mesh<kdtree3,3>::has_escaped(Particles *p_data, 	
				int part_id)
{
	//double *particles = p_data->pos;
	int dim = 3; //Macro
	
	//Check if inner is null, if not, check it
	//If it contains the thing, we ARE , and nearest_exit should be called
	//otherwise continue on
	if(inner != NULL) {
		if(inner->contains(p_data->pos+2*dim*part_id))
		{
//			int nearest_exit = find_point_id(p_data->pos+2*dim*part_id);
			int nearest_exit= p_data->p_id[part_id];
			return nearest_exit;
		}
	}

	//If outer contains the point
	//since the inner doesn't we're still inside
	//Problematic because we might have some generalized version
	if(outer->contains(p_data->pos+2*dim*part_id))
		return -1;
	//Check if nearest_exit is reflecting boundary
	int nearest_exit = find_point_id(p_data->pos+2*dim*part_id);
	//TODO: Implement is_reflecting and uncomment the following code
	//if(is_reflecting[nearest_exit])
	//{
	//	reflect(poly,p_data,part_id);
	//	return has_escaped(poly,p_data,part_id,mesh);
	//}
	return nearest_exit;
}
#ifndef MESH_CPP
#include "mesh.cpp"
#endif
#endif
