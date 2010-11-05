#ifndef PARTICLES_HPP
#define PARTICLES_HPP
#include <list>
#include <iostream> 
extern "C"{
#include "kdtree.h"
#include "kdtree3.h"
}
//Types for headers
template<class KD,int dim> class Mesh;

//Headers
#include "materials.hpp"

//Macros
#define POSITIONSTART(i) ((2*dim)*i)
#define MOMENTUMSTART(i) ((2*dim)*i+dim)
#define POSITIONX(i)  (4*(i))
#define POSITIONY(i)  (4*(i)+1)
#define MOMENTUMX(i)  (4*(i)+2)
#define MOMENTUMY(i)  (4*(i)+3)

//pnc expects 'dim' to be locally defined.
//pnc expects 'particles' to be locally defined
//Remember there is both momentum and position taken care of
//so for each dimension there are two variables. Hence the times 2
#define pnx(n,c) (particles[(dim*2)*n+c])
#define pknx(n,c) (particles[(dim*2)*n+dim+c])
//#define px(i)  particles[POSITIONX(i)]
//#define py(i)  particles[POSITIONY(i)]
//#define pkx(i)  particles[MOMENTUMX(i)]
//#define pky(i)  particles[MOMENTUMY(i)]
//#define pmass(i)  particles[MASS(i)]

using namespace std;

class Particles
{
public:
	double *pos; //Array of particle's positions, and momentums, TODO: RENAME
	int *p_charge; //array of particle's charges
	int *p_id; //array of particle's _Mesh_ Position ids, probaly worst possible name
	list<int>::iterator *local_id; //array of local_ids, which point into the mesh local list
	list<int>::iterator *live_id; //array of live_ids, which point into the p_live list
	double *p_mass; //array Particle masses, scaled to real units
	list<int> *p_live; //list of living particles, particle id's
	list<int> *p_dead; //List of available dead particles, particle ids
	int dim;//Dimension of particles
	Particles(double *_pos, int *_p_id, int *_p_charge,double *_p_mass,
			list<int> *_p_live, list<int> *_p_dead,int _dim)
	{
		//Initialize pointer blocks
		pos = _pos;
		p_charge = _p_charge;
		p_id = _p_id;
		p_live = _p_live; //should be empty
		p_dead = _p_dead;
		p_mass = _p_mass;
		int num = p_dead->size();
		this->local_id = new list<int>::iterator[num];
		this->live_id = new list<int>::iterator[num];
		dim = _dim;
	}
};

class Particle
{
public:
	int id;
	Particles *particles;
};

extern "C" int create_particleC(int mpos_id, Particles *p_data,int *density,
		        int charge, double mass,void *mesh);

template <class KD,int dim> 
void pick_up_particle(int part_id, 
		Particles *p_data, 
		int *density, 
		Mesh<KD,dim> *mesh);
template <class KD,int dim>
void put_down_particle(int part_id, 	
			Particles *p_data, 
			int *density,
			Mesh<KD,dim> *mesh);
list<int>::iterator destroy_particle(Particles *p_data, int part_id, list<int>::iterator pos);

/*templates*/
//TODO: This function should be made part of mesh
//Merge first two arguments into Particle &
template <class KD,int dim> 
void pick_up_particle(int part_id, Particles *p_data, 
			int *density, Mesh<KD,dim> *mesh)
{
	int mesh_pos_id = p_data->p_id[part_id];
	if(p_data->local_id[part_id] == NULL)
		cout << "We're about to die"<<endl;
	list<int>::iterator it = p_data->local_id[part_id];
	if(p_data->p_charge[part_id] < 0) //it's an electron
	{
		mesh->electrons_pos[mesh_pos_id].erase(it); //unlink particle from local list
	}
	else //it's a hole
	{
		mesh->holes_pos[mesh_pos_id].erase(it); //unlink particle from local list
	}
	density[mesh_pos_id] -= p_data->p_charge[part_id];  //Pick particle up from density
}

/*put down particle assumes that point is valid? If so, why do we bother with the kdtree_find_point_id call?*/
template <class KD,int dim>
void put_down_particle(int part_id, Particles *p_data, 
			int *density,Mesh<KD,dim> *mesh)
{	
	
	//Okay, this is now valid since 
	//We update immediately after move_particles
	int mesh_pos_id = p_data->p_id[part_id];
	//cout << "put_down_particle: exiting kd_tree"<<endl;
	density[mesh_pos_id] += p_data->p_charge[part_id];  //Put particle back down

	if(p_data->p_charge[part_id] < 0) //it's an electron
	{
		mesh->electrons_pos[mesh_pos_id].push_back(part_id); //add particle from local list
		p_data->local_id[part_id] = --(mesh->electrons_pos[mesh_pos_id].end()); //grab iterator position
	}
	else //it's a hole
	{
		mesh->holes_pos[mesh_pos_id].push_back(part_id); //add particle from local list
		p_data->local_id[part_id] = --(mesh->holes_pos[mesh_pos_id].end()); //grab iterator position
	}
	/*int holes = mesh->holes_pos[mpos_id].size();
	int electrons = mesh->electrons_pos[mpos_id].size();

	if((mesh->gen_num*empty_sign+(holes-electrons)) != density[mesh_pos_id])
	{	
		cout<<"INSUFFICIENT DENSITY"<<endl;
		exit(-1);	
	}*/
}

template<class KD,int dim> 
int create_particle(int mpos_id, Particles *p_data,int *density,
		        int charge, double mass,Mesh<KD,dim> *mesh)
{	
	int i = p_data->p_dead->back();
	p_data->p_dead->pop_back();	
	p_data->p_live->push_back(i);

	p_data->live_id[i] = --(p_data->p_live->end());  //Get inverse lookup of p_live
	if(charge < 0)
	{
		mesh->electrons_pos[mpos_id].push_back(i);
		p_data->local_id[i] = --(mesh->electrons_pos[mpos_id].end());
	}
	else
	{
		mesh->holes_pos[mpos_id].push_back(i);
		p_data->local_id[i] = --(mesh->holes_pos[mpos_id].end());
	}

	double *particles = p_data->pos;
	for(int c= 0; c < dim;c++)
	{
		pnx(i,c) = mesh->mpos[dim*mpos_id+c];
	}
	bool valid_trajectory = false;
	while(valid_trajectory != true){
		material_random_momentum(mesh->materials[mpos_id],
				 p_data->pos+MOMENTUMSTART(i),
				 mesh->particle_weight); //init pkx,pky
		if(mesh->valid_momentum(p_data->pos+MOMENTUMSTART(i),
					p_data->pos+POSITIONSTART(i)))
			valid_trajectory = true;
	}
	p_data->p_charge[i] = charge;
	p_data->p_mass[i] = mass;	
	p_data->p_id[i] = mpos_id;
	density[mpos_id] += charge;
	
	return i;
}

#ifndef PARTICLES_CPP
#include "particles.cpp"
#endif
#endif
