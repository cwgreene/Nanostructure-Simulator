#ifndef PARTICLES_HPP
#define PARTICLES_HPP
#include <list>
#include <iostream> 

#include "mesh.hpp"
#include "materials.hpp"

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
#define pnx(n,c) particles[(dim*2)*n+c]
#define pknx(n,c) particles[(dim*2)*n+c]
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
	int *p_id; //array of particle's mesh Position ids
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

extern "C" int create_particle(int mpos_id, Particles *p_data,int *density,
		        int charge, double mass,void *mesh);

template <class KD> 
	void pick_up_particle(int part_id, 
		Particles *p_data, 
		int *density, 
		Mesh<KD> *mesh);
template <class KD>
	void put_down_particle(int part_id, 	
				Particles *p_data, 
				int *density,
				Mesh<KD> *mesh);
list<int>::iterator destroy_particle(Particles *p_data, int part_id, list<int>::iterator pos);

#endif
