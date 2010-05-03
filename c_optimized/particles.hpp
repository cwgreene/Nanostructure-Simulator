#ifndef PARTICLES_HPP
#define PARTICLES_HPP
#include <list>
#include "mesh.hpp"
#include "materials.hpp"

#define POSITIONX(i)  4*i
#define POSITIONY(i)  4*i+1
#define MOMENTUMX(i)  4*i+2
#define MOMENTUMY(i)  4*i+3
#define MASS(i)  4*i+4

#define px(i)  particles[POSITIONX(i)]
#define py(i)  particles[POSITIONY(i)]
#define pkx(i)  particles[MOMENTUMX(i)]
#define pky(i)  particles[MOMENTUMY(i)]
//#define pmass(i)  particles[MASS(i)]

using namespace std;

class Particles
{
public:
	double *pos; //Array of particle's positions, and momentums, TODO: RENAME
	int *p_charge; //array of particle's charges
	int *p_id; //array of particle's mesh Position ids
	list<int>::iterator *local_id; //array of local_ids, which point into the mesh local list
	double *p_mass; //array Particle masses, scaled to real units
	list<int> *p_live; //list of living particles, particle id's
	list<int> *p_dead; //List of available dead particles, particle ids
	Particles(double *_pos, int *_p_id, int *_p_charge,double *_p_mass,
			list<int> *_p_live, list<int> *_p_dead,Mesh *mesh)
	{
		pos = _pos;
		p_charge = _p_charge;
		p_id = _p_id;
		p_live = _p_live;
		p_dead = _p_dead;
		p_mass = _p_mass;
		int num = p_dead->size();
		this->local_id = new list<int>::iterator[num];
		std::cout << "Imminent doom!" << p_dead->size() <<std::endl;
	}
};

extern "C" int create_particle(int mpos_id, Particles *p_data,int *density,
		        int charge, double mass,Mesh *mesh);
void pick_up_particle(int part_id, Particles *p_data, int *density, Mesh *mesh);
void put_down_particle(int part_id, Particles *p_data, int *density,Mesh *mesh);
list<int>::iterator destroy_particle(Particles *p_data, int part_id, list<int>::iterator pos);

#endif
