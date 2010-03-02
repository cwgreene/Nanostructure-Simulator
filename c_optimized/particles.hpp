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
	double *pos;
	int *p_charge;
	int *p_id;
	list<int> *p_live;
	list<int> *p_dead;
	Particles(double *_pos, int *_p_id, int *_p_charge,
			list<int> *_p_live, list<int> *_p_dead)
	{
		pos = _pos;
		p_charge = _p_charge;
		p_id = _p_id;
		p_live = _p_live;
		p_dead = _p_dead;
	}
};

int create_particle(int mpos_id, Particles *p_data,int *density,
		        double *mesh_pos, int charge, Mesh *mesh);
#endif
