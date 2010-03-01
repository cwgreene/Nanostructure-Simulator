#ifndef PARTICLES_HPP
#define PARTICLES_HPP
#include <list>
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
		        double *mesh_pos,
			double energy,int charge, Material *material);
#endif
