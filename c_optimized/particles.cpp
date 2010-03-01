#include <materials.h>
#include "particles.hpp"
extern "C" Particles *construct_particles(double *pos, int *p_id, int *p_charge,
			       list<int> *p_live,
			       list<int> *p_dead)
{
	return new Particles(pos,p_id,p_charge,p_live,p_dead);
}

void random_momentum(double energy,double *momentum//,Material *material)
{
	double theta = (((double)rand()/RAND_MAX)*2*PI;//0-1
	double speed = material->
	momentum[0] = cos(theta);
}

int create_particle(int mpos_id, Particles *p_data,int *density,
		        double *mesh_pos,
			double energy,int charge, Material *material)
{
	int i = p_data->p_dead.pop_back();
	p_data->p_live.push_back(i);
	double *particles = p_data->pos;
	px(i) = mesh_pos[2*mpos_id];
	py(i) = mesh_pos[2*mpos_id];
	
	//pkx(i) = material[mpos_id]->momentum
	//pky(u) = material[mpos_id]->
	material[mpos_id]->random_momentum(energy,(particles+4*i+2));
	p_charge[i] = charge;

	return i;
}
