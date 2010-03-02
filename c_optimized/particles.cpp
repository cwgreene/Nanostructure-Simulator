#include "materials.hpp"
#include "particles.hpp"
extern "C" Particles *construct_particles(double *pos, int *p_id, int *p_charge,
			       list<int> *p_live,
			       list<int> *p_dead)
{
	return new Particles(pos,p_id,p_charge,p_live,p_dead);
}


int create_particle(int mpos_id, Particles *p_data,int *density,
		        double *mesh_pos, int charge, Mesh *mesh)
{
	int i = p_data->p_dead->back();
	p_data->p_dead->pop_back();	
	p_data->p_live->push_back(i);
	double *particles = p_data->pos;
	px(i) = mesh_pos[2*mpos_id];
	py(i) = mesh_pos[2*mpos_id];
	
	//pkx(i) = material[mpos_id]->momentum
	//pky(u) = material[mpos_id]->
	material_random_momentum(mesh->materials[mpos_id],
				 p_data->pos+MOMENTUMX(i));
	p_data->p_charge[i] = charge;

	return i;
}
