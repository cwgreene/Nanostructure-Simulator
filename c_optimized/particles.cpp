#include <iostream>
#include "materials.hpp"
#include "particles.hpp"
using namespace std;
extern "C" Particles *init_particles(double *pos, int *p_id, 
				       int *p_charge,
				       double *p_mass,
				       list<int> *p_live,
				       list<int> *p_dead)
{
	return new Particles(pos,p_id,p_charge,p_mass,p_live,p_dead);
}


extern "C" int create_particle(int mpos_id, Particles *p_data,int *density,
		        int charge, double mass,Mesh *mesh)
{	
	//cout << "P_dead size" << p_data->p_dead->size() << endl;
	///cout << "P_live size" << p_data->p_live->size() << endl;
	//cout << "P_live loc" << p_data->p_live << endl;
	//cout << "P_dead loc" << p_data->p_dead << endl;
	int i = p_data->p_dead->back();
	p_data->p_dead->pop_back();	
	p_data->p_live->push_back(i);
	double *particles = p_data->pos;
	px(i) = mesh->mpos[2*mpos_id];
	py(i) = mesh->mpos[2*mpos_id+1];
	material_random_momentum(mesh->materials[mpos_id],
				 p_data->pos+MOMENTUMX(i));
	p_data->p_charge[i] = charge;
	p_data->p_mass[i] = mass;	
	density[mpos_id] += charge;
	//cout << "created: "<< px(i) << "/" << py(i) << " " << pkx(i) <<"/"
	//	<< pky(i) << endl;
	return i;
}
