#define PARTICLES_CPP
#include <iostream>
#include "particles.hpp"
#include "materials.hpp"
#include "mesh.hpp"
extern "C" {
#include "kdtree.h"
#include "kdtree3.h"
}
using namespace std;
extern "C" Particles *init_particles(double *pk_array, int *p_id, 
				       int *p_charge,
				       double *p_mass,
				       double *p_dist,
				       list<int> *p_live,
				       list<int> *p_dead,
				       int dim)
{
	cout << "dimension: "<<dim<<endl;
	return new Particles(pk_array,
			     p_id,
			     p_charge,p_mass,p_dist,
			     p_live,p_dead,
			     dim);
}

extern "C" int create_particleC(int mpos_id,Particles *p_data,int *density,
				int charge, double mass, void *mesh)
{
	if(p_data->dim == 3)
		return create_particle(mpos_id,p_data,density,charge,mass,
				(Mesh<kdtree3,3> *)mesh);
	if(p_data->dim == 2)
	{
		return create_particle(mpos_id,p_data,density,charge,mass,
				(Mesh<kdtree,2> *)mesh);
	}
	printf("Invalid dimension: %d\n",p_data->dim);
	exit(-3);
	return -3;
}



/**destroy_particle:
Remove particle from live list, add to dead list, 
return iterator to position in live list right afeter this particle.
Recycle particle id in array
**/
list<int>::iterator destroy_particle(Particles *p_data, int part_id, list<int>::iterator pos)
{
	list<int>::iterator it = p_data->p_live->erase(pos);
	p_data->p_dead->push_back(part_id);
	return it;
}
