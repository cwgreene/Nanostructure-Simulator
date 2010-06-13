#include <iostream>
#include "materials.hpp"
#include "particles.hpp"
extern "C" {
#include "kdtree.h"
}
using namespace std;
extern "C" Particles *init_particles(double *pos, int *p_id, 
				       int *p_charge,
				       double *p_mass,
				       list<int> *p_live,
				       list<int> *p_dead,
				       Mesh *mesh)
{
	return new Particles(pos,p_id,p_charge,p_mass,p_live,p_dead,mesh);
}


extern "C" int create_particle(int mpos_id, Particles *p_data,int *density,
		        int charge, double mass,Mesh *mesh)
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
	px(i) = mesh->mpos[2*mpos_id];
	py(i) = mesh->mpos[2*mpos_id+1];
	material_random_momentum(mesh->materials[mpos_id],
				 p_data->pos+MOMENTUMX(i)); //init pkx,pky
	p_data->p_charge[i] = charge;
	p_data->p_mass[i] = mass;	
	p_data->p_id[i] = mpos_id;
	density[mpos_id] += charge;
	
	
	return i;
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

void pick_up_particle(int part_id, Particles *p_data, int *density, Mesh *mesh)
{
	int mesh_pos_id = p_data->p_id[part_id];
	if(p_data->local_id[part_id] == NULL)
		cout << "We're about to die"<<endl;
	if(p_data->p_charge[part_id] < 0) //it's an electron
	{
		list<int>::iterator it = p_data->local_id[part_id];
		mesh->electrons_pos[mesh_pos_id].erase(it); //unlink particle from local list
	}
	else //it's a hole
	{
		list<int>::iterator it = p_data->local_id[part_id];
		mesh->holes_pos[mesh_pos_id].erase(it); //unlink particle from local list
	}
	density[mesh_pos_id] -= p_data->p_charge[part_id];  //Pick particle up from density
}

void put_down_particle(int part_id, Particles *p_data, int *density,Mesh *mesh)
{	
	int mesh_pos_id = p_data->p_id[part_id];
	double *particles = p_data->pos; //alias for macro
	vector2 pos;
	pos[0] = px(part_id);
	pos[1] = py(part_id);
	p_data->p_id[part_id] = kdtree_find_point_id(mesh->kdt,&pos); //find nearest spot
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

}
