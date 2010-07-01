//Standard Library Includes
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <list>
#include <iostream>

//HPP Includes
#include "particles.hpp"
#include "mesh.hpp"
#include "statistics.hpp"
#include "move_particles.hpp"
#include "Polytope.hpp"

extern "C" {
#include "kdtree.h"
#include "kdtree3.h"
}
//Macro Definitions
#define sq(x) (x*x)
#define EC (1.60217646e-19)
using namespace std;

extern "C" list<int> *init_particle_list(int n)
{
	list<int> *particle_list = new list<int>;
	srand(time(NULL));
	for(int i = 0; i < n;i++)
	{
		particle_list->push_back(i);
	}
	return particle_list;
}

extern "C" list<int> *init_dead_list(int start,int n)
{
	list<int> *particle_list = new list<int>;
	for(int i = start; i < n;i++)
	{
		particle_list->push_back(i);
	}
	return particle_list;
}

void randomElectronMovement(double *particles,
				double *p_mass,
				int *p_id,
				int *p_charge,
				int i,
				double *efield,
				double dt,
				double length_scale,
				double particle_weight,
				int dim)
{
	double dx[3] = {0.,0.,0.};
	//Move
	for(int c = 0; c < dim; c++)	
	{
		//if(i % 1000 == 0)
		//	cout<<"x_"<<c<<": "<< pnx(i,c) << endl;
		dx[c] = pknx(i,c)*dt/(length_scale*p_mass[i]*particle_weight);
		pnx(i,c) += dx[c];
		/*if(i % 1000 ==0)
		{
			cout<< "dx_"<<c<<": "<<dx[c]<<endl;
			cout<< "x_"<<c<<": "<<pnx(i,c) << endl;
		}*/
	}
	//Drift
	for(int c = 0; c < dim; c++)
	{
		pknx(i,c) += (efield[dim*p_id[i]]*p_charge[i]*dt)
				*EC*particle_weight;
		/*if(i %1000 == 0)
		{
		cout << "dir: "<<c<<" force: "<<
		 	((efield[dim*p_id[i]]*p_charge[i]*dt)
			*EC*particle_weight)<<endl;
		}*/
	}
	//Scatter... ugh... this is going to be complicated.
	//we're going to need to seperate drift and diffusion
	//we already knew that... but damn.
	double theta = rand()*(3.1415192)/RAND_MAX-3.141592/2;

	//handle x and y cooridinates.
	//this should be a bandstructure lookup
	double _pknx = pknx(i,0);
	double _pkny = pknx(i,1);
	pknx(i,0) = _pknx*cos(theta)-_pkny*sin(theta);
	pknx(i,1) = _pknx*sin(theta)+_pkny*cos(theta);
}

extern "C" void move_particlesC(Particles *pdata,
			double *efield,
			int *nextDensity,
			double dt,
			double length_scale,
			void *mesh)
{
	if(pdata->dim == 3)
	{
		move_particles(pdata,efield,
				nextDensity,dt,length_scale,
				(Mesh<kdtree3,3> *)mesh);
		return;
	}
	if(pdata->dim == 2)
	{
		move_particles(pdata,efield,
				nextDensity,dt,length_scale,
				(Mesh<kdtree,2> *)mesh);
		return;
	}
	cout << "Invalid dimension: "<<pdata->dim<<endl;
	exit(0);
}

template<class KD,int dim>
void move_particles(Particles *pdata,
			double *efield,
			int *nextDensity,
			double dt,
			double length_scale,
			Mesh<KD,dim> *mesh)
{
	int i;
	list<int>::iterator end = pdata->p_live->end();
	printf("moving particles now\n");
	
	//avg_momentum_grid(pdata,mesh); //Stats
	for(list<int>::iterator it = pdata->p_live->begin();
				it!= end;++it)
	{
		i = *it;
		pick_up_particle<KD>(i,pdata,nextDensity,mesh); //Lift particle
		randomElectronMovement(pdata->pos,pdata->p_mass,
					pdata->p_id,pdata->p_charge,i,
					efield,dt,length_scale,
					mesh->particle_weight,
					dim);
	}
	printf("Particles Moved\n");
}

void call_me()
{
	printf("I got called!\n");
}

typedef list<double *> Polygon;
//copies original points into new polygon array
//polygon will remain two dimensional but we need to handle
//polytopes and SOON.
extern "C" Polygon *new_polygon(double *points,int n)
{
	double *points_copy = new double[2*n];
	Polygon *polygon = new Polygon;
	for(int i = 0; i < n;i++)
	{
		points_copy[2*i] = points[2*i];
		points_copy[2*i+1] = points[2*i+1];
		polygon->push_back(points_copy+2*i);
	}
	return polygon;
}
/*
template<class KD>
double update_density(Particles *p_data,
			Mesh<KD> *mesh,
			int *nextDensity,
			Polygon *boundary,
			KD *kdt)
{
	cout << "update_density:"<<endl;
	cout << "Dimension ("<<p_data->dim<<") not supported"<<endl;
	exit(0);
	return 0;
}*/


//has_escaped
//returns id of nearest_exit
//or -1 if it has not escaped.
//This makes me unhappy, we would like 0 if not escaped so we can say
//if(has_escaped), but since the id can be zero we can't do it that way.
//the alternative would involve antoher call to nearest_exit, 
//which also makes me unhappy.
template<class KD,int _dim>
int has_escaped(Polytope<_dim> *poly, 
		Particles *p_data, 	
		int part_id, 
		Mesh<KD,_dim> *mesh)
{
	double *particles = p_data->pos;
	int dim = _dim;
	
	if(poly->contains(&pnx(part_id,0)))
		return -1;
	//Check if nearest_exit is reflecting boundary
	int nearest_exit = mesh->find_point_id(&pnx(part_id,0));
	//TODO: Implement is_reflecting and uncomment the following code
	//if(is_reflecting[nearest_exit])
	//{
	//	reflect(poly,p_data,part_id);
	//	return has_escaped(poly,p_data,part_id,mesh);
	//}
	return nearest_exit;
}

/*update_density<kdtree>:
Expected: All particles in nextdensity are 'lifted'
Output: Nextdensity has all particles 'put down'
NOTE: This is excessively complicated.
TODO: Keep track of only the particles.*/
template<class KD,int dim>
double update_density(Particles *ap,
			Mesh<KD,dim> *mesh,
			int *nextDensity,
			Polytope<dim> *boundary,
			KD *kdt)
{
	double *particles = ap->pos; //For Macros
	double current = 0;
	//int *p_id = ap->p_id;
	cout << "Updating Density" << endl;
	for(list<int>::iterator it = ap->p_live->begin();
				it != ap->p_live->end();++it)
	{
		int i = *it;
		int nearest_exit = has_escaped<KD,dim>(boundary,ap,i,mesh);
		if(nearest_exit != -1) //Outside particles are destroyed
		{
			//Nearest exit to determine which side 
			//int nearest_exit = kdtree_find_point_id(kdt,&pos);
			//Check if nearest_exit is reflecting boundary
 			/*Ntype is assumed to be high voltage, 
			so particles leaving it should have opposite their*/
			//current value
			//Ntype is assumed high voltage, so particles move
			int nearest_exit = mesh->find_point_id(&pnx(i,0));
			if(mesh->is_n_type[nearest_exit])	
				current -= (mesh->current_exit(ap,	
						i))
						*ap->p_charge[i];
			else
				current += (mesh->current_exit(ap,
						i))
						*ap->p_charge[i];
			it = destroy_particle(ap,i,it);
			--it;
		}
		else
		{
			put_down_particle(i,ap,nextDensity,mesh);
		}
	}
	cout << "Update current: " << current << endl;
	cout << "Density Updated" << endl;
	return current;
}

extern "C" double update_densityC(Particles *p_data,
			void *mesh,
			int *nextDensity,
			void *boundary,
			void *kdt)
{
	if(p_data->dim == 2)
	{
		return update_density(p_data,
					(Mesh<kdtree,2> *)mesh,
					nextDensity,
					(Polytope<2> *) boundary,
					(kdtree *)kdt);
	}
	if(p_data->dim == 3)
	{
		return update_density(p_data,
					(Mesh<kdtree3,3> *)mesh,
					nextDensity,
					(Polytope<3> *)boundary,
					(kdtree3 *)kdt);
	}
	printf("Invalid Dimension: %d\n",p_data->dim);
	exit(-7);
	return -7;
}

template<class KD, int dim>
double handle_region(int mpos_id,Mesh<KD,dim> *mesh, Particles *p_data,
			int *density, int sign)
{
	cout << "handle_region:"<<endl;
	cout << "Dimension ("<<dim<<") not supported"<<endl;
	exit(-1);
	return 0;
}

template<>
double handle_region(int mpos_id, Mesh<kdtree3,3> *mesh,
			Particles *p_data, int *density, int sign)
{
	double current = 0;
	while (density[mpos_id]*sign < 0) //not charge netural, need more 
	{	
		int i;
		i = create_particle(mpos_id,p_data,density,sign,
				    mesh->materials[mpos_id]->electron_mass,
				    mesh); 
		//ntype is higher voltage, so incoming particles
		//are going the 'right way'
		if(mesh->is_n_type[i])
		{
			//If you are leaving from ntype side
			Vector3f pos(p_data[i][0],p_data[i][1],p_data[i][2]);
			current += mesh->current_exit(p_data,i,
					mesh.nearest_face(pos))
					*sign; 
		}
		//Incoming particles on p side are going the 'wrong' way.
		if(mesh->is_p_type[i])
		{ 
			//leaving from ptype side
			current += mesh->current_exit(p_data,i,
					mesh.nearest_face(pos))
					*sign; 
		}
	}
	return current;
}

template<>
double handle_region(int mpos_id, Mesh<kdtree,2> *mesh, 
			Particles *p_data, int *density, int sign)
{
	double current = 0;
	//If we have too few carriers, inject them
	while (density[mpos_id]*sign < 0) //not charge netural, need more 
	{	
		int i;
		i = create_particle(mpos_id,p_data,density,sign,
				    mesh->materials[mpos_id]->electron_mass,
				    mesh); 
		//ntype is higher voltage, so incoming particles
		//are going the 'right way'
		if(mesh->is_n_type[i])
		{
			//If you are leaving from ntype side
			current -= mesh->current_exit(p_data,i)*sign; 
		}
		//Incoming particles on p side are going the 'wrong' way.
		if(mesh->is_p_type[i])
		{ 
			//leaving from ptype side
			current += mesh->current_exit(p_data,i)*sign; 
		}
	}
	//Not handling excess of charges, need to. 
	//Should probably figure out what failing to do this will result in.
/*
	while (density[mpos_id]*sign > 0) 
	{
		//Need Less, suck or inject opposite?
		//Pretty sure I'm supposed to suck...
		//or does it depend on which side we're on?
		//check bluebook
		if(sign < 0) //too many electrons
		{
			int i = mesh->electrons_pos[mpos_id];
			
		}
		if(sign > 0) //too many holes
		{
			int i = mesh->holes_pos[mpos_id];
			
		}
	}*/
	return current;
}

template<class KD,int dim>double replenish_boundary(Particles *p_data,
			  int *nextDensity, Mesh<KD,dim> *mesh)
{
	list<int>::iterator it = mesh->boundary.begin();
	list<int>::iterator end = mesh->boundary.end();
	double current = 0;
	int sign = 0;
	for(;it != end;++it)
	{
		int id = (*it);
		if(mesh->is_p_type[id])
		{
			sign = 1;//Holes get injected
			current += handle_region(id,mesh,p_data,
						 nextDensity,sign);
		}
		else if(mesh->is_n_type[id])
		{
			sign = -1;//Electrons get injected
			current += handle_region(id,mesh,p_data,
						nextDensity,sign);
		}else
		{
			cout << id << "Danger will robinson!" << endl;
		}
	}
	return current;
}


extern "C" double replenishC(Particles *p_data, 
			int *nextDensity, 
			void *mesh)

{
	if(p_data->dim == 2)
	{
		return replenish(p_data,nextDensity,(Mesh<kdtree,2> *) mesh);
	}
	if(p_data->dim == 3)
	{
		return replenish(p_data,nextDensity,(Mesh<kdtree3,3> *) mesh);
	}
	printf("Invalid Dimension: %d\n",p_data->dim);
	exit(-7);
}

template<class KD,int dim>
double replenish(Particles *p_data, int *nextDensity, Mesh<KD,dim> *mesh)
{	
	double current;
	cout << "Beginning replenish" << endl;
	current = replenish_boundary(p_data,nextDensity,mesh);
	cout << "replenish complete" << endl;
	return current;
}

template <class KD,int dim>
double recombinate(Particles *p_data, int *nextDensity, Mesh<KD,dim> *mesh)
{
	//For each cell, determine the number recombinations
	//delete particles from each
	//remember to remove both charge types from their position
	//Should be same
	
	for(int mesh_pos = 0; mesh_pos < mesh->npoints;mesh_pos++)
	{
	
		while( mesh->electrons_pos[mesh_pos].size() > 0 && 
		    mesh->holes_pos[mesh_pos].size() > 0)
		{
			//For now, obliterate the first two.
			int e_id = *(mesh->electrons_pos[mesh_pos].begin());
			int h_id = *(mesh->holes_pos[mesh_pos].begin());
			
			//electron
			pick_up_particle<KD>(e_id,p_data,nextDensity,mesh);
			destroy_particle(p_data,e_id,p_data->live_id[e_id]);
		
			//hole
			pick_up_particle<KD>(h_id,p_data,nextDensity,mesh);
			destroy_particle(p_data,h_id,p_data->live_id[h_id]);
		}
	}
	return 0;
}
/*recombinate:
Expected to be called after update_density, 
so need to pick up particles particles
before destruction*/
extern "C" double recombinateC(Particles *p_data, int *nextDensity, void *mesh)
{
	if(p_data->dim == 2)
	{
		return recombinate(p_data,nextDensity,(Mesh<kdtree,2> *) mesh);
	}
	if(p_data->dim == 3)
	{
		return recombinate(p_data,nextDensity,(Mesh<kdtree3,3> *) mesh);
	}
	printf("Invalid Dimension: %d\n",p_data->dim);
	exit(-7);
	return -7;
}
