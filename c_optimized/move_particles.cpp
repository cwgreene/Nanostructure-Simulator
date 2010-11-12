//Standard Library Includes
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <list>
#include <iostream>

//3rd Party
#include <Eigen/Core>

//HPP Includes
#include "mesh.hpp"
#include "particles.hpp"
#include "statistics.hpp"
#include "move_particles.hpp"
#include "Polytope.hpp"
#include "debug.hpp"
#include "myrand.hpp"


extern "C" {
#include "kdtree.h"
#include "kdtree3.h"
}
//Macro Definitions
#define sq(x) (x*x)
#define EC (1.60217646e-19)
using namespace std;
using namespace debug;

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

void scatter(double *particles, int dim,int i)
{
	double theta = rand()*(3.1415192)/RAND_MAX-3.141592/2;

	//handle x and y cooridinates.
	//this should be a bandstructure lookup
	double _pknx = pknx(i,0);
	double _pkny = pknx(i,1);
	pknx(i,0) = _pknx*cos(theta)-_pkny*sin(theta);
	pknx(i,1) = _pknx*sin(theta)+_pkny*cos(theta);
}

/**
 * randomElectronMovement
 *
 * movement and mommentum change for a single particle.
 * 
 * note: 02 optimization hopefully removes the overhead
 * of calling a function a few hundred thousand times
 * but if not, we can see what happens if we expliciyly
 * inline it. Possibly not a bad idea. As far as I know
 * only two functions actually call this.
 */
template<class KD,int dim>
void randomElectronMovement(Particles *p_data,
				int i,
				double dt,double length_scale,
				double *efield,
				Mesh<KD,dim> *mesh)
{
	double *particles = p_data->pos;
	double *p_mass = p_data->p_mass;
	int *p_id = p_data->p_id;
	int *p_charge = p_data->p_charge;
	double particle_weight=mesh->particle_weight;
	int mpos_id = p_data->p_id[i];

	double dx[3] = {0.,0.,0.};

	//Drift
	for(int c = 0; c < dim; c++)
	{
		pknx(i,c) += (efield[dim*p_id[i]]*p_charge[i]*dt)
				*EC*particle_weight;
	
	}
	//Move
	for(int c = 0; c < dim; c++)	
	{
#ifdef TRACK
		if(i == grabbed)
		{
			cerr<<"x_"<<c<<": "<< pnx(i,c) << endl;
			cerr<<"particle weight: "<< mesh->particle_weight<<endl;
		}
#endif 
		dx[c] = pknx(i,c)*dt/(length_scale*p_mass[i]*particle_weight);
		pnx(i,c) += dx[c];
		//update distance approximately
		//since dx and dy are small
		//If they aren't, probably dist > scattering length anyway
		p_data->p_dist[i] += (dx[c])*length_scale;
#ifdef TRACK
		if(i==grabbed)
		{
			cerr << "dx_"<<c<<": "<<dx[c]<<endl;
			//debug::dbg<< "x_"<<c<<": "<<pnx(i,c) << endl;
		}
#endif
	}
	
	//Scatter... ugh... this is going to be complicated.
	//we're going to need to seperate drift and diffusion
	//we already knew that... but damn.
	if(p_data->p_dist[i] > mesh->materials[mpos_id]->free_dist)
	{
		p_data->p_dist[i] = 0;
	}
}
/**
 * move_particles
 * 
 * External C Interface.
 * This times steps the entire simulated system by dt
 */
extern "C" void move_particlesC(Particles *p_data,
			double *efield,
			int *nextDensity,
			double dt,
			double length_scale,
			void *mesh)
{
	if(p_data->dim == 3)
	{
		/*move_particles(p_data,efield,
				nextDensity,dt,length_scale,
				(Mesh<kdtree3,3> *)mesh);*/
		std::cout << "Don't use 3 dimensions for now. Boundary code is stupid. Fix it.\n";
		exit(-1);
		return;
	}
	if(p_data->dim == 2)
	{
		move_particles(p_data,efield,
				nextDensity,dt,length_scale,
				(Mesh<kdtree,2> *)mesh);
		return;
	}
	//debug::dbg << "Invalid dimension: "<<pdata->dim<<endl;
	exit(0);
}

template<class KD,int dim>
void move_particles(Particles *p_data,
			double *efield,
			int *nextDensity,
			double dt,
			double length_scale,
			Mesh<KD,dim> *mesh)
{
	int i;
	list<int>::iterator end = p_data->p_live->end();
	printf("moving particles now\n");

	//avg_momentum_grid(pdata,mesh); //Stats
	for(list<int>::iterator it = p_data->p_live->begin();
				it!= end;++it)
	{
		double old_pos[dim];
		i = *it;

		for(int c = 0; c < dim;c++)
			old_pos[c] = p_data->pos[2*dim*i];

		//double *particles = pdata->pos;
		pick_up_particle<KD>(i,p_data,nextDensity,mesh); //Lift particle
		randomElectronMovement<KD,dim>(p_data,i,
					dt,length_scale,
					efield,
					mesh);
		p_data->p_id[i] = mesh->find_point_id(p_data->pos + 2*dim*i); //find nearest spot
		//After random movement, check to see if the nearest point
		//Is a reflecting boundary. If it is, we reflect
		int loc = mesh->find_point_id(p_data->pos+2*dim*i);
//		if(mesh->is_reflect[loc])
//		{
//			mesh->reflect(p_data,i,old_pos);
//		}
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
	int x;
		
	if(poly->contains(p_data+2*dim*part_id,0))
		return -1;
	//Check if nearest_exit is reflecting boundary
	int nearest_exit = p_data->p_id[part_id];
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
double update_density(Particles *p_data,
			Mesh<KD,dim> *mesh,
			int *nextDensity,
			KD *kdt)
{
	double current = 0;
	//int *p_id = ap->p_id;
	//debug::dbg << "Updating Density" << endl;
	for(list<int>::iterator it = p_data->p_live->begin();
				it != p_data->p_live->end();++it)
	{
		int i = *it;
		int nearest_exit = mesh->has_escaped(p_data,i);
		if(nearest_exit != -1) //Outside particles are destroyed
		{
			//Nearest exit to determine which side 
			//int nearest_exit = kdtree_find_point_id(kdt,&pos);
			//Check if nearest_exit is reflecting boundary
 			/*Ntype is assumed to be pos voltage, 
			so particles leaving it should have opposite their*/
			//charge value
			//Ntype is assumed high voltage, so particles move
			int nearest_exit = mesh->find_point_id(p_data->pos+2*dim*i);
			if(mesh->is_n_type[nearest_exit])	
				current -= (mesh->current_exit(p_data,	
						i))*EC
						*p_data->p_charge[i];
			else
				current += (mesh->current_exit(p_data,
						i))*EC
						*p_data->p_charge[i];
			it = destroy_particle(p_data,i,it);
			--it;
		}
		else
		{	
			//debug::dbg <<"Placing particle Down"<<endl;
			put_down_particle(i,p_data,nextDensity,mesh);
		}
	}
	//debug::dbg << "Update current: " << current << endl;
	//debug::dbg << "Density Updated" << endl;
	cerr << current << endl;

	//Invariants should be satisfied.
	return current;
}

extern "C" double update_densityC(Particles *p_data,
			void *mesh,
			int *nextDensity,
			void *kdt)
{
	if(p_data->dim == 2)
	{
		return update_density(p_data,
					(Mesh<kdtree,2> *)mesh,
					nextDensity,
					(kdtree *)kdt);
	}
	if(p_data->dim == 3)
	{
		return update_density(p_data,
					(Mesh<kdtree3,3> *)mesh,
					nextDensity,
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
	//debug::dbg << "handle_region:"<<endl;
	//debug::dbg << "Dimension ("<<dim<<") not supported"<<endl;
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
		i = create_particle(mpos_id,p_data,density,-sign,
				    mesh->materials[mpos_id]->electron_mass,
				    mesh); 
		//ntype is higher voltage, so incoming particles
		//are going the 'right way'
		if(mesh->is_n_type[i])
		{
			//If you are leaving from ntype side
			current += mesh->current_exit(p_data,i,
					mesh->nearest_edge(p_data->pos))
					*sign; 
		}
		//Incoming particles on p side are going the 'wrong' way.
		if(mesh->is_p_type[i])
		{ 
			//leaving from ptype side
			current += mesh->current_exit(p_data,i,
					mesh->nearest_edge(p_data->pos))
					*sign; 
		}
	}
	return current;
}

template<class KD,int dim>
void test_mesh_failure(Mesh<KD,dim> *mesh, int *density)
{
	for(int mpos_id = 0; mpos_id < mesh->npoints;mpos_id++)
	{
		if(mesh->is_n_type[mpos_id])
			test_sum_failure(mpos_id,mesh,density,1);
		else if(mesh->is_p_type[mpos_id])
			test_sum_failure(mpos_id,mesh,density,-1);
		else
		{
			cout << "POINT IS NEITHER NTYPE OR PTYPE"<<endl;
			exit(-1);
		}
	}
}

template<class KD, int dim>
void mesh_point_info(int mpos_id, Mesh<KD,dim> *mesh, int *density)
{
	int empty_sign;
	if(mesh->is_n_type[mpos_id])
	{
		empty_sign = 1;
	}else
	{
		empty_sign = -1;
	}
	
	cout << "holes, electrons, sign, density: "<<endl;
	cout << "h:" << mesh->holes_pos[mpos_id].size()<< " " ;
	cout << "e:" << mesh->electrons_pos[mpos_id].size() << " ";
	cout << "s:" << empty_sign << " ";
	cout << "d:"<< density[mpos_id] << endl;
}

template<class KD,int dim>
void test_sum_failure(int mpos_id, Mesh<KD,dim> *mesh, int *density, 
			int empty_sign)
{
	mesh_point_info(mpos_id,mesh,density);
	int holes = mesh->holes_pos[mpos_id].size();
	int electrons = mesh->electrons_pos[mpos_id].size();
	if((mesh->gen_num*empty_sign+(holes-electrons)) != density[mpos_id])
	{
		
		cout << "INSUFFICIENT DENSITY"<<endl;
		int x = mpos_id;
		do
		{
			mesh_point_info(x,mesh,density);
			cout << "Query Point"<<endl;
			cin >> x;
		}while(x >= 0 && x < mesh->npoints);
		exit(-1);
	}
}

template<>
double handle_region(int mpos_id, Mesh<kdtree,2> *mesh, 
			Particles *p_data, int *density, int empty_sign)
{
	double current = 0;

	//Things that can cause an imbalance
	//Too many carriers of one type
	//too few carriers of another

	//If the empty_sign of the region is the same as
	//sign as the density, then we're empty, and we need more particles
	while (density[mpos_id]*empty_sign > 0) //not charge netural, need more
	{
		int i = create_particle(mpos_id,p_data,density,-empty_sign,
				    mesh->materials[mpos_id]->electron_mass,
				    mesh); 
	
		//ntype is higher voltage, so incoming particles
		//are going the 'right way'
		if(mesh->is_n_type[mpos_id])
		{
			//Currently positive
			current -= mesh->current_exit(p_data,i)*empty_sign; 
		}
		//Incoming particles on p side are going the 'wrong' way.
		if(mesh->is_p_type[mpos_id])
		{ 
			//leaving from ptype side
			current += mesh->current_exit(p_data,i)*empty_sign; 
		}
	}
	//If the sign of the density is different than the empty_sign
	//of the region that means we have an excess of minority
	//carriers
	//In a p_region (empty means -100), this means we have
	while (density[mpos_id]*empty_sign < 0) 
	{	
		int i = create_particle(mpos_id,p_data,density,empty_sign,
				    mesh->materials[mpos_id]->electron_mass,
				    mesh); 
		//Need Less, suck or inject opposite?
		//Pretty sure I'm supposed to suck...
		//or does it depend on which side we're on?
		//check bluebook
		list<int>::iterator doomed;
		if(density[mpos_id] > 0) //Too many holes
		{
			//doomed = mesh->holes_pos[mpos_id].begin();
			current -= mesh->current_exit(p_data,i)*empty_sign; 
		}
		else if(density[mpos_id] < 0) //Too many electrons
		{
			//doomed = mesh->electrons_pos[mpos_id].begin();
			current += mesh->current_exit(p_data,i)*empty_sign; 
		}

//		int i = *doomed;
//		pick_up_particle<kdtree>(i,p_data,density,mesh);
//		destroy_particle(p_data,i,doomed);
	}
	return current;
}

template<class KD,int dim>double replenish_boundary(Particles *p_data,
			  int *nextDensity, Mesh<KD,dim> *mesh)
{
	vector<int>::iterator it = mesh->boundary.begin();
	vector<int>::iterator end = mesh->boundary.end();
	double current = 0;
	int sign = 0;
	for(;it != end;++it)
	{
		int id = (*it);
		if(mesh->is_p_type[id])
		{
			sign = -1; //empty is negative
			current += handle_region(id,mesh,p_data,
						 nextDensity,sign);
		}
		else if(mesh->is_n_type[id])
		{
			sign = 1;//empty is positive
			current += handle_region(id,mesh,p_data,
						nextDensity,sign);
		}else
		{
			//debug::dbg << id << "Danger will robinson!" << endl;
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
	//debug::dbg << "Beginning replenish" << endl;
	current = replenish_boundary(p_data,nextDensity,mesh);
	//debug::dbg << "replenish complete" << endl;
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
		int destroy_count = 0;
//		list<int> *electrons = &(mesh->electrons_pos[mesh_pos]);
//		list<int> *holes = &(mesh->holes_pos[mesh_pos]);
		for(list<int>::iterator it = mesh->electrons_pos[mesh_pos].begin();
		    it != mesh->electrons_pos[mesh_pos].end();
		    ++it)
		{
			if(mesh->holes_pos[mesh_pos].size() == 0)
				break;
			if(chance_recombinate(mesh->holes_pos[mesh_pos].size()))
			{
				destroy_count++;
			}
		}
		while(mesh->electrons_pos[mesh_pos].size() >0 && 
		      mesh->holes_pos[mesh_pos].size() > 0 && 
		      destroy_count > 0)
		{		
				//debug::dbg << mesh_pos << endl;
				//int x;
				//debug::dbg << "size:"<<mesh->electrons_pos[mesh_pos].size()<<endl;
				//debug::dbg << "size:"<<mesh->holes_pos[mesh_pos].size()<<endl;
				//For now, obliterate the last two.
				int e_id = *(mesh->electrons_pos[mesh_pos].begin());
				int h_id = *(mesh->holes_pos[mesh_pos].begin());
			
				//electron
				pick_up_particle<KD>(e_id,p_data,nextDensity,mesh);
				destroy_particle(p_data,e_id,p_data->live_id[e_id]);
			
				//hole
				pick_up_particle<KD>(h_id,p_data,nextDensity,mesh);
				destroy_particle(p_data,h_id,p_data->live_id[h_id]);
			
				//decrease needed destroys
				destroy_count--;
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

bool chance_recombinate(int other_type)
{
	if(randint(1,100)> 80)
		return true;
	return false;
}

bool chance_scatter()
{
	if(randint(1,100)>90)
		return true;
	return false;
}

template <class KD,int dim>
double photorecombinate(Particles *p_data, int *nextDensity, Mesh<KD,dim> *mesh)
{
	list<int>::iterator end = p_data->p_live->end();
	vector<list<int>::iterator> doomed;
	for(list<int>::iterator it = p_data->p_live->begin();
		it !=end;++it)
	{
		int p_id = *it;
		int i = p_id; //for uniformity of offset calc
		int m_id = mesh->find_point_id(p_data->pos+2*dim*i);

		if(p_data->p_charge[p_id] < 0)
		{
			if(chance_recombinate(mesh->holes_pos[m_id].size()))
			{
				doomed.push_back(it);
			}
		}
		else if (p_data->p_charge[p_id] > 0)
		{
			if(chance_recombinate(mesh->electrons_pos[m_id].size()))
			{
				doomed.push_back(it);
			}
		}
	}
	for(unsigned int i = 0; i < doomed.size();i++)
	{
		destroy_particle(p_data,i,doomed[i]);
	}
	return 0;
}
template<class KD,int dim>
void photo_move_particles(Particles *p_data,
			double *efield,
			int *nextDensity,
			double dt,
			double length_scale,
			Mesh<KD,dim> *mesh)
{
	int i;
	list<int>::iterator end = p_data->p_live->end();

	//avg_momentum_grid(pdata,mesh); //Stats
	for(list<int>::iterator it = p_data->p_live->begin();
				it!= end;++it)
	{
		i = *it;
		//only difference, don't deal with pick_up_put_down
		randomElectronMovement(p_data,i,
					dt,length_scale,
					efield,
					mesh);
	}
}

template<class KD, int dim>
double photo_exit_current(Particles *p_data, Mesh<KD,dim> *mesh)
{
	double total = 0;
	list<int>::iterator end = p_data->p_live->end();
	vector<list<int>::iterator > doomed;
	for(list<int>::iterator it = p_data->p_live->begin();
		it != end;++it)
	{
		if(mesh->has_escaped(p_data,*it))
		{
			doomed.push_back(it);
			total += 1;
		}
	}
	for(unsigned int i = 0; i < doomed.size();i++)
		destroy_particle(p_data,*(doomed[i]),doomed[i]);
	return total;
}

template<class KD,int dim>
double photocurrent(Particles *p_data,int *density, double *efield, Mesh<KD,dim> *mesh,
			double dt, double length_scale)
{
	//Mostly carbon copies of original
	//except none of that (pick up/put down) nonsense
	photo_move_particles(p_data,efield,density,dt,length_scale,mesh);
	photorecombinate(p_data,density,mesh);
	return photo_exit_current(p_data,mesh);
}

/*Photocurrent
takes electric field, density, and electric field, and mesh
*/
extern "C" double photocurrentC(Particles *p_data, int *density, double *efield,void *mesh, double dt,
				double length_scale)
{
	if(p_data->dim == 3)
	{
		//Move, check for recombination, check for exit.
		return photocurrent(p_data,density,efield,(Mesh<kdtree3,3> *)mesh,dt,length_scale);
	}
	if(p_data->dim == 2)
	{
		return photocurrent(p_data,density,efield,(Mesh<kdtree,2> *)mesh,dt,length_scale);
	}
	return 0;
}
