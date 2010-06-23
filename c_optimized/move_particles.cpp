#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <list>
#include <iostream>
#include "particles.hpp"
#include "mesh.hpp"
#include "statistics.hpp"

extern "C" {
#include "kdtree.h"
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
		dx[c] = pknx(i,c)*dt/(length_scale*p_mass[i]*particle_weight);
		pnx(i,c) += dx[c];
	}
	//Drift
	for(int c = 0; c < dim; c++)
		pknx(i,c) += (efield[dim*p_id[i]]*p_charge[i]*dt)
				*EC*particle_weight;

	//Scatter... ugh... this is going to be complicated.
	//we're going to need to seperate drift and diffusion
	//we already knew that... but damn.
	double theta = rand()*2*(3.1415192)/RAND_MAX;

	//handle x and y cooridinates.
	//this should be a bandstructure lookup
	pknx(i,0) = pknx(i,0)*cos(theta)-pknx(i,1)*sin(theta);
	pknx(i,1) = pknx(i,0)*sin(theta)+pknx(i,1)*cos(theta);
}

extern "C" void move_particles(Particles *pdata,
			double *efield,
			int *nextDensity,
			double dt,
			double length_scale,
			Mesh<kdtree> *mesh)
{
	int i;
	list<int>::iterator end = pdata->p_live->end();
	printf("moving particles now\n");
	
	//avg_momentum_grid(pdata,mesh); //Stats
	for(list<int>::iterator it = pdata->p_live->begin();
				it!= end;++it)
	{
		i = *it;
		pick_up_particle(i,pdata,nextDensity,mesh); //Lift particle
		randomElectronMovement(pdata->pos,pdata->p_mass,
					pdata->p_id,pdata->p_charge,i,
					efield,dt,length_scale,
					mesh->particle_weight,
					mesh->dim);
	}
	printf("Particles Moved\n");
}

void call_me()
{
	printf("I got called!\n");
}

double current_exit2(double *particles,int i,double mass)
{
	int dim = 2;
	double _pkx = pknx(i,0);
	double _pky = pknx(i,1);
	double current = sqrt(_pkx*_pkx+_pky*_pky)/mass;//velocity
	return current*EC;
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

double cross_segment_point(double seg1_x,double seg1_y,
			   double seg2_x,double seg2_y,
			   double point_x,double point_y)
{
	double segment_x = seg1_x-seg2_x;
	double segment_y = seg1_y-seg2_y;
	double cur_vec_x = seg1_x-point_x;
	double cur_vec_y = seg1_y-point_y;
	double cross_z = (segment_x)*(cur_vec_y)-
			 (cur_vec_x)*(segment_y);
	return cross_z;
}

bool point_in_polygon(Polygon *boundary, vector2 *pos)
{
	double *point;
	list<double *>::iterator it = boundary->begin();
	list<double *>::iterator end = boundary->end();

	double cur_x,cur_y;
	double prev_x = (*it)[0];
	double prev_y = (*it)[1];
	double first_x = prev_x;
	double first_y = prev_y;
	double p_x = (*pos)[0];
	double p_y = (*pos)[1];
	++it;
	double sign = 0;
	for(;it != end;++it)
	{
		point = *it;
		cur_x = point[0];
		cur_y = point[1];
		double cross_z = cross_segment_point(cur_x,cur_y,
						     prev_x,prev_y,
						     p_x,p_y);
		if(sign == 0)
		{
			sign = cross_z;
		}else if( cross_z*sign < 0)
		{
			return false;
		}
		prev_x = cur_x;
		prev_y = cur_y;
	}
	//Last with first
	if(cross_segment_point(first_x,first_y,
				prev_x,prev_y,
				p_x,p_y)*sign < 0)
		return false;
	
	return true;
}


/*update_density:
Expected: All particles in nextdensity are 'lifted'
Output: Nextdensity has all particles 'put down'
NOTE: This is excessively complicated.
TODO: Keep track of only the particles.*/
extern "C" double update_density2(Particles *ap,
			Mesh<kdtree> *mesh,
			int *nextDensity,
			Polygon *boundary,
			kdtree *kdt)
{
	int dim = 2; //For macros
	double *particles = ap->pos;
	double current = 0;
	//int *p_id = ap->p_id;
	cout << "Updating Density" << endl;
	for(list<int>::iterator it = ap->p_live->begin();
				it != ap->p_live->end();++it)
	{
		int i = *it;
		vector2 pos;
		pos[0] = pnx(i,0);
		pos[1] = pnx(i,1);
		if(!point_in_polygon(boundary,&pos)) //Outside particles are destroyed
		{	
			//Nearest exit to determine which side 
			int nearest_exit = kdtree_find_point_id(kdt,&pos);
 			/*Ntype is assumed to be high voltage, 
			so particles leaving it should have opposite their*/
			//current value
			//Ntype is assumed high voltage, so particles move
			if(mesh->is_n_type[nearest_exit])	
				current -= current_exit2(particles,i,ap->p_mass[i])*ap->p_charge[i];
			else
				current += current_exit2(particles,i,ap->p_mass[i])*ap->p_charge[i];
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

double handle_region(int mpos_id, Mesh<kdtree> *mesh, Particles *p_data, 
			int *density, int sign)
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
			current += current_exit2(p_data->pos,i,p_data->p_mass[i])*sign; //If you are leaving from ntype side
		//Incoming particles on p side are going the 'wrong' way.
		if(mesh->is_p_type[i])
			current += current_exit2(p_data->pos,i,p_data->p_mass[i])*sign; //leaving from ptype side
	}
	//Not handling too many, need to. Should probably figure out what failure to do this will resul tin.
	return current;
}

double replenish_boundary2(Particles *p_data,
			  int *nextDensity, Mesh<kdtree> *mesh)
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
			current -= handle_region(id,mesh,p_data,
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

extern "C" double replenish2(Particles *p_data, int *nextDensity, Mesh<kdtree> *mesh)
{	
	double current;
	cout << "Beginning replenish" << endl;
	current = replenish_boundary2(p_data,nextDensity,mesh);
	cout << "replenish complete" << endl;
	return current;
}

/*recombinate:
Expected to be called after update_density, so need to pick up particles particles
before destruction*/
extern "C" double recombinate2(Particles *p_data, int *nextDensity, Mesh<kdtree> *mesh)
{
	//For each cell, determine the number recombinations
	//delete particles from each
	//remember to remove both charge types from their position
	//Should be same
	
	for(int mesh_pos = 0; mesh_pos < mesh->npoints;mesh_pos++)
	{
		/*cout << "num particles:" << 
			mesh->electrons_pos[mesh_pos].size()+mesh->holes_pos[mesh_pos].size() << 
			"/" << mesh->electrons_pos[mesh_pos].size() <<" "<<
			" dens: "<<nextDensity[mesh_pos]<<endl;*/
		while( mesh->electrons_pos[mesh_pos].size() > 0 && 
		    mesh->holes_pos[mesh_pos].size() > 0)
		{
			//For now, obliterate the first two.
			int e_id = *(mesh->electrons_pos[mesh_pos].begin());
			int h_id = *(mesh->holes_pos[mesh_pos].begin());
			//cout << "Recombinate: "<<e_id <<" "<<h_id << endl;
			pick_up_particle(e_id,p_data,nextDensity,mesh);//electron
			destroy_particle(p_data,e_id,p_data->live_id[e_id]);
			pick_up_particle(h_id,p_data,nextDensity,mesh);//hole
			destroy_particle(p_data,h_id,p_data->live_id[h_id]);
		}
	}
	return 0;
}
