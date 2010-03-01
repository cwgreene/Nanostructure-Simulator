#include <stdio.h>
#include <math.h>
#include <list>
#include <iostream>
#include "particles.hpp"

extern "C" {
#include "kdtree.h"
}
//Macros
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

extern "C" list<int> *init_particle_list(int n)
{
	list<int> *particle_list = new list<int>;
	for(int i = 0; i < n;i++)
	{
		particle_list->push_back(i);
	}
	cout <<"Location " <<particle_list << endl;
	return particle_list;
}

extern "C" list<int> *init_dead_list()
{
	return new list<int>;
}

void randomElectronMovement(double *particles,
				double *p_mass,
				int *p_id,
				int *p_charge,
				int i,
				double *efield,
				double dt,
				double length_scale)
{
	//Move
	double dx = pkx(i)*dt/length_scale/p_mass[i];
	double dy = pky(i)*dt/length_scale;
	px(i) += dx;
	py(i) += dy;
	//Drift
	pkx(i) += efield[2*p_id[i]]*p_charge[i]*dt/length_scale;
	pky(i) += efield[2*p_id[i]+1]*p_charge[i]*dt/length_scale;
	//Scatter
}

extern "C" void move_particles(double *particles,
			double *p_mass,
			int *p_charge,
			int *p_id,
			list<int> *p_live,
			double *efield,
			double *nextDensity,
			double dt,
			double length_scale)
{
	int i;
	printf("moving particles now\n");
	cout <<p_live->size() << endl;
	for(list<int>::iterator it = p_live->begin();it!=p_live->end();++it)
	{
		i = *it;
		int id = p_id[i];
		int charge = p_charge[i];
		nextDensity[id] -= charge;
		randomElectronMovement(particles,p_mass,p_charge,p_id,i,
					efield,dt,length_scale);
	}
	printf("Particles Moved\n");
}

void call_me()
{
	printf("I got called!\n");
}

double current_exit(double *particles,int i)
{
	double _pkx = pkx(i);
	double _pky = pky(i);
	return sqrt(_pkx*_pkx+_pky*_pky);
}

typedef list<double *> Polygon;

//copies original points into new polygon array
extern "C" Polygon *new_polygon(double *points,int n)
{
	double *points_copy = new double[2*n];
	Polygon *polygon = new Polygon;
	for(int i = 0; i < n;i++)
	{
		points_copy[2*i] = points[2*i];
		points_copy[2*i+1] = points[2*i];
		polygon->push_back(points_copy+2*i);
	}
	return polygon;
}

bool point_in_polygon(Polygon *boundary, vector2 *pos)
{
	double *point;
	list<double *>::iterator it = boundary->begin();
	list<double *>::iterator end = boundary->end();

	double prev_x = (*it)[0];
	double prev_y = (*it)[1];
	double prev_vec_x = prev_x-(*pos)[0];
	double prev_vec_y = prev_y-(*pos)[1];
	++it;
	for(;it != end;++it)
	{
		point = *it;
		double cur_x = point[0];
		double cur_y = point[1];
	
		double cur_vec_x = cur_x-(*pos)[0];
		double cur_vec_y = cur_y-(*pos)[1];

		double cross_z = (prev_vec_x)*(cur_vec_y)-
				 (cur_vec_x)*(prev_vec_y);
		if( cross_z < 0)
		{
			return false;
		}
	}
	return true;
}

extern "C" double update_density(Particles *ap,
			int *p_id, int *p_charge, int *nextDensity,
			list<int> *p_live,list<int> *p_dead,
			Polygon *boundary,
			kdtree *kdt)
{
	double *particles = ap->pos;
	double current = 0;
	cout << "Updating Density" << endl;
	for(list<int>::iterator it = ap->p_live->begin();
				it != ap->p_live->end();++it)
	{
		int i = *it;
		vector2 pos;
		pos[0] = px(i);
		pos[1] = py(i);
		if(point_in_polygon(boundary,&pos))
		{
			current += current_exit(particles,i);
			it = ap->p_live->erase(it);
			ap->p_dead->push_back(i);
			--it;
		}
		else
		{
			p_id[i] = kdtree_find_point_id(kdt,&pos);
			nextDensity[p_id[i]] += ap->p_charge[i];
		}
	}
	cout << "Density Updated" << endl;
	return current;
}

double random_energy(double Temperature)
{
	
}

double handle_region(int mpos_id, Mesh *mesh, Particles *p_data, 
			int *density, int sign)
{
	double current;
	//If we have too few carriers, inject them
	while (denisty[id]*sign < 0) //not charge netural, need more 
	{	
		int i;
		energy = random_energy(300.);
		i = new_particle(mpos_id,p_data,density,sign);
		current += current_exit(p_data->p_pos,i);
	}
	return current;
}

double replenish_boundary(list<Point *> *boundary,Particles p_data,
			  double *nextDensity, Mesh *mesh)
{
	list<Point *>::iterator it = boundary->begin();
	list<Point *>::iterator it = boundary->end();
	double current = 0;
	for(;it != end;++it)
	{
		if(point->type == P_TYPE)
		{
			sign = +1;//Holes get injected
			current += handle_region(point->mpos_id,point->mpos);
		}else
		{
			sign = -1;//Electrons get injected
			current += handle_region(point->mpos_id,point->mpos);
		}
	}
	return current;
}

extern "C" replenish(list<int> *p_live, list<int> *p_dead)
{
	replenish_current(
}
