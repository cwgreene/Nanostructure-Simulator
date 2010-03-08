#include <stdio.h>
#include <math.h>
#include <list>
#include <iostream>
#include "particles.hpp"
#include "mesh.hpp"

extern "C" {
#include "kdtree.h"
}
//Macros
#define EC (1.60217646e-19)
using namespace std;

extern "C" list<int> *init_particle_list(int n)
{
	list<int> *particle_list = new list<int>;
	for(int i = 0; i < n;i++)
	{
		particle_list->push_back(i);
	}
	//cout <<"p_live Location " <<particle_list << endl;
	return particle_list;
}

extern "C" list<int> *init_dead_list(int start,int n)
{
	list<int> *particle_list = new list<int>;
	for(int i = start; i < n;i++)
	{
		particle_list->push_back(i);
	}
//	cout <<"pdead Location " <<particle_list << endl;
	return particle_list;
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
	double dx = pkx(i)*dt/p_mass[i];
	double dy = pky(i)*dt/p_mass[i];
	double old_pkx = pkx(i);
	double old_pky = pky(i);
	if(i % 1000 == 0)
	//cout << "x:" << px(i) << " y: " << py(i) <<
	//	" dx: " << dx <<  " dy: " << dy << " i: "<<i<<endl;
	//cout << "dx: "<< dx << " dy: " << dy;
	px(i) += dx;
	py(i) += dy;
	//Drift
	pkx(i) += (efield[2*p_id[i]]*p_charge[i]*dt/length_scale)
			*EC;
	pky(i) += (efield[2*p_id[i]+1]*p_charge[i]*dt/length_scale)
			*EC;
	//if(old_pkx != 0 && old_pky != 0 && i%1000 == 0)
	//	cout << p_id[i]<<" field:"<< efield[2*p_id[i]] <<" " <<efield[2*p_id[i]+1]<< " delta_p: "<< 
	//	(pkx(i)*dt/p_mass[i])<<"/"<<(old_pkx *dt/p_mass[i]) <<" "<<
	//	(pky(i)*dt/p_mass[i])<<"/"<<(old_pky *dt/p_mass[i])<< endl;

			/*	if(i %1000 == 0)
		cout << "after:" <<endl << 
		"x:" << px(i) << " y: " << px(i) <<
		" dx: " << dx <<  " dy: " << dy << 
		" pkx: " << pkx(i) << " pky: "<<pky(i)<<endl;*/

	//Scatter
}

extern "C" void move_particles(Particles *pdata,
			double *efield,
			int *nextDensity,
			double dt,
			double length_scale)
{
	int i;
	list<int>::iterator end = pdata->p_live->end();
	printf("moving particles now\n");
	
	for(list<int>::iterator it = pdata->p_live->begin();
				it!= end;++it)
	{
		//cout << "doom!" << endl;
		i = *it;
		int id = pdata->p_id[i];
		int charge = pdata->p_charge[i];
//		cout << "charge: "<<charge << endl;
		nextDensity[id] -= charge;
//		cout << "density = " << nextDensity[id];
		randomElectronMovement(pdata->pos,pdata->p_mass,
					pdata->p_id,pdata->p_charge,i,
					efield,dt,length_scale);
	//	nextDensity[id] = 3;
	//	cout << id << endl;
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
	double current = sqrt(_pkx*_pkx+_pky*_pky);
	return current;
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
		points_copy[2*i+1] = points[2*i+1];
		polygon->push_back(points_copy+2*i);
	}
	cout << "Copied Points" << endl;
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

extern "C" double update_density(Particles *ap,
			int *nextDensity,
			Polygon *boundary,
			kdtree *kdt)
{
	double *particles = ap->pos;
	double current = 0;
	int *p_id = ap->p_id;
	cout << "Updating Density" << endl;
	for(list<int>::iterator it = ap->p_live->begin();
				it != ap->p_live->end();++it)
	{
		int i = *it;
		vector2 pos;
		pos[0] = px(i);
		pos[1] = py(i);
		if(!point_in_polygon(boundary,&pos))
		{	
			/*cout <<i << " died " << endl
			  <<"x: "   <<pos[0] << " y: " <<pos[1] <<endl;*/
			current += current_exit(particles,i);
			it = ap->p_live->erase(it);
			ap->p_dead->push_back(i);
			--it;
		}
		else
		{
			p_id[i] = kdtree_find_point_id(kdt,&pos);
			//cout << "id" << p_id[i] << endl;
/*			cout << "x: "<<pos[0]<<" y: "
				<<pos[1]<<" id: "<<p_id[i] << endl;*/
			nextDensity[p_id[i]] += ap->p_charge[i];
		}
	}
	cout << "Update current: " << current << endl;
	cout << "Density Updated" << endl;
	return current;
}

double handle_region(int mpos_id, Mesh *mesh, Particles *p_data, 
			int *density, int sign)
{
	double current = 0;
	//If we have too few carriers, inject them
//	cout << "handle_region density(b):" << mpos_id<< " " <<
//			density[mpos_id]<<" sign: "<<sign<< endl;
	while (density[mpos_id]*sign < 0) //not charge netural, need more 
	{	
		int i;
	//	cout << "mpos_id before: " << density[mpos_id];
		i = create_particle(mpos_id,p_data,density,sign,
				    mesh->materials[mpos_id]->electron_mass,
				    mesh); 
	//	cout << "mpos_id after: "<<density[mpos_id];
		current -= current_exit(p_data->pos,i);
	}
	//cout << "handle_region density:" << density[mpos_id]<< endl;
	//cout << "handle_region current:" << current << endl;
	return current;
}

double replenish_boundary(Particles *p_data,
			  int *nextDensity, Mesh *mesh)
{
	list<int>::iterator it = mesh->boundary.begin();
	list<int>::iterator end = mesh->boundary.end();
	double current = 0;
	int sign = 0;
	for(;it != end;++it)
	{
		int id = (*it);
/*		cout << id << " mpos: "<<mesh->mpos[id]<<" "
			<< mesh->mpos[id+1]<<endl;*/
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

extern "C" double replenish(Particles *p_data, int *nextDensity, Mesh *mesh)
{	
	double current;
	cout << "Beginning replenish" << endl;
	current = replenish_boundary(p_data,nextDensity,mesh);
	cout << "replenish complete" << endl;
	return current;
}
