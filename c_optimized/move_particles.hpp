#ifndef MOVE_PARTICLES_HPP
#define MOVE_PARTICLES_HPP
#include "Polytope.hpp"
typedef list<double *> Polygon;

template<class KD, int _dim>
void move_particles(Particles *pdata,
			double *efield,
			int *nextDensity,
			double dt,
			double length_scale,
			Mesh<KD, _dim> *mesh);

template<class KD,int _dim>
double replenish(Particles *p_data, int *nextDensity, Mesh<KD,_dim> *mesh);

template<class KD,int _dim>
double handle_region(int mpos_id,Mesh<KD,_dim> *mesh, Particles *p_data,
			int *density, int sign);

template<class KD, int _dim>
double update_density(Particles *p_data,
			Mesh<KD,_dim> *mesh,
			int *nextDensity,
			//Polytope<_dim> *boundary,
			KD *kdt);

//forward declarations
bool chance_recombinate(int other);
bool chance_scatter();
#endif
