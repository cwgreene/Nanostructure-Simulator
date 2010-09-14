#include <stdlib.h>
#include <iostream>
#include "materials.hpp"


using namespace std;
Material::Material(double _electron_mass, double *_random_momentum,
		int _max_n,int _dim)
{
	random_momentum = new double[_dim*_max_n];
	dim = _dim;
	electron_mass = _electron_mass; 
	max_n = _max_n; 
	dim = _dim;
	for(int i = 0; i < dim*max_n; i++)
	{
		random_momentum[i] = _random_momentum[i];
	}
}

void material_random_momentum(	Material *material,double *momentum,
				double weight)
{
	int dim = material->dim;
	int max_n = material->max_n;
	int i = rand()%max_n;
	for(int c = 0; c< material->dim;c++)
	{
		momentum[c] = material->random_momentum[dim*i+c]*weight;
	}
}

extern "C" Material *new_material(double electron_mass,
				  double *random_momentum,
				  int max_n, int dim)
{
	return new Material(electron_mass,
				  random_momentum,max_n,dim);
}

extern "C" Material **material_pair(Material *p_type, Material *n_type)
{
	Material **pair = new Material*[2];
	pair[0] = p_type;
	pair[1] = n_type;
	return pair;
}
