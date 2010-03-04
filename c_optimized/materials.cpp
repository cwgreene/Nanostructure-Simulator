#include <stdlib.h>
#include <iostream>
#include "materials.hpp"

using namespace std;
Material::Material(double _electron_mass, double *_random_momentum,int _max_n)
{
	electron_mass = _electron_mass; 
	max_n = _max_n; 
	random_momentum = _random_momentum; 
}

void material_random_momentum(Material *material,double *momentum)
{
	int max_n = material->max_n;
	int i = rand()%max_n;
	momentum[0] = material->random_momentum[2*i];
	momentum[1] = material->random_momentum[2*i+1];
}

extern "C" Material *new_material(double electron_mass,
				  double *random_momentum,
				  int max_n)
{
	return new Material(electron_mass,
				  random_momentum,max_n);
}

extern "C" Material **material_pair(Material *p_type, Material *n_type)
{
	Material **pair = new Material*[2];
	pair[0] = p_type;
	pair[1] = n_type;
	return pair;
}
