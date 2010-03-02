#include <stdlib.h>
#include "materials.hpp"

Material::Material(double _electron_mass, int _max_n,
		 double *_random_momentum_x,double *_random_momentum_y)
{
	electron_mass = _electron_mass; 
	max_n = _max_n; 
	random_momentum_x = _random_momentum_x; 
	random_momentum_y = _random_momentum_y;
}

void material_random_momentum(Material *material,double *momentum)
{
	int max_n = material->max_n;
	momentum[0] = material->random_momentum_x[rand()%max_n];
	momentum[1] = material->random_momentum_y[rand()%max_n];
}

extern "C" Material *new_material(double electron_mass,
				  int max_n,
				  double *random_momentum_x,
				  double *random_momentum_y)
{
	return new Material(electron_mass,
				  max_n,
				  random_momentum_x,
				  random_momentum_y);
}

extern "C" Material **material_pair(Material *p_type, Material *n_type)
{
	Material **pair = new Material*[2];
	pair[0] = p_type;
	pair[1] = n_type;
	return pair;
}
