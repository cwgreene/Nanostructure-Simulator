#ifndef MATERIALS_H
#define MATERIALS_H
class Material
{
public:
	double electron_mass;

	//Bandstructure information.
	int max_n;
	double *random_momentum_x;
	double *random_momentum_y;

	Material(double _electron_mass, int _max_n,
		 double *_random_momentum_x,double *_random_momentum_y);
};
void material_random_momentum(Material *material,double *momentum);
#endif
