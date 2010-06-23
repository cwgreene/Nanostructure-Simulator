#ifndef MATERIALS_H
#define MATERIALS_H
class Material
{
public:
	double electron_mass;

	//Bandstructure information.
	int max_n;
	int dim;
	double *random_momentum;

	Material(double _electron_mass, double *_random_momentum,int max_n,int dim);
};
void material_random_momentum(Material *material,double *momentum);
#endif
