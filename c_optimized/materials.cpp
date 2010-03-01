class Material
{
public:
	double electron_mass;

	//Bandstructure information.
	int max_n;
	double *random_momentum_x;
	double *random_momentum_y;

	Material(double _electron_mass, int _max_n,
		 double *_random_momentum_x,double *_random_momementum_y)
	{
		electron_mass = _electron_mass; 
		max_n = _max_n; 
		random_momentum_x = _random_momentum_x; 
		random_momentum_y = _random_momentum_y;
	}
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
