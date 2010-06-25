#ifndef MOVE_PARTICLES_HPP
#define MOVE_PARTICLES_HPP
typedef list<double *> Polygon;

template<class KD>
void move_particles(Particles *pdata,
			double *efield,
			int *nextDensity,
			double dt,
			double length_scale,
			Mesh<KD> *mesh);

template<class KD>
double replenish(Particles *p_data, int *nextDensity, Mesh<KD> *mesh);

template<class KD>
double handle_region(int mpos_id,Mesh<KD> *mesh, Particles *p_data,
			int *density, int sign);

template<class KD>
double update_density(Particles *p_data,
			Mesh<kdtree> *mesh,
			int *nextDensity,
			Polygon *boundary,
			KD *kdt);
#endif
