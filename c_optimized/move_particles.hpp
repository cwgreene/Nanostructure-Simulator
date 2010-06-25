#ifndef MOVE_PARTICLES_HPP
#define MOVE_PARTICLES_HPP
template<class KD>
void move_particles(Particles *pdata,
			double *efield,
			int *nextDensity,
			double dt,
			double length_scale,
			Mesh<KD> *mesh);

#endif
