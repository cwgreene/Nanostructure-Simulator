#include <stdlib.h>

typedef double vector2[2];

typedef struct
{
	vector2 pos;
	vector2 momentum;
	vector2 dx;

	double mass_particle;
	int id;
	int part_id;
	int dead;
	double lifetime;
	double charge;
}Particle;

typedef struct
{
}Mesh;

vector2 *vector2_add(vector2 *vec1,vector2 *vec2)
{
	int i;
	for(i=0;i<2;i++)
		(*vec1)[i] = (*vec1)[i]+(*vec2)[i];
	return vec1;
}

vector2 *scale(vector2 *vec1, double scale)
{
	int i;
	for(i=0;i<2;i++)
		(*vec1)[i] *= scale;
	return vec1;
}

void drift(Mesh *mesh, vector2 *func, Particle *p,vector2 *force)
{
	int i;
	for(i = 0;i<2;i++)
	{
		(*force)[i] = func[p->id][i]*p->charge;
	}
}

void randomElectronMovement(Particle *p,
                            vector2 *electric_field,
                            vector2 *density_func,
                            int dt,
                            int *reaper, int *reaper_length)
{
	vector2 force;
	int i;

	drift(mesh,electric_field,p,&force);
	for(i=0;i < 2;i++)
	{
		p->momentum[i] += temp[i]*dt;
		p->pos[i] += p->momentum[i]*dt;
		p->dx[i] += p->momentum[i]*dt;
	}
}

void move_particles(	Particle *p, int num_particle
			vector2 *electric_field,
			vector2 *density_func,
			Mesh *mesh,
			int *reaper,
			int *reaper_length)
{
	
}
