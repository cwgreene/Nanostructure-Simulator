#include <stdio.h>

#define POSITIONX(i)  4*i
#define POSITIONY(i)  4*i+1
#define MOMENTUMX(i)  4*i+2
#define MOMENTUMY(i)  4*i+3
#define MASS(i)  4*i+4

#define px(i)  particle[POSITIONX(i)]
#define py(i)  particle[POSITIONY(i)]
#define pkx(i)  particle[MOMENTUMX(i)]
#define pky(i)  particle[MOMENTUMY(i)]
//#define pmass(i)  particle[MASS(i)]

void randomElectronMovement(double *particle,
				double *p_mass,
				int *p_id,
				int *p_charge,
				int i,
				double *efield,
				double dt,
				double length_scale)
{
	//Move
	double dx = pkx(i)*dt/length_scale/p_mass[i];
	double dy = pky(i)*dt/length_scale;
	px(i) += dx;
	py(i) += dy;
	//Drift
	pkx(i) += efield[2*p_id[i]]*p_charge[i]*dt/length_scale;
	pky(i) += efield[2*p_id[i]+1]*p_charge[i]*dt/length_scale;
	//Scatter
}

void move_particles(double *particles,
			double *p_mass,
			int *p_charge,
			int *p_id,
			int nparticles,
			double *efield,
			double *nextDensity,
			double dt,
			double length_scale)
{
	int i;
	printf("moving particles now\n");
	printf("particles: %d\n",nparticles);
	for(i = 0; i < nparticles;i++)
	{
	//	printf("started %d\n",i);
		int id = p_id[i];
		int charge = p_charge[i];
	//	printf("id:%d charge:%d\n",p_id[i],charge);
		nextDensity[id] -= charge;
	//	printf("next_density updated: %d\n",i);
		randomElectronMovement(particles,p_mass,p_charge,p_id,i,
					efield,dt,length_scale);
	//	printf("completed %d\n",i);
	}
}

void call_me()
{
	printf("I got called!\n");
}
