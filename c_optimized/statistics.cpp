#include <iostream>
#include <math.h>

/*#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>*/

#include "particles.hpp"
#include "mesh.hpp"

using namespace std;

inline double square(double x)
{
	return x*x;
}

void avg_momentum_grid(Particles *p, Mesh *mesh)
{
	double *particles = p->pos;
	for(int i = 0; i < mesh->npoints;i++)
	{
		double avg = 0;
//		cout << "npoints: "<< mesh->npoints;
		for(list<int>::iterator it = mesh->electrons_pos[i].begin();it != mesh->electrons_pos[i].end();++it)
		{
			int id = *it;
			avg += sqrt(square(pkx(id))+square(pky(id)));
		}
//		cout << "eMesh:" << i << " "<< mesh->mpos[2*i] <<","<< mesh->mpos[2*i+1] << " p: "<<avg << " n: " <<mesh->holes_pos[i].size()<<endl;
		double hole_avg  =0;
		for(list<int>::iterator it = mesh->holes_pos[i].begin();it != mesh->holes_pos[i].end();++it)
		{
			int id = *it;
			avg += sqrt(square(pkx(id))+square(pky(id)));
		}
//		cout << "hMesh:" << i << " "<< mesh->mpos[2*i] <<","<< mesh->mpos[2*i+1] << " p: " <<hole_avg << " p: "<<mesh->holes_pos[i].size() <<endl;
	
	}
}
