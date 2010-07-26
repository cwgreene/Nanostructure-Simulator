#include "Polytope.hpp"

#include <iostream>

using namespace std;

void testPoints(Polytope<2> poly)
{
	for(unsigned int i = 0; i < poly.points.size();i++)
		cout << "<"<<poly.points[i][0]<<","<<poly.points[i][1] << "> ";

	cout << endl;
	for(double y = -.2; y <= 1.2; y += .1)
	{
		for(double x = -.2;  x <=1.2; x+=.1)
		{
			double test[] = {x,y};
			cout <<poly.contains(test);
		}
		cout << endl;
	}
}

int main()
{
	//Square
	double points[] = {0.,0.,
			 0.,1.,
			 1.,1.,
			 1.,0.};
	Polytope<2> poly1(points,4);
	testPoints(poly1);

	double points2[] = {0.0,0.0,.5,.86,1.0,.0};
	Polytope<2> poly2(points2,3);
	testPoints(poly2);
}
