#include "Polytope.hpp"

#include <iostream>
#include <eigen/Geometry>

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

	//3D
	double face_bottomp[]={0.,0.,0.,
                              1.,0.,0.,
                              1.,1.,0.,
                              0.,1.,0.};
	double face_topp[] = {0.,0.,1.,
                             1.,0.,1.,
                             1.,1.,1.,
                             0.,0.,1.};
	double face_frontp[] = {0.,0.,0.,
                              1.,0.,0.,
                              1.,0.,1.,
                              0.,1.,1.};
	double face_backp[]  = {0.,1.,0.,
                              1.,1.,0.,
                              1.,1.,1.,
                              0.,1.,1.};
	double face_leftp[]  = {0.,0.,0.,
                              0.,1.,0.,
                              0.,1.,1.,
                              0.,0.,1.};
	double face_rightp[] = {1.,0.,0.,
                              1.,1.,0.,
	                      1.,1.,1.,
                              1.,0.,1.};
	Face face_top = Face(face_topp,4,3);
	Face face_bottom = Face(face_bottomp,4,3);
	Face face_left = Face(face_leftp,4,3);
	Face face_right = Face(face_rightp,4,3);
	Face face_back = Face(face_backp,4,3);
	Face face_front = Face(face_frontp,4,3);
	vector<Face> faces;
	faces.push_back(face_top);
	faces.push_back(face_bottom);
	faces.push_back(face_front);
	faces.push_back(face_back);
	faces.push_back(face_left);
	faces.push_back(face_right);
	Polytope<3> poly3(faces);

	for(double y=1.5;y>-.5;y-=.1)
	{
		for(double x = -.5;x<1.5;x+=.1)
		{
			double point3[] = {x,y,.5};
			cout <<poly3.contains(point3);
		}
		cout <<endl;
	}
}
