#ifndef FACE_HPP
#define FACE_HPP
/*Face is polygon in 3d*/
class Face
{
	void calculateNormalVector();
	
public:
	int num_points;
	double *points;
	bool cotainsPoint(double *point);
	Face(double points,int num_points,int dim);
	double *normalVector;
	Plane plane; //Array which points to two vectors;
};

class Plane
{
public:	
	int dim;
	double *spanVectors;
	Plane(double *spanvectors,int dim);
};
#endif
