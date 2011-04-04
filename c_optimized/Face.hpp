#ifndef FACE_HPP
#define FACE_HPP
#include <vector>
#include <list>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Geometry>

USING_PART_OF_NAMESPACE_EIGEN

using namespace std;

class Plane
{
public:	
	Vector3d v1;
	Vector3d v2;
	Vector3d point;
	Matrix3d projection;
	
	Vector3d project(const Vector3d &vec);
	Plane(Vector3d _v1,Vector3d _v2,Vector3d point);
};

/*Face is polygon in 3d*/
class Face
{
	void calculateNormalVector();
	
public:	
	//Constructor
	Face(double *_points,int num_points,int dim);

	//Fields, need to be initialized
	vector<Vector3d> points; //needs to be inited
	Vector3d normal; //Can be copied into
	Plane *plane; //Array which points to two vectors;

	//Methods	
	bool contains(const Vector3d &point);
	bool line_intersects(	const Vector3d &point,
				const Vector3d &vec,
				Vector3d &intersect);
	bool line_intersects_plane(const Vector3d &point,
				   const Vector3d &vec,
				   Vector3d &intersect);
//	bool line_intersects(const Vector3d &point,
//				const Vector3d &vec);//detection only
	Face *nearest_face(Vector3d point);
	Face *nearest_face(Vector3d point,double *dist);
	double distance(const Vector3d &point);
};

Plane::Plane(Vector3d _v1, Vector3d _v2,Vector3d _p)
{
	v1 = _v1/_v1.norm();
	v2 = _v2/_v2.norm();
	point = _p;
	projection << 	v1[0],v2[0],0,
			v1[1],v2[1],0,
			v1[2],v2[2],0;
}

void check_failure(bool test, string msg)
{
	if(test)
	{
		cout 	<< msg << endl;
		exit(-1);
	}
}

Face::Face(double *_points, int num_points, int dim)
{
	check_failure(dim!=3,"Serious error. Faces are three dimensional");
	check_failure(num_points <= 2, "Insufficient points in Face construction");

	for(int i = 0; i < num_points;i++)
	{
			Vector3d p(_points[i*dim],
				_points[i*dim+1],
				_points[i*dim+2]);
			points.push_back(p);
	}

	//Get Normal
	Vector3d v1((points)[1]-(points)[0]); //First vec
	int i =2;
	while(((points)[i]-(points)[i-1]).dot(v1) != 0 && i < num_points)
	{
		cout << "Point: "<< i << endl;
		cout << points[i] << endl << points[i-1] << endl;
		i++;
	}
	check_failure(i == num_points, "All points of face are collinear. Exiting");
	Vector3d v2((points)[i]-(points)[i-1]);
	normal = v1.cross(v2);
	plane = new Plane(v1,v2,points[0]);
}

bool Face::contains(const Vector3d &point)
{
	vector<Vector3d>::iterator bound= points.begin();
	vector<Vector3d>::iterator last = --(points.end());

	Vector3d vec1 = *bound-point;
	Vector3d vec2 = *(bound+1)-point;
	Vector3d prev = vec1.cross(vec2);
	++bound;
	
	for(;bound != last; ++bound)
	{
		vec1 = *bound-point;
		vec2 = *(bound+1)-point;
		Vector3d ncross = vec1.cross(vec2);
		if(ncross.dot(prev) < 0)
			return false;
		prev = ncross;
	}
	return true;
}

bool Face::line_intersects(const Vector3d &point,
			   const Vector3d &line,
			   Vector3d &intersect)
{
	if(line_intersects_plane(point,line,intersect))
		if(contains(intersect))
			return true;
	return false;
}


bool Face::line_intersects_plane(const Vector3d &point,
				 const Vector3d &line,
				 Vector3d &intersect)
{
	Matrix3d left;
	left <<	plane->v1[0],plane->v2[0],-line[0],
		plane->v1[1],plane->v2[1],-line[1],
		plane->v1[2],plane->v2[2],-line[2];
	if(left.determinant() == 0)
	{
		return false;
	}
	intersect = (left.inverse()*point)[2]*line+(point+plane->point);
	return true;
}

Vector3d Plane::project(const Vector3d &vec)
{
	return projection*vec;
}

double Face::distance(const Vector3d &point)
{
	Vector3d p_vec= plane->project(point);
	//We should figure out how far the projected point is if it's outside the hull,
	//But in the convex situation, we're guaranteed that the projection is _inside_ the face.
	//I think.
	return p_vec.norm();
}

extern "C" Face* create_face(double *x,int num_points, int dim)
{
	return new Face(x,num_points,dim);
}

#endif
