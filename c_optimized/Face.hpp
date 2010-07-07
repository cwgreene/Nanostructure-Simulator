#ifndef FACE_HPP
#define FACE_HPP
#include <vector>
#include <Eigen/Core>
#include <Eigen/LU>

USING_PART_OF_NAMESPACE_EIGEN

using namespace std;

class Plane
{
public:	
	Vector3f v1;
	Vector3f v2;
	Matrix3f projection;
	
	Vector3f Plane::project(const Vector3f &vec);
	Plane(double *spanvectors);
};

/*Face is polygon in 3d*/
class Face
{
	void calculateNormalVector();
	
public:	
	vector<Vector3f> points;
	bool contains(const Vector3f &point);
	Face(double *points,int num_points,int dim);
	Vector3f normal;
	Plane plane; //Array which points to two vectors;
	bool line_intersects(	const Vector3f &point,
				const Vector3f &vec,
				Vector3f &intersect);
	bool line_intersects_plane(const Vector3f &point,
				   const Vector3f &vec,
				   Vector3f &intersect);
	bool line_intersects(const Vector3f &point,
				const Vector3f &vec);//detection only
	Face *nearest_face(Vector3f point);
	Face *nearest_face(Vector3f point,double *dist);
	double Face::distance(const Vector3f &point);
};



bool Face::contains(const Vector3f &point)
{
	vector<Vector3f>::iterator bound= points.begin();
	vector<Vector3f>::iterator last = --(points.end());

	Vector3f vec1 = *bound-point;
	Vector3f vec2 = *(bound+1)-point;
	Vector3f prev = vec1.cross(vec2);
	++bound;
	
	for(;bound != last; ++bound)
	{
		vec1 = *bound-point;
		vec2 = *(bound+1)-point;
		Vector3f ncross = vec1.cross(vec2);
		if(ncross.dot(prev) < 0)
			return false;
		prev = ncross;
	}
	return true;
}

bool Face::line_intersects(const Vector3f &point,
			   const Vector3f &line,
			   Vector3f &intersect)
{
	if(line_intersects_plane(point,line,intersect))
		if(contains(point))
			return true;
	return false;
}


bool Face::line_intersects_plane(const Vector3f &point,
				 const Vector3f &line,
				 Vector3f &intersect)
{
	Matrix3f left;
	left <<	plane.v1[0],plane.v2[0],line[0],
		plane.v1[1],plane.v2[1],line[1],
		plane.v1[2],plane.v2[2],line[2];
	if(left.determinant() == 0)
		return false;
	intersect = left.inverse()*point;
	return true;
}

Plane::Plane(double *spanvectors)
{
	
}

Vector3f Plane::project(const Vector3f &vec)
{
	return projection*vec;
}

double Face::distance(const Vector3f &point)
{
	Vector3f p_vec= plane.project(point);
	//We should figure out how far the projected point is if it's outside the hull,
	//But in the convex situation, we're guaranteed that the projection is _inside_ the face.
	//I think.
	return p_vec.norm();
}

#endif
