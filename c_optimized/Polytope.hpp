#ifndef POLYTOPE_H
#define POLYTOPE_H

#include <math.h>
#include <vector>
#include <stdlib.h>

#include <Eigen/Geometry>

#include "Face.hpp"

//General Polytope Declaration.
template<int dim> class Polytope{};

//3d Polytope Class
template<>
class Polytope<3>
{
public:
	vector<Face> faces;
	bool contains(double point[3]);
	Face *nearest_face(Vector3d &point,double *dist);
	Face *nearest_face(Vector3d &point);
	Polytope(vector<Face> faces);
};

//2D Polytope Class
template<>
class Polytope<2>//Polygon
{
public:
	bool contains(double *point);
	vector<Vector2d> points;
	Polytope(double * _points,int num);
};

//Polytope<2> Constructor
Polytope<2>::Polytope(double *_points,int num)
{
	for(int i= 0; i < num;i++)
	{
		points.push_back(Vector2d(_points[2*i],_points[2*i+1]));
	}
}


//Polytope<3> Constructor
Polytope<3>::Polytope(vector<Face> _faces)
{
	faces = _faces;
}


//bool Polytope<3>::contains(double p[3])
//Checks to see if polytope contains the given point
//Presumably this is called by some function that has a view
//into an array that contains all points. This is probably
//why we pass in a double * instead of a Vector3d.
bool Polytope<3>::contains(double p[3])
{
	int dim =3;
	Vector3d ray_vector;
	Vector3d intersect;	
	Vector3d point = Vector3d(p[0],p[1],p[2]);
	
	for(int i = 0; i < dim; i++)
		ray_vector[i] = (1.*rand())/RAND_MAX;

	int intersection_count = 0;

	vector<Face >::iterator face = faces.begin();
	vector<Face >::iterator end = faces.end();
	for(;face != end;++face)
	{
			if(face->line_intersects(point,ray_vector,intersect))
				intersection_count++;
	}
	if(intersection_count >0 && intersection_count % 2 ==0) //enter in, enter out.
		return true;
	return false; 
}

//double cross_segment_point
//Helper function, given a point and a line, essentially
//determines the angle subtended by the point and the segment
//sin of the angle is calculated, since that's the cross product
//and we only care about the sign. (if it changes, we're outside)
double cross_segment_point(double seg1_x,double seg1_y,
			   double seg2_x,double seg2_y,
			   double point_x,double point_y)
{
/*	double segment_x = seg1_x-seg2_x;
	double segment_y = seg1_y-seg2_y;
	double cur_vec_x = seg1_x-point_x;
	double cur_vec_y = seg1_y-point_y;
	double cross_z = (segment_x)*(cur_vec_y)-
			 (cur_vec_x)*(segment_y);*/
	Vector3d v1(seg1_x-point_x,seg1_y-point_y,0);
	Vector3d v2(seg2_x-point_x,seg2_y-point_y,0);
	return v1.cross(v2)[2];
}

//bool Polytope<2>::contains(double *point)
//Determines if the polytope contains the point
//For two dimesions, we use the winding method.
//This was originally a Polygon function.
bool Polytope<2>::contains(double *point)
{
	double cur_x,cur_y;
	double prev_x = points[0][0];
	double prev_y = points[0][1];
	double p_x = point[0];
	double p_y = point[1];
	double sign = 0;
	
	//Check to see if the sign changes as we go around
	//the polygon and compute the cross product
	for(unsigned int i = 1;i < points.size()+1;i++)
	{
		Vector2d edge_point = points[i%points.size()];
		cur_x = edge_point[0];
		cur_y = edge_point[1];
		double cross_z = cross_segment_point(cur_x,cur_y,
						     prev_x,prev_y,
						     p_x,p_y);
		if(sign == 0)
		{
			sign = cross_z;
		}else if( cross_z*sign < -10e-10)
		{
			return false;
		}
		prev_x = cur_x;
		prev_y = cur_y;
	}
	return true;

}

//Face *Polytope<3>::nearest_face(Vector3d &point,double *dist)
//Finds the face of the polyope closest to a given point. 'dist'
//holds the distance.
Face *Polytope<3>::nearest_face(Vector3d &point,double *dist)
{
	vector<Face>::iterator face = faces.begin();
	vector<Face>::iterator end = faces.end();

	Face *closest = &(*face);
	double nearest_dist = face->distance(point);
	++face;
	for(;face!=end;++face)
	{
		double dist = face->distance(point);
		if(nearest_dist > dist)
		{
			closest = &(*face);
		}
	}
	*dist = nearest_dist;
	return closest;
}

//Face *Polytope<3>::nearest_face(Vector3d &point)
//Currying function to: nearest_face(point, &dist).
Face *Polytope<3>::nearest_face(Vector3d &point)
{
	double dist = 0;
	return nearest_face(point, &dist);
}

//C Interfaces
//Create 2 dimensional polytope, aka polygon
extern "C" Polytope<2> *create_polytope2(double *points, int num_points)
{
	Polytope<2> *res = new Polytope<2>(points,num_points);
	return res;
}

extern "C" Polytope<3> *create_polytope3(Face *faces, int num_faces)
{
	vector<Face> vfaces(faces,faces+num_faces);
	return new Polytope<3>(vfaces);
}

#endif
