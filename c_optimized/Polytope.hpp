#ifndef POLYTOPE_H
#define POLYTOPE_H

#include <math.h>
#include <vector>
#include <stdlib.h>

#include "Face.hpp"

template<int dim> class Polytope{};

template<>
class Polytope<3>
{
public:
	vector<Face> faces;
	bool contains(double point[3]);
	Face *nearest_face(Vector3f &point,double *dist);
	Face *nearest_face(Vector3f &point);
	Polytope(vector<Face> faces);
};

template<>
class Polytope<2>//Polygon
{
public:
	bool contains(double *point);
	vector<double *> points;
	Polytope(double * points);
};

/*To determine if you're inside a polytope, travel in one direction and count intersections.*/
bool Polytope<3>::contains(double p[3])
{
	int dim =3;
	Vector3f ray_vector;
	Vector3f intersect;	
	Vector3f point = Vector3f(p[0],p[1],p[2]);
	
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

/*bool point_in_polygon(Polygon *boundary, vector2 *pos)
{
}*/

double cross_segment_point(double seg1_x,double seg1_y,
			   double seg2_x,double seg2_y,
			   double point_x,double point_y)
{
	double segment_x = seg1_x-seg2_x;
	double segment_y = seg1_y-seg2_y;
	double cur_vec_x = seg1_x-point_x;
	double cur_vec_y = seg1_y-point_y;
	double cross_z = (segment_x)*(cur_vec_y)-
			 (cur_vec_x)*(segment_y);
	return cross_z;
}
bool Polytope<2>::contains(double *point)
{
	vector<double *>::iterator it = points.begin();
	vector<double *>::iterator end = points.end();

	double cur_x,cur_y;
	double prev_x = (*it)[0];
	double prev_y = (*it)[1];
	double first_x = prev_x;
	double first_y = prev_y;
	double p_x = point[0];
	double p_y = point[1];
	++it;
	double sign = 0;
	for(;it != end;++it)
	{
		point = *it;
		cur_x = point[0];
		cur_y = point[1];
		double cross_z = cross_segment_point(cur_x,cur_y,
						     prev_x,prev_y,
						     p_x,p_y);
		if(sign == 0)
		{
			sign = cross_z;
		}else if( cross_z*sign < 0)
		{
			return false;
		}
		prev_x = cur_x;
		prev_y = cur_y;
	}
	//Last with first
	if(cross_segment_point(first_x,first_y,
				prev_x,prev_y,
				p_x,p_y)*sign < 0)
		return false;
	
	return true;

}

Face *Polytope<3>::nearest_face(Vector3f &point,double *dist)
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

Face *Polytope<3>::nearest_face(Vector3f &point)
{
	double dist = 0;
	return nearest_face(point, &dist);
}
#endif
