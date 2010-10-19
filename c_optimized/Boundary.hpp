#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include <Eigen/Core>
#include <vector>

#include <iostream>

template<int dim>
struct Vector
{
	typedef Eigen::Matrix<double,dim,1> Type;
};

template<int dim>
class Line
{
public:
	typename Vector<dim>::Type start;
	typename Vector<dim>::Type end;

	double left;
	double right;
	double bottom;
	double top;
	Line(typename Vector<dim>::Type start,typename Vector<dim>::Type end);
	bool intersects(Line &line2);
	typename Vector<dim>::Type normal;
};

inline double min(double x,double y) { return x<y? x:y; }
inline double max(double x,double y) { return x>y? x:y; }

template<int dim>
Line<dim>::Line(typename Vector<dim>::Type start,typename Vector<dim>::Type end)
{
	this->start = start;
	this->end = end;

	//Thought these would be useful, but they aren't
	//Useful if we wanted to do a bounding rect thingy
	//But here we don't
	this->bottom = max(start[1],end[1]);
	this->top = min(start[1],end[1]);
	this->left = min(start[0],end[0]);
	this->right = min(start[0],end[0]);

}

template<int dim>
class Boundary
{
	std::vector<typename Vector<dim>::Type> *adjacent;//adjacent normals
	std::vector<Line<dim> > boundary_lines;//
public:
	Boundary(double *points,int n);

	bool reflect_trajectory(double *pos,typename Vector<dim>::Type old_pos);
	bool is_inward_facing(typename Vector<dim>::Type trajectory,int point);
	bool intersects_boundary(Line<dim> &line, int *id);

	void print_lines();
};

template<int dim>
void Boundary<dim>::print_lines()
{
	std::cout << "size is: "<<boundary_lines.size()<<std::endl;
	for(unsigned int i= 0; i < boundary_lines.size();i++)
	{
		for(int j = 0 ;j < dim;j++)
		{
			std::cout << boundary_lines[i].start[j] << " " ;
		}
		std::cout << "/";
		for(int j = 0 ; j < dim;j++)
		{
			std::cout << boundary_lines[i].end[j] << " " ;
		}
		std::cout << std::endl;
	}
}

template<int dim>
Boundary<dim>::Boundary(double *points,int n)
{
 	for(int i = 0 ;i < n;i++)
	{
		boundary_lines.push_back(Line<dim>
			  (typename Vector<dim>::Type(points+i*dim),
			   typename Vector<dim>::Type(points+((i*dim+dim)%(dim*n)))));
	}
}

/**/
/*I n t e r s e c t i o n*/
/*2 Dimensions*/
/*From Cormen et al*/
/**/

inline double direction(Eigen::Vector2d p1, Eigen::Vector2d p2, Eigen::Vector2d p3)
{
//	return ((p3-p1).cross(p2-p1))[2];
	Eigen::Vector2d a = p3-p1;
	Eigen::Vector2d b = p2-p1;

	return a[0]*b[1]-a[1]*b[0];
}

template<int dim>
bool Line<dim>::intersects(Line &line2)
{
	double d1 = direction(line2.end,line2.start,start);
	double d2 = direction(line2.end,line2.start,end);
	double d3 = direction(start,end,line2.start);
	double d4 = direction(start,end,line2.end);
	if (((d1>0 && d2<0) || (d1 < 0 && d2> 0)) &&
		(d3>0 && d4 < 0) || (d3 < 0 && d4>0))
		return true;
	//Going to ignore on segment concerns for now
	return false;
}
/*end I n t e r s e c t i o n*/

//is_inward_facing
//Takes a trajectory, and the nearest point,
//verifies that the two adjacent edges 
template<int dim>
bool Boundary<dim>::is_inward_facing(typename Vector<dim>::Type trajectory,int point)
{
	for(int i = 0 ; i < adjacent[point].size();i++)
	{
		double dot = 0;
		for(int d = 0; d < dim;d++)
			dot += adjacent[point][i][d]*trajectory[d]; //Dot with adjacent normals
		if(dot < 0)
			return false;
	}
	return true;
}

//Below 
template<int dim>
bool Boundary<dim>::intersects_boundary(Line<dim> &line, int *id)
{
	for(int i =0;i < boundary_lines.size();i++)
	{
		if(line.intersects(boundary_lines[i]))//2D only it seems
		{
			*id = i;
			return true; //Ignores possibility of two intersections
		}
	}
	return false;
}

//pos: direct reference to array location
template<int dim>
bool Boundary<dim>::reflect_trajectory(double *pos, typename Vector<dim>::Type oldpos)
{
	typename Vector<dim>::Type vpos(pos);
	Line<dim> trajectory(typename Vector<dim>::Type(pos), oldpos);
	int line_id;

	//Check for intersection
	if(!intersects_boundary(trajectory,&line_id))
		return false;

	//if intersects, get line	
	Line<dim> intersection_line = boundary_lines[line_id];

	//get intersection point
	typename Vector<dim>::Type ip = trajectory.intersection_point(intersection_line);
	typename Vector<dim>::Type remainder = vpos-ip;

	//Reflect "outside" line
	vpos = reflect_matrix(vpos);//Reflect remainder
	for(int i = 0; i < dim;i++)//Copy
		pos[i] = vpos[i];
	return true;
}
#endif
