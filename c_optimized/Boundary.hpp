#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include <Eigen/Core>
#include <vector>

#include <iostream>
#include <string>

template<int dim>
struct Vector
{
	typedef Eigen::Matrix<double,dim,1> Type;
};

template<int dim>
class Line
{
public:
	typedef typename Vector<dim>::Type VectorD;
	VectorD start;
	VectorD end;

	double left;
	double right;
	double bottom;
	double top;
	Line(VectorD start,VectorD end);

	bool intersects(Line &line2);
	VectorD intersection_point(Line<dim> &other_line);
	VectorD normal;

	std::string toString();
};

inline double min(double x,double y) { return x<y? x:y; }
inline double max(double x,double y) { return x>y? x:y; }

template<int dim>
Line<dim>::Line(VectorD start,VectorD end)
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
std::string Line<dim>::toString()
{
	std::string res = "";
	//start
	res += "(";
	for(int i= 0; i < dim;i++)
	{
		std::ostringstream num;
		num << start[i];
		res += num.str();
		res += ", ";
	}
	res = res.substr(0,res.length()-2);
	res += ") -> ";
	res += "(";
	for(int i= 0; i < dim;i++)
	{
		std::ostringstream num;
		num << end[i];
		res += num.str();
		res += ", ";
	}
	res = res.substr(0,res.length()-2);
	res += ")";
	return res;
}

//Intersection point
//Returns position of intersection
//Already have determined that the lines intersect
//Might be faster to just use this function and check t1 and t2 values
//
//Oh, right, the reason we probably DON'T do it that way is because
//we need start and end to be going in the right... no that's wrong
//No idea. Should probaby think about it more.
//
//Test status: Tested
template<>
Vector<2>::Type Line<2>::intersection_point(Line<2> &other)
{
	//Comments solve intersection equation
	//
	//To be honest, it'd probably just be faster to compute this
	//and then check that t2 and t1 are both inside [0,1], after checking that 
	//
	//
	//l1 = x1+(x2-x1)*t1
	//l2 = s1+(s2-s1)*t2
	//find t1,t2 such that l1=l2
	//return s1+(s2-s1)*t2
	//To solve for t1,t2
	//we note
	//l1-l2=0=x1+(x2-x1)*t1-s1-(s2-s1)*t2
	//resulting in
	//0=x1_x+(x2_x-x1_x)*t1-s1_x-(s2_x-s1_x)*t2
	//0=x1_y+(x2_y-x1_y)*t1-s1_y-(s2_t-s1)*t2
	//
	//  k1           a           b
	//s1_x-x1_x=(x2_x-x1_x)*t1+-(s2_x-s1_x)*t2
	//s1_y-x1_y=(x2_y-x1_y)*t1+-(s2_y-s1_y)*t2
	// k2            c           d
	//
	//k1 = a*t1+b*t2
	//k2 = c*t1+d*t2
	//
	//t1 = (k1-b*t2)/a
	//k2 = c*(k1-b*t2)/a+d*t2
	//k2 -c*k1/a=-c*b*t2/a+d*t2=(d-c*b/a)*t2
	//t2 = (k2-c*k1/a)/(d-c*b*/a)=(a*k2-c*k1)/(a*d-b*c) 
	//plug t2 back into t1
	//otherwise a=0 creates issues (which vanish when we plug back in, easiest way to think of this
	//is by thinking that t2(a) and we take the limit as a->0)
	//Cramer's proof will probably be better way to think so we need to look that back up.
	double a = end[0]-start[0];
	double b = -(other.end[0] - other.start[0]);
	double c = end[1]-start[1];
	double d = -(other.end[1] - other.start[1]);

	double k1 = other.start[0] - start[0];
	double k2 = other.start[1] - start[1];

	//This is where the magic happens!
	double t2 = (a*k2-c*k1)/(a*d-b*c);
	double t1 = (d*k1-b*k2)/(a*d-b*c);

#ifdef DEBUG
	std::cout << start[0] <<"/"<< start[0] << std::endl;
	std::cout << "times:"<<t1 << "/"<<t2<<std::endl;
#endif

	return (other.start+(other.end-other.start)*t2);
}

template<int dim>
class Boundary
{
public:
	std::vector<typename Vector<dim>::Type> *adjacent;//adjacent normals
	std::vector<Line<dim> > boundary_lines;//

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
	for(unsigned int i =0;i < boundary_lines.size();i++)
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
