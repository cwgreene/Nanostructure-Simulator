#include "Boundary.hpp"
#include <Eigen/Core>


bool test_Line()
{
	using namespace Eigen;
	Line<2> line0(Vector2d(0.,0.),Vector2d(1.,0.));
	Line<2> line1(Vector2d(.5,.8),Vector2d(.5,-.5));
	Line<2> line2(Vector2d(.25,.1),Vector2d(.6,1.0));
	
	std::cout<<"Intersects 0,1 (1):"<<line0.intersects(line1)<< std::endl;//True
	std::cout<<"Intersects 0,2 (0):"<<line0.intersects(line2)<< std::endl;//False
	std::cout<<"Intersects 1,2 (1):"<<line1.intersects(line2)<< std::endl;//False
	return true;
}

bool test_Boundary()
{
	double bounds[]={0.0,0.0, //>v
		         0.0,1.0, //^<
		         1.0,1.0,
		         1.0,0.0};
	Boundary<2> boundary(bounds,4);

	boundary.print_lines();

	double cur_pos[]={0.75,-0.25};
	Vector<2>::Type old_pos(.25,.25);
	//boundary.reflect_trajectory(cur_pos,old_pos);
	std::cout << "cur_pos[0]: "<<cur_pos[0]<< "cur_pos[1]: "<<cur_pos[1]<<std::endl;

	Line<2> line1(Vector<2>::Type(.5,.5),Vector<2>::Type(1.5,.5));

	int id;
	if(boundary.intersects_boundary(line1,&id))
	{
		std::cout << "Successfully detected intersection ("<<id<<")"<<":" << " ";
		std::cout << boundary.boundary_lines[id].toString()<<std::endl;
		Vector<2>::Type point = boundary.boundary_lines[id].intersection_point(line1);
		std::cout << point[0] << "," << point[1] << std::endl;
	}

	return true;
}

int main()
{
	test_Line();//Test line constructors and intersections
	test_Boundary();
}
