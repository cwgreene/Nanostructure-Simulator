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

bool test_inward_facing(Eigen::Vector2d x, Eigen::Vector2d y, Boundary<2> &b)
{
	std::cout << "is_inward_facing: ["<<vec_str<2>(x)<<"] ["<<vec_str<2>(y)<<"]\n";
	std::cout<< b.is_inward_facing(x,y)<<"\n";
	return true;
}

bool test_Boundary()
{
	using namespace Eigen;
	double bounds[]={0.0,0.0, //>v
		         0.0,1.0, //^<
		         1.0,1.0,
		         1.0,0.0};

	//Test Construction
	std::cout<< "\nTesting Construction of Boundary\n";
	Boundary<2> boundary(bounds,4,Vector2d(.1,.7));
	std::cout<< "\n Lines:\n";
	boundary.print_lines();
	std::cout<< "\n Line Normals:\n";
	boundary.print_normals();
	std::cout<< "\n points:\n";
	boundary.print_points();
	std::cout<< "\n adjacents:\n";
	boundary.print_adjacent();

	//Test Reflection
	std::cout<< "\nTesting Reflection\n";
	double cur_pos[]={.58,1.28};
	Vector<2>::Type old_pos(.45,.83);
	std::cout << " cur_pos[0]: " <<cur_pos[0]<< " cur_pos[1]: "<<cur_pos[1]<<std::endl;
	std::cout << " old_pos[0]: " <<old_pos[0]<< " old_pos[1]: "<<old_pos[1]<<std::endl;
	if(boundary.reflect_trajectory(cur_pos,old_pos))
	{
		std::cout<<" After Reflection:\n";
		std::cout << " cur_pos[0]: " <<cur_pos[0]<< " cur_pos[1]: "<<cur_pos[1]<<std::endl;
	}else
		std::cout<< " Error! Not Reflected!"<<std::cout;

	//Test Intersection	
	std::cout<< "\nTesting Intersection\n";
	Line<2> line1(Vector<2>::Type(.5,.5),Vector<2>::Type(1.5,.5));
	int id;
	if(boundary.intersects_boundary(line1,&id))
	{
		std::cout << "Successfully detected intersection ("<<id<<")"<<":" << " ";
		std::cout << boundary.boundary_lines[id].toString()<<std::endl;
		Vector<2>::Type point = boundary.boundary_lines[id].intersection_point(line1);
		std::cout << point[0] << "," << point[1] << std::endl;
	}

	//Test Line Normals
	std::cout<< "\nTesting Normals\n";
	for(unsigned int i = 0 ; i < boundary.normals.size();i++)
	{
		std::cout<< "Normal "<<i<<": "<< boundary.normals[i][0] << " " <<
			 boundary.normals[i][1] << std::endl;
	}

	//Test is-inward-facing
	std::cout<< "\nTesting Inward Facing\n";
	test_inward_facing(Vector2d(1,1),Vector2d(0,0),boundary);
	test_inward_facing(Vector2d(1,1),Vector2d(1,1),boundary);
	test_inward_facing(Vector2d(-1,1),Vector2d(1,1),boundary);
	test_inward_facing(Vector2d(-1,1),Vector2d(1,1),boundary);
	test_inward_facing(Vector2d(1,-1),Vector2d(0,0),boundary);

	return true;
}

int main()
{
	test_Line();//Test line constructors and intersections
	test_Boundary();
}
