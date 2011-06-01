#ifndef UTILS_HPP
#define UTILS_HPP
#include <Eigen/Core>
#include <iostream>
#include <string>

template<int dim>
struct Vector
{
	typedef Eigen::Matrix<double,dim,1,false> Type;
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template<int dim>
std::string vec_str(const typename Vector<dim>::Type &vec)
{
	std::string result = "";
	for(int i = 0 ; i  < dim;i++)
	{
		std::string s;
		std::stringstream out;
		out << vec[i];
		s = out.str();
		result += s+ " ";
	}
	return result.substr(0,result.length()-1);
}
#endif
