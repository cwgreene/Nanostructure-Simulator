#include "Face.hpp"

void Face::calculateNormalVector()
{
	for(int i = 0; i < num_points)
	{
		cross(
	}
}

//ReflectSurface
//The idea is to reduce the situation to
//   \ | y-c
//   |\|
//---|-*------ c
//   |/ \ | 
//   /   \|
//The above is done by taking the y distance and subtracting the
//y component of c. This is the same as calculating the projection
//onto the plane spanned by the face. 
//	We take that vector to the plane, add it once, getting us onto
//the plane, and then add it again, getting us to the reflected point.
//This function may need to be made part of the Material class. Because
//of this, we will not make it part of face.
void ReflectSurface(Face &face, Particles *p_data, int id)
{
	projectToPlane(face.plane,p_data->pos[id]);
}

Plane::Plane(double *spanvectors,int dim)
{
	span
}
