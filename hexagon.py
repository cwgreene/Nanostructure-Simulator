import numpy as np
import triangulate
from dolfin import *
import mesh_creator as mcreator

def deg(degrees):
	return (degrees/360.)*2*np.pi

upright = np.array([np.cos(deg(60)),np.sin(deg(60))])
upleft = np.array([np.cos(deg(120)),np.sin(deg(120))])
right = np.array([1.,0.])

point1 = np.array([0.,0.])
point2 = point1+upleft
point3 = point2+upright
point4 = point3+right
point5 = point4-upleft
point6 = point5-upright
point7 = point6+upleft

triangle1=[point1,point6,point7]
triangle2=[point2,point1,point7]
triangle3=[point3,point2,point7]
triangle4=[point4,point3,point7]
triangle5=[point5,point4,point7]
triangle6=[point6,point5,point7]

#x = triangulate.Triangle(triangle1)
#triangulated = triangulate.triangulate(x,10,11)
triangles = map(triangulate.Triangle,[triangle1,triangle2,triangle3,
				triangle4,triangle5,triangle6])
triangulated = map(lambda tri:triangulate.triangulate(tri,8,8),triangles)
triangulated = triangulate.flat(triangulated)

mesh = mcreator.mesh_creator(triangulated)
plot(mesh)
interactive()
