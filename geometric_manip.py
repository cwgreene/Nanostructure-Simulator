import numpy as np

def cross(v1,v2):
	"""Computes the z magnitude of the cross product of two vectors"""
	return v1[0]*v2[1]-v1[1]*v2[0]

def normal(triangle):
	vec1 = triangle[1]-triangle[0]
	vec2 = triangle[2]-triangle[1]
	vec3 = triangle[0]-triangle[2]
	across = [0,0,0]
	across[0]= cross(vec1,vec2)
	across[1]= cross(vec2,vec3)
	across[2] = cross(vec3,vec1)
	return across

def extrude(triangle,height):
	vertexList = []
	faceList = []
	dir = np.array(normal(triangle))
	#create other triangle
	dest = np.vectorize(lambda x: dir+x)(triangle)
	#create tetrahedrogens between the two
	for i in xrange(dest.size):
		vertexlist.append(triangle[i])
		vertexlist.append(dest[i])

def intersect(line,plane):
	#plane is represented as three points, numpy arrays
	v1 = plane[1]-plane[0]
	v2 = plane[2]-plane[0]
	#line represnted by two points, numpy arrays
	v3 = line[1]-line[0]

	#print v1,v2,v3

	#form matrix
	mat = np.transpose(np.array([v1,v2,-v3]))
	#print mat
	
	#return interesection, which is [a,b,t]
	result = np.linalg.solve(mat,-plane[0])
	point = result[0]*v1+result[1]*v2+plane[0]
	return (point,result[2])
