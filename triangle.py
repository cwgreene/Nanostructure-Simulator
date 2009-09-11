import itertools as it
import grid
import fdm

def vec(p1,p2):
	"""takes in two indexable objects with at least two dimensions and constructs a two dimensional vector"""
	return (p2[0]-p1[0],p2[1]-p1[1])

def dot(v1,v2):
	"""computes the dot product of two two-dimensional vectors. Does not handle higher dimensions"""
	return v1[0]*v2[0]+v1[1]*v2[1]

def cross(v1,v2):
	"""Computes the z magnitude of the cross product of two vectors"""
	return v1[0]*v2[1]-v1[1]*v2[0]

def point_in_triangle(point,triangle):
	"""Assumes that the triangle's vertices are ordered counterclockwise by angle,
	(probably should add a decorator to force this). Checks the z magnitudes of the cross products
	of the vectors formed by the point to the vertices. This should be > 0 for all (< 180 degs)"""

	"""Construct the vectors"""
	vec1 = vec(point,triangle[0])
	vec2 = vec(point,triangle[1])
	vec3 = vec(point,triangle[2])

	"""calculate their dot products. Simpler to be explicit, since the order is important
	we don't want to rely on it.combinations, because it would compute vec1,vec3, instead
	of vec3,vec1"""
	across = [0,0,0]
	across[0]= cross(vec1,vec2)
	across[1]= cross(vec2,vec3)
	across[2] = cross(vec3,vec1)

	for cp in across:
		if cp < 0:
			return False
	
	return True

def centroid(triangle):
	"""Calculates the centroid of the given triangle"""

	"""Starts with centroid being at 0,0"""
	cen = [0,0]
	for vertex in triangle:
		cen[0] += vertex[0]
		cen[1] += vertex[1]
	cen[0] /= 3.
	cen[1] /= 3.
	return cen

def str_grid(test,xmax,ymax):
	astr =""
	for y in range(ymax+1):
		for x in range(xmax+1):
			if test(x,(ymax-y)):
				astr+="x"
			else:
				astr+=" "
		astr +="\n"
	return astr

def points_in_triangle(triangle,increment,func=None):
	xvals = map(lambda ver: ver[0],triangle)
	yvals = map(lambda ver: ver[1],triangle)
	minx,miny,maxx,maxy = (min(xvals),min(yvals),max(xvals),max(yvals))

	agrid = grid.construct_grid(minx,
			    miny,
			    maxx,
			    maxy,
			    inc = increment,
			    test = lambda x,y: point_in_triangle([x,y],triangle),
			    value= lambda x,y: fdm.NodeData(0,1))
	#for node in agrid:
	#	agrid.label[node].node_id = agrid.nodes().index(node)
#	grid.number_grid(agrid)
	return agrid
	

def scale_triangle(triangle,scale):
	"""Scales given triangle around it's centroid"""

	"""Set up the new triangle (ntriangle) and the array for the
	distance vectors to the edges"""
	ntriangle = []
	vecs = []
	"""Calculate the centroid"""
	cen = centroid(triangle)
	"""Create the distance vectors from centroid to edge"""
	for edge in triangle:
		vecs.append(vec(cen,edge))
#	vecs.append(vec(cen,triangle[0]))
#	vecs.append(vec(cen,triangle[1]))
#	vecs.append(vec(cen,triangle[2]))

	"""For each of the distance vectors scale them. Note, they
	are at the origin, so we need to reoffset them once we're done
	back to the centroid"""
	for vector in vecs:
		ntriangle.append([vector[0]*scale,vector[1]*scale])
	for vector in ntriangle:
		vector[0] += cen[0]
		vector[1] += cen[1]
	return ntriangle

