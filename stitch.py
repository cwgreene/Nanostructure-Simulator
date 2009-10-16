import numpy as np

def flat(alist):
	flat_list = []
	if type(alist) != type([]):
		return [alist]
	for p in alist:
		flattened = flat(p)
		for x in flattend:
		
	return x

def stitch(triangles1,triangles2,small):
	"""Returns the collection of triangles
	formed by triangles 1 and triangles 2. Triangles
	that have nearby points get merged"""

	#get points from triangles
	points1,points2 = map(flat,[triangles1,triangles2])
	points_result = {}

	#merge points
	for x in points1+points2:
		points_result[tuple(x)] = x

	#Find the points in points1 that are (almost) in points2
	#since our new triangles need to use the same points
	#we keep track of 
	for x in points1:
		y = tree.find_point(x)
		if np.norm(x-y) < small:
			points_result[tuple(y)] = x
		else:
			points_result[tuple(y)] = y

	#create new triangles from old ones
	new_triangles = []
	for triangle in triangles1+triangles2:
		new_triangle
		for point in triangle:
			new_triangle.append(points_result[tuple(point)])
		new_triangles.append(new_triangle)
	return new_triangles
