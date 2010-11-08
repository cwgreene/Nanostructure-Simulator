import os
import itertools as it
import numpy as np
import inspect

threshold = 1.*10**-6
#Qhull communicator
def qhull_faces(vertices):
	#open connection
	child_in, child_out = os.popen2("qhull i")
	#print dimensions and number of points
	line = str(len(vertices[0]))+" "+str(len(vertices))
	child_in.write(line+"\n")
	
	#print vertices
	for vertex in vertices:
		line = ""
		for coordinate in vertex:
			line += str(coordinate)+" "
		line = line[:-1]#strip trailing space
		child_in.write(line+"\n")
	child_in.close() #all points are ready
	
	#grab output
	output = child_out.readlines()[1:]#first is count
	child_out.close()
	
	#get indices from output
	indices = map(lambda x: map(int,x.strip().split(" ")),output)
	print indices
	#construct faces from indices
	faces = []
	for face in indices:
		faces.append(map(lambda index: vertices[index],
				 face))
	#done
	return faces

#given a list of pairs of points that form 
#the line segments of a convex hull, return them as a set
#of ordered points.
def polygon(faces2d):
	x = []
	for head,tail in faces2d:
		if tuple(head) not in x:
			x.append(tuple(head))
		if tuple(tail) not in x:
			x.append(tuple(tail))
	return order_points(x)

#Given a set of points
#calculate the center of mass
#return as vector
def center_mass(points):
	x,y = transpose(points)
	x_m = sum(x)*1.0/len(x)
	y_m = sum(y)*1.0/len(y)
	return np.array((x_m,y_m))

#Given a set of points on a convex hull
#calculate the angle they are at from the
#center of mass
#does the order depend on the point?
#hmm... consider one of the points lies on the thingy. So yes.
def angles_center(points):
	v_m = center_mass(points)
	angles = []
	#get points to be array of vectors
	points = map(np.array, points)
	vectors = map(lambda x:(x-v_m)/np.linalg.norm(x-v_m),points)
	angles=map(lambda x:np.angle(x[0]+x[1]*1j)*360./(2*np.pi),vectors)
	return angles

#given a set of points on a convex hull
#return a list of those points in order, counter clockwise	
def order_points(points):
	def second(x):
		return x[1]
	angles = angles_center(points)
	ordered = map(second,sorted(zip(angles,points)))
	return ordered

#file reader
def qhull_file(filename):
	afile = open(filename) 
	afile.readline() #get rid of header
	values = [map(float,line.strip().split()) for line in afile]
	return values

#utils
def flatten(alist):
	res = []
	for obj in list(alist):
		if hasattr(obj,"__len__"):
			res += flatten(obj)
		else:
			res.append(obj)
	return res

#rotate
#takes list and returns a list of same length
#just counting from an offset
#if offset is greater than length of list it's
#equivalent to passing in (offset % len(alist))
def rotate(alist,offset):
	offset = offset % len(alist)
	return alist[offset:]+alist[:offset]

#inner_outer
#takes in list of points on a complex boundary with
#convex inner and outer boudaries
#finds the inner boundary and the outer boundary
#by taking the convex hull, which is the outer boundary
#then finding the convex hull of the remaining points
def inner_outer(points):
	outer = qhull_faces(points)
	print outer
	if len(points[0]) == 2:
		outer = polygon(outer)
	inner = []
	for point in points:
		if not point_on_curve(point,outer): #if on the outer curve, toss em
			inner.append(list(point))
	inner = order_points(inner)
	inner,outer = map(minimal_bound,(inner,outer))
	print "inner:",inner,"\n","outer:",outer
	return inner,outer

def test(func,val,expected,success="Success",failure="Failure:"):
	if inspect.argvalues(func)[0] == 1:
		result = func(val)
	else:
		result = func(*val)
	print "Testing:",func.__name__
	if result== expected:
		print success
		return False
	else:
		print failure,val
		print "Expected:",expected
		print "Result:",result
		return False

#minimal_bound
#takes a curve, returns the minimal number of points
#that specify that curve
#that is, removes unnecssary collinear points
#Assumes that curve contains at least three points.
#Should work on any polygon, but points must be
#ordered by placemnt on curv?
def test_minimal_bound():
	#simple triangle, should be unchanged
	triangle = [[0,0],[0,1],[1,1]]
	test(minimal_bound,triangle,triangle)
	triangle2 = [[0,0],[0,.5],[0,1],[1,1]]
	test(minimal_bound,triangle2,triangle)
	
def minimal_bound(curve):
	ring3 = zip(rotate(curve,-1),rotate(curve,0),rotate(curve,1))
	result_curve = []
	#compare vector of current with vector of next
	#if collinear, then we don't "need it
	for cpoint1, cpoint2,cpoint3 in ring3:
		vec1 = np.array(cpoint2)-np.array(cpoint1)
		vec2 = np.array(cpoint3)-np.array(cpoint2)
		if (np.abs(np.cross(vec1,vec2)) >= threshold).any(): #make sure we are not zero_vector
			result_curve.append(cpoint2)
	return result_curve
			

#geometric helper functions

#we claim that a point on a curve if is on the line
#between two points on that line
def test_between():
	point = [0,0]
	contains =  [[-1,0],[1,0]]
	does_not_contain = [[1,1],[1,2]]
	test(between,[point]+contains,True)
	test(between,[point]+does_not_contain,False)
	
def between(point,point1,point2):
	threshold = 1.*10**-6
	vec1 = np.array(point)-np.array(point1) #vec with point on curve as tail
	vec2 = np.array(point2)-np.array(point1) #vector passing through this and next point
	if (abs(np.cross(vec1,vec2)) <= threshold).all():
		return True
	return False

def point_on_curve(point,curve):
	ring = zip(curve,curve[1:]+[curve[0]])
	for cpoint1,cpoint2 in ring:
		if between(point,cpoint1,cpoint2):
			return True
	return False

def transpose(alist):
	a = []
	b = []
	for x,y in alist:
		a.append(x)
		b.append(y)
	return a,b
