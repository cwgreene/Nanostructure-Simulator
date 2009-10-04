import numpy as np
from numpy.linalg import norm

class Trapezoid():
	def __init__(self,origin,base,height1,height2):
		self.origin = np.array(origin)
		self.x=np.array(base)
		self.y=np.array(height1)
		self.y2=np.array(height2)

class BasicTrap(Trapezoid):
	def __init__(self):
		Trapezoid.__init__(self,[0.,0.],[1.,0.],[.2,.5],[-.2,.5])

class TriTrap(Trapezoid):
	def __init__(self):
		Trapezoid.__init__(self,[0.,0.],[1.,0.],
					[.5,.8660254],[-.5,.8660254])

def cull_triangles(triangles,cull):
	reduced = []
	for triangle in triangles:
		keep = True
		for point in triangle:
			if cull(point):
				keep = False
		if keep:
			reduced.append(triangle)
	return reduced

def triangulate(trap,nx,ny):
	unitx = trap.x/(1.*nx)
	unity = trap.y/(1.*ny)
	unity2 = trap.y2/(1.*ny)
	triangles = []
	for y in xrange(ny):
		#for each y we have the up triangles,
		#and the down triangles
		triangles += add_up_triangles(unitx,unity,unity2,y,trap)
		triangles += add_down_triangles(unitx,unity,unity2,y,trap)
	return triangles

def flat(list,max_depth=10,depth=0):
	"""word on the street says this method,
	using recursion, will start coughing
	at around 10 levels deep. Use iter methods for deeper
	or just a while loop."""
	flattened = []
	for obj in list:
		objtype = type(obj)
		if (objtype==type(list) or objtype==type(tuple))and depth < max_depth:
			flattened += flat(obj,depth+1,max_depth)
		else:
			flattened.append(obj)
	return flattened

def get_row(base,unitx,start,end):
	pos_x = start
	i = 1
	cur_row = []
	while np.dot((end+base)-pos_x,base)>=0:
		#print "dist to end",np.dot((end+base)-pos_x,base)/norm(unitx)
		cur_row.append(pos_x)
		pos_x = unitx*i+start
		i+=1
#	if norm((unity2+base)-pos_x) < norm(unitx)/100:
#		cur_row.append(pos_x)
	return cur_row

def add_up_triangles(unitx,unity,unity2,y,trap):
	#to add the up triangles,
	#we travel paralell to the base (trap.x)
	#take two of points at a time
	#and then take one point from the next row up
	#which will be this row 
	triangles = []
	cur_row = get_row(trap.x,unitx,unity*y,unity2*y)+trap.origin
	next_row = get_row(trap.x,unitx,unity*(y+1),unity2*(y+1))+trap.origin
	bottom_pairs = zip(cur_row[:-1],cur_row[1:])
	triangles = map(lambda x:flat(x,1),zip(bottom_pairs,next_row))
	return triangles


def add_down_triangles(unitx,unity,unity2,y,trap):
	#to add the up triangles,
	#we travel paralell to the base (trap.x)
	#take two of points at a time
	#and then take one point from the next row up
	#which will be this row 
	triangles = []
	top_row = get_row(trap.x,unitx,unity*(y+1),unity2*(y+1))+trap.origin	
	bottom_row = get_row(trap.x,unitx,unity*y,unity2*y)+trap.origin
	bottom_row = bottom_row[1:-1]
	top_pairs = zip(top_row[:-1],top_row[1:])
	triangles = map(flat,zip(top_pairs,bottom_row))
	return triangles

