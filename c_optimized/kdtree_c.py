import ctypes
#init
kdtree = ctypes.cdll.LoadLibrary("c_optimized/kdtree.so")
kdtree.call_this()

#types
vector2 = ctypes.c_double*2
vector2p = ctypes.c_double*3
kdtree_p = ctypes.pointer

#functions
kdtree.new_kdtree.restype = ctypes.POINTER(ctypes.c_int)
kdtree.find_point_r.restype = ctypes.POINTER(vector2)
kdtree.kdtree_find_point.restype = ctypes.POINTER(vector2)

def new_kdtree(points):
	for id in xrange(len(points)):
		points[id].append(id*1.)
	vec_points = vector2p*len(points)
	vecs = []
	for point in points:
		#print point
		vecs.append(vector2p(*point))
	vec_points = vec_points(*vecs)
	kd_p = kdtree.new_kdtree(vec_points,len(points),0)
	#print type(kd_p)
	for point in points:
		point.pop(2)
	return kd_p

def find_point(kd,point):
	best = vector2(100000,1000000)
	x = vector2(*point)
	bdist = ctypes.pointer(ctypes.c_double(1000000))
	
	#res =  kdtree.kdtree_find_point(kd,x)
	#res = kdtree.find_point_r(x,kd,best,bdist)
	res = kdtree.kdtree_find_point(kd,x)
	#print type(res)
	return res

def find_point_id(kd,point):
	best = vector2(100000,1000000)
	x = vector2(*point)
	bdist = ctypes.pointer(ctypes.c_double(1000000))
	
	#res =  kdtree.kdtree_find_point(kd,x)
	#res = kdtree.find_point_r(x,kd,best,bdist)
	id = kdtree.kdtree_find_point_id(kd,x)
	#print type(res)
	return id

def test():
	import time
	import random
	points = [[random.random()*10,random.random()*10] for x in range(100000)]
	start = time.time()
	kd = new_kdtree(points)
	print "construct:",time.time()-start
	start = time.time()
	for point in points:
		res =find_point(kd,point).contents
		if res[0] == point[0] and res[1] == point[1]:
			pass
			#print "success!"
		else:
			print "Unhappiness!"
			break
	print "exact:",time.time()-start
	points2 = []
	start = time.time()
	for point in points:
		point2 = [0,0]
		point2[0] = point[0]+.00000001*random.random()
		point2[1] = point[1]+.00000001*random.random()
		points2.append(point2)
	print "near_points construct:",time.time()-start
	start = time.time()
	for point2,point in zip(points2,points):
		res = find_point(kd,point2).contents	
		if res[0]==point[0] and res[1] == point[1]:
			pass
		else:
			print "unhappiness 2"
	print "near_points:",time.time()-start
	start = time.time()
	for index in xrange(len(points)):
		id = find_point_id(kd,points2[index])
		if id != index:
			print "Unhappiness 3"
	print "index:",time.time()-start
	
	print "all done"

if __name__ == "__main__":
	print "Testing"
	test()
