import ctypes
import time #test only
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
	points_copy = []
	for point in points:
		points_copy.append(list(point))
	for id in xrange(len(points_copy)):
		points_copy[id].append(id*1.)
	vec_points = vector2p*len(points_copy)
	vecs = []
	for point in points_copy:
		#print point
		vecs.append(vector2p(*point))
	vec_points = vec_points(*vecs)
	kd_p = kdtree.new_kdtree(vec_points,len(points),0)
	#print type(kd_p)
	return kd_p
acc = 0.
def find_point(kd,point):
	global acc
	best = vector2(100000,1000000)
	x = vector2(*point)
	bdist = ctypes.pointer(ctypes.c_double(1000000))
	
	#res =  kdtree.kdtree_find_point(kd,x)
	#res = kdtree.find_point_r(x,kd,best,bdist)
#	start = time.time()
	res = kdtree.kdtree_find_point(kd,x)
#	acc += time.time()-start
	#print type(res)
	return res

def find_point_id(kd,point):
	global acc
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
	global acc
	points = [[random.random()*10,random.random()*10] for x in range(10000)]
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
	z=[(random.random(),random.random()) for x in range(200000)]
	start = time.time()
	acc = 0
	for point in z:
		id = find_point_id(kd,point)
	print "random_lookup:",time.time()-start,acc
	start = time.time()
	acc = 0
	for point in z:
		kdtree.do_nothing(3)
		kdtree.do_nothing(3)
	print "do_nothing:",time.time()-start,acc

	print "all done"

if __name__ == "__main__":
	print "Testing"
	test()
