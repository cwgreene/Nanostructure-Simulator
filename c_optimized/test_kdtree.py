import kdtree_c as kd
import numpy as np
import numpy.random as random

def dist(point1,point2):
	vec1,vec2=np.array(point1),np.array(point2)
	return np.linalg.norm(vec1-vec2)

points = map(list,zip(random.rand(100000),random.rand(100000)))
tree = kd.new_kdtree(points)

for id,point in zip(range(len(points)),points):
	offset=list((point[0]+.01,point[1]+.01))
	res= kd.find_point_id(tree,offset)
	print res,id
	if id != res:
		better=dist(offset,points[res])
		expected = dist(offset,point)
		if better > expected:
			print "point:",point,expected
			print "offset:",offset
			print "believed:",points[res],better
			break
