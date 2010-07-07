import kdtree_c as kd
import kdtree3_c as kd3
import numpy as np
import numpy.random as random

def dist(point1,point2):
	vec1,vec2=np.array(point1),np.array(point2)
	return np.linalg.norm(vec1-vec2)

n=100000
no_fail = True
points = map(list,zip(random.rand(n),random.rand(n)))
tree = kd.new_kdtree(points)

for id,point in zip(range(len(points)),points):
	offset=list((point[0]+.01,point[1]+.01))
	res= kd.find_point_id(tree,offset)
	if id != res:
		better=dist(offset,points[res])
		expected = dist(offset,point)
		if better > expected:
			no_fail=False
			print "point:",point,expected
			print "offset:",offset
			print "believed:",points[res],better
			break
if no_fail==True:
	print "Success for 2d!",id
else:
	print "Failure for 2d",id

print "Longest chain 2D",kd.longest_chain(tree,0,0)

points = map(list,zip(random.rand(n),
			random.rand(n),
			random.rand(n)))
tree = kd3.new_kdtree3(points)

for id,point in zip(range(len(points)),points):
	offset=list((point[0]+.01,point[1]+.01,point[2]+.01))
	res= kd3.find_point3_id(tree,offset)
	#print res,id,point,points[res]
	if id != res:
		better=dist(offset,points[res])
		expected = dist(offset,point)
		if better > expected:
			no_fail=False
			print "point:",point,expected
			print "offset:",offset
			print "believed:",points[res],better
			break
if no_fail:
	print "Success for 3d!",id
else:
	print "Failure in 3d",id
