import numpy as np

class Node:pass
class kdTree:
	def __init__(self,pointList):
		self.root = self.kdtree(pointList)

	def find_point(self,point):
		return self.find_point_r(point, self.root,
				np.array([10000,10000]))

	def find_point_r(self,point,node,best):
		bdiff = point-best
		bdist = np.dot(bdiff,bdiff)
		diff = point-node.location
		dist = np.dot(diff,diff)

		#if our dist is smaller than the current best
		#we're the best
		if bdist > dist:
			best = node.location #we're the best!
		#are we to the right or left of this
		#points splitting plane?
		if point[node.split] < node.location[node.split]:
			if node.leftChild != None:
				best = self.find_point_r(point,
							 node.leftChild,best)
			other_branch = node.rightChild
		
		if point[node.split] >= node.location[node.split]:
			if node.rightChild != None:
				best = self.find_point_r(point,
							 node.rightChild,best)
			other_branch = node.leftChild

		#how close are we to the other side?
		bdiff = point-best
		bdist = np.dot(bdiff,bdiff)

		#if we are too close to split point
		#check other side
		dist = point[node.split]-node.location[node.split]
		dist = dist*dist
		if dist<bdist:
			#check other tree:
			if other_branch != None:
				best = self.find_point_r(point,other_branch,
							 best)
		return best

	#from wikipedia 
	def kdtree(self,pointList, depth=0):
		if len(pointList) == 0:
			return
	 
		# Select axis based on depth so that axis 
		#cycles through all valid values

		#assume all points have the same dimension
		k = len(pointList[0]) 
		axis = depth % k
	 
		# Sort point list and choose median as pivot element
		pointList.sort(key=lambda point: point[axis])
		median = len(pointList)/2 # choose median
	 
		# Create node and construct subtrees
		node = Node()
		node.location = np.array(pointList[median])
		node.split = axis
		node.leftChild = self.kdtree(pointList[0:median], depth+1)
		node.rightChild = self.kdtree(pointList[median+1:], depth+1)
		return node

#own function
