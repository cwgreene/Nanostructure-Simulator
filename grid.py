import numpy as np

def extremum(points,comp):
	ext_x = 0.
	ext_y = 0.
	bleft = array([ext_x,ext_y])
	for point in points:
		if comp(point[0],ext_x) < 0:
			ext_x = point[0]
		if comp(point[1], ext_y) < 0:
			ext_y = point[1]
	return bleft


def upperight(points)
	return extremum(points,lambda x:-cmp(x))

def bottomleft(points):
	return emtremum(points,cmp)
	
class Grid:
	def __init__(self,mesh,resolution):
		self.grid = {}
		self.resolution = resolution

		coord = mesh.coordinates()
		#box coordinates
		bottomleft(coord)
		upperright(coord)

		grid = self.grid
		for x in arange(bottomleft[0],upperright[0],resolution):
			for y in arange(bottomleft[1],upperright[1],resolution):
				grid[(x,y)] = []
		#now that grid is constructed, start dumping points in
		for point in coord:
			grid[tuple(point/resolution)].append(point)
		#now, at each grid point, dump in the nearest adjacent set of
		#points. If the list is still empty, extend radius
		for point in grid:
			offsetx = range(-1,1)
