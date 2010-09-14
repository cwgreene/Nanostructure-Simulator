import hexagon
import numpy as np

def print_points(test):
	for y in np.arange(1,-1,-.1):
		for x in np.arange(1,-1,-.1):
			p = np.array([x,y])
			if test(p):
				print "x",
			else:
				print ".",
		print ""
	print ""

def in_poly(poly):
	return lambda p: hexagon.point_in_polygon(p,poly)

def test():
	square = np.array([[0,0],[0,1],[1,1],[1,0]])+np.array([-.5,-.5])
	outerhex = hexagon.polygon(5)
	innerhex = hexagon.polygon(5,.5)
	print "outer:",outerhex
	print "inner:",innerhex
	in_inner = in_poly(innerhex)
	in_outer = in_poly(outerhex)
	print_points(in_poly(square))
	print_points(in_poly(innerhex))
	print_points(in_poly(outerhex))
	print_points(lambda p: (not in_inner(p)) and (in_outer(p)))

if __name__=="__main__":
	test()
