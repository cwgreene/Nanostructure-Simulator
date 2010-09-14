from dolfin import *
import numpy as np
import montecarlo_mockup as mc
import qhull

global mesh

#2D
#point: np array
#polygon: np array of 2 vecs
def point_in_polygon(point,polygon):
	sign = 0
	origin = sum(polygon)/len(polygon)
	point = point-origin
	polygon = polygon-origin
	polyvecs = list(polygon-point)
	polypair = zip(polyvecs,polyvecs[1:]+polyvecs[:1])
	for v1,v2 in polypair:
		nsign = np.cross(v1,v2)
		if sign == 0:
			sign = nsign
		elif sign*nsign < 0: #sign switch
			return False
	return True

def polygon(n,r=1,center=np.array([0,0])):
	degree = 2*np.pi/n
	points = []
	for i in range(n):
		points.append(np.array([r*np.cos(degree*i+np.pi/2),
		                       r*np.sin(degree*i+np.pi/2)])
		              +center)
	x,y = qhull.transpose(points)
	return points

def mesh_median(mesh):
	points = mesh.coordinates()
	median = sum(points)/len(points)
	return median

class HexagonMesh(mc.ParticleMesh):
	def __init__(self,options,n_material,p_material):
		global mesh
		mesh = Mesh("meshes/hexagon/Hexagon.xml")
		mesh.order() #why?

		#create geometric shadows
		median = mesh_median(mesh)
		hexagon = polygon(5,1,median) #shadows mesh
		innerhex = polygon(5,.75,median)

		#add points if necessary
		for x in range(options.size):
			mesh.refine()
		#create ParticleMesh
		mc.ParticleMesh.__init__(self,mesh,options.scale,options.length,
			options.dt,options.gen_num)
		
		#mark points as ntype or ptype
		#or reflecting
		self.populate_regions(lambda x: point_in_polygon(x,innerhex),
			0,0,
			n_material,
			p_material) 
		boundarymesh = BoundaryMesh(mesh)
		self.OuterBoundary = OuterHexagon()
		self.InnerBoundary = InnerHexagon()
		#why is this here, why don't I make it
		#associated with this object?
		mesh = self 
	def __del__(self):
		print "Hexagon Mesh Destroyed"
			

class InnerHexagon(SubDomain):
    def inside(self, x, on_boundary):
	if on_boundary:
		#if (x[0],x[1]) in repeated:
			#print "repeat:",(x[0],x[1])
		#repeated[(x[0],x[1])] = 0
		if mesh.in_p_region(x):
			print "inner:",x[0],x[1]
			return True
	return False
class OuterHexagon(SubDomain):
    def inside(self, x, on_boundary):
	if on_boundary:
		if not mesh.in_p_region(x):
			print "outer",x[0],x[1]
			return True
	return False

