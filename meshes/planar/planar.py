import dolfin
import numpy as np
from dolfin import *
import montecarlo_mockup as mc

global mesh

#Innerboundary Low voltage, ptype
#OuterBoundary High Voltage,ntpe
#Left: ndoped ->Inner
#right : pdoped -> Outer
class PlanarMesh(mc.ParticleMesh):
	def __init__(self,options,n_material,p_material):
		global mesh
		thetriangle = np.array([[0.,0.],[1.,0.],[.5,.8660254]])

		mesh = dolfin.UnitSquare(50,50)
		mc.ParticleMesh.__init__(self,mesh,
			options.scale,options.length,options.dt,options.gen_num)
		self.populate_regions(lambda x: x[0] < .5,
			0,0,
			n_material,
			p_material) 
		self.OuterBoundary = LeftRegion()
		self.InnerBoundary = RightRegion()
		mesh = self

class RightRegion(SubDomain):
    def inside(self, x, on_boundary):
	if on_boundary:
		if mesh.in_p_region(x):
			return True
	return False
class LeftRegion(SubDomain):
    def inside(self, x, on_boundary):
	if on_boundary:
		if not mesh.in_p_region(x):
			return True
	return False

