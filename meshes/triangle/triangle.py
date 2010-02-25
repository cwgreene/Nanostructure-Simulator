from dolfin import *
import meshtest
import trianglemesh as tm
import numpy as np
import montecarlo_mockup as mc
import triangle_util

global mesh

class TriangleMesh(mc.ParticleMesh):
	def __init__(self,options,n_material,p_material):
		global mesh #very bad, but innertriangle needs it.
		thetriangle = np.array([[0.,0.],[1.,0.],[.5,.8660254]])

		mesh = meshtest.TestMesh()
		for x in range(options.size):
			mesh.refine()
		mc.ParticleMesh.__init__(self,mesh,
			options.scale,options.length,options.dt,options.gen_num)

		inner = triangle_util.scale_triangle(thetriangle,.52)
		self.populate_regions(lambda x: 
			triangle_util.point_in_triangle(x,inner),0,0,
			n_material,
			p_material) 
		boundarymesh = BoundaryMesh(mesh)
		self.OuterBoundary = OuterTriangle()
		self.InnerBoundary = InnerTriangle()
		print len(boundarymesh.coordinates())
		mesh = self

class InnerTriangle(SubDomain):
	def in_p_region(x):
		pass
    	def inside(self, x, on_boundary):
		if on_boundary:
			if mesh.in_p_region(x):
				return True
		return False

class OuterTriangle(SubDomain):
	def inside(self, x, on_boundary):
		if on_boundary:
			if not mesh.in_p_region(x):
				return True
		return False
