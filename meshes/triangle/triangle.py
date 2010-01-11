from dolfin import *
import meshtest
import trianglemesh as tm
import numpy as np
import montecarlo_mockup as mc

class TriangleMesh(ParticleMesh):
	def __init__(self):
		mesh = self
		#thetriangle = np.array([[-.5,-.288675],[.5,-.288675],[0.,.57735]])
		thetriangle = np.array([[0.,0.],[1.,0.],[.5,.8660254]])

		mesh = meshtest.TestMesh()
		for x in range(options.size):
			mesh.refine()
		#mesh = mc.ParticleMesh(mesh,options.scale)
		mc.ParticleMesh.__init__(
		#mesh = mc.ParticleMesh(tm.innertriangle(5,.2,thetriangle))

		inner = triangle.scale_triangle(thetriangle,.52)
		self.populate_regions(
				lambda x: triangle.point_in_triangle(x,inner), 
				0,0)

		boundarymesh = BoundaryMesh(mesh)
		mesh.OuterBoundary = OuterTriangle
		mesh.InnerBoundary = InnerTriangle
		print len(boundarymesh.coordinates())
		return mesh

class InnerTriangle(SubDomain):
    def inside(self, x, on_boundary):
        if on_boundary:
		#if (x[0],x[1]) in repeated:
			#print "repeat:",(x[0],x[1])
		#repeated[(x[0],x[1])] = 0
		if mesh.in_p_region(x):
			#print "inner:",x[0],x[1]
			return True
	return False

class OuterTriangle(SubDomain):
    def inside(self, x, on_boundary):
        if on_boundary:
		if not mesh.in_p_region(x):
			#print "outer",x[0],x[1]
			return True
	return False

