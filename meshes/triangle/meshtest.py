from dolfin import *
import triangle_util, triangulate as trig, mesh_creator as mc

def point_in_triangle_checker(tri):
	tri = [[0.,0.],[1.,0.],[.5,.8660254]]
	tri = triangle_util.scale_triangle(tri,.12)
	def check(p):
		return triangle_util.point_in_triangle(p,tri)
	return check

def TestMesh():
	tri = trig.TriTrap()
	tris = trig.triangulate(tri,20,21)
	tris = trig.cull_triangles(tris,point_in_triangle_checker(tris))

	mesh = mc.mesh_creator(tris)
	#plot(mesh)
	return mesh
