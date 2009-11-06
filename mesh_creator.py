from dolfin import *
import triangulate as trig
from itertools import *

def mesh_creator(triangles):
	#flatten triangle list
	print "triangles",triangles[0:1]
	points = trig.flat(triangles)
	print "points",points[0:3]

	#create mesh and editor	
	mesh = Mesh()
	editor = MeshEditor()
	editor.open(mesh,"triangle",2,2)


	point_ids = {}

	#put points into hashtable,add them as vertices
	for point in points:
		try:
			point_ids[tuple(point)] = 0
		except:
			print point,type(point)
			raise
	print len(points),len(point_ids)

	#Init Points, now that we know how many
	editor.initCells(len(triangles))
	editor.initVertices(len(point_ids))
	for point,id in izip(point_ids,count()):
		point_ids[point] = id
		editor.addVertex(id,*point)
	#now add cells
	for tri,id in izip(triangles,count()):
		tri_id = map(lambda p: point_ids[tuple(p)],tri)
		editor.addCell(id,*tri_id)
	editor.close()
	print "Mesh triangles:points",len(triangles),":",len(mesh.coordinates())
	return mesh
