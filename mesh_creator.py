from dolfin import *
from triangulate import flat
from itertools import *

def mesh_creator(triangles):
	#flatten triangle list
	points = flat(triangles)

	#create mesh and editor	
	mesh = Mesh()
	editor = MeshEditor()
	editor.open(mesh,"triangle",2,2)

	editor.initCells(len(triangles))
	editor.initVertices(len(points))

	point_ids = {}

	#put points into hashtable,add them as vertices
	for point,id in izip(points,count()):
		point_ids[tuple(point)] = id
		editor.addVertex(id,*point)
	#now add cells
	for tri,id in izip(triangles,count()):
		tri_id = map(lambda p: point_ids[tuple(p)],tri)
		editor.addCell(id,*tri_id)
	editor.close()
	return mesh
