import numpy as np
import dolfin
import geometric_manip as geo
import itertools as it

class PreMesh:
	def __init__(self):
		self.faces = []
		self.vertices = []
	def add_vertex(self,vertex):
		self.vertices.append(vertex)
	def add_face(self,face):
		self.faces.append(face)

	def compile(self):
		mesh = dolfin.Mesh()
		editor = dolfin.MeshEditor(mesh,"tetrahedron",3,3)
		editor.open(mesh,"tetrahedron",3,3)
		editor.initVertices(len(self.vertices))
		editor.initCells(len(self.faces))
		for index,pos in it.izip(it.count(),self.vertices):
			editor.addVertex(index,pos)
		for index,face in it.izip(it.count(),self.faces):
			editor.addCell(index,face)
		editor.close()
