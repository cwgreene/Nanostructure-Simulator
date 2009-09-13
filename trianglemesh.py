from dolfin import *
import numpy as np

def refinemesh(mesh,refinements):
	for i in xrange(refinements):
		mesh.refine()


def tmesh(refinements=0):
	mesh = Mesh()
	editor = MeshEditor()
	
	editor.open(mesh,"triangle",2,2)
	#begin
	editor.addCell(0,2,5,0)
	editor.initVertices(3)
	editor.initCells(1)
	editor.addVertex(0,0.0,0.0)
	editor.addVertex(1,1.0,0.0)
	editor.addVertex(2,.5,.866)
	editor.addCell(0,0,1,2)
	#end
	editor.close()

	refinemesh(mesh,refinements)
	return mesh

def innertriangle(refinements=0):
	mesh = Mesh()
	editor = MeshEditor()
	
	editor.open(mesh,"triangle",2,2)
	#begin
	editor.initVertices(6)
	editor.initCells(6)

	#outer
	vert1 = np.array([-.5,-.288675])
	vert2 = np.array([.5,-.288675])
	vert3 = np.array([0.,.57735])
	editor.addVertex(0,*vert1)
	editor.addVertex(1,*vert2)
	editor.addVertex(2,*vert3)

	#inner
	scale = .25
	editor.addVertex(3,*(vert1*scale))
	editor.addVertex(4,*(vert2*scale))
	editor.addVertex(5,*(vert3*scale))

	editor.addCell(0,0,3,1)
	editor.addCell(1,1,3,4)
	editor.addCell(2,1,4,5)
	editor.addCell(3,2,5,1)
	editor.addCell(4,2,5,0)
	editor.addCell(5,0,3,5)
	#end
	editor.close()

	refinemesh(mesh,refinements)
	return mesh
#

def innerouter(divisions):
	vert1 = np.array([-.5,-.288675])
	vert2 = np.array([.5,-.288675])
	vert3 = np.array([0.,.57735])

	outer = np.array([vert1,vert2,vert3])
	
	#generate_points
	for x in range(ydivisions):
		
