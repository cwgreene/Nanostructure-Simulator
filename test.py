from dolfin import *

mesh = Mesh()
editor = MeshEditor()
topology = 2
dimension = 3
editor.open(mesh,"tetrahedron",3,3)
editor.initVertices(4)
editor.initCells(1)
editor.addVertex(0,0.,0.,0,)
editor.addVertex(1,1.,0.,0.)
editor.addVertex(2,0.,1.,0.)
editor.addVertex(3,0.,0.,1.)
editor.addCell(0,0,1,2,3)
editor.close()
plot(mesh)
