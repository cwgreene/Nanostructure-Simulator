import dolfin
import numpy as np
from dolfin import *
import montecarlo_mockup as mc

global mesh

def rectangular_mesh(x,y):
        mesh = Mesh()
        ed = MeshEditor()
        ed.open(mesh,"triangle",2,2)
        cells = 2*x*y
        verts = (x*y+1)**2
        ed.initVertices(verts)
        ed.initCells(cells)
        step = 1./(max(x,y))
        #create vertices
        for i in range(x+1):
                for j in range(y+1):
                        ed.addVertex(i*(y+1)+j,i*step,j*step)
                        print (i*(y+1)+j),i*step,j*step
        #create cells
        for i in range(x):
                for j in range(y):
                        cellid=i*y+j
                        vertid = i*(y+1)+j
                        ed.addCell(2*cellid,vertid,vertid+1,vertid+(y+1)+1)
                        ed.addCell(2*cellid+1,vertid,vertid+(y+1),vertid+(y+1)+1)
        ed.close()
        return mesh

#Innerboundary Low voltage, ptype
#OuterBoundary High Voltage,ntpe
#Left: ndoped ->Inner
#right : pdoped -> Outer
class PlanarMesh(mc.ParticleMesh):
	def __init__(self,options,n_material,p_material):
		global mesh
		thetriangle = np.array([[0.,0.],[1.,0.],[.5,.8660254]])

		mesh = rectangular_mesh(50,100)
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

