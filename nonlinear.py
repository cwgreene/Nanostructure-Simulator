from dolfin import *

class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return abs(x[0] - 1.0) < DOLFIN_EPS and on_boundary

mesh = UnitSquare(32,32)
V = FunctionSpace(mesh,"CG",1)

bc = DirichletBC(V,DirichletBoundary())

v = TestFunction(V)
u = TrialFunction(V)

a = dot(grad(v),U*U*grad(u))*dx
L = 
