import dolfin_util as du
from dolfin import *
import numpy as np

class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return ((x[0] < DOLFIN_EPS ) or (x[0] > 1.0 - DOLFIN_EPS) or
		(x[1] < DOLFIN_EPS) or (x[1] > 1.0 -DOLFIN_EPS))

mesh = UnitSquare(50,50)
V = FunctionSpace(mesh,"CG",1)
Vgrad = VectorFunctionSpace(mesh,"CR",1,2)
# Define boundary condition
u0 = Constant(mesh, 0.0)
bc = DirichletBC(V, u0, DirichletBoundary())

#fucntions
v = TestFunction(V)
u = TrialFunction(V)
f = Function(V)

du.set_cell(mesh,f,[.5,.5],10)
du.set_cell(mesh,f,[.58,.58],10)
#du.set_cell(mesh,f,[.49,.49],10000)
#du.set_cell(mesh,f,[.5,.51],10000)
#du.set_cell(mesh,f,[.51,.50],10000)
#du.set_cell(mesh,f,[.49,.5],10000)
#du.set_cell(mesh,f,[.49,.51],10000)
#du.set_cell(mesh,f,[.51,.49],10000)


#poisson variational problem
a = dot(grad(v), grad(u))*dx
L = v*f*dx

problem = VariationalProblem(a, L, bc)
res = problem.solve()
gradu = project(grad(-res),Vgrad)
print du.get_vec(mesh,gradu,np.array([.5,.5]))
print du.get_vec(mesh,gradu,np.array([.58,.58]))
plot(gradu)
interactive()
