"""This demo program solves Poisson's equation

    - div grad u(x, y) = f(x, y)

on the unit square with source f given by

    f(x, y) = 500*exp(-((x - 0.5)^2 + (y - 0.5)^2) / 0.02)

and boundary conditions given by

    u(x, y) = 0 for x = 0 or x = 1
"""

__author__ = "Anders Logg (logg@simula.no)/Chris Greene"
__date__ = "2007-08-16 -- 2008-12-13 / 2009-08-18"
__copyright__ = "Copyright (C) 2007-2008 Anders Logg"
__license__  = "GNU LGPL Version 2.1"

from dolfin import *
import montecarlo_mockup as mc
import numpy as np
import dolfin_util as du
import time
import mcoptions,sys
import trianglemesh as tm
options = mcoptions.get_options()
#dolfin_set("linear algebra backend","Epetra"

def custom_func(mesh,V,particles):
	f = Function(V)
#	for p in particles:
#		du.alter_cellid(mesh,f,p.id,p.charge)
	return f

# Create mesh and define function space
meshSizeX = 50
meshSizeY = 50
mesh = tm.innertriangle(5)#UnitSquare(meshSizeX,meshSizeY)
#plot(mesh)
V = FunctionSpace(mesh, "CG", 2)

# Define Dirichlet boundary (x = 0 or x = 1)
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        if on_boundary:
		if np.linalg.norm(x[0])+np.linalg.norm(x[1]) < .5:
			return True
	return False
# Define boundary condition
u0 = Constant(mesh, 0.0)
bc = DirichletBC(V, u0, DirichletBoundary())

# Define variational problem
v = TestFunction(V)
u = TrialFunction(V)

#function
particles = []

particles += mc.init_electrons(10,mesh.coordinates(),charge=-10,mesh=mesh)
#holes
particles += mc.init_electrons(10,mesh.coordinates(),charge=10,mesh=mesh)
f = custom_func(mesh,V,particles)
file = File("data/poisson_attract.pvd")
dfile = File("data/density_attract.pvd")
adfile = File("data/avg_density.pvd")
avfile = File("data/avg_voltage.pvd")
gradfile = File("data/grad_force.pvd")
avggradfile = File("data/avg_force.pvd")
avg_dens = mc.AverageFunc(f.vector().array())

def PoissonSolve(density):
	u = TrialFunction(V)
	v = TestFunction(V)

	a = dot(grad(v), grad(u))*dx
	L = v*density*dx
	# Compute solution
	problem = VariationalProblem(a, L, bc)
	sol = problem.solve()
	return sol

for x in range(options.num):
	sol = PoissonSolve(f)
	# Plot solution
	file << sol
	dfile << f
	print "Starting Step ",x
	start = time.time()
	electric_field = negGradient(mesh,potential_field)
	gradfile << electric_field
	mc.MonteCarlo(mesh,sol,electric_field,f,particles,avg_dens)
	print "Took: ",time.time()-start
#	plot(f)
file << sol
dfile << f

#dump average
f.vector().set(avg_dens.func)
adfile << f
avgE=mc.negGradient(mesh,PoissonSolve(f))
avggradfile << avgE
# Hold plot
plot(mc.negGradient(mesh,avgE))
interactive()
