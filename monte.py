"""This demo program solves Poisson's equation

    - div grad u(x, y) = f(x, y)

on the unit square with source f given by

    f(x, y) = 500*exp(-((x - 0.5)^2 + (y - 0.5)^2) / 0.02)

and boundary conditions given by

    u(x, y) = 0 for x = 0 or x = 1
"""

__author__ = "Chris Greene"
__date__ = "2009-08-18"
__copyright__ = "Copyright (C) 2007-2008 Anders Logg"
__license__  = "GNU LGPL Version 2.1"

from dolfin import *
import montecarlo_mockup as mc
import numpy as np
import dolfin_util as du
import time
import mcoptions,sys
import trianglemesh as tm
import triangle
import meshtest
options = mcoptions.get_options()
#dolfin_set("linear algebra backend","Epetra"

def custom_func(mesh,V):
	f = Function(V)
#	for p in particles:
#		du.alter_cellid(mesh,f,p.id,p.charge)
	return f

# Create mesh and define function space
thetriangle = np.array([[-.5,-.288675],[.5,-.288675],[0.,.57735]])
mesh = mc.ParticleMesh(meshtest.TestMesh())
mesh.populate_regions(lambda x: 
			triangle.point_in_triangle(x,thetriangle*.5),
		      0,0)
plot(mesh)
V = FunctionSpace(mesh, "CG", 2)

# Define Dirichlet boundary (x = 0 or x = 1)
class InnerTriangle(SubDomain):
    def inside(self, x, on_boundary):
        if on_boundary:
		if triangle.point_in_triangle(x,thetriangle*.5):
			return True
	return False

class OuterTriangle(SubDomain):
    def inside(self, x, on_boundary):
        if on_boundary:
		if not triangle.point_in_triangle(x,thetriangle*.5):
			return True
	return False

# Define boundary condition
u0 = Constant(mesh, 0.0)
u1 = Constant(mesh, options.V)
bc0 = DirichletBC(V, u0, InnerTriangle())
bc1 = DirichletBC(V, u1, OuterTriangle())
bcs = [bc0,bc1]

# Define variational problem
v = TestFunction(V)
u = TrialFunction(V)
#init particles
#electrons, holes
print "adding electrons to regions"
mc.init_electrons(10,mesh.n_region.keys(),charge=-10,mesh=mesh)
mc.init_electrons(10,mesh.p_region.keys(),charge=10,mesh=mesh)
print "Creating density functions"
f = custom_func(mesh,V)
g = Function(V)
g.vector().set(f.vector().array())

#init Files
print "Creating Files"
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
	problem = VariationalProblem(a, L, bcs)
	sol = problem.solve()
	return sol

current_values = []

print "Beginning Simulation"

for x in range(options.num):
	g.vector().set(avg_dens.func)
	sol = PoissonSolve(g)
	# Plot solution
	file << sol
	dfile << f
	adfile << g
	print "Starting Step ",x
	start = time.time()
	electric_field = mc.negGradient(mesh,sol)
	gradfile << electric_field
	mc.MonteCarlo(mesh,sol,electric_field,f,avg_dens,current_values)
	print "Took: ",time.time()-start
#	plot(f)
file << sol
dfile << f
print current_values

#dump average
f.vector().set(avg_dens.func)
adfile << f
avgE=mc.negGradient(mesh,PoissonSolve(f))
avggradfile << avgE
# Hold plot
plot(avgE)
interactive()
