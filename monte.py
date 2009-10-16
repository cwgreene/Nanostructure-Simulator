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
import mcoptions,sys,os
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
#thetriangle = np.array([[-.5,-.288675],[.5,-.288675],[0.,.57735]])
thetriangle = np.array([[0.,0.],[1.,0.],[.5,.8660254]])

mesh = meshtest.TestMesh()
mesh.refine()
mesh = mc.ParticleMesh(mesh)
#mesh = mc.ParticleMesh(tm.innertriangle(5,.2,thetriangle))

inner = triangle.scale_triangle(thetriangle,.52)
mesh.populate_regions(lambda x: triangle.point_in_triangle(x,inner), 0,0)
plot(mesh)

boundarymesh = BoundaryMesh(mesh)
print len(boundarymesh.coordinates())
plot(mesh)
V = FunctionSpace(mesh, "CG", 2)

# Define Dirichlet boundary (x = 0 or x = 1)
repeated = {}
class InnerTriangle(SubDomain):
    def inside(self, x, on_boundary):
        if on_boundary:
		#if (x[0],x[1]) in repeated:
			#print "repeat:",(x[0],x[1])
		repeated[(x[0],x[1])] = 0
		if mesh.in_p_region(x):
			#print "inner:",x[0],x[1]
			return True
	return False

class OuterTriangle(SubDomain):
    def inside(self, x, on_boundary):
        if on_boundary:
		if not mesh.in_p_region(x):
			#print "outer",x[0],x[1]
			return True
	return False

# Define boundary condition
pBoundary = Constant(mesh, options.V)
nBoundary = Constant(mesh, 0.0) 
bc0 = DirichletBC(V, pBoundary, InnerTriangle())
print "Doom@"
repeated_particles = {}
bc1 = DirichletBC(V, nBoundary, OuterTriangle())
print "Doom@"
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
datadir = options.datadir
file = File(datadir+"/poisson_attract.pvd")
dfile = File(datadir+"/density_attract.pvd")
adfile = File(datadir+"/avg_density.pvd")
avfile = File(datadir+"/avg_voltage.pvd")
gradfile = File(datadir+"/grad_force.pvd")
avggradfile = File(datadir+"/avg_force.pvd")
avg_dens = mc.AverageFunc(f.vector().array())

#other files
import time
import re
def new_results_file():
	files = os.listdir("results")
	num=max([0]+map(int,re.findall("([0-9]+)results"," ".join(files))))
	num += 1
	filename = ("results/"+str(num)+"results"+
			"_".join(map(str,time.gmtime())))
	print "Creating:",filename
	results_file = open(filename,"w")
	results_file.write(str(options.V)+"\n")
	results_file.write(str(options.num)+"\n")
	return results_file
		

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

print "Creating results File"
rf = new_results_file()

print "Beginning Simulation"

for x in range(options.num):
	#Solve equation
	start1 = time.time()
	g.vector().set(avg_dens.func)
	sol = PoissonSolve(g)

	#handle Monte Carlo
	print "Starting Step ",x
	start2 = time.time()
	electric_field = mc.negGradient(mesh,sol)
	gradfile << electric_field
	mc.MonteCarlo(mesh,sol,electric_field,f,avg_dens,current_values)

	#Report
	#Write Results
	file << sol
	dfile << f
	adfile << g
	rf.write(str(current_values[-1]));rf.write("\n");rf.flush()

	end = time.time()
	print "Monte Took: ",end-start2
	print "Loop Took:",end-start1
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
