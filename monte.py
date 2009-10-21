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

class Problem:
	pass

class DolfinFiles:
	pass
class ResultsFile:
	pass

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
#plot(mesh)

boundarymesh = BoundaryMesh(mesh)
print len(boundarymesh.coordinates())
#plot(mesh)

#these seem to need to be global
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
def init_problem(mesh,V):
	print "Initializing Probleming"
	problem = Problem()
	# Define boundary condition
	#for reasons I don't know, pBoundary needs to be 
	#kept globally
	pBoundary = Constant(mesh, options.V)
	nBoundary = Constant(mesh, 0.0) 
	problem.boundaryFuncs = [pBoundary,nBoundary]#prevent bad garbage?

	bc0 = DirichletBC(V, pBoundary, InnerTriangle())
	bc1 = DirichletBC(V, nBoundary, OuterTriangle())
	problem.bcs = [bc0,bc1]

	#init particles
	#electrons, holes
	print "adding electrons to regions"
	mc.init_electrons(10,mesh.n_region.keys(),charge=-10,mesh=mesh)
	mc.init_electrons(10,mesh.p_region.keys(),charge=10,mesh=mesh)

	print "Creating density functions"
	problem.f = custom_func(mesh,V)
	problem.g = Function(V)
	problem.g.vector().set(problem.f.vector().array())
	problem.avg_dens = mc.AverageFunc(problem.f.vector().array())
	return problem



def init_dolfin_files():
	#init Files
	print "Initializing Files"
	df = DolfinFiles()
	print "Creating Files"
	df.datadir = options.datadir
	df.file = File(df.datadir+"/poisson_attract.pvd")
	df.dfile = File(df.datadir+"/density_attract.pvd")
	df.adfile = File(df.datadir+"/avg_density.pvd")
	df.avfile = File(df.datadir+"/avg_voltage.pvd")
	df.gradfile = File(df.datadir+"/grad_force.pvd")
	df.avggradfile = File(df.datadir+"/avg_force.pvd")
	return df

#other files
import time
import re
def new_file(name):
	print "Creating results File"
	files = os.listdir("results")
	num=max([0]+map(int,re.findall("([0-9]+)"+name," ".join(files))))
	num += 1
	filename = ("results/"+str(num)+name+
			"_".join(map(str,time.gmtime())))
	print "Creating:",filename
	results_file = open(filename,"w")
	results_file.write(str(options.V)+"\n")
	results_file.write(str(options.num)+"\n")
	return results_file
		

def PoissonSolve(density,bcs):
	print "Solving Poisson Equation"
	u = TrialFunction(V)
	v = TestFunction(V)
	print "hi1"
	a = dot(grad(v), grad(u))*dx
	L = v*density*dx
	print "hi2"
	# Compute solution
	problem = VariationalProblem(a, L, bcs)
	print "hi3"
	sol = problem.solve()
	print "hi4"
	return sol


def mainloop(mesh,problem,df,rf):
	print "Beginning Simulation"
	current_values = []
	for x in range(options.num):
		#Solve equation
		start1 = time.time()
		problem.g.vector().set(problem.avg_dens.func)
		sol = PoissonSolve(problem.g,problem.bcs)

		#handle Monte Carlo
		print "Starting Step ",x
		start2 = time.time()
		electric_field = mc.negGradient(mesh,sol)
		df.gradfile << electric_field
		mc.MonteCarlo(mesh,sol,electric_field,
				problem.f,problem.avg_dens,current_values)

		#Report
		#Write Results
		df.file << sol
		df.dfile << problem.f
		df.adfile << problem.g

		#write current
		rf.current.write(str(current_values[-1]));
		rf.current.write("\n");rf.current.flush()

		end = time.time()
		print "Monte Took: ",end-start2
		print "Loop Took:",end-start1
	df.file << sol
	df.dfile << problem.f
	#dump average
	problem.f.vector().set(problem.avg_dens.func)
	for x in problem.avg_dens.func:
		rf.density.write(str(x)+" ")
	df.adfile << problem.f
	avgE=mc.negGradient(mesh,PoissonSolve(problem.f,problem.bcs))
	df.avggradfile << avgE

	print current_values

problem = init_problem(mesh,V)
dolfinFiles = init_dolfin_files()
rf = ResultsFile()
rf.current = new_file("current")
rf.density = new_file("density")
mainloop(mesh,problem,dolfinFiles,rf)
# Hold plot
#plot(avgE)
#interactive()
