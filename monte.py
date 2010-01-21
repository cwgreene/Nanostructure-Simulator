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
import re
import photocurrent as pc

import materials
import meshes

class Problem:
	pass

class DolfinFiles:
	pass
class ResultsFile:
	pass

options = mcoptions.get_options()

def custom_func(mesh,V):
	f = Function(V)
	return f

def init_problem(mesh,V,V2,options):
	print "Initializing Probleming"
	problem = Problem()
	# Define boundary condition
	#for reasons I don't know, pBoundary needs to be 
	#kept globally
	pBoundary = Constant(mesh, options.V)
	nBoundary = Constant(mesh, 0.0) 
	mesh.V = options.V
	problem.boundaryFuncs = [pBoundary,nBoundary]#prevent bad garbage?
	problem.V2 = V2

	bc0 = DirichletBC(V, pBoundary, mesh.InnerBoundary)
	bc1 = DirichletBC(V, nBoundary, mesh.OuterBoundary)
	problem.bcs = [bc0,bc1]

	#init particles
	#electrons, holes
	print "adding electrons to regions"
	mc.init_electrons(options.particles,mesh.n_region.keys(),
				charge=-10,mesh=mesh)
	mc.init_electrons(options.particles,mesh.p_region.keys(),
				charge=10,mesh=mesh)

	print "Creating density functions"
	problem.f = custom_func(mesh,V)
	problem.g = Function(V)
	problem.scaled_density = Function(V)
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
	results_file.write("c:"+str(options.num))
	results_file.write(" scale:"+str(options.scale))
	results_file.write(" particles:"+str(options.particles))
	results_file.write(" size:"+str(options.size))
	results_file.write(" tag:"+str(options.tag))
	results_file.write("\n")
	return results_file
		

def PoissonSolve(density,bcs):
	print "Solving Poisson Equation"
	u = TrialFunction(V)
	v = TestFunction(V)
	a = dot(grad(v), grad(u))*dx
	L = v*(density)*dx
	# Compute solution
	problem = VariationalProblem(a, L, bcs)
	sol = problem.solve()
	return sol


def mainloop(mesh,problem,df,rf,scale):
	print "Beginning Simulation"
	current_values = []
	for x in range(options.num):
		#Solve equation
		start1 = time.time()
		problem.g.vector().set(problem.avg_dens.func)
		print problem.g.vector().array()
#		problem.g.vector().set(problem.g.vector().array()*scale)
		sol = PoissonSolve(problem.g,problem.bcs)

		#handle Monte Carlo
		print "Starting Step ",x
		electric_field = mc.negGradient(mesh,sol,problem.V2)
		df.gradfile << electric_field
		start2 = time.time()
		mc.MonteCarlo(mesh,sol,electric_field,
				problem.f,problem.avg_dens,current_values)
		end2 = time.time()
		#Report
		#Write Results
		df.file << sol
		df.dfile << problem.f
		df.adfile << problem.g

		#write current
		rf.current.write(str(current_values[-1]));
		rf.current.write("\n");rf.current.flush()

		end = time.time()
		print "Monte Took: ",end2-start2
		print "Loop Took:",end-start1
		#del electric_field
	#pc.generate_photo_current(mesh,problem.avg_dens)
	df.file << sol
	df.dfile << problem.f
	#dump average
	problem.f.vector().set(problem.avg_dens.func)
	for x in problem.avg_dens.func:
		rf.density.write(str(x)+" ")
	df.adfile << problem.f
	avgE=mc.negGradient(mesh,PoissonSolve(problem.f,problem.bcs),problem.V2)
	df.avggradfile << avgE

	print current_values

mesh = meshes.TriangleMesh(options,materials.Silicon(),materials.Silicon())
#these seem to need to be global
V = FunctionSpace(mesh, "CG", 2)
V2 = VectorFunctionSpace(mesh,"CG",1,2)
problem = init_problem(mesh,V,V2,options)
dolfinFiles = init_dolfin_files()
rf = ResultsFile()
rf.current = new_file("current")
rf.density = new_file("density")
mainloop(mesh,problem,dolfinFiles,rf,options.scale)
# Hold plot
#plot(avgE)
#interactive()
