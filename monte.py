from dolfin import *
import montecarlo_mockup as mc
import move_particles_c as c_interface
import numpy as np
import dolfin_util as du
import time
import mcoptions,sys,os
import re
import photocurrent as pc

import density_funcs
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
	problem.space = V
	# Define boundary condition
	#for reasons I don't know, pBoundary needs to be 
	#kept globally
	pBoundary = Constant(mesh, 0.0)
	nBoundary = Constant(mesh, options.V) 

	bc0 = DirichletBC(V, pBoundary, mesh.InnerBoundary)
	bc1 = DirichletBC(V, nBoundary, mesh.OuterBoundary)

	mesh.V = options.V
	problem.bcs = [bc0,bc1]#prevent bad garbage?
	problem.boundaryFuncs = [pBoundary,nBoundary]
	problem.V = V
	problem.V2 = V2

	#init particles
	#electrons, holes
	print "adding electrons to regions"
#	mc.init_electrons(options.gen_num,mesh.n_region.keys(),
#				charge=-10,mesh=mesh)
#	mc.init_electrons(options.gen_num,mesh.p_region.keys(),
#				charge=10,mesh=mesh)

	print "Creating density functions"
	problem.density_funcs = density_funcs.DensityFuncs()
	problem.density_funcs.holes = Function(V)
	problem.density_funcs.electrons = Function(V)
	problem.density_funcs.combined_density = Function(V)
	problem.density_funcs.poisson_density = Function(V)
	problem.density_funcs.scaled_density = Function(V)
	problem.density_funcs.poisson_density.vector().set(problem.density_funcs.combined_density.vector().array())

	print "Creating Average Densities"
	problem.avg_dens = mc.AverageFunc(problem.density_funcs.combined_density.vector().array())
	problem.avg_holes = mc.AverageFunc(problem.density_funcs.combined_density.vector().array())
	problem.avg_electrons = mc.AverageFunc(problem.density_funcs.combined_density.vector().array())
	return problem

#this will be replaced with 
#database initialization and connection
def init_dolfin_files():
	#init Files
	print "Initializing Files"
	df = DolfinFiles()
	print "Creating Files"
	df.datadir = options.datadir
	if not os.path.exists(df.datadir):
		os.mkdir(df.datadir)
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
	results_file.write(" particles:"+str(options.gen_num))
	results_file.write(" size:"+str(options.size))
	results_file.write(" tag:"+str(options.tag))
	results_file.write("\n")
	return results_file

def init_database():
	#connect to local database
	if os.path.exists("database/runs.db"):
		raise IOError("database/runs.db does not exist. You'll need to initialize it")
	runs = sqlite3.connect("database/runs.db")
	if runs == None:
		raise IOError("Failure to open database/runs.db")
	return runs
		

def PoissonSolve(mesh,density,bcs,V):
	print "Solving Poisson Equation"
	lengthr = Constant(mesh,1./mesh.length_scale)
	length = Constant(mesh,mesh.length_scale)
	u = TrialFunction(V)
	v = TestFunction(V)
	a = dot(grad(v), grad(u))*lengthr*dx
	L = v*(density)*length*dx
	# Compute solution
	problem = VariationalProblem(a, L, bcs)
	sol = problem.solve()
	return sol

def write_results(df,rf,problem,sol,electric_field,current_values):
	#Write Results
	df.file << sol
	df.dfile << problem.density_funcs.combined_density
	df.adfile << problem.density_funcs.poisson_density
	df.gradfile << electric_field #moved, might have destabilized it.

	#write current
	rf.current.write(str(current_values[-1]));
	rf.current.write("\n");rf.current.flush()

def final_record_files(df,rf,sol,problem,mesh):
	df.file << sol
	df.dfile << problem.density_funcs.combined_density
	#dump average
	problem.density_funcs.combined_density.vector().set(problem.avg_dens.func)
	for x in problem.avg_dens.func:
		rf.density.write(str(x)+" ")
	df.adfile << problem.density_funcs.combined_density
	avgE=mc.negGradient(mesh,PoissonSolve(mesh,
					problem.density_funcs.combined_density,
					problem.bcs,problem.V),
				problem.V2)
	df.avggradfile << avgE

def mainloop(mesh,system,problem,df,rf,scale):
	print "Beginning Simulation"
	current_values = []
	for x in range(options.num):
		start1 = time.time()

		#Solve equation using avg_dens
		sol = PoissonSolve(mesh,
				problem.density_funcs.scaled_density,
				problem.bcs,problem.space)
		#handle Monte Carlo
		print "Starting Step ",x
		electric_field = (mc.negGradient(mesh,sol,problem.V2))
		start2 = time.time()
		mc.MonteCarlo(mesh,system,sol,electric_field,
				problem.density_funcs,problem.avg_dens,
				problem.avg_electrons,problem.avg_holes,
				current_values)
		end2 = time.time()
		#Report
		write_results(df,rf,problem,
				sol,electric_field,current_values)
		end = time.time()
		print "Monte Took: ",end2-start2
		print "Loop Took:",end-start1

	#photocurrent
	current= pc.generate_photo_current(mesh,electric_field,problem)
	rf.current.write("pc: "+str(current)+"\n")

	final_record_files(df,rf,sol,problem,mesh)
	avg_length = 0
	for particle in mesh.trajectories:
		avg_length += len(mesh.trajectories[particle])
		#rf.trajectory.write(str(mesh.trajectories[particle]))
		#rf.trajectory.write("\n")
	print current_values
#	avg_length /= 1.*len(mesh.trajectories)
	print "Average trajectory length:",avg_length

def init_files():
	dolfinFiles = init_dolfin_files() 
	rf = ResultsFile()
	rf.current = new_file("current")
	rf.density = new_file("density")
	rf.trajectory = new_file("trajectory")
	return (dolfinFiles,rf)

#main
def main():
	#init mesh
	mesh = options.mesh(options,
			    options.materials[0],
			    options.materials[1])
	#these seem to need to be global
	V = FunctionSpace(mesh, "CG", 1)
	V2 = VectorFunctionSpace(mesh,"CG",1,2)
	problem = init_problem(mesh,V,V2,options)
	system = c_interface.init_system(mesh,
		  problem.density_funcs.poisson_density.vector().array(),
		  options.gen_num, options.length)

	#init Files
	(dolfinFiles,rf)=init_files()
	
	#mainloop
	mainloop(mesh,system,problem,dolfinFiles,rf,options.scale)

if __name__=="__main__":
	main()
