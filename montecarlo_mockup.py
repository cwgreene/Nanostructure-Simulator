import random as rd
import dolfin_util as du
from dolfin import *
from numpy import *
import itertools as it
import time

class AverageFunc():
	def __init__(self,func):
		self.func = func
		self.count = 1
	def inc(self,func):
		self.count += 1
		self.func += func/self.count

class Particle():
	def __init__(self,pos,momentum,dx,lifetime,charge,mesh):
		self.pos = array(pos)
		self.momentum = array(momentum)
		self.dx = array(dx)
		self.lifetime = int(lifetime)
		self.charge = charge
		self.id = du.vert_index(mesh,self.pos)
		self.meshpos = mesh.coordinates()[self.id]

#globals
sim_total_force = array([0.,0.])
sim_force_count = 0
total_force = array([0.,0.])
force_count = 0
self_force = array([0.,0.])

#functions
def random_momentum():
	return (rd.random()-.5)*2,(rd.random()-.5)*2

def init_electrons(num,points,charge=-1,mesh=None):
	electrons = []
	for point in points:
		for i in xrange(num):
			dx = array([0,0])
			lifetime = 0
			electrons.append(Particle(point,random_momentum(),
					 dx,lifetime,charge,mesh))
	return electrons

def negGradient(mesh,field):
	V = VectorFunctionSpace(mesh,"CG",1,2)
	return project(grad(-field),V)

def reap_list(full,remove_ids):
	remove_ids.sort()
	count = 0
	for id in remove_ids:
		full.pop(id-count)
		count += 1

def replenish(mesh,density,boundary,particles):
	print "boundary",len(boundary)
	print "prior particles",len(particles)
	holes = []
	electrons = []
	for point in boundary:
		if point[0] < .5:
			holes.append(point)
		else:
			electrons.append(point)
	parts = init_electrons(1,electrons,-10,mesh)
	parts = init_electrons(1,holes,10,mesh)
	for p in parts:
		density[p.id] += p.charge
	for x in parts:
		particles.append(x)
	print "parts",len(parts)
	print "total particles",len(particles)

def MonteCarlo(mesh,potential_field,density_func,particles,avg_dens):
	#electrons = init_electrons()
	global total_force,force_count
	electric_field = negGradient(mesh,potential_field)
	#plot(electric_field)
	reaper = []

	total_momentum = array([0.,0.])
	total_force = array([0.,0.])
	count = 0
	force_count = 0

	#next_step density function array
	nextDensity = density_func.vector().array()
	start = time.time()
	bd = du.boundary_dict(mesh)
	for index in xrange(len(particles)):
		p = particles[index]
		nextDensity[p.id] -= p.charge #remove from old location
		randomElectronMovement(p,electric_field,
					density_func,mesh)
		total_momentum += p.momentum
		count += 1
		if(du.out_of_bounds(mesh,p.pos)): #need to figure out exit
			reaper.append(index)
			#print du.closest_exit(bd,p.pos)
		else:
			p.id = du.vert_index(mesh,p.pos) #get new p.id
			p.meshpos = mesh.coordinates()[p.id] #lock to grid
			nextDensity[p.id] += p.charge
	print (time.time()-start),len(reaper)
	start = time.time()
	reap_list(particles,reaper)
	print (time.time()-start)
	if count != 0:
		print "Avg momentum:",total_momentum/count,count
		print "Avg force:",total_force/force_count,force_count
		print "Avg sim force:",sim_total_force/sim_force_count,sim_force_count
	replenish(mesh,nextDensity,bd,particles)
	avg_dens.inc(nextDensity)
	density_func.vector().set(nextDensity)

def randomElectronMovement(particle,electric_field,density_func,mesh):
	meanpathlength = 1#getMeanPathLength(cell)
	lifetime = 1.#lifetime(cell)
	mass_particle = 1
	
	p = particle
	
	dt = lifetime/100.
	p.momentum += 100*drift(mesh,electric_field,p)*dt
	p.pos += p.momentum*dt/mass_particle
	p.dx += p.momentum*dt/mass_particle
	#check for out of bounds
	#	this needs to be changed, we should
	#	do all momentum changes BEFORE movement
	#	this guarantees that we are on the mesh.
	#if(dot(e.dx,e.dx) > meanpathlength**2):
	#	e.dx = array([0.,0.])
	#	e.momentum += array([0,0])#scatter(e.momentum,e.pos,mesh)
	#if(e.lifetime > lifetime):
	#	e.pos = array([.5,.5]) #Send back to middle
	#	e.lifetime = 0
	#print e.momentum
	p.lifetime += dt

#drift calculates drift force due to forces
#F = dp/dt
def drift(mesh,func,particle):
	global total_force,force_count
	global sim_total_force,sim_force_count
	global self_force
	p = particle
	force = du.get_vec(mesh,func,p.meshpos)*100*p.charge-self_force
	total_force += force
 	force_count += 1
	sim_force_count += 1
	sim_total_force += force
	return force

def scatter(momentum,pos,mesh):
	force = momentum*random_vec()
