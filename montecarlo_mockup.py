import random as rd
import dolfin_util as du
from dolfin import *
from numpy import *
import itertools as it

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
	V = VectorFunctionSpace(mesh,"Lagrange",1,2)
	return project(grad(-field),V)

def MonteCarlo(mesh,potential_field,density_func,particles):
	#electrons = init_electrons()
	global total_force,force_count
	electric_field = negGradient(mesh,potential_field)
	reaper = []

	total_momentum = array([0.,0.])
	total_force = array([0.,0.])
	count = 0
	force_count = 0

	#next_step density function array
	nextDensity = density_func.vector().array()
	for p in particles:
		nextDensity[p.id] -= p.charge #remove from old location
		randomElectronMovement(p,electric_field,
					density_func,mesh)
		total_momentum += p.momentum
		count += 1
		if(du.out_of_bounds(mesh,p.pos)): #need to figure out exit
			reaper.append(p)
		else:
			p.id = du.vert_index(mesh,p.pos) #get new p.id
			p.meshpos = mesh.coordinates()[p.id] #lock to grid
			nextDensity[p.id] += p.charge
	for reaped in reaper:
		particles.remove(reaped)#already removed cell
	if count != 0:
		print "Avg momentum:",total_momentum/count,count
		print "Avg force:",total_force/force_count,force_count
		print "Avg sim force:",sim_total_force/sim_force_count,sim_force_count
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
