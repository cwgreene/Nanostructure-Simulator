import random as rd
import dolfin_util as du
from dolfin import *
from numpy import *
import itertools as it
import time
#import triangle
import sys
import driftscatter
import bandstructure
import stats

#more path
sys.path.append("c_optimized/")
import kdtree_c
import materials
import constants

class ParticleMesh(Mesh):
	n_carrier_charge = -10
	p_carrier_charge = 10
	carrier_charge = 10
	def __init__(self,mesh,scale,length,time):
		Mesh.__init__(self,mesh)
		self.bd = du.boundary_dict(mesh)
		self.point_index = {}
		self.particles_point = {}
		self.p_region = {}
		self.n_region = {}
		self.particles = []
		self.scale = scale
		self.kdt = kdtree_c.new_kdtree(mesh.coordinates())
		self.bandstructure = bandstructure.ParabolicBand(self)
		self.material = {}
		self.length_scale = length
		self.dt = time
		self.charge_particle = constants.eC

		#init point->index map, point->particle map
		for x,id in it.izip(mesh.coordinates(),it.count()):
			self.point_index[tuple(x)] = id
			self.particles_point[id] = []
	def populate_regions(self,p_region_func,doping_p,doping_n,
				p_material,n_material):
		self.in_p_region = p_region_func
		def n_region_func(x):
			return not p_region_func
		self.in_n_region = n_region_func
		count = 0
		for x in self.coordinates():
			if(p_region_func(x)):
				self.material[tuple(x)] = n_material
				self.p_region[tuple(x)] = self.p_carrier_charge
			else:
				self.n_region[tuple(x)] = self.n_carrier_charge
				self.material[tuple(x)] = p_material

class AverageFunc():
	def __init__(self,func):
		self.func = func
		self.count = 1
	def inc(self,func):
		self.count += 1
		self.func = (self.func*(self.count-1)+func)/(self.count*1.)

class Particle():
	def __init__(self,pos,momentum,dx,lifetime,charge,mesh):
		self.pos = array(pos)
		self.momentum = array(momentum)*mesh.scale
		self.dx = array(dx)
		self.lifetime = 2*rd.random()#int(lifetime)
		self.charge = charge
		self.dead = False
		self.mass = mesh.material[tuple(pos)].electron_mass #mesh.mass[charge]

		#mesh_id data
		#self.id = du.vert_index(mesh,self.pos)
		self.id = kdtree_c.find_point_id(mesh.kdt,self.pos)
		self.meshpos = mesh.coordinates()[self.id]
		#particles_point must be update on move
		mesh.particles_point[mesh.point_index[tuple(self.meshpos)]].append(self)
		self.part_id = len(mesh.particles)
		mesh.particles.append(self)

#functions
def random_momentum(mesh,pos):
	theta = rd.random()*2*pi
#	energy = mesh.bandstructure.random_energy()
	material = mesh.material[tuple(pos)]
	bandgap = material.bandgap
	mass = material.electron_mass
	energy = (constants.kbT*
			random.exponential(bandgap/constants.kbT))
	magnitude = sqrt(energy*2*mass)
	#magnitude = random.exponential(V)
	return magnitude*cos(theta),magnitude*sin(theta)

def init_electrons(num,points,charge=-1,mesh=None):
	electrons = []
	for point in points:
		for i in xrange(num):
			if mesh.in_p_region(point):
				V = mesh.V
			else:
				V = 0
			dx = array([0.,0.])
			lifetime = 0
			electrons.append(Particle(point,
					 random_momentum(mesh,point),
					 dx,lifetime,charge,mesh))
	return electrons

def negGradient(mesh,field,V):
	return project(grad(-field),V)

def reap_list(full,remove_ids):
	#global avg_lifetime,lifetime_count
	remove_ids.sort()
	count = 0
	for id in remove_ids:
		p = full.pop(id-count)
		count += 1
	#	avg_lifetime += p.lifetime
	#	lifetime_count += 1
	for id in xrange(len(full)):
		p = full[id]
		p.part_id = id
def handle_region(mesh,density,point,add_list,reaper,sign,id):
	charge = mesh.carrier_charge*sign
	
	#TODO: This is wrong if density is not an integer.
	#need to fix this. seperate out holes from particles
	#add them together and scale appropriately.
	#Status: Fixed

	#create to balance
	#TODO: This is wrong. There are two different
	#phenomenon going on here. Imbalance due
	#to leaving particles, and imbalance
	#due to holes being present.
	#Status: The above TODO is believed to be wrong.
	if(density[id]*sign < 0):
		for i in xrange(int(-density[id]/charge)):
			add_list.append(array(point))
	#remove excess
	exit_current = 0
	if(density[id]*sign > 0):
		for i in xrange(int(density[id]/charge)):
			doom_particle = mesh.particles_point[id].pop()
			density[id] -= charge
			reaper.append(doom_particle.part_id)
			exit_current += current_exit(doom_particle,mesh)
	return exit_current

def replenish_boundary(mesh,density,holes,electrons,reaper):
	#these are the thermal equilibrium holes and
	#electrons provided by the contacts
	bmesh = BoundaryMesh(mesh)
	boundary = bmesh.coordinates()

	current = 0
	for point in boundary:
		id = mesh.point_index[tuple(point)]
		if mesh.in_p_region(point):
			current += handle_region(mesh,density,point,
					holes,reaper,1,id)
		else:
			current += handle_region(mesh,density,point,
					electrons,reaper,-1,id)
	return current

def photo_generate(mesh,density,holes,electrons):
	#these are the photogenerated electron hole pairs
	plot(mesh)
	coord = mesh.coordinates()
	for point in coord:
		holes.append(array(point))
		electrons.append(array(point))


def replenish(mesh,density,boundary,reaper):
	print "boundary",len(boundary)
	print "prior particles",len(mesh.particles)
	
	holes = []
	electrons = []	
	
	#currents
	exit_current = 0
	enter_current = 0
	
	exit_current=replenish_boundary(mesh,density,holes,electrons,reaper)
#	photo_generate(mesh,density,holes,electrons)

	new_particles = init_electrons(1,electrons,-10,mesh)
	new_particles += init_electrons(1,holes,10,mesh)
	for p in new_particles:
		density[p.id] += p.charge
		enter_current += -current_exit(p,mesh)
	print "parts",len(new_particles)
	print "total particles",len(mesh.particles)
	return enter_current+exit_current

def print_avg(name,value,count):
	print "Avg",name+":",value/count,count

def current_exit(particle,mesh):
	boundary = mesh.bd

	speed = sqrt(dot(particle.pos,particle.pos))
	#exit = du.closest_exit(boundary,particle.pos)
	exit = mesh.coordinates()[kdtree_c.find_point_id(mesh.kdt,particle.pos)]
	if mesh.in_p_region(exit):
		return particle.charge*-1*speed
	else:
		return particle.charge*speed

def recombinate(mesh,reaper,avg_dens,avg_electrons,avg_holes):
	#for point_id in mesh.point_index.values():
		#electron = None
		#hole = None
		#for p in mesh.particles_point[point_id]:
		#	if p.dead != True:
		#		if p.charge < 0:
		#			electron = p
		#		if p.charge > 0:
		#			hole = p
		#		if electron != None and hole != None:
		#			reaper.append(electron)
		#			reaper.append(hole)
		#			break
	avg = avg_dens.func
	for p in mesh.particles:
		#if we're in a  hole region,
		#recombinate
		if avg[p.id]*p.charge < 0: 
			p.dead = True
			reaper.append(p.part_id) 
				

def pre_compute_field(mesh,field):
	start = time.time()
	c_efield= []
	coord = mesh.coordinates()
	for index in xrange(len(coord)):
		pos = coord[index]
		vec = du.get_vec(mesh,field,pos)
		c_efield.append(vec)
	print "Forces Calculated:",time.time()-start
	return c_efield

def move_particles(mesh,c_efield,nextDensity,density_func,reaper):
	rem_time = 0.
	start=time.time()
	for index in xrange(len(mesh.particles)):
		p = mesh.particles[index]
		
		#remove from old locations
		nextDensity[p.id] -= p.charge
		mesh.particles_point[p.id].remove(p)

		#begin movement	
		start2 = time.time()
		driftscatter.randomElectronMovement(p,c_efield,
					density_func,mesh,reaper)
		rem_time += time.time()-start2
	print "drift scatter took:",rem_time
	print "Particle Movement took:",time.time()-start

def update_density(mesh,reaper,nextDensity,current):
	start = time.time()
	for index in xrange(len(mesh.particles)):
		p = mesh.particles[index]
		if p.dead == False: #if we didn't kill it.
			start2 = time.time()
			if(du.out_of_bounds(mesh,p.pos)): 
				#need to figure out exit
				reaper.append(index)
				p.dead = True
				current += current_exit(p,mesh)
			else:
				#get new p.id
				p.id = kdtree_c.find_point_id(mesh.kdt,p.pos)
				#lock to grid
				p.meshpos = mesh.coordinates()[p.id]
				mesh.particles_point[p.id].append(p)
				#associate charge with density func
				nextDensity[p.id] += p.charge
	print "Lookup time:",time.time()-start
	return current

def MonteCarlo(mesh,potential_field,electric_field,
		density_func,
		avg_dens,
		avg_electrons,
		avg_holes,
		current_values):
	reaper = []
	current = 0

	#next_step density function array
	nextDensity = density_func.vector().array()
	mesh_lookup_time = 0.

	#precompute electric field at all points
	c_efield = pre_compute_field(mesh,electric_field)
	move_particles(mesh,c_efield,nextDensity,density_func,reaper)
	current =update_density(mesh,reaper,nextDensity,current)
	recombinate(mesh,reaper,avg_dens,avg_electrons,avg_holes)

	#reap
	start = time.time()
	reap_list(mesh.particles,reaper)
	reap_time = time.time()-start
	reaper = []

	#replenish, with reaper		
	start = time.time()
	current += replenish(mesh,nextDensity,mesh.bd,reaper)
	replenish_time = time.time()-start

	#reap again
	start = time.time()
	reap_list(mesh.particles,reaper)
	reap_time += time.time()-start
	reaper = []
		
	current_values.append(current*constants.eC/mesh.dt)
	print "Current:",current
	
	print "Reaper:",reap_time
	print "Replenish took:",replenish_time
	print "Mesh lookup time:",mesh_lookup_time
	print "Average Momentum:",stats.avg_momentum
	print "Average Force:",stats.avg_force
	print "Average dx:",stats.avg_dx
	print "Average charge:",stats.avg_charge
	start = time.time()
	scaled_density = array(nextDensity)
	for point in mesh.coordinates():
		material = mesh.material[tuple(point)]
		#Q/eps=(particles*particle_charge*electrons_per_particle)*V/eps
		id = mesh.point_index[tuple(point)]
		scaled_density[id] = (nextDensity[id]*
				  constants.eC*
				  (material.doping3d*
					((mesh.length_scale)**2)
				  /material.epsilon))
		stats.avg_charge += abs(scaled_density[id])
	avg_dens.inc(scaled_density)
	density_func.vector().set(nextDensity)
	print "density increment time:",time.time()-start
