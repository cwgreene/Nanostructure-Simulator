import random as rd
import dolfin_util as du
from dolfin import *
from numpy import *
import itertools as it
import time
#import triangle
import sys

#import driftscatter
import stats
import convexhull
import materials
import constants

#more path
sys.path.append("c_optimized/")
import kdtree_c
import move_particles_c

p_count = it.count()

class ParticleMesh(Mesh):
	n_carrier_charge = -10
	p_carrier_charge = 10
	carrier_charge = 10
	def __init__(self,mesh,scale,length,time,gen_num,c_mesh=None):
		print "Hello!"
		Mesh.__init__(self,mesh) #this is why we needed to use it I think
		self.bd = du.boundary_dict(mesh)
		self.ibd = du.boundary_id_dict(mesh,self.bd)
		self.point_index = {}
#		self.particles_point = {}
		self.p_region = {}
		self.n_region = {}
		self.particles = []
		self.trajectories = {}
		self.scale = scale
		self.c_mesh = None
		self.super_particles_count = 0
		self.dim=self.topology().dim()
		if self.dim == 2:
			self.kdt = kdtree_c.new_kdtree(mesh.coordinates())
		else:
			self.kdt = kdtree_c.new_kdtree3(mesh.coordinates())
		#self.bandstructure = bandstructure.ParabolicBand(self) #should be part of the material
		self.material = {}
		self.length_scale = length
		self.dt = time
		self.charge_particle = constants.eC
		self.gen_num= gen_num

		#init point->index map, point->particle map
		for x,id in it.izip(mesh.coordinates(),it.count()):
			self.point_index[tuple(x)] = id
#			self.particles_point[id] = []
	def populate_regions(self,p_region_func,doping_p,doping_n,
				p_material,n_material):
		self.in_p_region = p_region_func
		def n_region_func(x):
			return not p_region_func
		self.in_n_region = n_region_func
		count = 0
		self.super_particles_count = (n_material.doping*
					     (self.length_scale**self.dim)/
						self.numCells())
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

particle_pos = []
class Particle():
	def __init__(self,pos,momentum,dx,lifetime,charge,mesh):
		self.uid = p_count.next()
		self.pos = array(pos)
		#mesh.trajectories[self.uid] = [tuple(self.pos)]
		#print mesh.trajectories[0]
		self.momentum = array(momentum)*mesh.scale
#		self.dx = array(dx)
#		self.lifetime = 2*rd.random()#int(lifetime)
		self.charge = charge
		self.dead = False
		self.mass = mesh.material[tuple(pos)].electron_mass #mesh.mass[charge]

		#mesh_id data
		#self.id = du.vert_index(mesh,self.pos)
		self.id = kdtree_c.find_point_id(mesh.kdt,self.pos)
		#self.meshpos = mesh.coordinates()[self.id]
		#particles_point must be update on move
#		mesh.particles_point[mesh.point_index[tuple(self.meshpos)]].append(self)
		self.part_id = len(mesh.particles)
		mesh.particles.append(self)

#functions
def random_momentum(mesh,pos):
	theta = rd.random()*2*pi
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
	start = time.time()
	remove_ids.sort()
	count = 0
	for id in remove_ids:
		p = full.pop(id-count)
		count += 1
	for id in xrange(len(full)):
		p = full[id]
		p.part_id = id
	remove_ids[:] = []
	stats.reap_time += time.time()-start
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
	if(density[id]*sign < 0):#too little
		for i in xrange(int(-density[id]/charge)):
			add_list.append(array(point))
	#remove excess
	exit_current = 0
	if(density[id]*sign > 0):#too much
		print "too much!"
	return exit_current


def photo_generate(mesh,density,holes,electrons):
	#these are the photogenerated electron hole pairs
	plot(mesh)
	coord = mesh.coordinates()
	for point in coord:
		holes.append(array(point))
		electrons.append(array(point))

def print_avg(name,value,count):
	print "Avg",name+":",value/count,count

def pre_compute_field(mesh,field):
	start = time.time()
	c_efield= []
	coord = mesh.coordinates()
	for index in xrange(len(coord)):
		pos = coord[index]
		vec = du.get_vec(mesh,field,pos)/mesh.length_scale
		c_efield.append(vec)
	print "Forces Calculated:",time.time()-start
	print max(c_efield,key=lambda x: x[0]**2+x[1]**2)
	return c_efield

def calculate_scaled_density(mesh,nextDensity):
	start = time.time()
	scaled_density = array(nextDensity.astype('double'))
	for point in mesh.coordinates():
		material = mesh.material[tuple(point)]
		#Q/eps=(particles*particle_charge*electrons_per_particle)*V/eps
		id = mesh.point_index[tuple(point)]
		spc = mesh.super_particles_count
		scaled_density[id] = spc*(nextDensity[id]*
				  constants.eC/mesh.gen_num*
				  (material.doping3d\
					*((mesh.length_scale)**3)
				  /material.epsilon))
		stats.avg_charge += abs(scaled_density[id])
	return scaled_density

def MonteCarlo(mesh,system,potential_field,electric_field,
		density_funcs,
		avg_dens,
		avg_electrons,
		avg_holes,
		current_values):
	reaper = []
	current = 0

	#next_step density function array

	nextDensity = density_funcs.combined_density.vector().array().astype('int')

	#holeDensity = density_funcs.holes.vector().array()
	#electronDensity = density_funcs.electrons.vector().array()

	#compute field, move_particles, calc current, recombinate, reap
	start = time.time()
	print "computing field"
	c_efield = pre_compute_field(mesh,electric_field)	#Evaluates field at each mesh point
	print "nextDensity max:",max(nextDensity)
	start =time.time()
	print "About to move particles"
	c_efield = array(c_efield)
	#print max(map(lambda x:max(abs(x[0]),abs(x[1])),c_efield)),\
	#	min(map(lambda x:min(abs(x[0]),abs(x[1])),c_efield))
	#raw_input()
	move_particles_c.move_particles(system,
					c_efield,nextDensity,
					 mesh.dt,mesh.length_scale)
	print "About to update"
	current = move_particles_c.update_density(system,
						  nextDensity,
						  mesh.kdt)
	print "Recombination"
#	move_particles_c.recombinate(system,nextDensity,mesh.kdt)
	print current
	print "time:",time.time()-start
	start = time.time()
#	current += move_particles_c.replenish(system,nextDensity)
	print "replenish:",time.time()-start
	#reap_list(mesh.particles,reaper)

	#update current values
	current_values.append(current*constants.eC/
				(mesh.dt*mesh.gen_num*mesh.numVertices()))

	#update avg_dens
	scaled_density = calculate_scaled_density(mesh,nextDensity)
	avg_dens.inc(scaled_density)
	density_funcs.combined_density.vector().set(
			nextDensity.astype('double'))
	density_funcs.scaled_density.vector().set(scaled_density)
	#print "doom",density_funcs.combined_density.vector().array()

	#print stats
	stats.print_stats(current)
