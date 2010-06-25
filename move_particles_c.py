import ctypes
import materials
import itertools as it
import numpy as np

import convexhull
lib = ctypes.cdll.LoadLibrary("c_optimized/move_particles.so")

lib.new_polygon.argtype = [ctypes.POINTER(ctypes.c_double),ctypes.c_int]
lib.new_polygon.restype = ctypes.POINTER(ctypes.c_int)
lib.init_particle_list.argtype = [ctypes.c_int]
lib.init_particle_list.restype = ctypes.POINTER(ctypes.c_int)
lib.init_dead_list.argtype = [ctypes.c_int]
lib.init_dead_list.restype = ctypes.POINTER(ctypes.c_int)
lib.update_densityC.argtype = [ctypes.POINTER(ctypes.c_int),
			      ctypes.POINTER(ctypes.c_int),
			      ctypes.POINTER(ctypes.c_int),
			      ctypes.POINTER(ctypes.c_int)]
lib.update_densityC.restype = ctypes.c_double
lib.replenishC.restype = ctypes.c_double

def print_list(string,list):
	print string,list,len(list)

def move_particles(system,
		   c_efield,nextDensity,
		   dt,length_scale):
	print "move_particles"
	lib.move_particlesC(system.particles.ptr,
			   c_efield.ctypes.data,
			   nextDensity.ctypes.data,
			   ctypes.c_double(dt),
			   ctypes.c_double(length_scale),
			   system.c_mesh);

def update_density(system,nextDensity,kdt):
	print "update_density"
	return lib.update_density2(system.particles.ptr,
			   system.c_mesh,
			   nextDensity.ctypes.data,
			   system.bounding_polygon,
			   kdt)
def new_polygon(points):
	print "new_polygon"
	return lib.new_polygon(points.ctypes.data,len(points))
def init_particle_list(n):
	print "init_particle_list"
	return lib.init_particle_list(n)
def init_dead_list(start,n):
	print "init_dead_list"
	return lib.init_dead_list(start,n)

class System():
	pass

class CParticles():
	def __init__(self,nparticles,c_mesh,dim):#particles,p_mass,p_charge,p_id,p_live):
		print "Creating CParticles:",nparticles,c_mesh,dim
		self.pos = np.zeros((nparticles,4))
		self.p_id = np.zeros((nparticles,1),'int')
		self.p_charge = np.zeros((nparticles,1),'int')
		self.p_mass = np.zeros((nparticles,1))
		self.p_live = lib.init_particle_list(0)
		print "p_live",self.p_live
		self.p_dead = lib.init_dead_list(0,nparticles)
		print "p_dead",self.p_dead
		self.ptr = lib.init_particles(self.pos.ctypes.data,
					      self.p_id.ctypes.data,
					      self.p_charge.ctypes.data,
					      self.p_mass.ctypes.data,
					      self.p_live,
					      self.p_dead,
					      dim)
		print "created particles"

def replenish(system,nextDensity):
	print "replenish"
	return lib.replenish(system.particles.ptr,	
			nextDensity.ctypes.data,
			system.c_mesh)

def init_system(mesh,nextDensity,particles_point):
	print "Creating system"
	dim = mesh.geometry().dim()
	system = System()
	ntypepy = materials.Silicon()
	ptypepy = materials.Silicon()

	#create materials
	system.n_dist = ntypepy.random_momentum("ntype")
	system.p_dist = ntypepy.random_momentum("ptype")
	print system.n_dist
	ntype = lib.new_material(ctypes.c_double(ntypepy.electron_mass),
				  system.n_dist.ctypes.data,
				  len(system.n_dist),dim)
	ptype = lib.new_material(ctypes.c_double(ptypepy.hole_mass),
				system.p_dist.ctypes.data,
				len(system.p_dist),dim)
	system.materials=lib.material_pair(ptype,ntype)

	#create mesh
	#create materials **
	index = it.count()
	boundary_ids = []
	ntype_ids = []
	ptype_ids = []
	
	mesh_coord = mesh.coordinates()
	for x in mesh_coord:
		i = index.next()
		if tuple(x) in mesh.p_region:
			ptype_ids.append(i)
		if tuple(x) in mesh.n_region:
			ntype_ids.append(i)
		if i in mesh.ibd:
			boundary_ids.append(i)
	ptype_ids = np.array(ptype_ids)
	ntype_ids = np.array(ntype_ids)
	boundary_ids = np.array(boundary_ids)
	print ""
	system.c_mesh = lib.create_mesh(mesh_coord.ctypes.data,
				    len(mesh_coord),
				    dim,
				    system.materials,
				    boundary_ids.ctypes.data,
				    len(boundary_ids),
				    ptype_ids.ctypes.data,
				    len(ptype_ids),
				    ntype_ids.ctypes.data,
				    len(ntype_ids),
				    mesh.kdt,
				    ctypes.c_double(1.*10**18))
				    
	#create bounding polygon
	convex_hull = list(convexhull.convexHull(map(tuple,list(mesh.bd))))
	convex_hull.reverse() 
	convex_hull = np.array(convex_hull)
	polygon = lib.new_polygon(convex_hull.ctypes.data,
					ctypes.c_int(len(convex_hull)))
	system.bounding_polygon = polygon
	print "Creating a ton of particles..."
	#create particles
	system.particles = CParticles(5*10**6,system.c_mesh,dim)
	#initialize particles
	for i in xrange(len(mesh_coord)):
		if tuple(mesh_coord[i]) in mesh.p_region:
			sign = 1
			mass = ptypepy.hole_mass
		if tuple(mesh_coord[i]) in mesh.n_region:
			sign = -1
			mass = ntypepy.electron_mass
		for delta in range(particles_point):
			lib.create_particleC(ctypes.c_int(i),
					system.particles.ptr,
					nextDensity.ctypes.data,
					ctypes.c_int(sign),
					ctypes.c_double(mass),
					system.c_mesh)
	for i in xrange(len(nextDensity)):
		nextDensity[i] = 0
	print "Created system"
	return system

def recombinate(system,nextDensity,mesh):
	print "doom"
	lib.recombinateC(system.particles.ptr,nextDensity.ctypes.data,system.c_mesh)
	print "don edoom"
