import ctypes
import materials
import itertools as it
import numpy as np

import convexhull
lib = ctypes.cdll.LoadLibrary("c_optimized/move_particles.so")

lib.new_polygon.restype = ctypes.POINTER(ctypes.c_int)
lib.init_particle_list.restype = ctypes.POINTER(ctypes.c_int)
lib.init_dead_list.restype = ctypes.POINTER(ctypes.c_int)
lib.update_density.restype = ctypes.c_double
lib.replenish.restype = ctypes.c_double

def print_list(string,list):
	print string,list,len(list)

def move_particles(particles_ptr,#particles_pos,p_mass,p_charge,p_id,p_live,
			 c_efield,nextDensity,
			 dt,length_scale):
	lib.move_particles(particles_ptr,
			   c_efield.ctypes.data,
			   nextDensity.ctypes.data,
			   ctypes.c_double(dt),
			   ctypes.c_double(length_scale))

def update_density(system,nextDensity,kdt):
	return lib.update_density(system.particles.ptr,
			   nextDensity.ctypes.data,
			   system.bounding_polygon,
			   kdt)
def new_polygon(points):
	return lib.new_polygon(points.ctypes.data,len(points))
def init_particle_list(n):
	return lib.init_particle_list(n)
def init_dead_list(start,n):
	return lib.init_dead_list(start,n)

class System():
	pass

class CParticles():
	def __init__(self,nparticles):#particles,p_mass,p_charge,p_id,p_live):
		nparticles = nparticles
		self.pos = np.zeros((nparticles,4))
		self.p_id = np.zeros((nparticles,1),'int')
		self.p_charge = np.zeros((nparticles,1),'int')
		self.p_mass = np.zeros((nparticles,1))
		self.p_live = lib.init_particle_list(0)
		print "p_live",self.p_live
		self.p_dead = lib.init_dead_list(0,nparticles)
		self.ptr = lib.init_particles(self.pos.ctypes.data,
					      self.p_id.ctypes.data,
					      self.p_charge.ctypes.data,
					      self.p_mass.ctypes.data,
					      self.p_live,
					      self.p_dead)

def replenish(system,nextDensity):
	return lib.replenish(system.particles.ptr,	
			nextDensity.ctypes.data,
			system.c_mesh)

def init_system(mesh,nextDensity,particles_point):
	system = System()
	ntypepy = materials.Silicon()
	ptypepy = materials.Silicon()

	#create materials
	system.n_dist = ntypepy.random_momentum("ntype")
	system.p_dist = ntypepy.random_momentum("ptype")
	print system.n_dist
	ntype = lib.new_material(ctypes.c_double(ntypepy.electron_mass),
				  system.n_dist.ctypes.data,
				  len(system.n_dist))
	ptype = lib.new_material(ctypes.c_double(ptypepy.hole_mass),
				system.p_dist.ctypes.data,
				len(system.p_dist))
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
	print boundary_ids
	system.c_mesh = lib.create_mesh(mesh_coord.ctypes.data,
				    len(mesh_coord),
				    system.materials,
				    boundary_ids.ctypes.data,
				    len(boundary_ids),
				    ptype_ids.ctypes.data,
				    len(boundary_ids),
				    ntype_ids.ctypes.data,
				    len(ntype_ids))
	print "c_mesh",hex(system.c_mesh)
				    
	#create bounding polygon
	convex_hull = list(convexhull.convexHull(map(tuple,list(mesh.bd))))
	convex_hull.reverse()
	polygon = lib.new_polygon(np.array(convex_hull).ctypes.data,
					len(convex_hull))
	system.bounding_polygon = polygon
	print "I'm on a boat!"
	#create particles
	system.particles = CParticles(10**6)
	#initialize particles
	for i in xrange(len(mesh_coord)):
		if tuple(mesh_coord[i]) in mesh.p_region:
			sign = 1
			mass = ptypepy.hole_mass
		if tuple(mesh_coord[i]) in mesh.n_region:
			sign = -1
			mass = ntypepy.electron_mass
		for delta in range(particles_point):
			lib.create_particle(ctypes.c_int(i),
					system.particles.ptr,
					nextDensity.ctypes.data,
					ctypes.c_int(sign),
					ctypes.c_double(mass),
					system.c_mesh)
	for i in xrange(len(nextDensity)):
		nextDensity[i] = 0
	return system
