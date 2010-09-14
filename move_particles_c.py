import ctypes
import materials
import itertools as it
import numpy as np
import dolfin

import qhull

import convexhull
lib = ctypes.cdll.LoadLibrary("c_optimized/move_particles.so")

lib.new_polygon.argtype = [ctypes.POINTER(ctypes.c_double),ctypes.c_int]
lib.new_polygon.restype = ctypes.POINTER(ctypes.c_int)
lib.create_polytope2.restype = ctypes.POINTER(ctypes.c_int)
lib.create_polytope3.restype = ctypes.POINTER(ctypes.c_int)
lib.init_particle_list.argtype = [ctypes.c_int]
lib.init_particle_list.restype = ctypes.POINTER(ctypes.c_int)
lib.init_dead_list.argtype = [ctypes.c_int]
lib.init_dead_list.restype = ctypes.POINTER(ctypes.c_int)
lib.update_densityC.argtype = [ctypes.POINTER(ctypes.c_int),
			      ctypes.POINTER(ctypes.c_int),
			      ctypes.POINTER(ctypes.c_int)]
lib.update_densityC.restype = ctypes.c_double
lib.replenishC.restype = ctypes.c_double
lib.create_particleC.restype = ctypes.c_int
lib.photocurrentC.restype = ctypes.c_double

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
	return lib.update_densityC(system.particles.ptr,
			   system.c_mesh,
			   nextDensity.ctypes.data,
			   kdt)
def new_polygon(points):
	print "new_polygon"
	return lib.new_polytope2(points.ctypes.data,len(points))
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
		self.pos = np.zeros((nparticles,2*dim)) #should be pk_array
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
	return lib.replenishC(system.particles.ptr,	
			nextDensity.ctypes.data,
			system.c_mesh)

#the below is WAY too big.
#should be broken into multiple functions
#in fact, I should make a system.py file
def init_system(mesh,nextDensity,particles_point,length_scale):
	print "Creating system"
	dim = mesh.geometry().dim()
	system = System()

	#this will need to be fixed.
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

	#Construct polytope from points
	if dim==2:
		#get inner and outer bounding arrays
		#at the moment, boundaries must be convex.
		bm = dolfin.BoundaryMesh(mesh);
		bmc = bm.coordinates().copy()
		inner,outer = 	map(np.array,
				 qhull.inner_outer(bmc))
		print "inner",inner
		print "outer",outer
		inner_p = lib.create_polytope2(inner.ctypes.data,len(inner))
		outer_p = lib.create_polytope2(outer.ctypes.data,len(outer))
	elif dim==3:
		qfaces = qhull.qhull_faces(mesh.coordinates())
		faces = []
		for face in qfaces:
			face = np.array(face)
			print face
			faces.append(lib.create_face(face.ctypes.data,
							len(face),3))
		raw_input()
		faces_p = (ctypes.void_p*len(faces))(*faces)
		inner_p = 0 #should be NULL
		outer_p = lib.create_polytope3(faces_p,len(faces))
	else:
		raise("Invalid dimension.")
	print "Polytopes",inner_p,outer_p

	super_particles_count = ntypepy.doping*(length_scale**dim)/mesh.numCells()
					    
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
			    mesh.gen_num,
			    ctypes.c_double(super_particles_count),
			    outer_p,inner_p)
	mesh.c_mesh = system.c_mesh
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
			index = lib.create_particleC(ctypes.c_int(i),
					system.particles.ptr,
					nextDensity.ctypes.data,
					ctypes.c_int(sign),
					ctypes.c_double(mass),
					system.c_mesh)
			#print hex(system.particles.pos.ctypes.data)
			#print system.particles.pos[index]
			#raw_input()
	for i in xrange(len(nextDensity)):
		nextDensity[i] = 0
	print "Created system"
	return system

def recombinate(system,nextDensity,mesh):
	print "doom"
	lib.recombinateC(system.particles.ptr,
			nextDensity.ctypes.data,
			system.c_mesh)
	print "don edoom"
