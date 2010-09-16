#import driftscatter as ds
import move_particles_c as mpc
import montecarlo_mockup as mc
import dolfin_util as du
import kdtree_c
import time
import numpy as np
import ctypes

#import bandstructure as bs

def reap_list(full,remove_ids):
	remove_ids.sort()
	count = 0
	for id in remove_ids:
		p = full.pop(id-count)
		count += 1
	for id in xrange(len(full)):
		p = full[id]
		p.part_id = id

def random_list(list,num):
	sublist = []
	for index in np.random.random_integers(0,len(list),num):
		sublist = list[index]
	return sublist

def run_simulation(mesh,e_field):
	reaper=[]
	current = 0
	for particle in mesh.particles:
		ds.randomElectronMovement(particle,e_field,
						None,mesh,reaper)
	for index in xrange(len(mesh.particles)):
		p = mesh.particles[index]
		if p.dead == False: #if we didn't kill it.
			if(du.out_of_bounds(mesh,p.pos)): 
				#need to figure out exit
				reaper.append(index)
				p.dead = True
				current += mc.current_exit(p,mesh)
			else:
				p.id = kdtree_c.find_point_id(mesh.kdt,p.pos)
	reap_list(mesh.particles,reaper)
	return current
	
def clean_mesh(mesh):
	mesh.particles = []
	for point in mesh.particles_point:
		point = []

def init_photo_pair(mesh,pos,energy):
	#create electron, and corresponding hole.
	momentum = mesh.bandstructure.random_momentum(mesh,pos,energy)
	
def new_particle(mesh_pos,particles,problem,charge_sign,mesh):
	coord = tuple(mesh.coordinates()[mesh_pos])
	material = mesh.material[coord]
	if charge_sign < 0:
		mass = material.electron_mass
	if charge_sign > 0:
		mass = material.hole_mass
	dim = mesh.geometry().dim()
	nextDensity = problem.density_funcs.combined_density.vector().array().astype('int')
	index =  mpc.lib.create_particleC(ctypes.c_int(mesh_pos),
					 particles.ptr,
					 nextDensity.ctypes.data,
					 ctypes.c_int(charge_sign),
					 ctypes.c_double(mass),
					 mesh.c_mesh)
	particles.pos[index][dim:2*dim]= [0.]*dim #bottom of gap

def photons_per_watt(wavelength_nm):
	#1 photon per X number of joules
	return 1/(constants.h*constants.c/wavelength_nm)

def generate_photo_current(mesh,e_field,problem):
	current = 0
	accumulated_charge = 0.
	total_photons = len(mesh.coordinates())
	particles = mpc.CParticles(2000,mesh.c_mesh,mesh.geometry().dim())
	e_field = np.array(mc.pre_compute_field(mesh,e_field))
	nextDensity = problem.density_funcs.combined_density.vector().array().astype('int')
	for point in xrange(len(mesh.coordinates())):
		new_particle(point,particles,problem,-1,mesh)
		new_particle(point,particles,problem,+1,mesh)
	for rep in xrange(300):
		accumulated_charge+=mpc.lib.photocurrentC(particles.ptr,
				nextDensity.ctypes.data,
				e_field.ctypes.data,
				mesh.c_mesh,
				ctypes.c_double(mesh.dt),
				ctypes.c_double(mesh.length_scale))
	power = 1000 #per meter squared
	photons_sec = power*photons_per_watt
	current = accumulated_charge*photons_sec
	return current
