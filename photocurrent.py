import driftscatter as ds
import montecarlo_mockup as mc
import dolfin_util as du
import kdtree_c
import time
import numpy

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
	for index in numpy.random.random_integers(0,len(list),num):
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
	

def generate_photo_current(mesh,e_field):
	current = 0
	e_field = e_field.func
	tot = 0.
	for point in random_list(mesh.coordinates(),100):
		clean_mesh(mesh)
		mc.init_electrons(100,[point],-10,mesh)
		start = time.time()
		for x in xrange(100):
			current += run_simulation(mesh,e_field)
		print "run",time.time()-start
	print "current:",current
