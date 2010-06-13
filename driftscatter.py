import ctypes
import random
import math
import numpy

import stats 
import constants
#energy = mesh.bandstructure.random_
ext = ctypes.cdll.LoadLibrary("c_optimized/driftscatter.so")

def randomElectronMovement(particle,electric_field,density_func,mesh,reaper):
#	global avg_lifetime
#	lifetime = .001#lifetime(cell)
	p = particle
	dt = mesh.dt
	p.momentum += drift(mesh,electric_field,p)*dt
	dx = (p.momentum*dt/mesh.length_scale)/p.mass
	#stats.avg_dx += numpy.linalg.norm(dx)
	#stats.avg_momentum += numpy.linalg.norm(p.momentum)
	#p.dx += dx
	p.pos += dx
#	scatter(mesh,particle)
	#check for out of bounds
	#	this needs to be changed, we should
	#	do all momentum changes BEFORE movement
	#	this guarantees that we are on the mesh.
	#if(dot(e.dx,e.dx) > meanpathlength**2):
	#	e.dx = array([0.,0.])
	#	e.momentum += array([0,0])#scatter(e.momentum,e.pos,mesh)
	#if(p.lifetime < 0):
	#	p.dead = True
	#	reaper.append(p.part_id)
	#print e.momentum
	#p.lifetime -= dt

def scatter(mesh,particle):
	scatter_type = random.random()
	if scatter_type < .01:
		print "scatter"
		mag = numpy.sqrt(numpy.dot(particle.momentum,particle.momentum))
		theta = random.random()*2*numpy.pi
		particle.momentum = mag*numpy.cos(theta),mag*numpy.sin(theta)

def drift(mesh,func,particle):
	p = particle
	force = func[p.id]*(p.charge/10)*mesh.charge_particle/mesh.length_scale #self force?
	#stats.avg_force += numpy.linalg.norm(force)
	return force
