import numpy as np
import distribution as dist
import constants
import itertools as it
import math
import pickle
import os

class Silicon:
	eV = 1.60217646*10**-19
	bandgap = 1.12 * eV
	free_mass = ( 9.10938188 *10**-31)

	electron_mass = 1.08*free_mass
#	hole_mass = 1.08*free_mass #test
	hole_mass = .58*free_mass
	dielectric = 11.7
	epsilon = dielectric*8.85418782*10**-12
	doping = 10.**24
	doping3d = 10.**24
	doping2d = 10.**22

	def random_momentum(self,mtype):
		filename="materials/momentum_tables/silicon_"+mtype+"array"
		if os.path.exists(filename):
			retarray= np.fromfile(open(filename))
			return retarray
		print "Generating",mtype,"momentum table"
		kbT = constants.kbT
		if mtype == "ptype":
			m = self.hole_mass
		else:
			m = self.electron_mass
		energy = dist.distribution(lambda x:(1/kbT)*math.exp(-x/(kbT)),
					0,kbT*10,boxes=100.,ddt=10**6.)
		velocity = map(lambda x: math.sqrt(x*2*m),energy)
		angles = np.arange(0,2*np.pi,.01)
		retarray = np.ndarray((len(velocity)*len(angles),2))
		index = it.count()
		for v in velocity:
			for theta in angles:
				i = index.next()
				retarray[i][0] = v*math.cos(theta)
				retarray[i][1] = v*math.sin(theta)
		retarray.tofile(open(filename,"w"))
		print "generated"
		return retarray
	
